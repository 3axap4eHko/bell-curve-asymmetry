"""
Spin Hypothesis - Full Angle Sweep on IBM Quantum
==================================================
Prepares a singlet Bell state, sweeps Bob's measurement angle from 0 to 180 degrees
in steps, measures P_disagree at each angle, then fits Model 2 and Model 3
to extract the asymmetry parameter alpha with confidence intervals.

Usage:
  pip install qiskit qiskit-ibm-runtime scipy numpy matplotlib
  python spin_angle_sweep.py --token YOUR_IBM_TOKEN [--backend ibm_fez] [--shots 8192] [--dry-run]
"""

import argparse
import json
import time
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2, norm

# ── CLI ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--token",   required=False, help="IBM Quantum API token")
parser.add_argument("--backend", default="ibm_sherbrooke", help="Backend name")
parser.add_argument("--shots",   type=int, default=8192)
parser.add_argument("--dry-run", action="store_true", help="Simulate locally, no IBM submission")
parser.add_argument("--results", default=None, help="Path to saved JSON results to skip to analysis")
args = parser.parse_args()

# ── Angle grid ───────────────────────────────────────────────────────────────

ANGLES_DEG = np.arange(0, 185, 5)   # 0, 5, 10, ..., 180  (37 points)
ANGLES_RAD = np.deg2rad(ANGLES_DEG)

# ── Model definitions ─────────────────────────────────────────────────────────

def qm_singlet(d):
    """Standard QM prediction for singlet P_disagree."""
    return np.cos(d / 2) ** 2

def model2(d, a, V):
    """
    Perturbative asymmetric pole model with visibility correction.
    P = V * [cos^2(d/2) + (2a-1)*sin^2(d)/4] + (1-V)*0.5
    QM recovered at a=0.5. Max deviation at d=90deg.
    """
    qm = np.cos(d / 2) ** 2 + (2 * a - 1) * np.sin(d) ** 2 / 4
    return V * qm + (1 - V) * 0.5

def model3(d, a, V):
    """
    Power-law asymmetric pole model with visibility correction.
    P = [cos^2(d/2)]^(1/(2a))
    QM recovered at a=0.5.
    """
    base = np.clip(np.cos(d / 2) ** 2, 1e-12, 1)
    qm = base ** (1.0 / (2 * a))
    return V * qm + (1 - V) * 0.5

def qm_with_visibility(d, V):
    """QM singlet with free visibility (alpha fixed at 0.5)."""
    return V * np.cos(d / 2) ** 2 + (1 - V) * 0.5

# ── Circuit builder ───────────────────────────────────────────────────────────

def build_circuits(angles_rad):
    """
    Build one circuit per angle.
    State: singlet |psi-> = (|01> - |10>) / sqrt(2)
    Alice: always measures in Z basis
    Bob:   rotates by Ry(delta) before measuring Z
    Outcome disagreement = one of them gets 0, other gets 1.
    """
    from qiskit import QuantumCircuit

    circuits = []
    for delta in angles_rad:
        qc = QuantumCircuit(2, 2)

        # Prepare singlet |Ψ-> = (|01> - |10>)/sqrt(2)
        # X(q1) -> H(q0) -> CNOT(q0,q1) -> Z(q0)
        qc.x(1)
        qc.h(0)
        qc.cx(0, 1)
        qc.z(0)
        # Verified: X(q1)|00>=|01>, H(q0)->(|01>+|11>)/√2,
        # CNOT->(|01>+|10>)/√2=|Ψ+>, Z(q0)->(|01>-|10>)/√2=|Ψ->

        # Bob rotates by delta (Ry on qubit 1)
        qc.ry(delta, 1)

        qc.measure([0, 1], [0, 1])
        circuits.append(qc)

    return circuits

# ── P_disagree from counts ────────────────────────────────────────────────────

def p_disagree_from_counts(counts, shots):
    disagree = counts.get("01", 0) + counts.get("10", 0)
    return disagree / shots

# ── Readout error correction ──────────────────────────────────────────────────

def fetch_readout_cal(backend, qubit0=0, qubit1=1):
    """
    Fetch per-qubit readout error rates from backend calibration.
    Returns dict with e01 (prep0→meas1) and e10 (prep1→meas0) for each qubit.
    Qiskit convention: prob_meas1_prep0 = e01, prob_meas0_prep1 = e10
    """
    props = backend.properties()
    cal = {}
    for q, label in [(qubit0, "q0"), (qubit1, "q1")]:
        qp = props.qubit_property(q)
        cal[label] = {
            "e01": float(qp.get("prob_meas1_prep0", (0.0,))[0]),  # 0→1 error
            "e10": float(qp.get("prob_meas0_prep1", (0.0,))[0]),  # 1→0 error
        }
    print(f"  Readout cal q0: e01={cal['q0']['e01']:.4f}  e10={cal['q0']['e10']:.4f}")
    print(f"  Readout cal q1: e01={cal['q1']['e01']:.4f}  e10={cal['q1']['e10']:.4f}")
    return cal

def correct_counts(counts_dict, shots, cal):
    """
    Apply 2-qubit readout correction via matrix inversion.

    Qiskit bit strings are little-endian: "q1q0"
      "00" → q0=0,q1=0  idx 0
      "01" → q0=1,q1=0  idx 1
      "10" → q0=0,q1=1  idx 2
      "11" → q0=1,q1=1  idx 3

    Assignment matrix A = M_q1 ⊗ M_q0  where
      M_qi = [[1-e01_i,   e10_i ],
              [  e01_i, 1-e10_i ]]
    p_obs = A @ p_true  →  p_true = A^{-1} @ p_obs
    """
    def assignment_matrix(e01, e10):
        return np.array([[1 - e01, e10],
                         [e01,     1 - e10]])

    M0 = assignment_matrix(cal["q0"]["e01"], cal["q0"]["e10"])
    M1 = assignment_matrix(cal["q1"]["e01"], cal["q1"]["e10"])
    A  = np.kron(M1, M0)  # 4×4

    p_obs = np.array([
        counts_dict.get("00", 0),
        counts_dict.get("01", 0),
        counts_dict.get("10", 0),
        counts_dict.get("11", 0),
    ], dtype=float) / shots

    try:
        p_true = np.linalg.solve(A, p_obs)
        p_true = np.clip(p_true, 0, 1)          # numerical safety
        p_true /= p_true.sum()                   # renormalize
    except np.linalg.LinAlgError:
        p_true = p_obs                           # fallback: uncorrected

    # disagreement = "01" + "10" → indices 1 and 2
    p_disagree_corr = float(p_true[1] + p_true[2])
    err = np.sqrt(p_disagree_corr * (1 - p_disagree_corr) / shots)
    return p_disagree_corr, err

# ── Dry-run simulator ─────────────────────────────────────────────────────────

def simulate_locally(angles_rad, shots, alpha=0.5, V=1.0, noise=0.02):
    """
    Simulates ideal QM + Gaussian shot noise.
    Used for testing without IBM access.
    """
    results = []
    for d in angles_rad:
        p_true = model2(d, alpha, V)
        p_true = float(np.clip(p_true, 0, 1))
        n_disagree = np.random.binomial(shots, p_true)
        p_obs = n_disagree / shots
        err = np.sqrt(p_obs * (1 - p_obs) / shots)
        results.append({
            "angle_deg": float(np.degrees(d)),
            "p_disagree": p_obs,
            "err": err,
            "shots": shots,
        })
    return results

# ── IBM runner ────────────────────────────────────────────────────────────────

def run_on_ibm(angles_rad, shots, token, backend_name):
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    print(f"Connecting to IBM Quantum (backend: {backend_name})...")
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(backend_name)

    circuits = build_circuits(angles_rad)

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuits = pm.run(circuits)

    # Extract physical qubits from transpiled layout
    layout = isa_circuits[0].layout
    phys_q0 = layout.initial_layout[circuits[0].qubits[0]]
    phys_q1 = layout.initial_layout[circuits[0].qubits[1]]
    print(f"Physical qubit mapping: q0 -> {phys_q0}, q1 -> {phys_q1}")

    print("Fetching readout calibration...")
    cal = fetch_readout_cal(backend, qubit0=phys_q0, qubit1=phys_q1)

    print(f"Submitting {len(isa_circuits)} circuits ({shots} shots each)...")
    sampler = Sampler(backend)
    job = sampler.run(isa_circuits, shots=shots)
    print(f"Job ID: {job.job_id()}")
    print("Waiting for results (this may take several minutes)...")

    result = job.result()

    data = []
    for i, (delta, pub_result) in enumerate(zip(angles_rad, result)):
        counts = pub_result.data.c.get_counts()
        p_dis = p_disagree_from_counts(counts, shots)
        err = np.sqrt(p_dis * (1 - p_dis) / shots)
        p_dis_corr, err_corr = correct_counts(counts, shots, cal)
        data.append({
            "angle_deg": float(np.degrees(delta)),
            "p_disagree": float(p_dis),
            "err": float(err),
            "p_disagree_corr": float(p_dis_corr),
            "err_corr": float(err_corr),
            "shots": shots,
            "job_id": job.job_id(),
            "cal": cal,
        })

    return data

# ── Fitter ────────────────────────────────────────────────────────────────────

def fit_model(model_fn, angles_rad, p_obs, p_err, label):
    # Initial guess: a=0.5, V=0.9
    p0 = [0.5, 0.9]
    bounds = ([0.3, 0.5], [1.0, 1.0])

    try:
        popt, pcov = curve_fit(
            model_fn, angles_rad, p_obs,
            p0=p0, sigma=p_err, absolute_sigma=True,
            bounds=bounds, maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov))
        a_fit, V_fit = popt
        a_err, V_err = perr

        # chi2
        residuals = p_obs - model_fn(angles_rad, *popt)
        chi2_val = np.sum((residuals / p_err) ** 2)
        dof = len(angles_rad) - 2
        chi2_dof = chi2_val / dof

        print(f"\n  {label}:")
        print(f"    alpha = {a_fit:.4f} ± {a_err:.4f}")
        print(f"    V     = {V_fit:.4f} ± {V_err:.4f}")
        print(f"    chi2/dof = {chi2_dof:.3f}  ({chi2_val:.2f}/{dof})")

        return {
            "label": label,
            "a": a_fit, "a_err": a_err,
            "V": V_fit, "V_err": V_err,
            "chi2": chi2_val, "dof": dof, "chi2_dof": chi2_dof,
            "popt": popt.tolist(),
        }

    except Exception as e:
        print(f"  {label}: fit failed — {e}")
        return None

def fit_qm_v(angles_rad, p_obs, p_err):
    """Fit QM + visibility (1 free param: V). Alpha fixed at 0.5."""
    try:
        popt, pcov = curve_fit(
            qm_with_visibility, angles_rad, p_obs,
            p0=[0.9], sigma=p_err, absolute_sigma=True,
            bounds=([0.5], [1.0]), maxfev=10000,
        )
        V_fit = popt[0]
        V_err = float(np.sqrt(pcov[0, 0]))
        residuals = p_obs - qm_with_visibility(angles_rad, V_fit)
        chi2_val = np.sum((residuals / p_err) ** 2)
        dof = len(angles_rad) - 1
        chi2_dof = chi2_val / dof

        print(f"\n  QM + V:")
        print(f"    V     = {V_fit:.4f} ± {V_err:.4f}")
        print(f"    chi2/dof = {chi2_dof:.3f}  ({chi2_val:.2f}/{dof})")

        return {"label": "QM+V", "V": V_fit, "V_err": V_err,
                "chi2": float(chi2_val), "dof": dof, "chi2_dof": float(chi2_dof)}
    except Exception as e:
        print(f"  QM+V: fit failed — {e}")
        return None

# ── Plotter ───────────────────────────────────────────────────────────────────

def plot_results(angles_rad, data, fits, fits_corr=None, out_path="spin_sweep_results.png"):
    angles_deg = np.degrees(angles_rad)
    p_obs      = np.array([d["p_disagree"] for d in data])
    p_err      = np.array([d["err"] for d in data])
    has_corr   = "p_disagree_corr" in data[0]
    if has_corr:
        p_corr     = np.array([d["p_disagree_corr"] for d in data])
        p_corr_err = np.array([d["err_corr"] for d in data])

    d_fine = np.linspace(0, np.pi, 500)
    d_fine_deg = np.degrees(d_fine)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: full curves
    ax = axes[0]
    ax.errorbar(angles_deg, p_obs, yerr=p_err, fmt="o", color="red",
                label="IBM data (raw)", zorder=5, ms=4)
    if has_corr:
        ax.errorbar(angles_deg, p_corr, yerr=p_corr_err, fmt="s", color="blue",
                    label="IBM data (corrected)", zorder=5, ms=4, alpha=0.7)
    ax.plot(d_fine_deg, qm_singlet(d_fine), "k-", lw=2, label="QM (singlet)")

    colors = {"Model 2": "steelblue", "Model 3": "darkorange"}
    for fit in fits:
        if fit is None:
            continue
        fn = model2 if "2" in fit["label"] else model3
        y = fn(d_fine, fit["a"], fit["V"])
        lbl = f"{fit['label']} α={fit['a']:.3f}±{fit['a_err']:.3f}"
        ax.plot(d_fine_deg, y, "--", lw=1.5, label=lbl, color=colors.get(fit["label"], "gray"))

    ax.set_xlabel("δ (degrees)")
    ax.set_ylabel("P_disagree")
    ax.set_title("Correlation curve — IBM hardware vs models")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Right: residuals from QM
    ax2 = axes[1]
    qm_at_angles = qm_singlet(angles_rad)
    residuals = p_obs - qm_at_angles
    ax2.errorbar(angles_deg, residuals, yerr=p_err, fmt="o", color="red",
                 label="Raw - QM", ms=4)
    if has_corr:
        res_corr = p_corr - qm_at_angles
        ax2.errorbar(angles_deg, res_corr, yerr=p_corr_err, fmt="s", color="blue",
                     label="Corrected - QM", ms=4, alpha=0.7)
    ax2.axhline(0, color="black", lw=1.5, label="QM baseline")

    # raw fits
    for fit in fits:
        if fit is None:
            continue
        fn = model2 if "2" in fit["label"] else model3
        y_model = fn(angles_rad, fit["a"], fit["V"])
        ax2.plot(angles_deg, y_model - qm_at_angles, "--", lw=1.5,
                 label=f"{fit['label']} raw residual", color=colors.get(fit["label"], "gray"))

    # corrected fits
    if fits_corr:
        corr_colors = {"Model 2": "navy", "Model 3": "darkorange"}
        for fit in fits_corr:
            if fit is None:
                continue
            fn = model2 if "2" in fit["label"] else model3
            y_model = fn(angles_rad, fit["a"], fit["V"])
            ax2.plot(angles_deg, y_model - qm_at_angles, "-.", lw=1.5,
                     label=f"{fit['label']} corr residual", color=corr_colors.get(fit["label"], "gray"))

    ax2.set_xlabel("δ (degrees)")
    ax2.set_ylabel("P_disagree - QM")
    ax2.set_title("Deviation from QM prediction")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"\nPlot saved: {out_path}")

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if args.results:
        print(f"Loading saved results from {args.results}")
        with open(args.results) as f:
            data = json.load(f)
    elif args.dry_run:
        print("DRY RUN — simulating locally (QM + shot noise, alpha=0.5, V=1.0)")
        data = simulate_locally(ANGLES_RAD, args.shots)
    else:
        if not args.token:
            raise ValueError("Provide --token or use --dry-run")
        data = run_on_ibm(ANGLES_RAD, args.shots, args.token, args.backend)

    # Save raw results
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_json = f"spin_sweep_raw_{ts}.json"
    with open(out_json, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Raw results saved: {out_json}")

    # Extract arrays
    angles_deg_data = np.array([d["angle_deg"] for d in data])
    angles_rad_data = np.deg2rad(angles_deg_data)
    p_obs = np.array([d["p_disagree"] for d in data])
    p_err = np.array([d["err"] for d in data])
    # Avoid zero errors
    p_err = np.maximum(p_err, 1.0 / args.shots)

    print("\n=== FIT RESULTS (RAW) ===")
    print(f"  Data points: {len(data)}")
    print(f"  Shots per angle: {args.shots}")

    fit_qmv = fit_qm_v(angles_rad_data, p_obs, p_err)
    fit2 = fit_model(model2, angles_rad_data, p_obs, p_err, "Model 2")
    fit3 = fit_model(model3, angles_rad_data, p_obs, p_err, "Model 3")

    fits_corr = None
    has_corr = "p_disagree_corr" in data[0]
    if has_corr:
        p_corr     = np.array([d["p_disagree_corr"] for d in data])
        p_corr_err = np.array([d["err_corr"] for d in data])
        p_corr_err = np.maximum(p_corr_err, 1.0 / args.shots)
        print("\n=== FIT RESULTS (READOUT-CORRECTED) ===")
        fit_qmv_c = fit_qm_v(angles_rad_data, p_corr, p_corr_err)
        fit2c = fit_model(model2, angles_rad_data, p_corr, p_corr_err, "Model 2")
        fit3c = fit_model(model3, angles_rad_data, p_corr, p_corr_err, "Model 3")
        fits_corr = [fit2c, fit3c]

    print("\n=== MODEL COMPARISON (Likelihood Ratio Test) ===")
    pairs = [("raw", fit_qmv, fit2, fit3)]
    if has_corr:
        pairs.append(("corrected", fit_qmv_c, fit2c, fit3c))
    for suffix, qmv, m2, m3 in pairs:
        if qmv is None:
            continue
        print(f"\n  --- {suffix} ---")
        print(f"  QM+V:    chi2/dof = {qmv['chi2_dof']:.2f}  (V = {qmv['V']:.4f})")
        for m, mlabel in [(m2, "Model 2"), (m3, "Model 3")]:
            if m is None:
                continue
            delta_chi2 = qmv["chi2"] - m["chi2"]
            p_lrt = 1 - chi2.cdf(delta_chi2, 1)
            sigma_lrt = norm.isf(p_lrt / 2) if p_lrt > 0 else float('inf')
            print(f"  {mlabel}: chi2/dof = {m['chi2_dof']:.2f}  (alpha = {m['a']:.4f} +/- {m['a_err']:.4f})")
            print(f"    LRT vs QM+V: dchi2 = {delta_chi2:.1f}, p = {p_lrt:.2e} ({sigma_lrt:.1f}sigma)")

    plot_results(angles_rad_data, data, [fit2, fit3], fits_corr=fits_corr, out_path=f"spin_sweep_results_{ts}.png")

if __name__ == "__main__":
    main()
