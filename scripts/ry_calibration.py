"""
Ry Gate Calibration Check
==========================
Prepares |0>, applies Ry(delta) to qubit 1 (Bob's qubit), measures P(|1>).
Compares to ideal sin^2(delta/2) to detect any systematic gate offset.

If best-fit epsilon != 0, then the Bell sweep was measuring cos^2((delta+eps)/2)
instead of cos^2(delta/2), which would explain alpha < 0.5 without any physics.

Usage:
  python ry_calibration.py --token YOUR_TOKEN --backend ibm_marrakesh
  python ry_calibration.py --dry-run
  python ry_calibration.py --results ry_cal_raw_TIMESTAMP.json
"""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--token",   required=False)
parser.add_argument("--backend", default="ibm_marrakesh")
parser.add_argument("--shots",   type=int, default=8192)
parser.add_argument("--dry-run", action="store_true")
parser.add_argument("--results", default=None)
args = parser.parse_args()

ANGLES_DEG = np.arange(0, 185, 5)
ANGLES_RAD = np.deg2rad(ANGLES_DEG)
QUBIT = 1  # Bob's qubit — the one Ry acts on in the Bell sweep

# ── Models ────────────────────────────────────────────────────────────────────

def ideal(d):
    return np.sin(d / 2) ** 2

def with_offset(d, eps):
    """Ry gate fires at d+eps instead of d."""
    return np.sin((d + eps) / 2) ** 2

def with_offset_and_scale(d, eps, V):
    """Offset + visibility (decoherence compresses toward 0.5)."""
    return V * np.sin((d + eps) / 2) ** 2 + (1 - V) * 0.5

# ── Circuits ──────────────────────────────────────────────────────────────────

def build_circuits(angles_rad):
    from qiskit import QuantumCircuit
    circuits = []
    for delta in angles_rad:
        qc = QuantumCircuit(2, 1)   # 2 qubits to match Bell sweep layout
        qc.ry(delta, QUBIT)         # same qubit as Bell sweep
        qc.measure(QUBIT, 0)
        circuits.append(qc)
    return circuits

# ── IBM runner ────────────────────────────────────────────────────────────────

def run_on_ibm(angles_rad, shots, token, backend_name):
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(backend_name)
    circuits = build_circuits(angles_rad)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuits = pm.run(circuits)

    print(f"Submitting {len(isa_circuits)} circuits ({shots} shots each)...")
    sampler = Sampler(backend)
    job = sampler.run(isa_circuits, shots=shots)
    print(f"Job ID: {job.job_id()}")
    result = job.result()

    data = []
    for delta, pub_result in zip(angles_rad, result):
        counts = pub_result.data.c.get_counts()
        n1 = counts.get("1", 0)
        p1 = n1 / shots
        err = np.sqrt(p1 * (1 - p1) / shots)
        data.append({
            "angle_deg": float(np.degrees(delta)),
            "p1": float(p1),
            "err": float(err),
            "shots": shots,
            "job_id": job.job_id(),
        })
    return data

# ── Dry-run ───────────────────────────────────────────────────────────────────

def simulate_locally(angles_rad, shots, eps=0.0, V=1.0):
    data = []
    for d in angles_rad:
        p_true = with_offset_and_scale(d, eps, V)
        n1 = np.random.binomial(shots, float(np.clip(p_true, 0, 1)))
        p1 = n1 / shots
        err = np.sqrt(p1 * (1 - p1) / shots)
        data.append({
            "angle_deg": float(np.degrees(d)),
            "p1": p1, "err": err, "shots": shots,
        })
    return data

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if args.results:
        with open(args.results) as f:
            data = json.load(f)
        print(f"Loaded {args.results}")
    elif args.dry_run:
        print("DRY RUN — ideal Ry, eps=0, V=1")
        data = simulate_locally(ANGLES_RAD, args.shots)
    else:
        if not args.token:
            raise ValueError("Provide --token or use --dry-run")
        data = run_on_ibm(ANGLES_RAD, args.shots, args.token, args.backend)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_json = f"ry_cal_raw_{ts}.json"
    with open(out_json, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Raw results saved: {out_json}")

    angles_rad = np.deg2rad([d["angle_deg"] for d in data])
    angles_deg = np.degrees(angles_rad)
    p1  = np.array([d["p1"]  for d in data])
    err = np.array([d["err"] for d in data])
    err = np.maximum(err, 1.0 / args.shots)

    # Fit 1: offset only
    try:
        (eps_fit,), pcov1 = curve_fit(with_offset, angles_rad, p1,
                                      p0=[0.0], sigma=err, absolute_sigma=True,
                                      bounds=([-np.pi/4], [np.pi/4]))
        eps_err = float(np.sqrt(pcov1[0, 0]))
        eps_deg = float(np.degrees(eps_fit))
        eps_deg_err = float(np.degrees(eps_err))
        res1 = p1 - with_offset(angles_rad, eps_fit)
        chi2_1 = float(np.sum((res1 / err) ** 2) / (len(p1) - 1))
    except Exception as e:
        print(f"Offset-only fit failed: {e}")
        eps_fit, eps_deg, eps_deg_err, chi2_1 = 0.0, 0.0, 0.0, 0.0

    # Fit 2: offset + visibility
    try:
        (eps_fit2, V_fit), pcov2 = curve_fit(
            with_offset_and_scale, angles_rad, p1,
            p0=[0.0, 0.9], sigma=err, absolute_sigma=True,
            bounds=([-np.pi/4, 0.5], [np.pi/4, 1.0]))
        perr2 = np.sqrt(np.diag(pcov2))
        eps_deg2     = float(np.degrees(eps_fit2))
        eps_deg2_err = float(np.degrees(perr2[0]))
        chi2_2 = float(np.sum((p1 - with_offset_and_scale(angles_rad, eps_fit2, V_fit))**2 / err**2) / (len(p1) - 2))
    except Exception as e:
        print(f"Offset+V fit failed: {e}")
        eps_fit2, V_fit, eps_deg2, eps_deg2_err, chi2_2 = 0.0, 1.0, 0.0, 0.0, 0.0

    print("\n=== Ry CALIBRATION RESULTS ===")
    print(f"  Offset-only fit:    eps = {eps_deg:.3f} +/- {eps_deg_err:.3f} deg   chi2/dof = {chi2_1:.2f}")
    print(f"  Offset + V fit:     eps = {eps_deg2:.3f} +/- {eps_deg2_err:.3f} deg   V = {V_fit:.4f}   chi2/dof = {chi2_2:.2f}")

    # How much alpha shift would this epsilon cause in the Bell sweep?
    # cos^2((d+eps)/2) vs cos^2(d/2): at d=90deg, shift ~ eps/2 * sin(90) = eps/2
    # alpha deviation ~ eps/(2*pi) * (crossover sensitivity)
    # Rough estimate: delta_alpha ~ -eps_rad / (2*pi) * 2 = -eps_rad/pi
    if abs(eps_fit2) > 0:
        delta_alpha_estimate = -float(eps_fit2) / np.pi
        observed_deviation = 0.030
        fraction = abs(delta_alpha_estimate) / observed_deviation
        print(f"\n  Implied alpha shift from gate offset: {delta_alpha_estimate:+.4f}")
        print(f"  Observed alpha deviation in Bell sweep: ~-0.030")
        print(f"  Gate offset accounts for ~{fraction*100:.0f}% of the observed deviation")
        if fraction >= 0.8:
            print("  *** Gate offset accounts for most of the Bell sweep alpha deviation ***")
        elif fraction >= 0.1:
            print("  *** Gate offset is insufficient to explain the full deviation ***")
        else:
            print("  *** Gate offset contribution is negligible ***")

    # Plot
    d_fine = np.linspace(0, np.pi, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.errorbar(angles_deg, p1, yerr=err, fmt="o", color="red", ms=4, label="IBM data")
    ax.plot(np.degrees(d_fine), ideal(d_fine), "k-", lw=2, label="Ideal sin²(δ/2)")
    ax.plot(np.degrees(d_fine), with_offset(d_fine, eps_fit),
            "b--", lw=1.5, label=f"Offset fit ε={eps_deg:.2f}°")
    ax.plot(np.degrees(d_fine), with_offset_and_scale(d_fine, eps_fit2, V_fit),
            "g-.", lw=1.5, label=f"Offset+V fit ε={eps_deg2:.2f}°, V={V_fit:.3f}")
    ax.set_xlabel("δ (degrees)")
    ax.set_ylabel("P(|1⟩)")
    ax.set_title(f"Ry calibration — qubit {QUBIT} on {args.backend}")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    ax2.errorbar(angles_deg, p1 - ideal(angles_rad), yerr=err,
                 fmt="o", color="red", ms=4, label="Data - ideal")
    ax2.plot(np.degrees(d_fine),
             with_offset(d_fine, eps_fit) - ideal(d_fine),
             "b--", lw=1.5, label=f"Offset fit residual (ε={eps_deg:.2f}°)")
    ax2.axhline(0, color="black", lw=1.5)
    ax2.set_xlabel("δ (degrees)")
    ax2.set_ylabel("P(|1⟩) - ideal")
    ax2.set_title("Deviation from ideal Ry gate")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = f"ry_cal_results_{ts}.png"
    plt.savefig(out_png, dpi=150)
    print(f"Plot saved: {out_png}")

if __name__ == "__main__":
    main()
