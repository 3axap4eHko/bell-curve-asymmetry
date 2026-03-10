"""
T2 Decoherence Model for Bell Sweep
=====================================
Simulates the singlet Bell sweep under realistic T1/T2 thermal relaxation noise
using IBM marrakesh published hardware parameters (or token-fetched values).

Fits the resulting P_disagree curve with Model 2/3 and reports whether the
thermal-relaxation channel produces any systematic alpha shift. In the
fixed-duration model used here, alpha should remain near 0.5; any single-run
deviation is expected to be finite-shot scatter rather than channel bias.

Usage:
  python t2_noise_model.py                          # default marrakesh params
  python t2_noise_model.py --token TOKEN            # fetch live T1/T2 from IBM
  python t2_noise_model.py --t1 200 --t2 120        # override T1/T2 in microseconds
  python t2_noise_model.py --sweep-t2               # sweep T2 from 50 to 300 us
"""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--token",    default=None)
parser.add_argument("--backend",  default="ibm_marrakesh")
parser.add_argument("--shots",    type=int, default=8192)
# T1/T2 in microseconds -- marrakesh typical published values
parser.add_argument("--t1-q0",   type=float, default=250.0)
parser.add_argument("--t2-q0",   type=float, default=130.0)
parser.add_argument("--t1-q1",   type=float, default=220.0)
parser.add_argument("--t2-q1",   type=float, default=110.0)
# Gate times in nanoseconds -- marrakesh typical
parser.add_argument("--t-sq",    type=float, default=36.0,   help="Single-qubit (SX) gate time (ns)")
parser.add_argument("--t-cx",    type=float, default=68.0,   help="Two-qubit (CZ) gate time (ns)")
parser.add_argument("--sweep-t2", action="store_true")
args = parser.parse_args()

ANGLES_DEG = np.arange(0, 185, 5)
ANGLES_RAD = np.deg2rad(ANGLES_DEG)

# -- Models (same as Bell sweep script) ----------------------------------------

def qm_singlet(d):
    return np.cos(d / 2) ** 2

def model2(d, alpha, V):
    return V * (np.cos(d/2)**2 + (2*alpha - 1) * np.sin(d)**2 / 4) + (1 - V) * 0.5

def model3(d, alpha, V):
    base = np.cos(d / 2) ** 2
    exp = 1.0 / (2 * alpha)
    return V * np.power(np.clip(base, 0, 1), exp) + (1 - V) * 0.5

def fit_model(angles, p_data, err, model_fn, label):
    try:
        (alpha, V), pcov = curve_fit(
            model_fn, angles, p_data,
            p0=[0.48, 0.95], sigma=err, absolute_sigma=True,
            bounds=([0.3, 0.5], [1.0, 1.0]))
        perr = np.sqrt(np.diag(pcov))
        res = p_data - model_fn(angles, alpha, V)
        chi2 = float(np.sum((res / err)**2) / (len(p_data) - 2))
        z = (alpha - 0.5) / perr[0]
        print(f"  {label}: alpha={alpha:.4f}+/-{perr[0]:.4f}  V={V:.4f}  chi2/dof={chi2:.2f}  z={(z):.1f}sigma")
        return alpha, V, perr[0], chi2
    except Exception as e:
        print(f"  {label}: fit failed -- {e}")
        return 0.5, 1.0, 0.0, 0.0

# -- Fetch live T1/T2 from IBM -------------------------------------------------

def fetch_hardware_params(token, backend_name):
    from qiskit_ibm_runtime import QiskitRuntimeService
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(backend_name)
    props = backend.properties()
    params = {}
    for q in [0, 1]:
        params[f"t1_q{q}"] = props.t1(q) * 1e6   # seconds -> microseconds
        params[f"t2_q{q}"] = props.t2(q) * 1e6
    try:
        params["t_sq"] = props.gate_length("sx", [0]) * 1e9
        params["t_cx"] = props.gate_length("cx", [0, 1]) * 1e9
    except Exception:
        params["t_sq"] = args.t_sq
        params["t_cx"] = args.t_cx
    print(f"Live hardware params from {backend_name}:")
    for k, v in params.items():
        print(f"  {k}: {v:.2f}")
    return params

# -- Build singlet circuit ------------------------------------------------------

def build_bell_circuits(angles_rad):
    from qiskit import QuantumCircuit
    circuits = []
    for delta in angles_rad:
        qc = QuantumCircuit(2, 2)
        # Singlet |Psi-> prep
        qc.x(1)
        qc.h(0)
        qc.cx(0, 1)
        qc.z(0)
        # Bob rotates
        qc.ry(delta, 1)
        qc.measure([0, 1], [0, 1])
        circuits.append(qc)
    return circuits

# -- Build noise model from T1/T2 ----------------------------------------------

def build_noise_model(t1_q0, t2_q0, t1_q1, t2_q1, t_sq_ns, t_cx_ns):
    from qiskit_aer.noise import NoiseModel, thermal_relaxation_error

    # Convert to seconds
    t1 = [t1_q0 * 1e-6, t1_q1 * 1e-6]
    t2 = [t2_q0 * 1e-6, t2_q1 * 1e-6]
    t_sq = t_sq_ns * 1e-9
    t_cx = t_cx_ns * 1e-9

    # Clamp T2 <= 2*T1 (physical constraint)
    for q in range(2):
        if t2[q] > 2 * t1[q]:
            t2[q] = 2 * t1[q]

    noise_model = NoiseModel()

    # Single-qubit gate errors on each qubit
    for q in range(2):
        err_sq = thermal_relaxation_error(t1[q], t2[q], t_sq)
        noise_model.add_quantum_error(err_sq, ["rx", "ry", "rz", "sx", "x", "h", "z"], [q])

    # 2-qubit CX error -- apply to both qubits during gate
    err_cx_q0 = thermal_relaxation_error(t1[0], t2[0], t_cx)
    err_cx_q1 = thermal_relaxation_error(t1[1], t2[1], t_cx)
    err_cx = err_cx_q0.expand(err_cx_q1)
    noise_model.add_quantum_error(err_cx, ["cx"], [0, 1])

    return noise_model

# -- Run simulation -------------------------------------------------------------

def simulate_sweep(angles_rad, shots, noise_model, label=""):
    from qiskit_aer import AerSimulator
    from qiskit import transpile

    simulator = AerSimulator(noise_model=noise_model)
    circuits = build_bell_circuits(angles_rad)
    t_circuits = transpile(circuits, simulator)

    p_disagree = []
    for qc in t_circuits:
        job = simulator.run(qc, shots=shots)
        counts = job.result().get_counts()
        n_dis = counts.get("01", 0) + counts.get("10", 0)
        p_disagree.append(n_dis / shots)

    p_disagree = np.array(p_disagree)
    err = np.sqrt(p_disagree * (1 - p_disagree) / shots)
    err = np.maximum(err, 1.0 / shots)
    return p_disagree, err

# -- T2 sweep ------------------------------------------------------------------

def run_t2_sweep(base_params):
    t2_values = np.arange(50, 320, 25)  # microseconds
    alphas_m2, alphas_m3 = [], []

    print(f"\n{'T2 (us)':>10} {'alpha_M2':>10} {'alpha_M3':>10}")
    print("-" * 35)

    for t2_val in t2_values:
        nm = build_noise_model(
            base_params["t1_q0"], t2_val,
            base_params["t1_q1"], t2_val,
            base_params["t_sq"], base_params["t_cx"])
        p, err = simulate_sweep(ANGLES_RAD, args.shots, nm)
        try:
            (a2, _), _ = curve_fit(model2, ANGLES_RAD, p, p0=[0.48, 0.95],
                                   bounds=([0.3,0.5],[1.0,1.0]))
            (a3, _), _ = curve_fit(model3, ANGLES_RAD, p, p0=[0.48, 0.95],
                                   bounds=([0.3,0.5],[1.0,1.0]))
        except Exception:
            a2, a3 = 0.5, 0.5
        alphas_m2.append(a2)
        alphas_m3.append(a3)
        print(f"{t2_val:>10.0f} {a2:>10.4f} {a3:>10.4f}")

    return t2_values, np.array(alphas_m2), np.array(alphas_m3)

# -- Main -----------------------------------------------------------------------

def main():
    if args.token:
        hw = fetch_hardware_params(args.token, args.backend)
        t1_q0, t2_q0 = hw["t1_q0"], hw["t2_q0"]
        t1_q1, t2_q1 = hw["t1_q1"], hw["t2_q1"]
        t_sq, t_cx   = hw["t_sq"],  hw["t_cx"]
    else:
        t1_q0, t2_q0 = args.t1_q0, args.t2_q0
        t1_q1, t2_q1 = args.t1_q1, args.t2_q1
        t_sq,  t_cx  = args.t_sq,  args.t_cx

    print(f"\nNoise parameters:")
    print(f"  Q0: T1={t1_q0:.1f}us  T2={t2_q0:.1f}us")
    print(f"  Q1: T1={t1_q1:.1f}us  T2={t2_q1:.1f}us")
    print(f"  Gate times: sq={t_sq:.0f}ns  cx={t_cx:.0f}ns")

    base_params = dict(t1_q0=t1_q0, t2_q0=t2_q0, t1_q1=t1_q1, t2_q1=t2_q1,
                       t_sq=t_sq, t_cx=t_cx)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    if args.sweep_t2:
        print("\n=== T2 SWEEP ===")
        t2_vals, a2, a3 = run_t2_sweep(base_params)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.axhline(0.5,   color="black", lw=1.5, ls="--", label="QM (alpha=0.5)")
        ax.axhline(0.467, color="red",   lw=1.5, ls=":",  label="Observed alpha=0.467 (marrakesh)")
        ax.plot(t2_vals, a2, "b-o", ms=5, label="Model 2")
        ax.plot(t2_vals, a3, "g-s", ms=5, label="Model 3")
        ax.set_xlabel("T2 (us)")
        ax.set_ylabel("Fitted alpha")
        ax.set_title("Alpha vs T2 -- does decoherence fake the signal?")
        ax.legend()
        ax.grid(True, alpha=0.3)
        out_png = f"t2_sweep_alpha_{ts}.png"
        plt.tight_layout()
        plt.savefig(out_png, dpi=150)
        print(f"\nPlot saved: {out_png}")
        return

    # Single noise model run
    print("\nRunning noisy simulation...")
    noise_model = build_noise_model(t1_q0, t2_q0, t1_q1, t2_q1, t_sq, t_cx)
    p_noisy, err_noisy = simulate_sweep(ANGLES_RAD, args.shots, noise_model)

    # Also run ideal (no noise) for reference
    print("Running ideal simulation...")
    from qiskit_aer.noise import NoiseModel
    ideal_nm = NoiseModel()  # empty
    p_ideal, err_ideal = simulate_sweep(ANGLES_RAD, args.shots, ideal_nm)

    print("\n=== FIT RESULTS ===")
    print("Ideal (shot noise only):")
    fit_model(ANGLES_RAD, p_ideal, err_ideal, model2, "Model 2")
    fit_model(ANGLES_RAD, p_ideal, err_ideal, model3, "Model 3")
    print(f"\nWith T2 noise (T2_q0={t2_q0}us, T2_q1={t2_q1}us):")
    a2, V2, a2_err, chi2_2 = fit_model(ANGLES_RAD, p_noisy, err_noisy, model2, "Model 2")
    a3, V3, a3_err, chi2_3 = fit_model(ANGLES_RAD, p_noisy, err_noisy, model3, "Model 3")

    print(f"\n=== VERDICT ===")
    print("  Single finite-shot runs can fluctuate around alpha=0.5.")
    print("  The relevant test is whether the channel produces a systematic shift across T2 values.")
    for a, label in [(a2, "Model 2"), (a3, "Model 3")]:
        print(f"  {label}: noise produces alpha={a:.4f} in this run")

    # Save results
    out = {
        "params": base_params,
        "angles_deg": ANGLES_DEG.tolist(),
        "p_ideal": p_ideal.tolist(),
        "p_noisy": p_noisy.tolist(),
        "err_noisy": err_noisy.tolist(),
        "fits": {"m2_alpha": a2, "m2_V": V2, "m3_alpha": a3, "m3_V": V3},
    }
    out_json = f"t2_noise_results_{ts}.json"
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)

    # Plot
    d_fine = np.linspace(0, np.pi, 500)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(np.degrees(d_fine), qm_singlet(d_fine), "k-", lw=2, label="QM ideal")
    ax.plot(ANGLES_DEG, p_ideal, "g.", ms=6, alpha=0.6, label="Sim: no noise")
    ax.errorbar(ANGLES_DEG, p_noisy, yerr=err_noisy, fmt="ro", ms=4,
                label=f"Sim: T2 noise (q0={t2_q0}us, q1={t2_q1}us)")
    ax.plot(np.degrees(d_fine), model2(d_fine, a2, V2), "b--", lw=1.5,
            label=f"M2 fit alpha={a2:.4f}, V={V2:.3f}")
    ax.set_xlabel("delta (degrees)")
    ax.set_ylabel("P_disagree")
    ax.set_title("T2 noise effect on Bell sweep curve")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    ax2.errorbar(ANGLES_DEG, p_noisy - qm_singlet(ANGLES_RAD), yerr=err_noisy,
                 fmt="ro", ms=4, label="T2 noisy sim - QM")
    ax2.plot(np.degrees(d_fine), model2(d_fine, a2, V2) - qm_singlet(d_fine),
             "b--", lw=1.5, label=f"M2 fit residual (alpha={a2:.4f})")
    ax2.axhline(0, color="black", lw=1.5)
    ax2.set_xlabel("delta (degrees)")
    ax2.set_ylabel("P_disagree - QM")
    ax2.set_title("Deviation from QM under T2 noise")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = f"t2_noise_results_{ts}.png"
    plt.savefig(out_png, dpi=150)
    print(f"\nResults: {out_json}")
    print(f"Plot:    {out_png}")

if __name__ == "__main__":
    main()
