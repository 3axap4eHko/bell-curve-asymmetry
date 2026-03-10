"""
Crosstalk Check
===============
Applies Ry(delta) to qubit 1, measures qubit 0 only.
No entanglement. If P(|1>) on q0 varies with delta -> Ry on q1 is leaking
into q0, and crosstalk is a viable explanation for the Bell sweep alpha deviation.

Usage:
  python crosstalk_check.py --token TOKEN --backend ibm_marrakesh
  python crosstalk_check.py --results crosstalk_raw_TIMESTAMP.json
  python crosstalk_check.py --dry-run
"""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--token",   required=False)
parser.add_argument("--backend", default="ibm_marrakesh")
parser.add_argument("--shots",   type=int, default=4096)
parser.add_argument("--dry-run", action="store_true")
parser.add_argument("--results", default=None)
args = parser.parse_args()

ANGLES_DEG = np.arange(0, 190, 10)   # 19 points, 0-180 in 10-degree steps
ANGLES_RAD = np.deg2rad(ANGLES_DEG)

# ── Circuit ───────────────────────────────────────────────────────────────────

def build_circuits(angles_rad):
    from qiskit import QuantumCircuit
    circuits = []
    for delta in angles_rad:
        qc = QuantumCircuit(2, 1)
        qc.ry(delta, 1)       # Ry on Bob's qubit
        qc.measure(0, 0)      # measure Alice's qubit only
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
            "p1_q0": float(p1),
            "err": float(err),
            "shots": shots,
            "job_id": job.job_id(),
        })
    return data

# ── Dry-run ───────────────────────────────────────────────────────────────────

def simulate_locally(angles_rad, shots):
    """No crosstalk: q0 should stay at ~0 (just readout error ~0.01)."""
    data = []
    for d in angles_rad:
        p_true = 0.01  # readout error floor only
        n1 = np.random.binomial(shots, p_true)
        p1 = n1 / shots
        err = np.sqrt(p1 * (1 - p1) / shots)
        data.append({
            "angle_deg": float(np.degrees(d)),
            "p1_q0": p1, "err": err, "shots": shots,
        })
    return data

# ── Analysis ──────────────────────────────────────────────────────────────────

def analyze(data):
    angles_deg = np.array([d["angle_deg"] for d in data])
    p1  = np.array([d["p1_q0"] for d in data])
    err = np.array([d["err"]   for d in data])
    err = np.maximum(err, 1.0 / args.shots)

    mean_p1 = float(np.mean(p1))
    std_p1  = float(np.std(p1))
    # Peak-to-peak variation relative to baseline
    ptp     = float(np.ptp(p1))

    # Fit a sinusoidal leakage model: p = a + b*sin^2(delta/2)
    # If b is nonzero -> Ry on q1 is rotating q0
    def leakage_model(d, a, b):
        return a + b * np.sin(d / 2) ** 2

    from scipy.optimize import curve_fit
    try:
        (a_fit, b_fit), pcov = curve_fit(
            leakage_model, np.deg2rad(angles_deg), p1,
            p0=[mean_p1, 0.0], sigma=err, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        z_b  = abs(b_fit) / perr[1]
        chi2 = float(np.sum(((p1 - leakage_model(np.deg2rad(angles_deg), a_fit, b_fit)) / err)**2) / (len(p1) - 2))
    except Exception as e:
        a_fit, b_fit, perr, z_b, chi2 = mean_p1, 0.0, [0, 0], 0.0, 0.0
        print(f"Fit failed: {e}")

    print("\n=== CROSSTALK RESULTS ===")
    print(f"  Q0 baseline (mean P(|1>)):  {mean_p1:.5f}")
    print(f"  Q0 std across angles:       {std_p1:.5f}")
    print(f"  Q0 peak-to-peak:            {ptp:.5f}")
    print(f"  Leakage fit: a={a_fit:.5f}  b={b_fit:.5f}±{perr[1]:.5f}  z={z_b:.1f}σ  chi2/dof={chi2:.2f}")

    print("\n  VERDICT:")
    if z_b > 5 and abs(b_fit) > 0.01:
        print(f"  *** CROSSTALK DETECTED: b={b_fit:.4f} at {z_b:.0f}σ ***")
        print( "  Ry on q1 is leaking into q0. This is a viable explanation for alpha<0.5.")
    elif z_b > 3:
        print(f"  *** MARGINAL: b={b_fit:.4f} at {z_b:.1f}σ — weak crosstalk, probably not sufficient ***")
    else:
        print(f"  *** NO CROSSTALK: b consistent with zero at {z_b:.1f}σ ***")
        print( "  Ry on q1 does not measurably affect q0.")
        print( "  No evidence for crosstalk as explanation for Bell sweep alpha deviation.")

    return angles_deg, p1, err, a_fit, b_fit

# ── Plot ──────────────────────────────────────────────────────────────────────

def plot(angles_deg, p1, err, a_fit, b_fit, ts):
    d_fine = np.linspace(0, np.pi, 300)
    leakage_curve = a_fit + b_fit * np.sin(d_fine / 2) ** 2

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.errorbar(angles_deg, p1, yerr=err, fmt="ro", ms=5, label="P(|1⟩) on q0")
    ax.plot(np.degrees(d_fine), leakage_curve, "b--", lw=1.5,
            label=f"Leakage fit: {a_fit:.4f} + {b_fit:.4f}·sin²(δ/2)")
    ax.axhline(a_fit, color="gray", lw=1, ls=":", label=f"Baseline a={a_fit:.4f}")
    ax.set_xlabel("Ry angle on q1 (degrees)")
    ax.set_ylabel("P(|1⟩) on q0")
    ax.set_title("Crosstalk check: Ry(q1) → q0 leakage")
    ax.set_ylim(-0.01, max(0.08, float(np.max(p1)) * 1.4))
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = f"crosstalk_results_{ts}.png"
    plt.savefig(out, dpi=150)
    print(f"  Plot: {out}")

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if args.results:
        with open(args.results) as f:
            data = json.load(f)
        print(f"Loaded {args.results}")
    elif args.dry_run:
        print("DRY RUN — no crosstalk simulation")
        data = simulate_locally(ANGLES_RAD, args.shots)
    else:
        if not args.token:
            raise ValueError("Provide --token or use --dry-run")
        data = run_on_ibm(ANGLES_RAD, args.shots, args.token, args.backend)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_json = f"crosstalk_raw_{ts}.json"
    with open(out_json, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Raw results saved: {out_json}")

    angles_deg, p1, err, a_fit, b_fit = analyze(data)
    plot(angles_deg, p1, err, a_fit, b_fit, ts)

if __name__ == "__main__":
    main()
