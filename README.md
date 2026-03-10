# Bell Curve Asymmetry

Anomalous curve-shape deviation from the quantum singlet correlation in full angle-sweep Bell measurements on IBM superconducting qubits.

Measurements on **ibm_fez** and **ibm_marrakesh** yield a systematic crossover shift: the quantum model with free visibility gives χ²/dof = **4.83** on ibm_marrakesh; a two-parameter deformation model gives χ²/dof = **1.41**, a Δχ² = 124.5 improvement for one additional degree of freedom. The crossover-shift parameter α = 0.470 ± 0.003 is consistent across both chips after readout error correction (Δχ² = 124.5, 1 dof, p < 10⁻²⁸, likelihood ratio test). Five candidate instrumental explanations are individually tested and none found sufficient to account for the anomaly.

## Key results

| Backend | Correction | QM+V χ²/dof | Model 2 α | Model 2 χ²/dof | Model 3 α | Model 3 χ²/dof |
|---|---|---|---|---|---|---|
| ibm_fez | Raw | 4.97 | 0.468 ± 0.003 | 2.11 | 0.480 ± 0.002 | 2.17 |
| ibm_fez | Corrected | 7.93 | 0.467 ± 0.003 | 4.10 | 0.480 ± 0.002 | 4.49 |
| ibm_marrakesh | Raw | 4.52 | 0.470 ± 0.003 | 1.22 | 0.481 ± 0.002 | 1.41 |
| ibm_marrakesh | Corrected | **4.83** | **0.470 ± 0.003** | **1.41** | 0.482 ± 0.002 | 1.62 |

Δχ² = 124.5 for 1 extra degree of freedom on ibm_marrakesh (corrected). p < 10^-28.

## Systematic error candidates

| Candidate | Result | Status |
|---|---|---|
| Readout asymmetry | Δα < 0.001 after 4×4 correction | No significant effect detected |
| Ry gate offset | ε = 1.12° ± 0.12°, explains ~20% | Insufficient |
| T2 decoherence | No systematic α shift; visibility-only effect | Not supported |
| Angle-dependent gate duration | 0.03% effect vs ~2.7% observed | No significant effect detected* |
| Qubit-qubit crosstalk | b = 0.0000 ± 0.0002, absent | No significant effect detected |

\* FakeMarrakesh approximation; pulse-level schedule not independently verified.

## Repository structure

```
bell-curve-asymmetry/
  README.md                  # this file
  paper/
    paper_draft.md           # full paper draft
  data/
    bell_sweep_fez_raw.json
    bell_sweep_fez_corrected.json
    bell_sweep_marrakesh_raw.json
    bell_sweep_marrakesh_corrected.json
    fakemarrakesh_transpile_counts.json
    ry_cal_raw_20260310_181931.json
    crosstalk_raw_20260310_183555.json
    t2_exact_verification.json
  scripts/
    spin_angle_sweep.py      # Bell sweep experiment
    ry_calibration.py        # Ry gate offset check
    crosstalk_check.py       # qubit-qubit crosstalk check
    t2_noise_model.py        # T2 decoherence simulation
    verify_nonqpu.py         # local non-QPU verification
  figures/
    *.png                    # all plots
  jobs.md                    # IBM job IDs mapped to experiments
```

## How to reproduce

### Requirements

```
pip install qiskit qiskit-ibm-runtime qiskit-aer scipy numpy matplotlib
```

### Bell sweep on real hardware

```bash
python scripts/spin_angle_sweep.py --token YOUR_IBM_TOKEN --backend ibm_fez --shots 8192
```

### Analyze existing data

```bash
python scripts/spin_angle_sweep.py --results data/bell_sweep_fez_raw.json
```

### Run systematic checks

```bash
# Ry gate calibration
python scripts/ry_calibration.py --token YOUR_IBM_TOKEN --backend ibm_marrakesh

# Crosstalk isolation
python scripts/crosstalk_check.py --token YOUR_IBM_TOKEN --backend ibm_marrakesh

# T2 noise simulation (no token needed)
python scripts/t2_noise_model.py
python scripts/t2_noise_model.py --sweep-t2

# Full local verification (no QPU)
python scripts/verify_nonqpu.py
```

### Dry-run (no IBM account needed)

```bash
python scripts/spin_angle_sweep.py --dry-run
python scripts/ry_calibration.py --dry-run
python scripts/crosstalk_check.py --dry-run
```

## Hardware

- **ibm_fez** and **ibm_marrakesh**: IBM superconducting transmon processors
- Accessed via IBM Quantum Platform (ibm_quantum_platform channel)
- Qiskit IBM Runtime SDK

## IBM Quantum job IDs

See [jobs.md](jobs.md) for full mapping of job IDs to experiments and data files.
