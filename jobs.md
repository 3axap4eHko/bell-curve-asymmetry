# IBM Quantum Job IDs

All experiments run on 2026-03-10 via the `ibm_quantum_platform` channel.

## ibm_fez

| Job ID | Experiment | Script | Data file |
|---|---|---|---|
| d6o91869td6c73ap2mig | Bell sweep (raw only, no readout correction) | spin_angle_sweep.py | spin_sweep_raw_20260310_174546.json |
| d6o99u8fh9oc73eq2pk0 | Bell sweep (with readout correction) | spin_angle_sweep.py | spin_sweep_raw_20260310_180445.json |

## ibm_marrakesh

| Job ID | Experiment | Script | Data file |
|---|---|---|---|
| d6o9cqobfi7c73a5kq0g | Bell sweep (with readout correction) | spin_angle_sweep.py | spin_sweep_raw_20260310_181014.json |
| d6o9h68bfi7c73a5kvi0 | Ry gate calibration | ry_calibration.py | ry_cal_raw_20260310_181931.json |
| d6o9pb8fh9oc73eq3dfg | Crosstalk check | crosstalk_check.py | crosstalk_raw_20260310_183555.json |

## Local simulations (no job ID)

| Timestamp | Experiment | Script | Data file |
|---|---|---|---|
| 20260310_174142 | Bell sweep dry-run (ideal QM + shot noise) | spin_angle_sweep.py --dry-run | spin_sweep_raw_20260310_174142.json |
| 20260310_180027 | Bell sweep dry-run (ideal QM + shot noise) | spin_angle_sweep.py --dry-run | spin_sweep_raw_20260310_180027.json |

## Local verification artifacts

- `data/fakemarrakesh_transpile_counts.json` -- FakeMarrakesh transpilation audit backing Appendix B / Section 5.4
- `data/t2_exact_verification.json` -- exact Aer density-matrix verification backing Section 5.3

## Data file mapping in this repo

The raw source files above (from the `quantum/` working directory) were processed into the `data/` files:

- `bell_sweep_fez_raw.json` -- raw P_disagree extracted from d6o99u8fh9oc73eq2pk0 (same job as corrected)
- `bell_sweep_fez_corrected.json` -- readout-corrected P_disagree from d6o99u8fh9oc73eq2pk0
- `bell_sweep_marrakesh_raw.json` -- raw P_disagree from d6o9cqobfi7c73a5kq0g
- `bell_sweep_marrakesh_corrected.json` -- readout-corrected P_disagree from d6o9cqobfi7c73a5kq0g
- `fakemarrakesh_transpile_counts.json` -- generated locally from FakeMarrakesh at optimization level 1
- `ry_cal_raw_20260310_181931.json` -- copied directly
- `crosstalk_raw_20260310_183555.json` -- copied directly
- `t2_exact_verification.json` -- generated locally from exact Aer density-matrix simulation
