#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import DensityMatrix
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, thermal_relaxation_error
from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
from scipy.optimize import curve_fit
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
ANGLES_DEG = np.arange(0, 185, 5)
ANGLES_RAD = np.deg2rad(ANGLES_DEG)


BELL_EXPECTED = {
    "fez_raw": {
        "file": "bell_sweep_fez_raw.json",
        "p_key": "p_disagree",
        "err_key": "err",
        "qmv_chi2_dof": 4.965796954820693,
        "model2_alpha": 0.4679715006836887,
        "model2_chi2_dof": 2.107225970981562,
        "model3_alpha": 0.47978808094033787,
        "model3_chi2_dof": 2.171565739547133,
    },
    "fez_corrected": {
        "file": "bell_sweep_fez_corrected.json",
        "p_key": "p_disagree_corr",
        "err_key": "err_corr",
        "qmv_chi2_dof": 7.934592716855572,
        "model2_alpha": 0.46709494965426324,
        "model2_chi2_dof": 4.095874452561847,
        "model3_alpha": 0.4804383957271304,
        "model3_chi2_dof": 4.485075751068294,
    },
    "marrakesh_raw": {
        "file": "bell_sweep_marrakesh_raw.json",
        "p_key": "p_disagree",
        "err_key": "err",
        "qmv_chi2_dof": 4.521727942947301,
        "model2_alpha": 0.46951290671019863,
        "model2_chi2_dof": 1.220398725693709,
        "model3_alpha": 0.48136919466998707,
        "model3_chi2_dof": 1.4103418875528946,
    },
    "marrakesh_corrected": {
        "file": "bell_sweep_marrakesh_corrected.json",
        "p_key": "p_disagree_corr",
        "err_key": "err_corr",
        "qmv_chi2_dof": 4.825223894800445,
        "model2_alpha": 0.46966504979633084,
        "model2_chi2_dof": 1.4062822266909532,
        "model3_alpha": 0.48162188935780526,
        "model3_chi2_dof": 1.6173902532790014,
        "lrt_dchi2": 124.48818227863265,
        "lrt_p": 6.587027514556295e-29,
    },
}


RY_EXPECTED = {
    "epsilon_deg": 1.1194014235144487,
    "epsilon_err_deg": 0.11654028171856091,
    "visibility": 0.9820697860056928,
    "delta_alpha": -0.006218896797302493,
    "fraction_of_observed": 0.2072965599100831,
}


CROSSTALK_EXPECTED = {
    "b": 4.474691279671446e-05,
    "b_err": 0.00015440808951447816,
    "z": 0.2897964280072175,
    "stray_counts": 5,
    "total_shots": 77824,
}


READOUT_EXPECTED = {
    "q0_e01": 0.00244140625,
    "q0_e10": 0.0068359375,
    "q1_e01": 0.0009765625,
    "q1_e10": 0.005126953125,
}


def load_json(name):
    with open(DATA_DIR / name) as handle:
        return json.load(handle)


def qm_with_visibility(delta, visibility):
    return visibility * np.cos(delta / 2) ** 2 + (1 - visibility) * 0.5


def model2(delta, alpha, visibility):
    return visibility * (np.cos(delta / 2) ** 2 + (2 * alpha - 1) * np.sin(delta) ** 2 / 4) + (1 - visibility) * 0.5


def model3(delta, alpha, visibility):
    base = np.clip(np.cos(delta / 2) ** 2, 1e-12, 1.0)
    return visibility * np.power(base, 1.0 / (2 * alpha)) + (1 - visibility) * 0.5


def fit_qmv(angles, values, errors):
    popt, pcov = curve_fit(
        qm_with_visibility,
        angles,
        values,
        p0=[0.9],
        sigma=errors,
        absolute_sigma=True,
        bounds=([0.5], [1.0]),
        maxfev=10000,
    )
    visibility = float(popt[0])
    residuals = values - qm_with_visibility(angles, visibility)
    chi2_value = float(np.sum((residuals / errors) ** 2))
    dof = len(angles) - 1
    return {
        "visibility": visibility,
        "visibility_err": float(np.sqrt(pcov[0, 0])),
        "chi2": chi2_value,
        "chi2_dof": chi2_value / dof,
    }


def fit_model(model_fn, angles, values, errors):
    popt, pcov = curve_fit(
        model_fn,
        angles,
        values,
        p0=[0.5, 0.9],
        sigma=errors,
        absolute_sigma=True,
        bounds=([0.3, 0.5], [1.0, 1.0]),
        maxfev=10000,
    )
    alpha = float(popt[0])
    visibility = float(popt[1])
    residuals = values - model_fn(angles, alpha, visibility)
    chi2_value = float(np.sum((residuals / errors) ** 2))
    dof = len(angles) - 2
    return {
        "alpha": alpha,
        "alpha_err": float(np.sqrt(pcov[0, 0])),
        "visibility": visibility,
        "visibility_err": float(np.sqrt(pcov[1, 1])),
        "chi2": chi2_value,
        "chi2_dof": chi2_value / dof,
    }


def fit_bell_dataset(spec):
    data = load_json(spec["file"])
    angles = np.deg2rad(np.array([row["angle_deg"] for row in data], dtype=float))
    values = np.array([row[spec["p_key"]] for row in data], dtype=float)
    errors = np.maximum(np.array([row[spec["err_key"]] for row in data], dtype=float), 1.0 / data[0]["shots"])
    qmv = fit_qmv(angles, values, errors)
    m2 = fit_model(model2, angles, values, errors)
    m3 = fit_model(model3, angles, values, errors)
    result = {"qmv": qmv, "model2": m2, "model3": m3}
    if "lrt_dchi2" in spec:
        dchi2 = qmv["chi2"] - m2["chi2"]
        result["lrt"] = {"dchi2": dchi2, "p": float(chi2.sf(dchi2, 1))}
    return result


def fit_ry():
    data = load_json("ry_cal_raw_20260310_181931.json")
    angles = np.deg2rad(np.array([row["angle_deg"] for row in data], dtype=float))
    values = np.array([row["p1"] for row in data], dtype=float)
    errors = np.maximum(np.array([row["err"] for row in data], dtype=float), 1.0 / data[0]["shots"])

    def with_offset_and_scale(delta, epsilon, visibility):
        return visibility * np.sin((delta + epsilon) / 2) ** 2 + (1 - visibility) * 0.5

    (epsilon, visibility), pcov = curve_fit(
        with_offset_and_scale,
        angles,
        values,
        p0=[0.0, 0.9],
        sigma=errors,
        absolute_sigma=True,
        bounds=([-np.pi / 4, 0.5], [np.pi / 4, 1.0]),
    )
    epsilon = float(epsilon)
    return {
        "epsilon_deg": math.degrees(epsilon),
        "epsilon_err_deg": math.degrees(float(np.sqrt(pcov[0, 0]))),
        "visibility": float(visibility),
        "delta_alpha": -epsilon / math.pi,
        "fraction_of_observed": abs(epsilon / math.pi) / 0.03,
    }


def fit_crosstalk():
    data = load_json("crosstalk_raw_20260310_183555.json")
    angles = np.deg2rad(np.array([row["angle_deg"] for row in data], dtype=float))
    values = np.array([row["p1_q0"] for row in data], dtype=float)
    errors = np.maximum(np.array([row["err"] for row in data], dtype=float), 1.0 / data[0]["shots"])

    def leakage_model(delta, offset, leak):
        return offset + leak * np.sin(delta / 2) ** 2

    (offset, leak), pcov = curve_fit(
        leakage_model,
        angles,
        values,
        p0=[float(np.mean(values)), 0.0],
        sigma=errors,
        absolute_sigma=True,
    )
    leak_err = float(np.sqrt(pcov[1, 1]))
    return {
        "offset": float(offset),
        "b": float(leak),
        "b_err": leak_err,
        "z": abs(float(leak)) / leak_err,
        "stray_counts": int(round(float(np.sum(values)) * data[0]["shots"])),
        "total_shots": data[0]["shots"] * len(data),
    }


def build_transpile_artifact():
    backend = FakeMarrakesh()
    pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=1)
    rows = []
    for angle in [0, 30, 60, 90, 120, 150, 180]:
        circuit = QuantumCircuit(2, 2)
        circuit.x(1)
        circuit.h(0)
        circuit.cx(0, 1)
        circuit.z(0)
        circuit.ry(math.radians(angle), 1)
        circuit.measure([0, 1], [0, 1])
        transpiled = pass_manager.run(circuit)
        counts = transpiled.count_ops()
        sx_count = int(counts.get("sx", 0))
        x_count = int(counts.get("x", 0))
        rows.append(
            {
                "angle_deg": angle,
                "sx_count": sx_count,
                "x_count": x_count,
                "rz_count": int(counts.get("rz", 0)),
                "cz_count": int(counts.get("cz", 0)),
                "single_qubit_duration_ns": 36 * (sx_count + x_count),
            }
        )
    return {
        "backend": "FakeMarrakesh",
        "optimization_level": 1,
        "basis_gates": backend.configuration().basis_gates,
        "durations_ns": {
            "sx": backend.target["sx"][(0,)].duration * 1e9,
            "x": backend.target["x"][(0,)].duration * 1e9,
            "rz": backend.target["rz"][(0,)].duration * 1e9,
            "cz": backend.target["cz"][(0, 1)].duration * 1e9,
        },
        "angles_deg": [0, 30, 60, 90, 120, 150, 180],
        "counts": rows,
        "max_single_qubit_duration_variation_ns": max(row["single_qubit_duration_ns"] for row in rows) - min(row["single_qubit_duration_ns"] for row in rows),
        "differential_t2_factor_at_110us": math.exp(-(36e-9) / (110e-6)),
    }


def build_noise_model(t1_q0, t2_q0, t1_q1, t2_q1, sx_ns, cz_ns):
    t1 = [t1_q0 * 1e-6, t1_q1 * 1e-6]
    t2 = [min(t2_q0 * 1e-6, 2 * t1[0]), min(t2_q1 * 1e-6, 2 * t1[1])]
    noise_model = NoiseModel()
    sx_duration = sx_ns * 1e-9
    cz_duration = cz_ns * 1e-9
    for qubit in range(2):
        single = thermal_relaxation_error(t1[qubit], t2[qubit], sx_duration)
        noise_model.add_quantum_error(single, ["rx", "ry", "rz", "sx", "x", "h", "z"], [qubit])
    pair = thermal_relaxation_error(t1[0], t2[0], cz_duration).expand(thermal_relaxation_error(t1[1], t2[1], cz_duration))
    noise_model.add_quantum_error(pair, ["cx"], [0, 1])
    return noise_model


def build_density_matrix_circuit(angle):
    circuit = QuantumCircuit(2)
    circuit.x(1)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.z(0)
    circuit.ry(angle, 1)
    circuit.save_density_matrix()
    return circuit


def exact_t2_curve(noise_model):
    simulator = AerSimulator(method="density_matrix", noise_model=noise_model)
    circuits = transpile([build_density_matrix_circuit(angle) for angle in ANGLES_RAD], simulator)
    result = simulator.run(circuits).result()
    values = []
    for index in range(len(circuits)):
        density_matrix = DensityMatrix(result.data(index)["density_matrix"])
        probabilities = density_matrix.probabilities_dict()
        values.append(float(probabilities.get("01", 0.0) + probabilities.get("10", 0.0)))
    return np.array(values)


def build_t2_artifact():
    params = {
        "t1_q0_us": 250.0,
        "t2_q0_us": 130.0,
        "t1_q1_us": 220.0,
        "t2_q1_us": 110.0,
        "sx_ns": 36.0,
        "cz_ns": 68.0,
    }
    base_curve = exact_t2_curve(build_noise_model(250.0, 130.0, 220.0, 110.0, 36.0, 68.0))
    base_fit = fit_model(model2, ANGLES_RAD, base_curve, np.maximum(np.sqrt(base_curve * (1 - base_curve) / 8192), 1.0 / 8192))
    base_fit_model3 = fit_model(model3, ANGLES_RAD, base_curve, np.maximum(np.sqrt(base_curve * (1 - base_curve) / 8192), 1.0 / 8192))
    sweep = []
    for t2_value in range(50, 301, 25):
        curve = exact_t2_curve(build_noise_model(250.0, t2_value, 220.0, t2_value, 36.0, 68.0))
        fit2 = fit_model(model2, ANGLES_RAD, curve, np.maximum(np.sqrt(curve * (1 - curve) / 8192), 1.0 / 8192))
        fit3 = fit_model(model3, ANGLES_RAD, curve, np.maximum(np.sqrt(curve * (1 - curve) / 8192), 1.0 / 8192))
        sweep.append({"t2_us": t2_value, "model2_alpha": fit2["alpha"], "model3_alpha": fit3["alpha"]})
    return {
        "method": "Aer density_matrix exact evolution",
        "params": params,
        "angles_deg": ANGLES_DEG.tolist(),
        "base_fit": {
            "model2_alpha": base_fit["alpha"],
            "model2_visibility": base_fit["visibility"],
            "model3_alpha": base_fit_model3["alpha"],
            "model3_visibility": base_fit_model3["visibility"],
            "max_abs_deviation_from_qm": float(np.max(np.abs(base_curve - np.cos(ANGLES_RAD / 2) ** 2))),
        },
        "sweep": sweep,
    }


def compare(actual, expected, tolerance, label, failures):
    if abs(actual - expected) > tolerance:
        failures.append(f"{label}: expected {expected}, got {actual}")


def compare_nested(actual, expected, tolerance, prefix, failures):
    for key, expected_value in expected.items():
        actual_value = actual[key]
        if isinstance(expected_value, dict):
            compare_nested(actual_value, expected_value, tolerance, f"{prefix}.{key}", failures)
        elif isinstance(expected_value, list):
            if len(actual_value) != len(expected_value):
                failures.append(f"{prefix}.{key}: expected length {len(expected_value)}, got {len(actual_value)}")
                continue
            for index, item in enumerate(expected_value):
                if isinstance(item, dict):
                    compare_nested(actual_value[index], item, tolerance, f"{prefix}.{key}[{index}]", failures)
                elif isinstance(item, str):
                    if actual_value[index] != item:
                        failures.append(f"{prefix}.{key}[{index}]: expected {item}, got {actual_value[index]}")
                else:
                    compare(actual_value[index], item, tolerance, f"{prefix}.{key}[{index}]", failures)
        elif isinstance(expected_value, str):
            if actual_value != expected_value:
                failures.append(f"{prefix}.{key}: expected {expected_value}, got {actual_value}")
        else:
            compare(actual_value, expected_value, tolerance, f"{prefix}.{key}", failures)


def build_report():
    bell = {}
    failures = []
    for label, spec in BELL_EXPECTED.items():
        fit = fit_bell_dataset(spec)
        bell[label] = fit
        compare(fit["qmv"]["chi2_dof"], spec["qmv_chi2_dof"], 1e-9, f"{label}.qmv_chi2_dof", failures)
        compare(fit["model2"]["alpha"], spec["model2_alpha"], 1e-9, f"{label}.model2_alpha", failures)
        compare(fit["model2"]["chi2_dof"], spec["model2_chi2_dof"], 1e-9, f"{label}.model2_chi2_dof", failures)
        compare(fit["model3"]["alpha"], spec["model3_alpha"], 1e-9, f"{label}.model3_alpha", failures)
        compare(fit["model3"]["chi2_dof"], spec["model3_chi2_dof"], 1e-9, f"{label}.model3_chi2_dof", failures)
        if "lrt_dchi2" in spec:
            compare(fit["lrt"]["dchi2"], spec["lrt_dchi2"], 1e-9, f"{label}.lrt_dchi2", failures)
            compare(fit["lrt"]["p"], spec["lrt_p"], 1e-36, f"{label}.lrt_p", failures)

    marrakesh_corrected = load_json("bell_sweep_marrakesh_corrected.json")
    cal = marrakesh_corrected[0]["cal"]
    compare(cal["q0"]["e01"], READOUT_EXPECTED["q0_e01"], 1e-12, "readout.q0_e01", failures)
    compare(cal["q0"]["e10"], READOUT_EXPECTED["q0_e10"], 1e-12, "readout.q0_e10", failures)
    compare(cal["q1"]["e01"], READOUT_EXPECTED["q1_e01"], 1e-12, "readout.q1_e01", failures)
    compare(cal["q1"]["e10"], READOUT_EXPECTED["q1_e10"], 1e-12, "readout.q1_e10", failures)

    ry = fit_ry()
    compare(ry["epsilon_deg"], RY_EXPECTED["epsilon_deg"], 1e-9, "ry.epsilon_deg", failures)
    compare(ry["epsilon_err_deg"], RY_EXPECTED["epsilon_err_deg"], 1e-9, "ry.epsilon_err_deg", failures)
    compare(ry["visibility"], RY_EXPECTED["visibility"], 1e-9, "ry.visibility", failures)
    compare(ry["delta_alpha"], RY_EXPECTED["delta_alpha"], 1e-12, "ry.delta_alpha", failures)
    compare(ry["fraction_of_observed"], RY_EXPECTED["fraction_of_observed"], 1e-12, "ry.fraction_of_observed", failures)

    crosstalk = fit_crosstalk()
    compare(crosstalk["b"], CROSSTALK_EXPECTED["b"], 1e-12, "crosstalk.b", failures)
    compare(crosstalk["b_err"], CROSSTALK_EXPECTED["b_err"], 1e-12, "crosstalk.b_err", failures)
    compare(crosstalk["z"], CROSSTALK_EXPECTED["z"], 1e-12, "crosstalk.z", failures)
    compare(crosstalk["stray_counts"], CROSSTALK_EXPECTED["stray_counts"], 0, "crosstalk.stray_counts", failures)
    compare(crosstalk["total_shots"], CROSSTALK_EXPECTED["total_shots"], 0, "crosstalk.total_shots", failures)

    transpile_actual = build_transpile_artifact()
    transpile_expected = load_json("fakemarrakesh_transpile_counts.json")
    compare_nested(transpile_actual, transpile_expected, 1e-9, "transpile", failures)

    t2_actual = build_t2_artifact()
    t2_expected = load_json("t2_exact_verification.json")
    compare_nested(t2_actual, t2_expected, 1e-9, "t2", failures)

    return {
        "bell": bell,
        "readout": cal,
        "ry": ry,
        "crosstalk": crosstalk,
        "transpile": transpile_actual,
        "t2": t2_actual,
        "failures": failures,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    report = build_report()
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("Bell fits: OK")
        print("Readout calibration: OK")
        print("Ry calibration: OK")
        print("Crosstalk isolation: OK")
        print("FakeMarrakesh transpilation: OK")
        print("T2 exact verification: OK")
        if report["failures"]:
            print("\nMismatches:")
            for failure in report["failures"]:
                print(f"  - {failure}")
        else:
            print("\nAll non-QPU checks passed.")
    raise SystemExit(1 if report["failures"] else 0)


if __name__ == "__main__":
    main()
