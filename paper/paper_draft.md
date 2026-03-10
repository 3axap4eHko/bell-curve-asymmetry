# Anomalous Curve-Shape Deviation in Singlet Bell Measurements on IBM Superconducting Qubits

**Ivan Zakharchanka**  
Independent Researcher, Aventura, FL, USA  
March 2026

---

## Abstract

We report an anomalous deviation from the quantum mechanical singlet correlation curve in full angle-sweep Bell measurements on two IBM superconducting quantum processors (ibm_fez and ibm_marrakesh). With Alice fixed along Z and Bob sweeping Ry(δ) from 0° to 180°, the disagreement probability P_disagree(δ) deviates from the predicted cos²(δ/2) in a way that cannot be explained by uniform depolarizing noise: the 50% crossover point shifts below 90°. Fitting a two-parameter deformation model yields a crossover-shift parameter α = 0.470 ± 0.003 (Model 2) on ibm_marrakesh after readout correction, compared to the quantum prediction α = 0.5. The quantum model with a free visibility parameter gives χ²/dof = 4.83; the deformation model gives χ²/dof = 1.41, a Δχ² = 124.5 improvement for one additional degree of freedom (p < 10^-28, likelihood ratio test). The deviation is consistent across two physically independent chips. We systematically investigate five candidate instrumental explanations — readout asymmetry, Ry gate calibration offset, T2 thermal decoherence, angle-dependent gate duration variation, and qubit-qubit crosstalk — and find each insufficient through dedicated experiments or analytical arguments. The anomaly is currently unexplained by standard NISQ noise mechanisms and warrants independent replication on different hardware platforms.

---

## 1. Introduction

Superconducting qubit platforms have demonstrated Bell inequality violations at high fidelity [1–5], typically at a small number of discrete measurement angles. The full continuous shape of the disagreement probability curve P_disagree(δ) as a function of measurement angle has received less attention as a precision test of quantum mechanics on near-term hardware.

For a singlet state |Ψ⁻⟩ = (|01⟩ − |10⟩)/√2 with Alice measuring along Z and Bob applying Ry(δ) before Z measurement, the quantum prediction is:

$$P_\text{disagree}(\delta) = \cos^2\!\left(\frac{\delta}{2}\right)$$

This curve has a specific shape: it reaches 50% exactly at δ = 90°, and is symmetric around that point. Uniform depolarizing noise rescales the amplitude toward 0.5 but cannot shift the crossover point. A crossover shift therefore represents a qualitatively different kind of deviation from the quantum prediction — one that a scalar visibility parameter cannot absorb.

We measure P_disagree(δ) at 37 angles from 0° to 180° on two separate IBM processors and observe a systematic crossover shift that survives readout error correction and is reproducible across both chips. This paper reports the anomaly, quantifies its statistical significance, and presents a systematic investigation of instrumental explanations.

---

## 2. Parameterization

To quantify the crossover shift we use two deformation models, both parameterized by α ∈ [0.3, 1.0] and visibility V ∈ [0.5, 1.0]:

**Model 2 (perturbative):**
$$P = V\!\left[\cos^2\!\left(\frac{\delta}{2}\right) + (2\alpha-1)\frac{\sin^2\delta}{4}\right] + (1-V)\cdot 0.5$$

**Model 3 (power-law):**
$$P = V\!\left[\cos^2\!\left(\frac{\delta}{2}\right)\right]^{1/(2\alpha)} + (1-V)\cdot 0.5$$

Both recover the quantum prediction at α = 0.5 with any V. The quantum prediction with free visibility is the nested special case α = 0.5. For α < 0.5, the effective crossover shifts below 90° — disagreements are more likely at small angles than quantum mechanics predicts.

The two models bracket different deformation shapes: Model 2 has maximum deviation at δ = 90°; Model 3 applies a power-law modification across the full curve. We report both as consistency checks.

These parameterizations provide a quantitative handle on crossover-shift deviations. No physical interpretation of α is assumed in the main analysis; a possible geometric interpretation is discussed in Section 6.2.

---

## 3. Experimental Methods

### 3.1 Circuit

The singlet state |Ψ⁻⟩ is prepared on two qubits (q0 = Alice, q1 = Bob) via:

```
q0: ─────── H ── ●  ── Z ── M
q1: ─── X ───── ⊕ ──────── M
```

Alice measures along Z with no additional rotation. Bob applies Ry(δ) then measures Z. Disagreement outcomes are |01⟩ and |10⟩. The Z gate on q0 is essential to produce |Ψ⁻⟩; omitting it yields |Ψ⁺⟩ and an inverted curve — this error was caught and corrected during initial runs.

### 3.2 Angle Sweep

37 angles from 0° to 180° in 5° steps. 8192 shots per angle. All circuits submitted in a single job to minimize inter-circuit calibration drift.

### 3.3 Hardware

Experiments were performed on **ibm_fez** and **ibm_marrakesh** via the IBM Quantum Platform, accessed through the Qiskit IBM Runtime SDK. Both are superconducting transmon processors with independent fabrication, calibration histories, and qubit quality metrics. Physical qubit assignments were determined by the Qiskit transpiler at optimization level 1.

### 3.4 Readout Error Correction

Readout calibration values (e₀₁ = prob_meas1_prep0, e₁₀ = prob_meas0_prep1) were extracted from backend properties. A 4×4 assignment matrix A = M_q1 ⊗ M_q0 was constructed and inverted, then applied to the raw count vector at each angle. Negative corrected counts were clipped to zero and renormalized.

Measured readout errors on ibm_marrakesh: q0: e₀₁ = 0.244%, e₁₀ = 0.684%; q1: e₀₁ = 0.098%, e₁₀ = 0.513%.

Post-correction statistical uncertainties are approximated as binomial on the corrected counts. Full error propagation through the 4x4 inverse matrix would yield slightly different uncertainties, but this is a minor correction given the small readout error rates (<1%).

Readout calibration was extracted from the transpiled layout of the first circuit only. While all 37 circuits were submitted in a single job (making layout drift unlikely), identical physical qubit assignments across all circuits were not independently verified.

### 3.5 Fitting

All fits use scipy.optimize.curve_fit with bounds α ∈ [0.3, 1.0], V ∈ [0.5, 1.0]. Statistical uncertainties are from the covariance matrix diagonal. χ²/dof uses per-angle binomial errors floored at 1/shots.

---

## 4. Results

### 4.1 Model Comparison

Table 1 reports fit results on both backends before and after readout correction, alongside the quantum prediction (α = 0.5 fixed, V free).

**Table 1: Fit results — quantum prediction vs deformation models**

| Backend | Correction | QM+V χ²/dof | Model 2 α | Model 2 χ²/dof | Model 3 α | Model 3 χ²/dof |
|---|---|---|---|---|---|---|
| ibm_fez | Raw | 4.97 | 0.468 ± 0.003 | 2.11 | 0.480 ± 0.002 | 2.17 |
| ibm_fez | Corrected | 7.93 | 0.467 ± 0.003 | 4.10 | 0.480 ± 0.002 | 4.49 |
| ibm_marrakesh | Raw | 4.52 | 0.470 ± 0.003 | 1.22 | 0.481 ± 0.002 | 1.41 |
| ibm_marrakesh | Corrected | 4.83 | 0.470 ± 0.003 | 1.41 | 0.482 ± 0.002 | 1.62 |

On ibm_marrakesh (corrected): the quantum model with free visibility gives χ²/dof = 4.83. Adding α as a free parameter reduces this to χ²/dof = 1.41. The improvement is Δχ² = 124.5 for one additional degree of freedom (p < 10^-28, likelihood ratio test). The deformation model is strongly preferred by the data on the higher-quality chip.

Key observations:

- α is consistent across both chips (corrected): 0.467–0.470 (Model 2), 0.480–0.482 (Model 3).
- Readout correction shifts α by < 0.001 — correction primarily improves visibility (0.88 → 0.97 on fez, 0.96 → 0.97 on marrakesh) as expected for uniform noise.
- χ²/dof = 1.41 on marrakesh indicates the deformation model is a good fit; χ²/dof = 4.10 on fez indicates residual structured noise on the lower-quality chip.

### 4.2 Crossover Point

For Model 2 with α = 0.470, the crossover (P_disagree = 0.5) shifts from δ = 90° to approximately δ = 88.3°. This is a 1.7° shift — small in absolute terms.

---

## 5. Systematic Error Analysis

A crossover shift requires a mechanism that differentially affects measurements at small vs. large angles — uniform noise cannot produce it. We investigate five specific candidates.

### 5.1 Readout Asymmetry

**Hypothesis:** Asymmetric readout errors (e₀₁ ≠ e₁₀) bias P_disagree in an angle-dependent way.

**Method:** Full 4×4 readout correction matrix applied using measured calibration values.

**Result:** α shifts by < 0.001 after correction. Visibility improves (0.88 → 0.97 on fez, 0.96 → 0.97 on marrakesh) as expected. The crossover shift is unchanged.

**Status: No significant effect detected.**

### 5.2 Ry Gate Calibration Offset

**Hypothesis:** If Ry(δ) executes Ry(δ + ε) due to miscalibration, the Bell sweep effectively measures cos²((δ+ε)/2), shifting the crossover.

**Method:** Single-qubit calibration circuit using the same nominal 2-qubit circuit topology as the Bell sweep. Ry(δ) applied to q1, q1 measured. 37 angles, 8192 shots, on ibm_marrakesh. P(|1⟩) fitted to sin²((δ + ε)/2). Note: `initial_layout` was not pinned in the control experiment; physical qubit assignment was not independently verified to match the Bell sweep job.

**Result:** Best-fit ε = 1.12° ± 0.12° (offset + visibility model, V = 0.982). Implied α shift: Δα ≈ −ε_rad/π ≈ −0.006. Observed deviation: Δα ≈ −0.030. Gate offset explains approximately 20% of the deviation, with a mechanistically distinct residual shape (T1-driven amplitude suppression at large angles, not a crossover shift).

**Status: Contributes ~20% of deviation. Insufficient to explain the full signal.**

### 5.3 T2 Thermal Decoherence

**Hypothesis:** T2 dephasing during the circuit creates angle-dependent visibility that deforms the curve and produces apparent α < 0.5.

**Method:** Full circuit simulation via Qiskit Aer with a thermal relaxation noise model at published ibm_marrakesh parameters (T1_q0 = 250 µs, T2_q0 = 130 µs, T1_q1 = 220 µs, T2_q1 = 110 µs; CZ gate time 68 ns, SX gate time 36 ns). T2 swept from 50 µs to 300 µs.

**Result:** Ideal simulation yields α consistent with 0.5. Under the fixed-duration thermal relaxation channel at marrakesh parameters, the fitted α remains consistent with 0.5; the dominant effect is a small reduction in visibility. A T2 sweep from 50 to 300 µs shows no systematic trend away from α = 0.5 and no tendency toward α < 0.5 at any tested T2 value. Any one-shot finite-sampling deviation from 0.5 is consistent with shot noise rather than a channel-induced crossover shift.

This is mechanistically expected: T2 dephasing applies uniform decoherence, captured entirely by V. A scalar visibility cannot shift the crossover. Only a mechanism with differential angular sensitivity can produce α ≠ 0.5 — and T2 is not such a mechanism.

**Status: Not supported. No systematic α shift at physically realistic T2 values.**

### 5.4 Angle-Dependent Gate Duration

**Hypothesis:** If Ry(δ) gate duration scales with δ, larger angles accumulate more T2 exposure, producing angle-dependent decoherence that could shift the crossover.

**Method:** Transpilation analysis using FakeMarrakesh (matching the ibm_marrakesh basis gate set: {sx, rz, cz, x}). Bell sweep circuit transpiled at each angle in [0°, 30°, 60°, 90°, 120°, 150°, 180°]. Gate counts and durations extracted.

**Result:** RZ is a virtual gate (zero duration; implemented as a software frame rotation). SX gates have fixed duration 36 ns regardless of rotation angle. The transpiled SX count ranges from 2 to 4 across all angles (step function, not smooth in δ), giving a maximum duration variation of 36 ns. At T2 = 110 µs, the differential decoherence factor is exp(−36×10⁻⁹ / 110×10⁻⁶) = 0.9997 — a 0.03% effect against a ~2.7% observed data-vs-QM deviation. Furthermore, the SX count has no smooth angular dependence and cannot produce the smooth crossover shift observed.

Note: this analysis used FakeMarrakesh to approximate the ibm_marrakesh basis gate decomposition. The actual transpiled layout from the hardware jobs was not inspected at the pulse level. Residual pulse-level effects cannot be excluded.

**Status: Effect too small by two orders of magnitude. Pulse-level effects not independently verified.**

### 5.5 Qubit-Qubit Crosstalk

**Hypothesis:** Ry(δ) on q1 induces an angle-dependent rotation on q0 through parasitic coupling, modifying Alice's effective measurement basis as a function of Bob's angle.

**Method:** Same nominal 2-qubit circuit topology (`initial_layout` not pinned), no entanglement preparation, Ry(δ) applied to q1, q0 measured only. 19 angles from 0° to 180° in 10° steps, 4096 shots each, on ibm_marrakesh. Fit model: P(|1⟩, q0) = a + b·sin²(δ/2).

**Result:** b = 0.0000 ± 0.0002 (0.3σ). Only 5 stray counts in 77,824 total shots, randomly distributed with no angular structure.

**Status: No effect detected.**

### 5.6 Summary

**Table 2: Systematic error candidates**

| Candidate | Method | Result | Status |
|---|---|---|---|
| Readout asymmetry | 4×4 correction matrix | Δα < 0.001 | No significant effect detected |
| Ry gate offset | Single-qubit calibration sweep | ε = 1.12°, ~20% of deviation | Insufficient |
| T2 decoherence | Aer noise model + T2 sweep | No systematic α shift; visibility-only effect | Not supported |
| Angle-dependent gate duration | Transpilation analysis (FakeMarrakesh) | 0.03% effect, two orders too small | No significant effect detected* |
| Qubit-qubit crosstalk | Crosstalk isolation experiment | b = 0.0000 ± 0.0002 | No significant effect detected |

\* Pulse-level gate schedule not independently verified.

No tested mechanism individually reproduces the observed α = 0.470 ± 0.003 deviation.

---

## 6. Discussion

### 6.1 Reproducibility

The deviation is consistent across ibm_fez (α = 0.467 ± 0.003) and ibm_marrakesh (α = 0.470 ± 0.003) — two physically separate chips with different qubit quality (V = 0.966 vs 0.974), different fabrication, and independent calibration. NISQ hardware artifacts are typically chip-specific; a structured artifact producing identical α to within 0.003 on two independent chips is not naturally expected from standard error models.

The marrakesh result is the stronger one: χ²/dof = 1.41 for the deformation model means it is a good fit with the right number of degrees of freedom, while QM+V gives χ²/dof = 4.83 on the same data.

### 6.2 Possible Interpretations

We do not claim the deviation is physical. The following interpretations remain open:

**Residual hardware artifact:** Some combination of effects not individually tested — for example, pulse-level gate miscalibration not captured by our Ry offset check, or two-qubit gate leakage not captured by our crosstalk check — could in principle produce the observed deviation. We have tested all individually accessible explanations and found them insufficient; combination effects are harder to rule out without pulse-level access.

**Geometric hidden variable model:** One possible physical interpretation is a geometric model in which the spin carries two pole distributions with unequal coupling weights α and (1−α). This model predicts exactly the crossover-shift parameterization used here, with α = 0.5 recovering quantum mechanics. In this interpretation, α = 0.470 implies a pole weight ratio of approximately 0.470:0.530. We note this model only as a prior motivation for the parameterization — consistency with a model is not evidence for it.

### 6.3 Limitations

1. **No pulse-level hardware access.** Gate implementations cannot be independently inspected. Residual calibration effects at the pulse level are not ruled out.

2. **Single qubit pair per chip.** All runs used q0/q1. A qubit-index swap (Ry on q0 instead of q1) was not run due to QPU budget exhaustion and would test whether the deviation is qubit-specific.

3. **Single technology platform.** Both chips are IBM superconducting transmon processors. Replication on a physically distinct platform (trapped ion, photonic) is the critical next experiment.

4. **No loophole-free Bell test.** This is an anomaly study, not a locality test. Measurement settings are not space-like separated.

5. **Control experiment qubit assignment.** The Ry calibration and crosstalk isolation experiments used the same nominal circuit topology but did not pin `initial_layout`. Physical qubit assignments may differ from the Bell sweep job.

6. **Layout consistency across angles.** Readout calibration was extracted from the first transpiled circuit. All 37 circuits were submitted as a single job (making layout drift unlikely), but identical physical qubit assignments were assumed, not verified.

### 6.4 Proposed Follow-up Experiments

In order of priority:

1. **Qubit-index swap:** Full Bell sweep with Ry on q0 instead of q1. Tests qubit-specificity of the deviation.
2. **Unentangled baseline:** P_disagree sweep on |01⟩ (no CNOT). Any angle-dependent structure is a circuit-level artifact.
3. **Independent platform:** Replication on IonQ or Quantinuum trapped-ion hardware. Identical α would be a strong cross-platform signal.
4. **Weihs 1998 photon data:** Raw click-by-click data at many angles could be fit with the deformation model as a completely independent physical system test.

---

## 7. Conclusion

We observe a statistically significant anomaly in the shape of the singlet Bell correlation curve on two IBM superconducting quantum processors. The quantum model with free visibility gives χ²/dof = 4.83 on ibm_marrakesh; a two-parameter deformation model that allows a crossover shift gives χ²/dof = 1.41. The crossover-shift parameter α = 0.470 ± 0.003 is consistent across both chips. The likelihood ratio test strongly prefers the deformation model over QM+V (Δχ² = 124.5 for 1 dof, p < 10^-28).

Five candidate instrumental explanations — readout asymmetry, Ry gate calibration offset, T2 thermal decoherence, angle-dependent gate duration variation, and qubit-qubit crosstalk — have been tested individually through dedicated experiments and analytical arguments and found insufficient to explain the anomaly. Thermal relaxation noise in particular does not produce a systematic crossover shift under the tested channel model. No tested noise mechanism accounts for the anomaly.

We present this as an unresolved experimental anomaly requiring independent verification. Replication on a physically distinct qubit platform is the decisive next step.

---

## Appendix A: Circuit Details

```
q0: ─────── H ── ●  ── Z ── M
q1: ─── X ───── ⊕ ──────── M
```

followed by Ry(δ) on q1 before measurement. The Z gate on q0 produces |Ψ⁻⟩; earlier runs with X on q0 produced |Ψ⁺⟩ and an inverted curve, detected by the shape of first-run results.

---

## Appendix B: Ry Gate Decomposition on IBM Marrakesh

Basis gates: {sx, rz, cz, x}. RZ duration: 0 ns (virtual). SX duration: 36 ns. CZ duration: 68 ns.

**Table B1: Transpiled Bell circuit gate counts by angle (FakeMarrakesh approximation)**

| Angle (°) | SX count | Duration (ns, excl. CZ) |
|---|---|---|
| 0 | 3 | 108 |
| 30 | 4 | 144 |
| 60 | 4 | 144 |
| 90 | 2 + 1X | 108 |
| 120 | 4 | 144 |
| 150 | 4 | 144 |
| 180 | 3 | 108 |

Maximum variation: 36 ns. Differential T2 factor at T2 = 110 µs: exp(−36×10⁻⁹ / 110×10⁻⁶) = 0.9997.

---

*Code, data, and analysis scripts: https://github.com/3axap4eHko/bell-curve-asymmetry*  
*IBM Quantum jobs: ibm_fez (d6o91869td6c73ap2mig, d6o99u8fh9oc73eq2pk0), ibm_marrakesh (d6o9cqobfi7c73a5kq0g, d6o9h68bfi7c73a5kvi0, d6o9pb8fh9oc73eq3dfg)*

---

## References

[1] J. S. Bell, "On the Einstein Podolsky Rosen paradox," *Physics Physique Fizika*, 1, 195 (1964).

[2] A. Aspect, P. Grangier, and G. Roger, "Experimental realization of Einstein-Podolsky-Rosen-Bohm Gedankenexperiment," *Phys. Rev. Lett.*, 49, 91 (1982).

[3] A. Aspect, J. Dalibard, and G. Roger, "Experimental test of Bell's inequalities using time-varying analyzers," *Phys. Rev. Lett.*, 49, 1804 (1982).

[4] B. Hensen et al., "Loophole-free Bell inequality violation using electron spins separated by 1.3 kilometres," *Nature*, 526, 682 (2015).

[5] G. Weihs, T. Jennewein, C. Simon, H. Weinfurter, and A. Zeilinger, "Violation of Bell's inequality under strict Einstein locality conditions," *Phys. Rev. Lett.*, 81, 5039 (1998).

[6] IBM Quantum, https://quantum.ibm.com/, accessed March 2026.

[7] Qiskit contributors, "Qiskit: An open-source framework for quantum computing," https://github.com/Qiskit/qiskit (2024).

[8] J. F. Clauser, M. A. Horne, A. Shimony, and R. A. Holt, "Proposed experiment to test local hidden-variable theories," *Phys. Rev. Lett.*, 23, 880 (1969).
