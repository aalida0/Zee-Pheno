# Step 0 validation ledger — 2026-11-18

## Scope and release

**MODEL_STRUCTURAL_VALIDATION_PASS (restricted model).** This covers scalar
tree reconstruction, field/vertex conventions, native input propagation,
charged-lepton basis, and the one-card one-loop Zee mass **calculator**.
It does **not** establish physical viability of the benchmark, neutrino
oscillation acceptance, exact Z₂ protection, RG closure, global vacuum
minimality, LFV, flavor, Higgs bounds or collider sensitivity.

The source-checked root-cause investigation and its precise Fortran line
locators are retained separately:
[NEUTRINO_ROOT_CAUSE.md](NEUTRINO_ROOT_CAUSE.md).

## Evidence matrix

| Check | Evidence | Status |
| --- | --- | --- |
| Scalar masses, tadpoles, charged Goldstones | Model.test_model; symbolic 10×10 real-field audit, both gamma signs | PASS (tree) |
| Random scalar closure stress | 2,000 algebraic points, max squared-mass relative residual 1.24e−13 | PASS (algebraic) |
| SARAH 4.15.4 CheckModel | Completed; zero tree neutrino mass expected; missing allowed terms documented | PASS (structural with scope) |
| Native SPheno 4.0.5 generation + execution | Forced build attestation and one-card input propagation | PASS (one card) |
| Complex F/Ye2 native inputs | F relative error 7.76e−10; Ye2 4.69e−10; F antisym 0 | PASS |
| Tree charged-lepton Yukawa basis | Relative offdiag 2.285e−10; imaginary diagonal 4.019e−11 | PASS |
| Physical charged scalar Zee vertex | Explicit full F, charged mixing sign | PASS (structural) |
| Native 1-loop Majorana matrix | Analytic vs generated relative residual 1.721e−14 at matched running inputs | PASS (one card) |
| Two charged loop orientations | Independent vertex-sum residual 2.22e−15 | PASS (one card) |
| Native Takagi matrix | Residual 2.45e−15; arbitrary tree UV independence | PASS (one card) |
| F=0 / F doubled | Zero masses / doubled masses and matrix | PASS (one card) |
| Original pre-Step-1 neutrino oscillation fit | Historic baseline does not survive loop-enabled matching | **FAIL (historical card)** |
| LFV, full flavor, Higgs, vacuum-globality, RG, scan | Outside structural model sign-off | UNKNOWN |

Original reference SHA-256:
- Source input: f2efe759aa3196978d466739fc788b1dd7fb01d71cbc6e883c8a87830da454bf
- Loop-enabled input: e83f315928a3906bc14e48aec39b2bd842abf66949f53b6cd78824d4c45a8e07
- Production executable: 2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52
- Original loop spectrum: 2d72ae9f6d7e87b90a7fbb9db7825a6e2cbe919a8331d3f0b064e6a722e10ff9

The old mismatch arose because loop-enabled automatic matching changes Ye
and running masses, especially cancellation-sensitive D33. The generated
calculation at matched inputs gives masses (0.012785582, 0.041449883,
0.062565579) eV, faithfully reproduced by the analytic amplitude.
These **do not** reproduce NuFIT's solar mass-squared splitting. This is a
Step 1 fit problem, not a failure of the one-loop mass formula.

## Remaining model-scope limitations

The conventional Type-II quark Z₂ is hard-broken by generic lepton Ye2
and the charged-singlet k12 quartic. λ6/λ7 and wrong-Higgs quark operators
are absent only by the declared low-scale texture, **not** by proven
radiative protection. The scalar potential assumes a restricted slice
and does not yet have a global vacuum proof. The quark–singlet possible-term
warning from SARAH has not been promoted to a genuine interaction without
independent Lorentz/color/Fermi-statistics confirmation.

**Step 1 update (2026-11-18):** Two separately generated NO/IO native SPheno
benchmarks now reproduce their NuFIT 6.1 mass splittings and PMNS targets;
each passes all six one-dimensional 3σ screens and matrix-level self-energy
checks. These replace the failed *historical* card for the purpose of the
initial neutrino pilot. See [Neutrino/README.md](../Neutrino/README.md) for
observables, original proposal hashes and remaining scope limitations.
Neither matching-scale stability nor a joint NuFIT likelihood is certified;
LFV/Higgs/B/vacuum/scan viability remains outside the Step 0 sign-off.
