# Zee-Pheno

Minimal validation-first phenomenology for the **restricted CP-conserving Zee
model**: Type-II **quark texture**, two general complex lepton Yukawa matrices,
an antisymmetric charged-singlet coupling, and the charged-scalar cascade
H₂± → H₁± Z. The conventional Z₂ is **not** exact; see [Model](Model/README.md).

| Module | Scope |
| --- | --- |
| [Model](Model/README.md) | Scalar and SARAH/SPheno structural validation; [evidence](Model/VALIDATION.md) |
| [Neutrino](Neutrino/README.md) | Step 1: native NO/IO 1D NuFIT 3σ matched controls; not global neutrino scan |
| [LFV](LFV/README.md) | Step 2: general per-point scalar one-loop, source-matched photon Barr–Zee and tree LFV subsets; full rates UNKNOWN |
| [Higgs](Higgs/README.md) | Not yet promoted |
| [Flavor](Flavor/README.md) | Not yet promoted |
| [Collider](Collider/README.md) | Not yet promoted |
| [Scan](Scan/README.md) | Scalar-first proposals, multi-family neutrino inverse and native-card LFV gates; full allowed-region scan NOT implemented |

Run regression tests from the repository root:

~~~bash
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit LFV.test_pilot LFV.test_constraints Scan.test_native_pipeline Scan.test_physical_proposals -v
~~~

NumPy and mpmath are needed for neutrino/native-loop calculations; SciPy for
source-matched LFV Barr-Zee quadrature. The historical SARAH/SPheno build and
large simulation products are **not** copied into this repository. The saved
native reference is one fixed card, not evidence of global parameter viability.
