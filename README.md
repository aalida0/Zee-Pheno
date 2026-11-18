# Zee-Pheno

Minimal validation-first phenomenology for the **restricted CP-conserving Zee
model**: Type-II **quark texture**, two general complex lepton Yukawa matrices,
an antisymmetric charged-singlet coupling, and the charged-scalar cascade
H₂± → H₁± Z. The conventional Z₂ is **not** exact; see [Model](Model/README.md).

| Module | Scope |
| --- | --- |
| [Model](Model/README.md) | Scalar and SARAH/SPheno structural validation; [evidence](Model/VALIDATION.md) |
| [Neutrino](Neutrino/README.md) | Step 1: NuFIT 6.1 NO/IO pilot, native oscillation revalidation outstanding |
| [LFV](LFV/README.md) | Not yet promoted |
| [Higgs](Higgs/README.md) | Not yet promoted |
| [Flavor](Flavor/README.md) | Not yet promoted |
| [Collider](Collider/README.md) | Not yet promoted |
| [Scan](Scan/README.md) | Not yet promoted |

Run regression tests from the repository root:

~~~bash
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit -v
~~~

NumPy is needed for neutrino calculations; mpmath for the independent
high-precision native-loop regression. The historical SARAH/SPheno build and
large simulation products are **not** copied into this repository. The saved
native reference is one fixed card, not evidence of global parameter viability.
