# Zee-Pheno

Minimal, validation-first code for the **restricted CP-conserving Zee model with Type-II quark Yukawas and general lepton Yukawas**, focusing on the charged cascade `H2+ -> H1+ Z`.

Only promote a calculation when its inputs, conventions, and scoped tests are documented. A passing unit test is **not** a global phenomenological exclusion or complete LFV validation.

| Directory | Status |
| --- | --- |
| [Model](Model/README.md) | Tree-level scalar inversion and complex antisymmetric Zee coupling; [validation record](Model/VALIDATION.md) |
| [Neutrino](Neutrino/README.md) | Reserved; no promoted code |
| [LFV](LFV/README.md) | Reserved; native SPheno/FZEE inputs still need independent rerun |
| [Higgs](Higgs/README.md) | Reserved; no promoted code |
| [Flavor](Flavor/README.md) | Reserved; no promoted code |
| [Collider](Collider/README.md) | Reserved; no promoted code |
| [Scan](Scan/README.md) | Reserved; no promoted code |

Run the model checks (Python standard library only):

```bash
python3 -m unittest Model.test_model -v
```

No generated SPheno/SARAH, old scan artifacts, or unvalidated global acceptance claims are included.
