# Model

**Promoted scope (tree level only):** `model.py` contains the independently cross-checkable scalar inversion and the antisymmetric complex Zee Yukawa constructor `zee_f(f12,f13,f23)`.

Source provenance: adapted from the audited scalar formulas in [the earlier Zee repository's Gate 6](https://github.com/aalida0/zee/blob/research/zee-full-constraints/scripts/gate6_benchmark.py), with the tadpole parameters from `module1/scalar_constraints.py`. The present port's independent two-by-two neutral/charged spectrum closure, pseudoscalar mass, signed mixing, and antisymmetry checks pass in `test_model.py`.

**Model convention:** H1 couples to down-type quarks, H2 to up-type quarks; both doublets have general lepton Yukawas, with an extra charged SU(2)-singlet antisymmetric coupling `f=-f.T`. Scalar potential is real and CP conserving; `lambda6=lambda7=0`. Thus 'Type-II' describes the quark sector, **not** a restrictive Type-II lepton texture.

**Not included/validated here:** a canonical SARAH model copy, SPheno-generated code, charged-lepton basis matching, loop-induced neutrino masses, one-/two-loop LFV, RG running, physical widths, vacuum-globality, or full scans. The older main-repo SARAH files and the model used for local SPheno generation differ in their neutrino-rotation definition. Do not silently promote either snapshot until the difference is reconciled and the corrected SLHA inputs are validated natively.

Validation evidence and outstanding native checks: [VALIDATION.md](VALIDATION.md).

Run: `python3 -m unittest Model.test_model -v`.
