# Model validation record — 2026-11-18

Scope: **restricted, CP-conserving Zee scalar potential at tree level**, Type-II quark Yukawas, general charged-lepton Yukawas and a complex antisymmetric singlet coupling `f=-f^T`. This is **not** a globally accepted Zee benchmark.

## Executed checks

- `python3 -m unittest Model.test_model -v`: **4/4 passed**, Python 3.13.5. Tests cover exact antisymmetry, original scalar spectrum, tadpoles, full charged Goldstone/invariants, and signed charged mixing.
- Independent potential audit (SymPy differentiation of the 10 real scalar degrees of freedom, NumPy diagonalization): benchmark `mh=125, mH=580.7, mA=633.7, mH1+=300, mH2+=700 GeV`, for `gamma=+1.089` and `gamma=-1.089`. Three massless Goldstones; largest |Goldstone squared mass| `6.8e-11 GeV^2`, maximum physical squared-mass relative residual `2.8e-15`, maximum tree-level tadpole `6.3e-10 GeV^3`.
- Independent deterministic stress test of 2,000 unconstrained parameter points (seed `20261118`), covering both charged-mixing signs: maximum relative mass-squared reconstruction residual `1.24e-13`; maximum normalized tadpole residual `6.23e-16`; maximum absolute charged-Goldstone mass-squared residual `1.93e-9 GeV^2`. **These are algebraic points, not phenomenologically accepted scan points.**

The 10x10 symbolic audit and stress test were executed independently of the GitHub production model tests; only the four small standard-library regression tests are retained in this minimal repository.

## Unvalidated / do not promote

- **SARAH/SPheno canonical model:** not established. The older repository's `models/SARAH/ZeeTypeII/ZeeTypeII.m` and the generated-model snapshot at `/home/amine/zee-conditional-cascade-scan/validation/module14/local/models/ZeeTypeII/` differ in their neutrino rotation definitions. Do not silently select either snapshot.
- **Native input propagation:** the older repo's SLHA serialization fix restores a fourth field for generated Fortran matrix readers, but a corrected one-card native execution has not yet been independently inspected here.
- **Radiative neutrino masses, LFV, EDM, Yukawa basis transformations, loop normalization and global vacuum:** not validated by these tree-level tests.
- **BFB/unitarity/EWPO/Higgs/flavor/collider constraints:** no acceptance or exclusion claim from this module.

## Required native-model acceptance gate

On the user's machine with SARAH 4.15.4 and existing compiled SPheno 4.0.5:

1. Reconcile the two historical SARAH model snapshots, preserving neutrino conventions; run `CheckModel[]` and independently inspect `Fzee` antisymmetry and charged scalar vertices.
2. On one corrected four-field SLHA benchmark, verify independent changes to `FZEEIN` and `YEE2IN` affect only their intended output parameters, reconstruct `f+f^T=0`, and match the scalar spectrum and physical Yukawa basis.
3. Validate the radiative neutrino mass and source-matched LFV amplitudes separately before promoting any generated model, code, or global exclusions.

Result: **TREE_SCALAR_VALIDATED; NATIVE_SARAH_SPHENO_PENDING**.
