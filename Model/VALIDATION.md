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

## 2026-11-18 SARAH/SPheno validation preparation

### Inspectable evidence

- The older repository's `models/SARAH/ZeeTypeII/parameters.m` declares `Fzee` complex and `Form -> AntiSymmetric`. Its `ZeeTypeII.m` uses doublet Yukawas `Ye, Ye2`, quark `Yd, Yu`, and `Fzee/2 l.l.Splus`. This is **static source evidence only**.
- The November 14 native source audit attributed the anomalous output Yukawas to three-field rows of `FZEEIN/IMFZEEIN/YEE2IN/IMYEE2IN` being read with the generated Fortran's four-field `Read(io,*) i1,i2,wert,read_line`. A separate local gfortran reader check confirms `#` is accepted as the fourth character token.
- `Model/check_native.py` generates a corrected baseline and independent antisymmetric-real-`F` and imaginary-`Ye2` probes, keeps the original benchmark untouched, verifies the executable against the previously recorded forced-build SHA, and writes the diagnostics in a new directory. `Model/test_native.py` tests serialization, antisymmetry and a **mock** executable. Mock success is not physical validation.
- `Model/check_sarah.wls` runs the local SARAH `Start["ZeeTypeII"]` and `CheckModel[]`. Its output **requires human review**; no automatic PASS is inferred from the exit status.

### Execute the native gate on the user's machine

From `~/Zee-Pheno` after pulling the current `main`:

```bash
python3 -m unittest Model.test_model Model.test_native -v

# The two historical model snapshots are not identical; compare first.
diff -u \
  ~/zee/models/SARAH/ZeeTypeII/ZeeTypeII.m \
  ~/zee-conditional-cascade-scan/validation/module14/local/models/ZeeTypeII/ZeeTypeII.m || true

# Supply the real SARAH installation's SARAH.m path.
SARAH_M="$(find "$HOME" -name SARAH.m -type f -print -quit)"
wolframscript -file Model/check_sarah.wls "$SARAH_M" \
  "$HOME/zee-conditional-cascade-scan/validation/module14/local/models" \
  2>&1 | tee /tmp/zee_sarah_checkmodel.log

WORK="$HOME/zee-lfv64-serial-20261112-094549"
CARD="$WORK/EDM3ZERO_SOURCE64_20261114_172527/native-source-64-results/Zee_2264_NO_20261111_newEDM3zero/input.in"
python3 -m Model.check_native \
  --input "$CARD" \
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
  --attestation "$WORK/model-rebuild-attestation.json" \
  --output-dir "$WORK/model-validation-native-20261118"
```

The output directory must not exist beforehand. Inspect `native-input-report.json` and the three `.spc`/`.log` files. If `SARAH.m` cannot be found, provide its explicit path rather than using a substitute version.

### Interpretation and release gate

- `NATIVE_INPUT_GATE_PASS` means the **specific one-card input/output propagation checks** passed; other loop-related correctness remains UNKNOWN.
- `NATIVE_INPUT_GATE_FAIL` or an exception requires investigating SPheno reading, RG/output scale, source linkage, or missing matrices before proceeding. The diagnostic tolerances (`5e-4` relative input matching and `5e-9` absolute antisymmetry) may be too restrictive if output and input are evaluated at different scales. Do **not** loosen them simply to obtain a pass.
- `F_probe_unintended_Ye2_delta` and `Ye2_probe_induced_F_delta` are measurements, not zero-coupling requirements: RGE mixing can be physical.
- Keep **native model validation UNKNOWN** until the SARAH model-snapshot/rotation discrepancy, `CheckModel[]` warnings, source-to-binary provenance and one-card native results are reviewed. Then check the charged vertices, charged-lepton basis `D` and radiative neutrino-mass normalization independently.

**Current release decision: TREE_SCALAR_VALIDATED / NATIVE_SARAH_SPHENO_PENDING / PHYSICAL_LFV_UNKNOWN.**
