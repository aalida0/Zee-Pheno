# Model

**Promoted scope (tree level only):** `model.py` contains the independently cross-checkable scalar inversion and the antisymmetric complex Zee Yukawa constructor `zee_f(f12,f13,f23)`.

Source provenance: adapted from the audited scalar formulas in [the earlier Zee repository's Gate 6](https://github.com/aalida0/zee/blob/research/zee-full-constraints/scripts/gate6_benchmark.py), with the tadpole parameters from `module1/scalar_constraints.py`. The present port's independent two-by-two neutral/charged spectrum closure, pseudoscalar mass, signed mixing, and antisymmetry checks pass in `test_model.py`.

**Model convention:** H1 couples to down-type quarks, H2 to up-type quarks; both doublets have general lepton Yukawas, with an extra charged SU(2)-singlet antisymmetric coupling `f=-f.T`. Scalar potential is real and CP conserving; `lambda6=lambda7=0`. Thus 'Type-II' describes the quark sector, **not** a restrictive Type-II lepton texture.

**Not included/validated here:** a canonical SARAH model copy, SPheno-generated code, charged-lepton basis matching, loop-induced neutrino masses, one-/two-loop LFV, RG running, physical widths, vacuum-globality, or full scans. The older main-repo SARAH files and the model used for local SPheno generation differ in their neutrino-rotation definition. Do not silently promote either snapshot until the difference is reconciled and the corrected SLHA inputs are validated natively.

Validation evidence and outstanding native checks: [VALIDATION.md](VALIDATION.md).

Run: `python3 -m unittest Model.test_model -v`.

## SARAH/SPheno validation (pending on the local installation)

`check_sarah.wls` runs SARAH 4.15.4 `CheckModel[]` on an explicitly supplied historical model snapshot; its warnings must be reviewed. `check_native.py` checks real/imaginary `FZEEIN` and `YEE2IN` propagation in an already compiled, hash-attested SPheno binary, using one baseline card and two tiny probes. No generation, recompilation or full scan is performed.

```bash
python3 -m unittest Model.test_model Model.test_native -v
```

For the exact original paths and native commands, see [VALIDATION.md](VALIDATION.md). **A unit-test pass does not certify SPheno physics**; a native model pass does not certify LFV.

## Lepton basis and published one-loop formula

`check_basis.py` inspects the existing corrected `baseline.spc` without rebuilding SPheno. `test_conventions.py` cross-checks the leading one-loop neutrino mass convention against 2017 and 2023 Zee formulas, including the SARAH `Fzee/2` normalization.

```bash
python3 -m unittest Model.test_conventions Model.test_basis -v
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.check_basis --spectrum "$WORK/model-validation-native-20261118/baseline.spc"
```

The tree-level neutrino `UV` is not the loop-corrected PMNS matrix. Symmetry protection for every gauge-allowed omitted coupling and the native one-loop mass are still **UNKNOWN**; see the warning register in [VALIDATION.md](VALIDATION.md).

## Controlled one-loop neutrino test

The 2026-11-18 native baseline had two `SPhenoInput[55]=0` controls, so its zero neutrino masses are **not a one-loop test**. Run `python3 -m unittest Model.test_loop_probe -v`, then follow the attested one-card `Model.probe_neutrino_loop` command in [VALIDATION.md](VALIDATION.md). It enables both loop-mass switches in a **new** card without modifying the baseline. A nonzero output neutrino pole mass is still **not** a validated Majorana mass matrix.

## Native Zee mass comparison

The loop-enabled one-card run now returns nonzero PDG 12/14/16 masses. They have **not** been matched against the published Zee neutrino mass prediction. With NumPy available, run:

```bash
python3 -m unittest Model.test_compare_neutrino_masses -v
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.compare_neutrino_masses \
  --input "$WORK/model-validation-native-20261118/baseline.in" \
  --baseline "$WORK/model-validation-native-20261118/baseline.spc" \
  --loop "$WORK/model-neutrino-loop-20261118/loop_enabled.spc"
```

The output compares invariant **neutrino singular values only**, not the complex Majorana matrix or PMNS angles. See [VALIDATION.md](VALIDATION.md) for provenance and the mismatch gate.

### November 18 numerical mismatch — follow-up controls

The one-card native result gives `[0.01278558, 0.04144988, 0.06256558]` eV whereas the published-formula reconstruction gives `[0.00930183, 0.01272101, 0.05095806]` eV: **do not use native neutrino pole masses as a physical constraint yet**. `compare_neutrino_masses.py` now also prints tree-to-loop charged-scalar/lepton mass shifts. `probe_neutrino_scaling.py` checks that the same attested native binary produces zero neutrino masses at `Fzee=0` and responds approximately linearly at `Fzee*2`.

```bash
python3 -m unittest Model.test_compare_neutrino_masses Model.test_neutrino_scaling -v
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.probe_neutrino_scaling \
  --input "$WORK/model-neutrino-loop-20261118/loop_enabled.in" \
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
  --baseline "$WORK/model-neutrino-loop-20261118/loop_enabled.spc" \
  --output-dir "$WORK/model-neutrino-F-controls-20261118"
```

