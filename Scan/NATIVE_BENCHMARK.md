# Native Higgs calibration and neutrino refit

Continuation of `08da284`, on `investigation/spheno-scalar-inputs`.

**Two fitted s00017 benchmarks, NO and IO, pass the native Higgs-mass,
oscillation and cascade necessary screens.** They preserve the pilot's three
established scalar preselection checks. They are **not globally allowed
points**: full LFV, flavor, Higgs rates/searches, unitarity, vacuum and EWPO
remain UNKNOWN. This establishes reproducible native benchmark controls
before expanding the scan.

The [runner](calibrate_native_benchmark.py) uses the unchanged analytical scalar
inversion, complex neutrino inverse, matched one-loop audit, input editor and
SPheno installation. No generated physics source or original benchmark changes.
The [ledger](reference/native_benchmark/ledger.json) retains all 30 native
runs and their failures, inputs, parameters, observables and artifact hashes.
Final cards, original-executable spectra and matched snapshots are committed:

* [NO card](reference/native_benchmark/s00017_NO_benchmark/input.in),
  [spectrum](reference/native_benchmark/s00017_NO_benchmark/spectrum.slha),
  [native snapshot](reference/native_benchmark/s00017_NO_benchmark/native.json).
* [IO card](reference/native_benchmark/s00017_IO_benchmark/input.in),
  [spectrum](reference/native_benchmark/s00017_IO_benchmark/spectrum.slha),
  [native snapshot](reference/native_benchmark/s00017_IO_benchmark/native.json).

## Quantitative diagnosis of the four configurations

All values below come from rerunning the mapping-corrected cards, not the
old ignored-block inputs. The first two families per scalar are tabulated;
all four families and all 16 failures are in the ledger.

| Scalar | First native tree mh | Final matching tree mh | Native full one-loop mh, flag 67=0 | Default native mh, flag 67=1 |
|---|---:|---:|---:|---:|
| s00001 NO | 125.001427 | 125.698126 | 146.702936 | 146.700634 |
| s00010 NO | 125.004527 | 125.560331 | 139.561802 | 139.561224 |
| s00017 NO | 125.004991 | 125.709877 | 153.759019 | 153.755294 |
| s00020 NO | 125.010479 | 125.688894 | 136.194759 | 136.195028 |

Masses are GeV. The first tree shift from the requested 125.25 GeV follows
from the reference-to-native vev change. Matching and tadpoles shift the
recorded tree masses by less than 1 GeV. The dominant remaining difference
comes from the native loop pole-mass calculation: approximately 21.00, 14.00,
28.05 and 10.51 GeV above the final matching tree values. The tree invocation
is identified by its `Match_and_Run|OneLoopMasses|` call stack, excluding
later decay/coupling calls. Removing the final
effective-SM replacement changes these numbers by only 2.302, 0.578, 3.725
and −0.269 MeV. Thus changing flag 67 does not solve the Higgs failure.
These four flag-67 runs are diagnostics only; fitted cards retain flag 67=1.

Relevant unchanged Fortran: `Boundaries_ZeeTypeII.f90:1757–1763` sets the first
zero-loop tadpoles and tree spectrum; `:1514` evaluates masses during matching;
`:1608–1621` computes final loop masses and the configured effective-SM Higgs
replacement. The previous report gives the full reader/boundary conventions.

| Scalar | Order | sin²θ12 | sin²θ13 | sin²θ23 | Δm²21, eV² | Δm²3l, eV² | δ, degrees |
|---|---|---:|---:|---:|---:|---:|---:|
| s00001 | NO | .005775 | .030155 | .473881 | .000339697 | .00464498 | 341.706 |
| s00001 | IO | .512386 | .354084 | .999019 | .0740514 | −.572257 | 73.764 |
| s00010 | NO | .317945 | .239313 | .341596 | 9.07449 | 9.33506 | 359.243 |
| s00010 | IO | .526652 | .546713 | .998898 | .0630619 | −.511557 | 340.543 |
| s00017 | NO | .020108 | .005940 | .486972 | .0582762 | .0777486 | 295.243 |
| s00017 | IO | .343183 | .044309 | .527092 | .000755290 | −.00299292 | 295.073 |
| s00020 | NO | .302390 | .055310 | .403174 | .00118086 | .00177464 | 44.599 |
| s00020 | IO | .389901 | .027385 | .562791 | .000153996 | −.00264814 | 287.192 |

For IO, Δm²3l is Δm²32. The existing NuFIT 6.1 IC24+SK central targets and
six projected 3σ ranges are reused unchanged from `Neutrino/fit.py`; this is
not a new data update or a joint likelihood.

Define the SM-direction Yukawa `SY=cosβ Ye+sinβ Ye2` and the orthogonal one
`Dg=−sinβ Ye+cosβ Ye2`. The generator used the archived `SYref`, so its card
contained `Ye2=sinβ SYref+cosβ Dg_proposal`. Native matching instead solves
Ye with its own `SYnative`. Since this LOW model directly sets Ye2 from its
input, at the same β the exact gauge-basis error is

```text
Dg_native − Dg_proposal = tanβ (SYref − SYnative).
```

Transporting both sides with the actual ZEL/ZER rotations explains the
observed D shifts to below 2e-17 absolute Yukawa units. This is the previously
validated running-input effect, quantified at the new scalar points.

| Scalar | Native kernel-only target error, NO | Full target error, NO | Full target error, IO | Matching-induced real ΔDττ, NO |
|---|---:|---:|---:|---:|
| s00001 | .003635 | .367041 | 12.4584 | 1.21555e-4 |
| s00010 | .002863 | 83.4435 | 11.5107 | 6.45413e-5 |
| s00017 | .003827 | 4.52006 | .120272 | 1.74706e-4 |
| s00020 | .003757 | 1.50001 | .038624 | 5.56973e-5 |

Errors are `||M−Mtarget||F/||Mtarget||F`. “Kernel-only” freezes the proposed
mass-basis F,D and substitutes native charged masses, charged mixing and
running lepton masses. “Full” also uses native running F,D and rotations.
Kernel changes are only about 0.3–0.4%; the mismatch of SY dominates the
failure. For example s00017 NO Dττ changes from approximately
`2.35355e-5+7.26739e-6 i` to `1.98241e-4+7.26736e-6 i`.
The matrix error can greatly exceed one because some proposed D entries are
small cancellations. Matched analytical/native calculator agreement remains
at approximately 1e-13 or better for the original 16 cases.

## Minimal scalar calibration

Selection was fixed by the smallest original maximum absolute quartic among
the four configurations, then scalar ID. This selects s00017 with maximum
quartic 1.3136072905. There was no selection using LFV or a wider scan.

Only the **tree light CP-even mass coordinate** is calibrated. β, α, γ, v,
the heavy CP-even, CP-odd and both charged tree masses, soft M12Sq, LambdaC
and K1/K2/K12 stay fixed. Reapplying `Model/model.py` changes Lambda1/2/3
and dependent M112/M222; it leaves Lambda4/5, Mu4 and MCS2 fixed. The
charged tree mass matrix therefore stays fixed at the proposal vev.

A prespecified `[70,125.25]` GeV bracket and Brent root solve in mh_tree²
minimize the physical residual `mh_native²−125.25²`. Every new potential
must preserve all three preselection PASS results: necessary doublet BFB,
the conservative sufficient full-Zee **quartic** BFB condition, and
`max|quartic|<4π`. Every executed candidate must keep a positive native
spectrum and the on-shell cascade. A failed constraint or absent root bracket
ends calibration with an explicit blocker, not a changed physics gate.

| Tree mh coordinate, GeV | Native mh, GeV | Outcome |
|---:|---:|---|
| 125.250000000 | 153.755294 | Baseline reused as high endpoint |
| 70.000000000 | 107.630933 | New low endpoint; theory and cascade PASS |
| 92.581645580 | 125.215854 | Intermediate root evaluation |
| 92.623006114 | 125.250000 | Calibrated |

The root-coordinate tolerance is 0.1 GeV²; required native numerical accuracy
is 0.02 GeV. These solver criteria do not replace the existing ±3 GeV Higgs
necessary screen, and do not represent a Higgs-mass theory uncertainty.
The final independent quartics are:

```text
Lambda1 = 0.9913961366090698     Lambda2 = 0.07131547405679571
Lambda3 = 0.5427544843746859     Lambda4 = 1.313607290470687
Lambda5 = −0.9242430772423869    LambdaC = 0.18029868498871982
K1 = 0.0874750535214805         K2 = 0.5648002187234421
K12 = −0.028930108119468036
M12Sq = 98976.16533692155 GeV²; MCS2 = 1338909.6115773146 GeV²
Mu4 = 814.638153103011 GeV.
```

These are deliberate changes to the model inputs to fit a physical pole
mass, not suppression of a warning. Original cards remain unchanged. The
fitted tree mh coordinate departs from the old generator's fixed 125.25 GeV
coordinate; these fitted benchmarks are not unbiased draws from that prior.

## Native-consistent complex neutrino refit

For each of the selected scalar's original NO f00 and IO f01 families:

1. Run the calibrated card and identify the OneLoopFv invocation that uniquely
   matches the written neutrino MASS entries. Later Higgs-coupling invocations
   are not used for fitting.
2. Read its MFe, charged tree MHm, charged ZP mixing, vd/vu, ZEL/ZER, Ye, Ye2
   and F. Verify the LOW-boundary F/Ye2 input-to-running identity before using
   them as input coordinates. All matching flags and SM inputs remain fixed.
3. Use the existing finite-lepton-mass analytical inverse with these running
   inputs, retaining the family's original lightest mass, Majorana phases,
   F scale/phase, null-cone chart, t, root and homogeneous coefficients.
4. Rotate back with `Fg=ZEL.T Fmass ZEL`, `Dg=ZER.T Dmass ZEL`; set
   `Ye2=sinβ SYnative+cosβ Dg`. Native matching continues to solve Ye.
5. Rerun SPheno. Repeat with newly measured running inputs until the full
   target matrix closes to relative 1e-8 and all six oscillation screens,
   Higgs accuracy, scalar propagation, cascade, calculator and charged-lepton
   identities pass. Require running and proposed Yukawas below 1. Maximum
   six updates; unavailable inverses, failed matching, constraints or exhaustion
   retain explicit failures. No optimizer changes the target or screen ranges.

Both orderings require two updates. First updates already pass all six
oscillation intervals, but the stronger full-matrix criterion intentionally
requires a second update:

| Order | After scalar calibration | After refit 1 | After refit 2 |
|---|---:|---:|---:|
| NO target matrix relative error | 5.09465 | 1.32482e-6 | 5.75237e-14 |
| IO target matrix relative error | .135558 | 3.38237e-7 | 2.09767e-12 |

Maximum final running Yukawa magnitudes are 0.00977346 (NO) and 0.01197638
(IO). Native charged-lepton mass-direction matrix residuals are below
2.1e-15. The matched analytical/native neutrino calculator residuals are
2.02e-14 (NO) and 4.98e-15 (IO). This is numerical closure at the retained
perturbative order; no higher-order accuracy is inferred.

## Final native benchmark values and gates

The unchanged original executable, SHA256
`2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52`,
independently reproduces each entire final SLHA spectrum **byte for byte**.
This checks masses and mixing outputs as well as the instrumented snapshot
calculation. The archived native files retain all spectrum-associated complex
matrices and generated vertex arrays for independent calculator regression.

| Observable | NO | IO |
|---|---:|---:|
| mh pole, GeV | 125.250000 | 125.250135 |
| mH pole, GeV | 689.790923 | 689.791504 |
| mA pole, GeV | 735.458874 | 735.459604 |
| mH1± pole, GeV | 666.305092 | 666.305622 |
| mH2± pole, GeV | 1184.553650 | 1184.553680 |
| mH2±−mH1±−mZ, GeV | 427.060958 | 427.060458 |
| sin²θ12 | .308800000000009 | .308800000008933 |
| sin²θ13 | .022480000000001 | .022619999999728 |
| sin²θ23 | .470000000000014 | .549999999999405 |
| Δm²21, eV² | 7.53699999998e-5 | 7.53699999873e-5 |
| Δm²3l, eV² | .002511000000000 | −.002482999999994 |
| δ, degrees | 211.999999999983 | 273.999999999258 |
| Neutrino masses (physical label order), eV | .02854870, .02983954, .05767173 | .05416538, .05485670, .02294031 |
| Sum neutrino masses, eV | .11605997 | .13196239 |

Pole e,mu,tau remain `(0.00051099893,0.105658372,1.77686)` GeV. GF and
mZ inputs are unchanged; mZ=91.1876 GeV. Native running v is approximately
247.745019 GeV. Charged tree masses used in the validated neutrino loop are
approximately `(671.646401,1174.085121)` GeV: they are distinct from the
charged pole masses used for the physical cascade threshold.

| Gate | NO | IO | Evidence / limitation |
|---|---|---|---|
| Native positive spectrum | PASS | PASS | Original executable and diagnostic copy agree |
| Native h125 necessary screen | PASS | PASS | Also within 0.02 GeV numerical target accuracy |
| Native six NuFIT projected 3σ screens | PASS | PASS | Central target reproduced, not a correlated likelihood |
| Matched analytical/native one-loop neutrino calculator | PASS | PASS | Actual running parameters and vertex arrays |
| Charged-lepton mass identity and EW inputs | PASS | PASS | Native mass direction and unchanged SM inputs |
| On-shell charged cascade | PASS | PASS | Native pole masses and nonzero mixing; no BR assertion |
| Three declared scalar preselection conditions | PASS | PASS | Rechecked after calibration; no global vacuum certificate |
| Yukawa magnitude below 1 | PASS | PASS | Full complex matrices, not only diagonals |
| Full LFV / flavor / Higgs rates and direct searches | UNKNOWN | UNKNOWN | These screens alone cannot establish viability |
| Full tree unitarity, global vacuum, EWPO, cosmology | UNKNOWN | UNKNOWN | Must evaluate before allowed-region claims |

## Downstream pipeline and approximations

The final s00017 NO SCALARMIX has determinant −1, an allowed orthogonal
matrix with an arbitrary real eigenstate sign. `LFV/pilot.py:neutral_alpha`
previously required only the determinant +1 convention. It now accepts either
sign of the second row while keeping the same normalization and orthogonality
checks and tolerance. The extracted canonical alpha and internal scalar
amplitudes are unchanged by that row rephasing. Nonorthogonal matrices are
still rejected; no spectrum is edited.

The existing source-bound partial pipeline now evaluates both final cards,
reports `UNKNOWN_PENDING_GATES`, and keeps `globally_allowed=null`. Its
**scalar one-loop subset** BR values are:

| Process | NO | IO |
|---|---:|---:|
| mu→e gamma | 3.6400e-19 | 2.8396e-13 |
| tau→e gamma | 4.1441e-19 | 2.2837e-15 |
| tau→mu gamma | 4.8094e-17 | 1.9982e-13 |

These are diagnostics, not LFV acceptance/exclusion: omitted two-loop,
other LFV and meson contributions remain unbounded. This subset mixes native
charged tree masses with neutral pole masses as already documented by that
pipeline. No arbitrary quark running masses were supplied to enable extra
Barr–Zee terms. The neutrino fit targets do not impose cosmological mass-sum
constraints; the displayed sums are not a cosmology PASS.

The scalar pole calibration uses the existing native matching (`19=-2`,
resolving to 2), one-loop masses (both flag-55 entries remain 1), and the
existing effective-SM Higgs treatment (`67=1`). A sizeable tree-to-pole Higgs
shift remains: the 92.623 GeV coordinate is a running/tree parameter, not a
prediction of a second light pole state. Higher-order scalar uncertainties
and charge-breaking/global-vacuum stability are UNKNOWN. No precision claim
or demand for a two-loop calculation is needed to describe the present scope.

## Complete attempt statistics

| Phase | Native runs | Outcome |
|---|---:|---|
| Mapping-corrected original pilot | 16 | 16 Higgs and 16 oscillation screen FAIL; spectra complete |
| Flag-67 diagnostic controls | 4 | 4 spectra complete; never eligible as fitted exports |
| Scalar root evaluations beyond reused baseline | 3 | All theory/cascade PASS; neutrino targets still FAIL |
| IO start at the NO-calibrated potential | 1 | Higgs PASS, oscillation FAIL |
| Native inverse updates | 4 | All six oscillation screens PASS after each; first update per ordering fails the stronger 1e-8 matrix criterion |
| Independent original-executable final controls | 2 | Full SLHA identical to diagnostic spectra |
| **Total** | **30** | **0 native failures, 0 rejected theory trials, 2 fitted necessary-screen benchmarks** |

No original failure was removed or renamed as a success. The 30 initial
scalar proposals / four independent scalar configurations / 16 conditional
families remain the original sampling denominators. Two deliberately fitted
benchmarks do not define a scan survival fraction. Initial exploratory
experiments in `/tmp/zee-native-fit-experiment` and the development reruns
are retained locally; the reproducible 30-call experiment is recorded in
`/tmp/zee-native-benchmark-release` and the committed ledger.

The [complete execution index](reference/native_benchmark/execution_statistics.json)
also records the development repeats: 104 native executions in total
(74 development/control repetitions plus the 30-call release run), zero
native spectrum failures. It records each card/log/spectrum hash and each
available Higgs/oscillation result; unevaluated diagnostics remain UNKNOWN.
The repeated original pilot failures are retained. Two downstream parser
rejections encountered during development were traced to the SCALARMIX
row-sign check above and resolved before the release run. They were not
native SPheno failures.

## Reproduction and regression

```bash
cd /home/amine/Zee-Pheno
WORK=/home/amine/zee-lfv64-serial-20261112-094549
python3 -m Scan.calibrate_native_benchmark \
  --spheno "$WORK/SPheno-4.0.5" \
  --pilot "$WORK/scalar_neutrino_pilot_01" \
  --output-dir /tmp/zee-native-benchmark-reproduce
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit LFV.test_pilot LFV.test_constraints \
  LFV.test_strict_2hdm_reference LFV.test_gauge_vertices \
  LFV.test_approximation_study LFV.test_nonzero_mixing_study \
  LFV.test_doublet_proxy_study Scan.test_native_pipeline \
  Scan.test_physical_proposals Scan.test_spectrum_preflight \
  Scan.test_diagnose_scalar_spectrum Scan.test_calibrate_native_benchmark -v
```

Use a fresh output directory; existing outputs are never overwritten.
Versions used: Python 3.10.12, numpy 2.2.6, scipy 1.15.3, GNU Fortran 11.4.0,
existing SARAH-4.15.4-generated ZeeTypeII/SPheno-4.0.5 sources and libraries.
Their hashes and the original reference-snapshot/card provenance are in the
ledger. The original binary, seven relevant original Fortran sources and
libraries are checked unchanged. Original reference-neutrino regressions
remain PASS. **88 regression tests PASS** (17 Model + 71 others), including
six new native-benchmark regressions and the reflected-row mixing regression.
CI runs archived numerical evidence and synthetic contracts; it does not
pretend to build the local SPheno installation or rerun external-author LFV CI.

**Release decision:** PASS for using these two source-bound controls to test
the native workflow; UNKNOWN for full phenomenological viability. The next
step is evaluating the outstanding constraints on these fixed cards before
using them to justify expansion into an allowed-region scan. No wider scan
was launched.
