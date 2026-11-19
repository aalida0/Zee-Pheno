# Native scalar-spectrum investigation — 19 November 2026

**Root cause: incorrect scalar input serialization, demonstrated at tree level.**
The old generator wrote tanbeta to `MINPAR[3]`, which this generated model
reads as **Lambda3Input**. It left `MINPAR[10]` (the actual tanbeta) and the
other scalar MINPAR values at the neutrino-reference values. Its newly appended
`ZEEPOT` and `ZEEQUARTIC` blocks were ignored. The intended analytical point
and the native point were different scalar potentials.

**Correction:** write all 13 scalar MINPAR entries using the existing card
editor. Remove synthesis of the ignored blocks. Keep `M112,M222` as tadpole
solutions. No change to the analytical inversion, matching flags, generated
model, executable, loop neutrino calculation, priors or Yukawa proposals.

**Decision:** PASS for repairing the interface and resuming exploratory
candidate generation with native rejection gates. FAIL for promoting any of
these 16 pilot points to a viable h125/neutrino benchmark. Full allowed-region
and publication claims remain UNKNOWN. No tachyon check was bypassed.

## Frozen provenance and isolation

Installation: `/home/amine/zee-lfv64-serial-20261112-094549/SPheno-4.0.5`.
All Fortran references below are relative to its **`ZeeTypeII/`** directory,
with line numbers in the original sources, before diagnostic writes.

* Original binary SHA256:
  `2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52`.
* Original `model-neutrino-loop-20261118/loop_enabled.in` SHA256:
  `e83f315928a3906bc14e48aec39b2bd842abf66949f53b6cd78824d4c45a8e07`.
* Existing pilot and preflight directories were read only. Existing untracked
  outputs in the repository were left untouched. The two identical flag-55
  entries remain 1; matching remains `19=-2`.
* Raw output-only diagnostic build/runs from the final reproducibility check:
  `/tmp/zee-scalar-integration-check`; initial ordered controls:
  `/tmp/zee-scalar-first-controls`. Corrected regenerated cards and original
  executable preflight: `/tmp/zee-scalar-corrected-pilot-v1` and
  `/tmp/zee-scalar-corrected-preflight-v1`.
* [Archived numerical evidence](reference/native_scalar_investigation.json)
  records all six relevant original source hashes, library/binary/card hashes,
  first tree matrices and call stacks, 16 corrected matched-neutrino checks,
  matching controls, and the original-executable preflight. Raw artifact hashes
  permit later comparison; `/tmp` files are not durable publication storage.

The diagnostic links copies of the main program, tree-mass and loop-mass
sources against the existing libraries. Changes in those copies only write
matrices, parameters, call stacks and the previously validated neutrino
snapshots. The original production binary independently reproduces the control
and all 16 regenerated spectra. No SARAH regeneration or library rebuild.

## Exact reader and boundary chain

| Original source | Lines | Demonstrated behavior |
|---|---:|---|
| `InputOutput_ZeeTypeII.f90` | 128–134, 276–338 | MINPAR reader; entry 3 is Lambda3Input, entry 10 is TanBeta |
| same | 213–247, 10891–10984 | Recognizes `ZEEQUARTICIN`, `ZEEPOTIN`; no corresponding bare output-block input reader |
| same | 1663–1682 | Writes bare `ZEEQUARTIC`/`ZEEPOT` at output Q |
| `SPhenoZeeTypeII.f90` | 169, 176–199 | Initializes Q=160 GeV; automatic matching with loops enabled selects order 2 |
| same | 201–282 | Even explicit matching −1 resets scalar couplings from MINPAR inputs |
| `Boundaries_ZeeTypeII.f90` | 1667–1718 | FirstGuess constructs electroweak inputs and matching conditions |
| same | 1724–1752 | Assigns `...IN` values, then **overwrites** independent scalar values with MINPAR inputs |
| same | 1757–1763 | Solves zero-loop tadpoles, then calls TreeMasses |
| same | 1465–1507 | Repeats the LOW scalar boundary assignments during iteration |
| `Model_Data_ZeeTypeII.f90` | 2074–2106 | vd=vSM/sqrt(1+tanbeta²), vu=tanbeta vd; Ye=vSM YeSM/vd−vu Ye2/vd |
| `TadpoleEquations_ZeeTypeII.f90` | 23–69 | M112/M222 solved, including supplied loop tadpoles on later calls |
| `TreeLevelMasses_ZeeTypeII.f90` | 57–88 | TreeMasses calls CalculateMhh and terminates on negative eigenvalues |
| same | 399–493 | Builds and diagonalizes the **tree** CP-even Hessian; sets SignOfMassChanged |
| same | 495–596, 598–699 | CP-odd and charged matrices, including gauge-fixing Goldstones |
| `Boundaries_ZeeTypeII.f90` | 1607–1621 | Later computes loop masses and replaces light Higgs mass by effective-SM result when configured |

The source model is `/home/amine/zee/models/SARAH/ZeeTypeII/ZeeTypeII.m`;
the matching configuration agrees with
`/home/amine/zee/module14/physics_adapters/SPheno.m`. These files explain
the conventions; the actual compiled installation's Fortran determines the
accepted interface. SARAH permits model-specific MINPAR assignments; see its
[input-parameter documentation](https://gitlab.in2p3.fr/goodsell/sarah/-/wikis/Input_parameters_of_SPheno).
`parameters.m` Les Houches metadata alone is insufficient to infer inputs.

| MINPAR index | Parameter | Dimension |
|---:|---|---|
| 1–5 | Lambda1, Lambda2, Lambda3, Lambda4, Lambda5 | 1 |
| 6 | LambdaC | 1 |
| 7–9 | K1, K2, K12 | 1 |
| 10 | tanbeta=vu/vd | 1 |
| 11 | M12Sq (Fortran M12S) | GeV² |
| 12 | MCS2 (Fortran MC2) | GeV² |
| 13 | Mu4 | GeV |

M112/M222 are dependent GeV² parameters. Although `ZEEPOTIN[1,2]` can be
read, the tadpole solver replaces them. Writing the independent parameters
only through `...IN` also fails because of the explicit MINPAR overrides.
The accepted correction uses MINPAR, with default Q=160 GeV and the existing
SM-to-BSM matching. Input tree masses are proposal coordinates, not pole masses.

## Ordered controls and first divergence

1. Original reference: original binary and instrumented copy both give
   mh=484.284443 GeV. The initial tree eigenvalues are
   `(15625.0, 862495.5580124758)` GeV². Its matched neutrino masses reproduce
   `(0.0127855820, 0.0414498832, 0.0625655793)` eV; matrix residual
   `1.72e-14`. This preserves the previous calculator validation, not an
   assertion that the reference is a physical h125/neutrino fit.
2. Original s00001 with unchanged Yukawas: first call stack is exactly
   `FirstGuess|TreeMassesZeeTypeII|CalculateMhh|`. At
   `(vd,vu)=(169.50340582127222,178.58530583760682)` GeV,
   lambda3 is 3.3993935040359005, the sampled **tanbeta**; actual tanbeta
   remains 1.0535794544795800. Other independent scalar inputs remain those
   of the reference, including M12S=351743.6879578036 GeV².
3. Replacing only all four complex Yukawa blocks with the reference values
   gives **identical first tree matrices** and the identical failure.
4. Direct recorded CP-even matrices, in GeV²:

   ```text
   reference: [[447124.547325785, -423358.480585379],
               [-423358.480585379, 430996.010686691]]
   s00001:    [[447124.547325785, -513434.007851738],
               [-513434.007851738, 430996.010686691]]
   ```

   Only lambda3 changes among independent parameters. Tadpoles then alter
   M112/M222, canceling lambda3 from the diagonal entries while changing
   the off-diagonal. The smaller eigenvalue becomes −74437.05578940437 GeV².
   The failure precedes RG/matching iteration and scalar self-energies.
5. Other original scalar inputs reproduce all four reported negative values
   in the same first-stage call. The existing 16-case preflight is consistent
   with this deterministic scalar-only error.

Thus the first divergence is already **card-to-boundary propagation**. The
unintended native potentials really are tachyonic at their first tree guess;
this does not classify the intended analytical potentials as tachyonic.

## Independent scalar reconstruction

Take H1,H2 of hypercharge 1/2 and S+ of hypercharge 1, real scalar couplings,
and neutral fields `(vi+phi_i+i sigma_i)/sqrt(2)`. Define
`a=H1†H1`, `b=H2†H2`, `z=H1†H2`, `c=|S|²`. The implemented potential is

```text
V = M112 a + M222 b − 2 M12Sq Re(z) + MCS2 c
  + Lambda1 a²/2 + Lambda2 b²/2 + Lambda3 ab
  + Lambda4 |z|² + Lambda5 Re(z²) + LambdaC c²
  + K1 ca + K2 cb − 2 K12 c Re(z)
  + 2 Mu4 Re[S− (H1+ H20 − H10 H2+)].
```

With L345=Lambda3+Lambda4+Lambda5, stationarity gives

```text
M112 = M12Sq vu/vd − (Lambda1 vd²+L345 vu²)/2
M222 = M12Sq vd/vu − (Lambda2 vu²+L345 vd²)/2.
M_even² = [[M12Sq tanbeta+Lambda1 vd², −M12Sq+L345 vd vu],
           [−M12Sq+L345 vd vu, M12Sq/tanbeta+Lambda2 vu²]].
mA² = M12Sq/(sinbeta cosbeta) − Lambda5 v².
```

Before stationarity, the CP-even 11 component is
`M112+3 Lambda1 vd²/2+L345 vu²/2`, exactly as generated. The CP-odd
matrix has diagonal `M112+Lambda1 vd²/2+(L345−2 Lambda5)vu²/2` and its
2↔1 counterpart, off-diagonal `−M12Sq+Lambda5 vd vu`, plus the Z
gauge-fixing projector onto `(cosbeta,sinbeta)`.

In charged basis `(H1+,H2+,S+)`, the ungauge-fixed entries are

```text
11 = M112+(Lambda1 vd²+Lambda3 vu²)/2
22 = M222+(Lambda3 vd²+Lambda2 vu²)/2
12 = −M12Sq+(Lambda4+Lambda5)vd vu/2
13 = +Mu4 vu/sqrt(2); 23 = −Mu4 vd/sqrt(2)
33 = MCS2+(K1 vd²+K2 vu²)/2−K12 vd vu.
```

Project onto `HD+=−sinbeta H1++cosbeta H2+` and S+:

```text
M_charged² = [[M12Sq/(sinbeta cosbeta)−(Lambda4+Lambda5)v²/2, −Mu4 v/sqrt(2)],
              [−Mu4 v/sqrt(2), MCS2+(K1 vd²+K2 vu²)/2−K12 vd vu]].
```

Conventions in `Model/model.py`: `h=cosalpha phi1+sinalpha phi2`,
`H=−sinalpha phi1+cosalpha phi2`, and
`H1+=cosgamma HD++singamma S+`, `H2+=−singamma HD++cosgamma S+`.
Masses are ordered light to heavy. Native mixing rows can change overall
sign; compare projectors or rephase consistently. Goldstones occupy the first
native odd/charged row; physical charged rows are 2 and 3. This gives
`Mu4=−sqrt(2)(mH1²−mH2²)cosgamma singamma/v`, matching both the model
and the generated off-diagonal signs. The soft input is positive M12Sq in
the **negative** potential bilinear; it is not the rescaled soft mass
`M12Sq/(sinbeta cosbeta)`. No inversion sign/normalization bug was found.

`potential_hessian` evaluates the polynomial itself and differentiates it
with symmetric differences and Richardson cancellation, exact through quartic
order in exact arithmetic. This avoids using the inversion's mass formulas
as their own check. A step v/8 avoids small-step cancellation. Across the
reference and all four first native matrices, maximum physical-matrix error
is 2.87e-8 GeV² (roundoff); the numerical test bound 1e-6 GeV² is for this
independent derivative comparison only and never changes a physics gate.
The odd/charged comparisons project out gauge-fixing Goldstones.

At each proposal's own v=246.70928224668677 GeV, the same independent Hessian
reconstructs requested mh=125.25 GeV, mH, mA, both charged masses and the
specified rotations. Native FirstGuess instead uses v=246.21965 GeV.
It re-solves tadpoles, giving first mh approximately 125.00143, 125.00453,
125.00499 and 125.01048 GeV for the four corrected points. This small shift
is a vev-scheme effect, distinct from the original large negative eigenvalues.
Later matching, loop tadpoles and the configured loop/effective-SM Higgs
calculation change the spectrum further. The analytical tree target must
never be substituted for that native Higgs pole result.

## Before/after native results

Values below use each scalar's unchanged NO f00 Yukawas, automatic matching,
flag 55=1 and the original effective-SM Higgs setting. All four families per
scalar were evaluated; full family results are in the evidence JSON.

| Scalar | Before smallest tree eigenvalue, GeV² | After mh pole | After mH1± | After mH2± | mH2±−mH1±−mZ, GeV |
|---|---:|---:|---:|---:|---:|
| s00001 | −74437.05578940437 | 146.700634 | 562.896999 | 846.402526 | 192.317927 |
| s00010 | −99521.37184410309 | 139.561224 | 438.451330 | 910.844904 | 381.205974 |
| s00017 | −35290.08577607828 | 153.755294 | 666.861520 | 1184.549060 | 426.499940 |
| s00020 | −130954.31698470516 | 136.195028 | 564.874014 | 1000.846300 | 344.784686 |

Original reference remains mh=484.284443 GeV, mH1±=696.792942 GeV,
mH2±=1272.125990 GeV. Its card and spectrum calculation are unchanged.

The regenerated 30-attempt pilot retains the same 4 analytical scalar
preselections and 16 families. All four complex Yukawa blocks are numerically
**identical** to those in the corresponding old cards. Original-binary native
preflight changes from 0/16 spectra to 16/16 spectra. The corrected independent
scalar couplings propagate to `ZEEQUARTIC`/`ZEEPOT` output at Q=160 GeV;
largest scaled difference is below 5e-9, consistent with SLHA's nine significant
digits. Output tanbeta agrees with entry 10. M112/M222 differ from proposal
values as expected from native vevs and loop tadpoles; they are not ignored
independent physics inputs.

## Controlled matching comparison

Only flag 19 changes in this table, using corrected s00001 NO f00. Flag 55
and the Yukawas remain fixed. Automatic −2 resolves to 2 in this LOW model.

| Flag 19 | Spectrum | mh, GeV | Neutrino execution | Oscillation screen |
|---:|---|---:|---|---|
| −2 | PASS | 146.700634 | PASS | FAIL |
| 0 | PASS | 145.267052 | PASS | FAIL |
| 1 | PASS | 146.626376 | PASS | FAIL |
| 2 | PASS | 146.700634 | PASS | FAIL |
| −1 | FAIL | unavailable | UNKNOWN | UNKNOWN |

All configurations start with the same corrected, positive first tree matrix.
Explicit −1 terminates with `NaN appearing in solution of tadpole equations
for M112`; the existing card lacks a validated complete explicit running-input
definition for that path. It is not an accepted alternative. No matching
configuration was promoted. These comparisons support the input-interface
diagnosis and document scheme sensitivity, not a preferred scheme selection.

## Physics evidence gates

| Claim / gate | Evidence / method and tolerance | Result | Severity / next step |
|---|---|---|---|
| Original native reference reproduced | Original binary, frozen card; same printed masses; matched matrix residual 1.72e-14 | PASS | Regression control only |
| Scalar input propagation | Reader/boundary trace; all 12 independent scalar values in output; nine-digit serialization residual <5e-9 | PASS | Keep per-card provenance |
| Analytical scalar inversion | Independent potential Hessian, all four points, matrix absolute residual <1e-6 GeV² | PASS | Tree-level assertion only |
| Native scalar spectra | Original executable, all 16 regenerated cards | PASS | Does not imply physical acceptance |
| Charged-lepton masses | Native e,mu,tau MASS=(0.00051099893,0.105658372,1.77686) GeV; running mass-direction matrix residual <2e-15 | PASS | Native matching definition preserved |
| Electroweak inputs preserved | Same GF and mZ in cards; output mZ=91.1876 GeV; native running vev/matching retained | PASS | EW precision likelihood remains UNKNOWN |
| Native one-loop neutrino execution | Unique spectrum-associated OneLoopFv snapshot for every corrected card | PASS | Use the spectrum-associated snapshot, not later coupling calls |
| Matched analytic/native neutrino matrices | Existing validated finite-lepton kernel and actual running masses/rotations; maximum residual 2.24e-13 across 16 | PASS | Calculator consistency, not oscillation compatibility |
| h125 necessary mass screen | Existing abs(mh−125.25)<=3 GeV screen | FAIL (16/16) | STOP for accepting these points |
| H2± > H1±+mZ | Native pole masses, all 16 positive margins | PASS | Kinematics only, no rate/BR claim |
| NO/IO oscillation screens | Existing source-versioned NuFIT 6.1 one-dimensional screens with native matched results | FAIL (16/16) | STOP for accepting these points; no refit performed |
| Full unitarity/global vacuum/EWPO/Higgs rates/searches/flavor/full LFV | No new evidence from this spectrum investigation | UNKNOWN | STOP for global allowed-region claims |

The previous neutrino running-input resolution is preserved. For example,
corrected s00001 NO gives `(0.02783708738,0.03338562508,0.07361984279)` eV,
analytical/native matrix residual 1.08e-14, but dm21=3.39697e-4 eV² and
dm31=4.64498e-3 eV²: a clear necessary-screen FAIL. This is a matched
calculation of an incompatible proposal, not a renewed calculator discrepancy.

## Reproduction and regression commands

Use fresh output paths; all scripts refuse overwriting existing directories.

```bash
cd /home/amine/Zee-Pheno
WORK=/home/amine/zee-lfv64-serial-20261112-094549
python3 -m Scan.diagnose_scalar_spectrum \
  --spheno "$WORK/SPheno-4.0.5" \
  --reference "$WORK/model-neutrino-loop-20261118/loop_enabled.in" \
  --pilot "$WORK/scalar_neutrino_pilot_01" \
  --output-dir /tmp/zee-scalar-investigation-reproduce
python3 -m Scan.physical_proposals \
  --config Scan/physical_pilot_config.json \
  --template "$WORK/model-neutrino-loop-20261118/loop_enabled.in" \
  --snapshot Model/reference/native_neutrino_snapshot.json \
  --output-dir /tmp/zee-scalar-pilot-reproduce
python3 -m Scan.spectrum_preflight \
  --proposals /tmp/zee-scalar-pilot-reproduce/proposals.json \
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
  --output-dir /tmp/zee-scalar-preflight-reproduce
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit LFV.test_pilot LFV.test_constraints \
  LFV.test_strict_2hdm_reference LFV.test_gauge_vertices \
  LFV.test_approximation_study LFV.test_nonzero_mixing_study \
  LFV.test_doublet_proxy_study Scan.test_native_pipeline \
  Scan.test_physical_proposals Scan.test_spectrum_preflight \
  Scan.test_diagnose_scalar_spectrum -v
```

The diagnostic needs the installed gfortran/modules/libraries, numpy and the
existing Python dependencies. CI uses the archived native numerical evidence
and synthetic contract tests; CI does not claim to rerun local SPheno.
Observed regressions: 17 Model tests, 60 existing/updated remaining tests,
5 new native scalar tests all PASS. Public external-author integration CI
was not rerun; this change does not alter the LFV amplitudes.

## Release assessment and remaining decisive work

The systematic interface blocker is resolved with a smaller generator. The
production change touches only scalar MINPAR serialization and its required
template entries. `Model/model.py`, `Model/check_native.py`,
`Scan/spectrum_preflight.py` and the neutrino code need no physics change.
Review of the call sites and recorded native parameters confirms the mapping.

Exploratory scans can proceed with fresh cards and mandatory native gates.
The current tree-mass and reference-running-Yukawa proposals do not impose
h125 or oscillation data in the native scheme. A larger scan alone does not
repair that mismatch. The smallest next scientific experiment is a single
predeclared benchmark whose scalar mass definition and lepton/Yukawa inputs
are solved in the **same native renormalized scheme**, then checked against
the existing pole Higgs and oscillation screens. Such a physical fit must
document input changes and its objective; it is separate from making SPheno
converge. No tuning, matching change, new priors, broad scan or higher-order
precision claim was introduced here.
