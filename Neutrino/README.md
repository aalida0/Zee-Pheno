# Step 1 — neutrino oscillations

Minimal **pilot**, not a physical scan. **fit.py** uses the Zee one-loop
finite-mass formula checked in [Model/NEUTRINO_ROOT_CAUSE.md](../Model/NEUTRINO_ROOT_CAUSE.md)
and the saved **native running parameters** to invert an NO or IO NuFIT target.

**Experimental input:** [NuFIT 6.1 (2025), IC24 *with SK atmospheric*](https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf),
lower block, separate NO/IO best fits and one-dimensional 3-sigma intervals.
The convention is Δm²3ℓ = Δm²31 for NO and Δm²32 for IO. Selected best fits are
**targets**, not an official correlated NuFIT likelihood. Majorana phases and
the lightest neutrino mass are user inputs.

## Generate NO and IO proposals

Requires Python 3 and NumPy; native validation needs the prebuilt SPheno.
Run from the repository root:

~~~bash
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit -v

WORK="$HOME/zee-lfv64-serial-20261112-094549"
CARD="$WORK/model-neutrino-loop-20261118/loop_enabled.in"

python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --reference-card "$CARD" --write-card "$WORK/step1_NO_proposal.in" \
  --output "$WORK/step1_NO_proposal.json"

python3 -m Neutrino.fit --ordering IO --mlight-ev 0.001 \
  --reference-card "$CARD" --write-card "$WORK/step1_IO_proposal.in" \
  --output "$WORK/step1_IO_proposal.json"
~~~

The **single-target** `Neutrino.fit` CLI defaults to a minimum-norm
rank-five inverse. For the scalar-first scan, its `inverse()` function
now also supports both null-cone roots, alternate projective charts,
complex F phases, and arbitrary four-complex-dimensional homogeneous
D coefficients. See [Scan/README.md](../Scan/README.md) for the
deterministic multi-family proposal generator. These solutions are
algebraic at an archived charged-lepton running reference; none
is a native SPheno fit until revalidated on its own card.
All outputs remain **ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN** until native
revalidation. Files use exclusive creation: use new filenames when rerunning.

The proposal uses the *reference* running masses and mixing in the saved
[one-card native snapshot](../Model/reference/native_neutrino_snapshot.json).
It changes **only** the four SLHA real/imag blocks for Fzee and Ye2. SPheno
reconstructs Ye, running leptons and scalar inputs during matching; the first
proposal **must not** be assumed to pass oscillation constraints natively.

## Check native masses and PMNS

Run instrumented SPheno in a fresh directory without changing the original
installation:

~~~bash
python3 -m Model.diagnose_neutrino_native \
  --spheno "$WORK/SPheno-4.0.5" \
  --card "$WORK/step1_NO_proposal.in" \
  --output-dir "$WORK/step1_NO_native"

python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --native-result "$WORK/step1_NO_native/result.json" \
  --output "$WORK/step1_NO_native_screen.json"
~~~

Repeat for IO with new output directories. The screen uses the *one-loop*
native Takagi rotation and charged-lepton rotation to extract mixing angles,
δCP and mass splittings, checking published one-dimensional 3-sigma ranges.
**NATIVE_3SIGMA_SCREEN_PASS** is not a correlated χ² or full Zee acceptance.
A failed first proposal requires retuning in the native matching scheme.

To make the next **analytic** iteration from the actual new native running
parameters, pass its result as the snapshot while retaining the original
hash-pinned loop-enabled card as the template:

~~~bash
python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --snapshot "$WORK/step1_NO_native/result.json" \
  --reference-card "$CARD" --write-card "$WORK/step1_NO_iteration2.in" \
  --output "$WORK/step1_NO_iteration2.json"
~~~

This is a small fixed-point proposal, **not a certified convergent algorithm**;
rerun native diagnostics after each change and stop on mismatched quantities,
tachyons, tadpole failure, or nonperturbative couplings. Do not overwrite cards.


**Step 1 benchmark gate: PASS for the two native reference cards.**
Those cards may now be used as input to the **Step 2 LFV pilot**, with their
one-loop Yukawa scheme and input hashes preserved. This is not permission to
treat an unconstrained scan as neutrino-allowed: new scalar points, lightest
masses, or Majorana phases require fresh native fit validation.
Matching-scale stability is an optional **theory-systematic study** for a
full parameter scan, not an observed failure of the existing two native cards.
A full correlated NuFIT likelihood is not implemented; only published
one-dimensional 3σ screens are claimed. Do not promote the full physical model
before LFV, flavor, Higgs and vacuum gates.


## User-run native benchmark results (2026-11-18)

The original user-run SPheno 4.0.5 diagnostic produced **native one-loop
oscillation agreement for both NO and IO** on the proposed cards. These are
benchmarks, not an accepted global Zee parameter scan.

| Observable | NO native | IO native |
| --- | ---: | ---: |
| m1, m2, m3 (eV) | .00929997206, .01272241882, .05096558374 | .04907780206, .04983977318, .001000018592 |
| sin²θ12 | .3088016078 | .3088013655 |
| sin²θ13 | .0224798411 | .0226197952 |
| sin²θ23 | .4699999626 | .5499995722 |
| δCP (deg) | 212.0000835 | 273.9999020 |
| Δm21² (eV²) | 7.53704604e−5 | 7.53723353e−5 |
| Δm3ℓ² (eV²) | 2.51100125e−3 | −2.48300295e−3 |
| Native vs analytic complex matrix residual | 5.892e−14 | 1.092e−15 |
| Individual NuFIT 6.1 3σ screen | **PASS (6/6)** | **PASS (6/6)** |

Proposed-card SHA-256 values recorded in the generating reports:
NO `71bb94b67a4db185a6d44fa56c96da751799211f5da23761408de75a5acbbdad`;
IO `26e705a9115cafd0d20833ac5420803f998f48013f97b98649526fb6c35b6e4d`.
Native run paths: `$WORK/step1_NO_native/result.json` and
`$WORK/step1_IO_native/result.json`; screen outputs:
`$WORK/step1_NO_native_screen.json` and
`$WORK/step1_IO_native_screen.json`. These locally generated results were
provided in the project conversation; they are not bundled into Git.

**Interpretation:** source scheme and chosen reference card are native-matched
for these points; no retuning was needed. The generic top-level
`ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN` in the CLI is the status of the
initial *proposal object*, not of the nested `native_screen`. Do not use the
top-level string to override the actual user-run
`NATIVE_3SIGMA_SCREEN_PASS` evidence.

**Scope:** central-point and separate 1D NuFIT 3σ agreement demonstrated at
the reference matching conditions. Matching-scale sensitivity and correlated
NuFIT likelihood are not evaluated; neither is the full LFV/Higgs/B/vacuum
acceptance. The preliminary IO `dmsq21` printed by
`diagnose_neutrino_native` refers to *ascending* singular-value ordering and
must not be interpreted as physical IO Δm21². The **screen** explicitly maps
IO to (m1,m2,m3)=(middle,largest,lightest) before computing splittings.



## Absolute-mass benchmarks (conditional, not extra oscillation constraints)

For the specified **zero Majorana phases** and the two chosen lightest
masses, the standard three-neutrino expressions give:

| Derived parameter | NO | IO |
| --- | ---: | ---: |
| Σmν (eV) | 0.072987998 | 0.099917539 |
| mβ (eV) | 0.012871998 | 0.048753638 |
| mββ (eV) | 0.010676045 | 0.048175223 |

These numbers are *derived predictions* at the chosen fixed lightest mass
and phase inputs. They are **not** exclusions from cosmological, beta-decay
or neutrinoless-double-beta-decay data. NuFIT does not fix these two
Majorana phases or the absolute neutrino mass. The new fit tests cover
additional lightest masses and phase pairs *algebraically*; those variants
have **not** been rerun natively.

## Native result provenance check

For final benchmark acceptance, bind the native output to its exact proposal
card, rather than trusting a free-standing report:

~~~bash
python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --native-result "$WORK/step1_NO_native/result.json" \
  --checked-card "$WORK/step1_NO_proposal.in" \
  --output "$WORK/step1_NO_cardbound_report.json"

python3 -m Neutrino.fit --ordering IO --mlight-ev 0.001 \
  --native-result "$WORK/step1_IO_native/result.json" \
  --checked-card "$WORK/step1_IO_proposal.in" \
  --output "$WORK/step1_IO_cardbound_report.json"
~~~

A PASS is labeled `NATIVE_3SIGMA_PASS_CARD_BOUND` only when the native
input SHA-256 matches the original candidate input card. The native SPheno
output is still a **one-card benchmark**, not a full NuFIT joint likelihood.
