# Charged-Higgs production -> HiggsBounds: source-locked 13 TeV LO bridge

**Status (21 Nov 2026):** implemented and CI-tested **input bridge**;
no *new* MG5 production run for the two `s00017` SPheno cards has been
executed in this repository. The actual 13 TeV source-matched rates
are UNKNOWN. This code does not compute them and does not certify full
HiggsConstraints or a globally allowed Zee point.

The existing HiggsTools v1.3 HiggsSignals workflow already runs with
full native Zee `hgg` and `hγγ` effective loop couplings preserved,
`calcggH=False, calcHgamgam=False`, and the preferred HXSWG-normalized
h125 effective-width prescription. That established scoped h125 screen
is not rewritten here.

## Preferred 13 TeV central production prescription: official HiggsTools tables

The official HiggsTools v1.3 reference API contains
`Higgs.predictions.EffectiveCouplingCxns.ppHpmtb`, documented at
https://higgsbounds.gitlab.io/higgstools/EffC_Refs.html
and based on HXSWG 4FS/5FS matched charged-Higgs cross section
tables. This is a *better leading central reference* than requiring
a new LO 5FS MG integration on every scan point. It is defined for
13 TeV and charged masses 145--2000 GeV. We restrict this adapter to
mass > max(200 GeV, mt+mb) so on-shell top exotic decay is closed,
and set its input BR(t->H+b)=0. Outside the declared domain the
adapter raises rather than interpreting the backend's 0 as a physical
cross-section bound.

The restricted Zee's charged singlet does not couple to quarks. In the
native SPheno basis `CHARGEMIX` rows are (G+,H1+,H2+); columns are
(Hd+,Hu+,charged singlet), with Goldstone doublet (vd,vu)/v. The
physical charged-doublet direction is (vu,-vd,0)/v, so

`d_i = ZP[i,Hd]*vu/v - ZP[i,Hu]*vd/v`,

`kappa_t(H_i) = d_i/tan(beta)`,
`kappa_b(H_i) = d_i*tan(beta)`.

The overall sign of each charged eigenstate cancels from this rate.
The code checks the Goldstone direction and
`d_1^2+d_2^2 ~= 1`, source mass consistency and the
quadratic rescaling of the official rate versus its same-mass
pure-doublet Type-II reference. Both states and the
source `tan(beta)=vu/vd` are re-derived for each spectrum; no fixed
s00017 mixing is reused.

Execute the official HiggsBounds evaluation with this central
tabulated production prescription:

```bash
~/zee-higgstools/venv/bin/python -m ZeeHiggs.higgsbounds_preferred \
  --hbdataset ~/zee-higgstools/hbdataset \
  --charged-reference \
  --output /tmp/zee-hb-tabulated-typeII-charged-reference.json
```

The branch's `higgstools-pilot.yml` also executes and archives this
numerical calculation in GitHub Actions when triggered. This is a
**scoped parametric reference**: assumptions include the Type-II
Yukawa structure, doublet projection factorization, table's
higher-order/MSSM normalization, and the NWA of the experimental
searches. It is not an independent MG5/UFO calculation, and
not every 8/13 TeV or non-tb mode is covered. Even if HiggsBounds
returns `allowed=True`, full Zee Higgs viability remains UNKNOWN.

## Executed official Type-II Zee mixed-charged production and HiggsBounds

**Actual official numerical run 21 Nov 2026:** GitHub Actions run
[44901675064](https://github.com/aalida0/Zee-Pheno/actions/runs/44901675064),
archived artifact ID **14706401577**, official HiggsTools v1.3,
HiggsBounds data v1.7. The run compiled the backend, evaluated both
native NO/IO s00017 spectra, and passed all output assertions.
These are not fake unit-test production numbers.

| Observable (13 TeV) | s00017 NO | s00017 IO |
|---|---:|---:|
| Native mH1 [GeV] | 666.305092 | ~666.305 |
| Native mH2 [GeV] | 1184.55365 | ~1184.554 |
| Doublet H1 amplitude | +0.9887401 | ~same |
| Doublet H2 amplitude | +0.1496431 | ~same |
| Official sigma(tbH1), pb | **0.0118937708214** | **0.0118937360337** |
| Official sigma(tbH2), pb | **2.0415759190e-5** | **2.0415755972e-5** |
| Native BR(H1->tb) | 0.973016169 | 0.970631419 |
| Native BR(H2->H1 Z) | 0.278192521 | 0.278190181 |
| ATLAS 2102.10076 H1->tb selected observed ratio | **0.1021373869** | **0.1018868918** |

H1 and H2 production are inclusive in charge, and are **not** combined
with 13.6 TeV rates from the former repository. The selected search
was ATLAS 139 fb^-1, `2102.10076`,
`LHC13 [Hpmtb>tb]`. H1->tb predicted signal at 13 TeV is
11.57283 fb (NO) or 11.54443 fb (IO). The associated H2 cascade
`sigma(tbH2)*BR(H2->H1Z)` is 0.0056795 fb (NO/IO),
about **17 produced H2 cascade events before daughter BRs and
detector acceptance at 3000 fb^-1**. No inference of signal
sensitivity, full Higgs 95% CL, or global acceptance is made.
H2 does not appear in the selected-limit list: **absence of selected
limit is not an independently validated absence of relevant searches**.

These s00017 points separately fail the already established
**exploratory** h125 HiggsSignals SM-referenced screen. Therefore these
charged results are methodology controls, *not* allowed benchmarks.

## Reusable native Higgs scan gate (not s00017-hardcoded)

The new `ZeeHiggs/native_higgs_scan.py` accepts the **same**
`schema_version=1` manifest of exact `point_id`, `ordering`,
`card`, `native` and `spectrum` used by
`Scan.native_pipeline.py` and `Scan.theory_gate_native.py`.
It source-locks each SPheno output, independently checks physical
scalar masses and native branching fractions, uses its own
`CHARGEMIX`, top mass, scalar masses and `HIGGSLHC13` neutral rates,
then evaluates the official HiggsTools HiggsSignals chi2 and selected
HiggsBounds analyses.

```bash
~/zee-higgstools/venv/bin/python -m ZeeHiggs.native_higgs_scan \
  --manifest /absolute/path/to/actual_native_manifest.json \
  --hbdataset ~/zee-higgstools/hbdataset \
  --hsdataset ~/zee-higgstools/hsdataset \
  --output /tmp/zee-native-higgs-screen.json
```

Every point is evaluated separately; failed source/backend checks are
`NOT_EVALUATED_SOURCE_OR_BACKEND_ERROR`, not a physical Higgs FAIL.
The HiggsSignals Δchi2 against the same-mass SM is **diagnostic**, not a
profiled model likelihood or claimed 95% CL. If charged masses fall
outside the declared tabulated domain, only the charged-production
subset is unavailable. `globally_allowed` stays `null`. A full Zee
HiggsBounds exclusion awaits search topology and production-coverage
validation, nonstandard heavy width closure and uncertainty policy.

The official HiggsTools CI workflow also runs a two-entry NO/IO
manifest through this generic gate, to check that the reusable scan
path reproduces the reference points.

## Existing verified production machinery (different benchmark set)

The previous [aalida0/zee module8](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module8/D6_PRODUCTION_ROBUSTNESS.md)
has real MG5 3.6.7 5FS `gb -> tH_i^-` + conjugate LO runs,
generator/UFO graph checks, and scale variations. Its documented rates
are at **13.6 TeV**, with native points DIFFERENT from `s00017`:
**they cannot be inserted into 13 TeV HiggsBounds**, reused for newly
generated scalars, or treated as 13 TeV rates. The
[old module13 HiggsTools exercise](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module13/D13_HIGGSTOOLS_AUDIT.md)
also contains 13 TeV charged searches for its own D7/F13 samples.
It provides independent evidence for the HiggsTools `Hpmtb` production
mode and the ATLAS/CMS selected charged-search coverage, **not**
source-matched rates for our new cards.

The direct-search input mode used here is official
HiggsPredictions `Hpmtb`: associated `pp -> tb H_i^+-` production,
inclusive in both charge conjugates, in **pb**, in `LHC13` at
`sqrt(s)=13 TeV`. It is distinct from the `brtHpb` entry
produced by SLHA translation: the latter is a *top decay*, not
high-mass heavy-charged-Higgs associated production.
Do **not** double count charge conjugates, sum 4FS+5FS, infer
an associated-production rate solely from a charged mass, or reuse
13.6 TeV results at 13 TeV. No Type-II one-state mass exclusion can
be transferred directly to both mixed Zee charged states.

## How to prepare the two exact SPheno benchmark receipts

Generate **unfilled** (deliberately invalid until real MG integration)
input templates directly from the SHA-checked native NO/IO files:

```bash
cd ~/Zee-Pheno
git pull --ff-only origin research/pre-scan-theory-completion

python3 -m unittest ZeeHiggs.test_charged_production_bridge \
  ZeeHiggs.test_higgsbounds_preferred -v

python3 -m ZeeHiggs.charged_production_bridge \
  --order both --output-dir /tmp/zee-s00017-charged13-unrun
```

For future points, use `template_for_source(source)` where `source`
is the independently audited dictionary for **that specific** newly
generated spectrum. Never substitute s00017 hashes or masses.

To fill each row, use an authenticated Zee UFO sourced from the
same SARAH model and a source-checked MG5 card. The old D5 generator
process syntax for the relevant native 2021-style UFO particle labels
is (repeat in separate output folders for `hm2/hm2c` and `hm3/hm3c`):

```text
# H2-, charge conjugate, and both beam orders (5FS)
generate g d3 > u3 hm3 QCD=1 QED=1
add process g d3bar > u3bar hm3c QCD=1 QED=1
add process d3 g > u3 hm3 QCD=1 QED=1
add process d3bar g > u3bar hm3c QCD=1 QED=1
```

**Do not run this syntax with an unaudited UFO**: confirm the actual
MG5 names/charge assignments and mass eigenstate ordering and preserve
source CKM/Yukawa phases. Set *both* beam energies to **6500 GeV**,
a validated LO five-flavour PDF, consistent running Yukawa and factorization/
renormalization scales, and a source-locked scalar parameter card.
Confirm the banner and result totals, statistical integration error,
diagram set and charge inclusivity. This does not yet automate the
physical SARAH/SPheno-to-UFO mapping.

The receipt JSON rows require, for each PDG 37 and 9000037:

```text
pdg, mass_GeV, mode="Hpmtb", collider="LHC13",
sqrt_s_TeV=13, flavor_scheme="5FS",
perturbative_order="LO", charge_inclusive=true,
sigma_pb, integration_error_pb, generator,
param_card_sha256, run_card_sha256,
process_card_sha256, generator_result_sha256
```

The source `card_sha256`, `spectrum_sha256`, and
`spheno_binary_sha256` must exactly match the audited *same* point.
The receipt hashes identify the external files but are **not**
a proof that the reported rate agrees with an MG5 `results.dat` file.
Such direct numerical agreement, PDF/scale variations and diagram/width
consistency still need an independent generator validation before
assigning a physical charged production bound.

## Apply only when real, checked numerical receipts exist

Store completed files under an external path as
`s00017_NO.json` and `s00017_IO.json` (do not commit private raw
models and multi-GB outputs). Then run:

```bash
~/zee-higgstools/venv/bin/python -m ZeeHiggs.higgsbounds_preferred \
  --hbdataset ~/zee-higgstools/hbdataset \
  --charged-rates-dir /absolute/path/to/real/charged13-rates \
  --output /tmp/zee-hb-source-matched-charged13.json
```

The optional `--charged-rates-dir` injects `cxn_37_Hpmtb_LHC13`
and `cxn_9000037_Hpmtb_LHC13` into both the preferred HXSWG
and historic-SPheno h125 ggF choices. It leaves top-decay
`brtHpb`, neutral H/A, BRs, and all loop-induced couplings untouched.
Malformed, mixed-card, mixed-energy, missing or duplicate receipts
fail closed. Default command (without this option) keeps all charged
production **missing** and reports its physical scope UNKNOWN.

This adapter does **not** make HiggsBounds fully complete. Missing
charged production channels at 8 TeV and other topologies, heavy
H/A production and exotic-width validation, search-specific narrow
width/topology assumptions and calibrated HiggsSignals statistical
interpretation still constrain the final claim. Apply only
searches with demonstrably populated rates and validated topology,
and report all unsupported analyses as **UNEVALUATED**.

## Release decision

- `H125_HIGGSSIGNALS_REFERENCE_SCREEN_READY`:
  already implemented and executed for s00017. Both NO/IO reference
  points are disfavoured by the robust **exploratory** subset screen.
- `CHARGED_13TEV_HPMTB_INPUT_ADAPTER_READY`: yes, with CI tests.
- `CHARGED_13TEV_NATIVE_PRODUCTION_CALCULATED`: **NO** for s00017.
- `FULL_HIGGSBOUNDS_ZEE_ALLOWED_OR_EXCLUDED`: **UNKNOWN**.
- `GLOBAL_ZEE_ALLOWED`: **UNKNOWN**.

**Updated priority:** the HiggsTools tabulated Type-II Zee projection is
the economical **central charged-production approximation** once
official CI execution verifies the API and results. The more explicit
MG5 receipt path described above remains valuable as a complementary
independent test and for mass/topologies outside table validity.
Do not spend large generator resources on the known
HiggsSignals-disfavoured s00017 controls except as a brief interface
validation. For new scan candidates, first require scalar theory
and necessary h125/neutrino screens before costly production work.
