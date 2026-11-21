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

**Priority:** do not spend large generator resources on the known
HiggsSignals-disfavoured s00017 controls except as a brief interface
validation. For new scan candidates, first require scalar theory
and necessary h125/neutrino screens before costly production work.
