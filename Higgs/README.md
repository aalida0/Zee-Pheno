# Higgs constraints in Zee-Pheno: SPheno -> official HiggsTools

**Status: native NO/IO SPheno inputs audited and official HiggsTools v1.3
executed in CI; full Zee-model Higgs acceptance is still UNKNOWN.**

[Executed numerical pilot report](HIGGSTOOLS_PILOT_RESULT.md) documents the
actual HiggsBounds selected analyses, matched-mass Standard Model
HiggsSignals reference and large h125 decay-precision sensitivity.

The successful official numerical run and complete JSON evidence are available
in [GitHub Actions run 44732975159](https://github.com/aalida0/Zee-Pheno/actions/runs/44732975159).

No 2HDMC or ScannerS is used. Our Python modules live in the **ZeeHiggs** namespace to avoid shadowing the upstream HiggsTools package, whose own import name is **Higgs**. Native SARAH 4.15.4 / SPheno 4.0.5 is the
source of the full Zee mass spectrum, all physical charged mixing matrices,
tree/loop decays and effective neutral-Higgs couplings. Official
**HiggsTools v1.3** supplies experimental HiggsBounds 6.1 / HiggsSignals 3.1.
Its numerical engine and published experimental datasets are NOT copied,
changed, or reimplemented in Zee-Pheno.

Upstream references:
- HiggsTools v1.3: https://gitlab.com/higgsbounds/higgstools/-/tags/v1.3
- HiggsTools input: https://higgsbounds.gitlab.io/higgstools/alternative_input.html
- HiggsBounds: https://higgsbounds.gitlab.io/higgstools/HiggsBoundsAPI.html
- HiggsSignals: https://higgsbounds.gitlab.io/higgstools/HiggsSignalsAPI.html
- HiggsBounds dataset v1.7: https://gitlab.com/higgsbounds/hbdataset/-/tags/v1.7
- HiggsSignals dataset v1.1: https://gitlab.com/higgsbounds/hsdataset/-/tags/v1.1

## Implementation

- [ZeeHiggs/spheno_audit.py](../ZeeHiggs/spheno_audit.py) verifies original SPheno binary, exact
  native card association, spectrum mass/total width/decay BR closure,
  HiggsCouplingsFermions, HiggsCouplingsBosons, and SPheno HiggsLHC13 block.
  Only the two frozen fitted s00017 NO/IO cards are initially supported.
- [ZeeHiggs/higgstools_pilot.py](../ZeeHiggs/higgstools_pilot.py) calls the official
  Higgs.tools.Input.readHB5SLHA with neutralPDGs [25,35,36] and
  chargedPDGs [37,9000037], detects the imported particle IDs by their
  audited masses, and calls predictionsFromDict, HiggsBounds and
  HiggsSignals. It reports both all-explicit SPheno decay BR input and the
  HiggsTools-recommended SM-like-h125 effective-coupling BR prescription.
  Both computations use the **same official experimental datasets**.
- Official HiggsTools supports multiple singly charged Higgs bosons and
  nonstandard particle IDs. Neither charged mass eigenstate is discarded.
- SPheno HiggsLHC13 gluon-fusion predictions are copied as **explicit
  approximate cross sections** for h/H/A, replacing any inferred ggH 13 TeV
  cross section. They are based on historical SusHi-1.5.0 rescaling and
  MUST NOT be silently called modern QCD-accurate production rates.
- SPheno's h/H/A Zgamma effective coupling is explicitly marked
  "not yet calculated"; such entries are removed from the imported
  effective-coupling dictionary where present. Their physical prediction
  is NOT zero and its uncertainty/coverage is NOT validated.
- Charged H1/H2 production cross sections are **absent** from native
  SPheno data. Do NOT assume their values are zero, approximate them by
  masses alone, or claim full charged-Higgs direct-search exclusions.
- The h125 SPheno total width and BRs are not certified at the precision
  required by signal-strength tests. Both all-explicit BR and effective-
  coupling h125 predictions are retained as a sensitivity check. The
  HiggsTools input documentation specifically recommends effective
  coupling input for SM-like h when precision BRs are insufficient.
- H, A and charged widths/decay channels are currently native SPheno
  predictions, NOT source-matched to HDECAY/SusHi for all relevant modes.
  Charged scalar EW/topology models and non-standard lepton Yukawas are
  not guaranteed to be represented by naive Type-II scaling.

**CRITICAL**: A HiggsBounds numerical allowed/excluded output is a result for
the **provided inputs**. As long as charged-production and other important
predictions are missing, it is NOT a globally valid Zee exclusion/pass. An
HiggsSignals chi2 is likewise not a complete likelihood or 95% CL label until
the dataset, prediction precision, degrees of freedom and relevant
uncertainties are validated.

The adapter therefore always writes:

- physical_higgs = UNKNOWN
- higgssignals_confidence = UNKNOWN
- higgsbounds_full_model_constraint = UNKNOWN
- globally_allowed = null

No missing values are converted into a physical PASS or FAIL. The routine
repository unit tests exercise source-bound SPheno inputs and mocked APIs,
while the separate official HiggsTools workflow has **completed a real
HiggsBounds/HiggsSignals calculation** on the two native benchmarks.
The experimental constraint is not physics-certified merely because the
external program executes.

## Installation on the user's local machine (persistent and isolated)

A compiler with C++17 support and CMake >=3.17 is required. The external
HiggsTools source may download its C++ dependencies at build time. This
script installs only in a user-owned Python virtual environment and
user-owned source directories; it does not install 2HDMC/ScannerS.

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/spheno-higgstools-pilot
    python3 -m unittest ZeeHiggs.test_spheno_audit -v
    bash Higgs/install_higgstools.sh

The reproducible tags in the install script are HiggsTools v1.3
(2026-09-17), HBDataSet v1.7 and HSDataSet v1.1. The installer prints the
resolved Git SHA of each tag. These datasets have different experimental
update histories, so record the tag/commit for the publication.

Run the **actual official HiggsTools evaluation**:

    cd ~/Zee-Pheno
    ~/zee-higgstools/venv/bin/python -m ZeeHiggs.higgstools_pilot \
      --hbdataset ~/zee-higgstools/hbdataset \
      --hsdataset ~/zee-higgstools/hsdataset \
      --order both \
      --output /tmp/zee-higgstools-s00017-pilot.json

And separately generate the source-bound SLHA audit:

    python3 -m ZeeHiggs.spheno_audit \
      --output /tmp/zee-higgs-s00017-input-audit.json

The outputs refuse overwrite. The backend evaluation intentionally fails
if HiggsTools is uninstalled, required HB/HS datasets are absent, native
masses cannot be uniquely identified, or the supported official API is
incompatible. **Do not silence those errors or set missing production to
zero to make it run.** Fix the actual interface with version-matched
source/documentation and tests.

## Remaining scientific work before full Higgs constraint

1. Numerically execute the official HiggsBounds/HiggsSignals pilot and
   inspect selected analyses and the sensitivity to h125 branching ratio
   treatment. Verify the datasets actually load and applicable searches
   are present; archive raw outputs, version hashes and selected channels.
2. Validate the SPheno h125 total width and h->gamma gamma/Zgamma
   calculations (including both mixed charged scalars). Source-match
   ggF/bbH/VBF/VH/ttH production and QCD precision at the relevant energies.
3. Supply both charged scalar production rates (not only doublet weights),
   check cascade final-state acceptance/topology coverage, and audit H/A
   exotic widths, multibody/off-shell decays and nonstandard lepton couplings.
4. Establish consistent collider likelihood application with theoretical
   uncertainties. An unsupported search topology is NOT an exclusion and
   cannot establish that a point is globally viable.

No LO-independent new Higgs likelihood or detector recast is attempted.
The first genuinely executed HiggsTools NO/IO diagnostic is complete.
The next task is understanding the ~15-point HiggsSignals chi2 spread
between width treatments and supplying source-validated charged-Higgs
production. Numerical results and validation limits are archived in
[HIGGSTOOLS_PILOT_RESULT.md](HIGGSTOOLS_PILOT_RESULT.md).
