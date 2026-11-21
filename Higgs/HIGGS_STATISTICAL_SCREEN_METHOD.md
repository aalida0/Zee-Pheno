# 2026 HiggsSignals screening: statistical prescription and boundaries

**Result:** The NO and IO fitted s00017 Zee benchmarks are robustly
DISFAVORED by a conditional, SM-referenced two-dimensional HiggsSignals
screen. They have NOT been excluded at a calibrated 95% confidence level.
Full-model Higgs status remains UNKNOWN.

## Authority and statistical choice

The HiggsTools v1.3 authors explicitly discuss this question in
[HiggsTools for LHC Run 3 and Beyond, arXiv:2608.05401,
Section 4.3](https://arxiv.org/abs/2608.05401), especially Eq. (28),
footnote 12 and the discussion of HiggsSignals chi2 interpretation.

They discourage testing the total chi2 against an assumed distribution
with degrees of freedom equal to number of observables minus number of
free model parameters, because correlations and weakly constrained
parameters can make that misleading. Instead they recommend relative
Delta chi2 = chi2(model) - chi2(best fit).
An SM reference chi2 is acceptable as a proxy if the SM fits well and
is a good approximation to the model best fit.

For a genuinely profiled two-parameter inference under Wilks assumptions,
the standard 95% CL chi2_2 threshold is 5.991464547. In a broader scan
without profiling, applying the same threshold is at best a consistent
relative screening criterion, not necessarily a literal 95% CL exclusion.
HiggsTools leaves the model-specific statistical interpretation to users.

In our current Zee study:
- The number of true HiggsSignals observables reported by the pinned
  backend is 159, and SM chi2 is about 151.556 at mh=125.25.
  This is compatible with the SM behaving as a reasonable reference,
  but neither proves that it is the best fit nor that 159 residuals are
  independent Gaussian observations.
- We have NOT profiled all other Zee parameters, fixed a physical
  two-parameter scan, calibrated missing Zee theory nuisances, or
  established a minimum model chi2.
- Therefore the physical 95% CL status is UNKNOWN even if the
  conditional SM-referenced two-dimensional screen says FAIL.
- In particular chi2/159 is NOT used as an exclusion test, nor is
  a one-parameter 3.84 threshold imposed on this many-parameter model.

## Source-bound input

The official HiggsTools v1.3 / HiggsSignals 3.1 was executed with
HiggsSignals dataset v1.1, the same h125 mass and the same NO/IO
Zee effective couplings for all four production/decay prescriptions.
We keep the original full Zee loop inputs by explicitly setting
calcggH=False and calcHgamgam=False.

- [Official numerical source run](https://github.com/aalida0/Zee-Pheno/actions/runs/44895059643).
- [Reference four-corner JSON](reference/s00017_HXSWG_ggF_width_validation_summary.json),
  with original numerical and archive hashes.
- [Zee full-loop closure JSON](reference/s00017_h125_loop_preserved_summary.json).
- [Earlier HiggsBounds numerical diagnostic](reference/s00017_higgstools_v1_3_executed.json).
- [Production and decay validation](H125_PRODUCTION_DECAY_RESULT.md).

## Actual comparison

| Source of h125 ggF / BR | NO chi2 | NO Delta chi2 vs SM | IO chi2 | IO Delta chi2 vs SM |
|---|---:|---:|---:|---:|
| Same-mass SM | 151.555560 | -- | 151.555561 | -- |
| HXSWG ggF / effective BR (PREFERRED) | 179.240823 | **27.685263** | 182.728126 | **31.172565** |
| SPheno ggF / effective BR | 167.450431 | 15.894871 | 171.219382 | 19.663821 |
| HXSWG ggF / explicit BR | 192.118116 | 40.562556 | 195.197379 | 43.641819 |
| SPheno ggF / explicit BR | 179.416807 | 27.861247 | 182.843994 | 31.288433 |

Both fitted points exceed 5.99146 in **all four** source-qualified
calculations. It is therefore not a fragile result driven by a
particular one of these two production or two decay prescriptions.

    s00017_NO   HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL   physical=UNKNOWN
    s00017_IO   HIGGS_TESTED_SUBSET_REFERENCE_SCREEN_FAIL   physical=UNKNOWN

The four-corner spread is NOT a measured or calibrated Gaussian
theory uncertainty. The two-dimensional threshold has NOT been shown to
have 95% coverage for the full Zee parameter space.

## HiggsBounds neutral subset

The earlier official HiggsBounds run on the **historical SPheno ggF**
input yielded allowed=true for the supplied rate inputs; its selected
limits were only for neutral states PDG 25, 35, 36, with the largest
observed-limit ratio about 0.73.

This is a historic numerical diagnostic only. Because the old rate
inputs do not match the preferred HXSWG ggF prescription, the result
is labelled HB_HISTORIC_NEUTRAL_SELECTED_NOT_EXCLUDED, NOT a validated
synchronized neutral limit, and certainly NOT a full Zee Higgs PASS.
Neither charged state has validated complete production/topology coverage.

## Implementation and reproduction

[ZeeHiggs/higgs_subset_screen.py](../ZeeHiggs/higgs_subset_screen.py)
source-verifies the archived run, checks four chi2 differences, uses
the exact two-dimensional chi2 threshold and refuses to convert the
result into a full-model experimental pass or exclusion.

Run:

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/higgs-statistical-screen
    python3 -m unittest ZeeHiggs.test_higgs_subset_screen -v
    python3 -m ZeeHiggs.higgs_subset_screen \
      --output /tmp/zee-s00017-higgs-subset-screen.json

The JSON records experimental versions, source hashes, decision margins,
all four chi2 values, historic HiggsBounds input caveats, and the fields:

    physical_higgs: UNKNOWN
    higgsSignals_95CL_exclusion: null
    higgsBounds_full_model: UNKNOWN
    globally_allowed: null

The reusable pure function screen_point(point_id, sm_chi2, corner_chi2)
may be applied to new source-qualified scan predictions under the same
dataset and conventions. Never hard-veto globally allowed parameters
based only on this exploratory indicator.

## Before a physical 95% CL statement

1. Decide the scan dimensionality and evaluate a Zee HiggsSignals
   best-fit, profiling nuisance/non-interest parameters as necessary.
2. Validate remaining model theory and total-width systematics;
   an SM-only HDECAY benchmark is not a complete Zee NLO computation.
3. Re-run official HiggsBounds with the **same preferred HXSWG
   production normalization** and validate H/A and both charged
   production rates and search topology coverage.
4. Address the age of HiggsSignals v1.1 measurements versus 2026 data.

We can therefore prioritize scanning now while reporting full physical
uncertainty faithfully.
