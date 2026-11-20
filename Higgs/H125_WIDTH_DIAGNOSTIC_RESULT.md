# h125 HiggsSignals BR discrepancy: source-level diagnosis and fixed-loop rerun

20 November 2026. Scope: fixed source-attested s00017 NO and IO spectra.
**Official HiggsTools v1.3 was executed twice, including an independent
SM reference at each point's own mh. The experimental Higgs constraint
is still UNKNOWN.** No 2HDMC or ScannerS was used.

## Finding 1: Source-level error in the original bridge

At upstream official HiggsTools v1.3 git commit
`8732eafee86597c696fca959df52c771bb29bcf2`, the public
`python/Higgs/tools/Input.py` signature of `predictionsFromDict` has:

    calcggH=True
    calcHgamgam=True

It forwards both booleans into `HP.effectiveCouplingInput`.
The original Zee bridge left these settings at their defaults. Therefore
HiggsTools recomputed the loop-induced effective `gg` and `gamgam`
couplings instead of preserving SPheno's full-Zee charged-scalar loop inputs.

**Code fix:** In `ZeeHiggs/higgstools_pilot.py` and
`ZeeHiggs/diagnose_h125.py`, explicitly specify
`calcggH=False, calcHgamgam=False` for all SPheno-based predictions,
and demand all three neutral scalar native `effc_*_gg` and
`effc_*_gamgam` entries be present and finite. Neither SPheno nor
HiggsTools' physics code was modified. The SM reference uses official
SM-like effective couplings. The old default-loop numerical outputs
are retained for an honest before/after comparison, NOT promoted
as the final Zee prediction.

**Independent closure witness:** For NO at mh=125.25 GeV:

    (Gamma_eff(h->gamma gamma)/Gamma_SM(h->gamma gamma)) = 1.0803674937169983
    (kappa_gamma_from_SPheno)^2 = 1.03940728^2 = 1.08036749...
    (Gamma_eff(h->gg)/Gamma_SM(h->gg)) = 1.0638323941347003
    (kappa_gg_from_SPheno)^2 = 1.03142251^2 = 1.06383239...

The ratios agree at printed precision, demonstrating the *corrected
effective branch* consumes the SPheno loop inputs, rather than
regenerating amplitudes from ordinary SM-like fields. These checks
are now regression-tested. They validate the *interface contract*,
not the precision of the original SPheno loop calculation itself.

## Official numerical execution and receipts

- Corrected pinned upstream HiggsTools actual workflow:
  https://github.com/aalida0/Zee-Pheno/actions/runs/44737211306
- Original full native/measurement JSON workflow artifact:
  https://github.com/aalida0/Zee-Pheno/actions/runs/44737211306/artifacts/14630125099
- The actual two-file artifact contains original HiggsTools pilot
  and the entire 22-measurement, channel-resolved diagnostic JSON.
- Hashes of the actual raw JSON are in
  `Higgs/reference/s00017_h125_loop_preserved_summary.json`.
- The original, uncorrected auto-loop reference is
  `Higgs/reference/s00017_higgstools_v1_3_executed.json`.
  Its effective-branch `chi2` is NOT the correct Zee charged-loop
  prescription; it is a transparent historical diagnostic.

### HiggsSignals numerical comparison

| Quantity | s00017 NO | s00017 IO |
|---|---:|---:|
| SM chi2, same Higgs mass and dataset | 151.55556002 | 151.55556065 |
| **Corrected effective BR chi2, with SPheno loops kept** | **167.45043147** | **171.21938194** |
| Delta chi2 corrected effective vs SM | **15.89487145** | **19.66382129** |
| Corrected explicit native SPheno BR chi2 | 179.41680734 | 182.84399375 |
| Delta chi2 explicit vs corrected effective | **11.96637587** | **11.62461181** |
| Historical (auto-loop default) effective BR chi2 | 163.01017406 | 166.75111446 |

The earlier baseline differences +11.45 and +15.20 were missing
the source-supplied loop effects. The corrected effective branch makes
the difference from the SM LARGER, not smaller.

### Higgs h125 widths in NO (MeV)

| Decay | SPheno explicit | HiggsTools corrected effective | SPheno/eff - 1 |
|---|---:|---:|---:|
| total | 3.4521311 | 3.6499125 | -5.42% |
| bb | 1.7725419 | 1.8968073 | -6.55% |
| cc | 0.1061009 | 0.1195143 | -11.22% |
| tau tau | 0.2574170 | 0.2556966 | +0.67% |
| WW* | 0.8622989 | 0.8975259 | -3.92% |
| ZZ* | 0.0939469 | 0.1104188 | -14.92% |
| gg | 0.3478590 | 0.3582457 | -2.90% |
| gamma gamma | 0.0104092 | 0.0101115 | +2.94% |

The old default HiggsTools gamma gamma partial width had been
0.0093201 MeV, i.e. a misleading 11.7% difference from native SPheno;
after preserving the Zee charged-loop input the difference is only 2.94%.

The approximately 6.6% bb discrepancy and 14.9% ZZ* discrepancy
remain **not source-resolved**. They likely involve differences in
SM reference partial widths, perturbative QCD and VV off-shell
treatment, but this is a hypothesis pending a matched external
SM HDECAY/Prophecy4f/HXSWG comparison — not a demonstrated mechanism.
No ad hoc correction factors are applied.

### Dominant experimental sources in NO

Each entry is the official *correlated per-measurement* chi2, not
an independent single-bin sum:

| Experiment/observable | arXiv | Delta chi2 explicit - corrected eff | Delta chi2 corrected eff - SM |
|---|---|---:|---:|
| ATLAS h->gamma gamma 139 fb-1 STXS | 2207.00348 | +6.168 | +6.180 |
| CMS h->tau tau CP 137 fb-1 | 2110.04836 | +2.819 | +3.116 |
| CMS h->tau tau STXS 138 fb-1 | 2204.12957 | +1.973 | +2.144 |
| CMS h->ZZ->4l | 2103.04956 | -1.664 | +1.794 |

One-channel replacement with native SPheno partial widths in the
corrected effective-branch baseline (not additive):
- bb: +10.042 chi2
- gamma gamma: +2.540 chi2
- ZZ*: -2.356 chi2

**Major additional sensitivity:** native SPheno explicitly provided
13 TeV ggF cross section h125 of 48.379 pb. If we remove that explicit
override, HiggsTools with retained SPheno gg effective coupling predicts
51.495 pb and NO HiggsSignals chi2 rises from 167.450 to 179.241
(an additional +11.790 chi2). Thus production input order and source
must also be validated. Do not silently switch to the recomputed
ggH rate to force a stronger or weaker exclusion.

## Working physics prescription for a non-precision Zee scan

Use the corrected HiggsTools **effective-coupling BR** h125 branch
as the *primary diagnostic*, retaining the complete source-computed
Zee gg/gamma gamma amplitudes. Keep the SPheno explicit BR branch as
a theory-prescription sensitivity check. Do not average or take the
difference as a Gaussian error bar. Freeze the h125 production
cross-section convention for the scan and retain provenance.

SPheno native H/A and both H+ physical widths, exotic decays and
masses remain inputs for diagnostic HiggsBounds. Native charged
production cross sections and full HB search-topology coverage
remain missing; no full-model Higgs exclusion/pass follows from
HiggsBounds allowed=true. The HiggsSignals dataset itself is v1.1
(older than 2026 measurements).

A HiggsSignals Delta chi2 relative to a reference Standard Model
is not by itself a valid 95%-CL exclusion without a documented
statistical interpretation and prediction uncertainty treatment.
For both s00017 points:

    physical_higgs = UNKNOWN
    higgsSignals_exclusion = null
    globally_allowed = null

Next validations before promoting a Higgs criterion: independent
h125 SM partial widths (bb and ZZ* first), matching QCD orders,
official h->gamma gamma calculation with both Zee charged scalars,
production ggF/VBF/VH/ttH schemes and full experimental coverage.
