# h125 production and decay validation — final executed pilot results

**Date:** 21 November 2026.
**Scientific verdict:** official numerical HXSWG HiggsTools production/BR
validation **executed**, independent pinned public SM HDECAY reference
**executed**, globally acceptable Zee model **UNKNOWN**. No 2HDMC or ScannerS
and no modifications to HiggsTools/HDECAY physics code.

Source evidence and immutable numerical references:

- [Successful source-bound 4-corner HiggsTools v1.3 actual run](https://github.com/aalida0/Zee-Pheno/actions/runs/44895059643)
  (JSON artifact 14703179437, complete NO/IO widths and all four chi2).
- [Successful HDECAY 6.62 SM-only gfortran compilation and execution](https://github.com/aalida0/Zee-Pheno/actions/runs/44895491646)
  (JSON artifact 14703608040).
- [HXSWG/SPheno numerical summary](reference/s00017_HXSWG_ggF_width_validation_summary.json).
- [Independent HDECAY SM output summary](reference/s00017_independent_HDECAY_SM_control.json).
- [Original physics and source-contract specification](H125_PRODUCTION_DECAY_VALIDATION.md).

Both reference calculations evaluate the same h125 native masses: NO 125.250000
and IO 125.250135 GeV. HDECAY prints only 125.250 to 3 decimal places,
but both input files are distinct and source-hashed.

## A. The 13 TeV h125 ggF normalization (pb)

| Observable | NO | IO |
|---|---:|---:|
| Native source-supplied kappa_g^2 | 1.0638323941 | 1.0638323529 |
| SPheno native Zee ggF | 48.3793449 | 48.3792470 |
| Implied historical SPheno SM ggF | 45.47647277 | 45.47638250 |
| HiggsTools SM HXSWG ggF | 48.4050000 | 48.4049055 |
| HXSWG SM times full source-Zee kappa_g^2 | **51.49480704** | **51.49470451** |
| Change in Zee ggF rate, HXSWG vs SPheno | **+6.43965%** | **+6.43966%** |

This difference is a historical SM normalization/order input effect,
not an inferred Zee NP matching discrepancy or calibrated theory error.
The full *source-supplied Zee* gg and gamma-gamma loop couplings are
retained in HiggsTools with calcggH=False, calcHgamgam=False; the
charged-scalar loop contribution must not silently disappear.

### Official HiggsSignals v3.1 chi2: full crossed experiment

Both choices of h125 production and both choices of decay prescription
are evaluated in the **same** official HiggsSignals v1.1 measurement
dataset, with identical h125 mass and all other SPheno input rates.

| h125 ggF normalization | h125 BR prescription | NO chi2 | IO chi2 |
|---|---|---:|---:|
| SPheno old absolute | HXSWG effective-coupling BR | 167.450431 | 171.219382 |
| **HXSWG SM × native Zee kappa_g^2** | **HXSWG effective-coupling BR** | **179.240823** | **182.728126** |
| SPheno old absolute | SPheno explicit BR | 179.416807 | 182.843994 |
| HXSWG SM × native Zee kappa_g^2 | SPheno explicit BR | 192.118116 | 195.197379 |
| Same-mass SM, original HS dataset | native SM | **151.555560** | **151.555561** |

The preferred **exploratory** prescription yields delta chi2 = +27.685263
(NO), +31.172565 (IO) relative to the same-mass SM. These are meaningful
tensions but **not** 95%-CL p-values, model exclusion decisions, or
theory-sigma measurements. The official dataset is old (HSDataSet v1.1);
no new 2026 experimental likelihood was inferred.

H/A production remains source-native SPheno in every corner. The test
does NOT extend an h125 kappa approximation to heavy scalar ggF or
either of the two charged-Higgs production modes.

## B. Independent official HDECAY SM-only calculation

The third-party official HDECAY source is pinned to PSI repository
commit 4884213286e1edbdf3b2bb9fe30095b29e336dff, with source file
SHA256 f83c9e2cbaa42b8590e111a227218055cb703d15c99a5fc23ff0f588642c8a8e.

Its correct **Standard Model** input mode is HIGGS=0, SLHAOUT=0,
COUPVAR=0, 2HDM=0, SM4=0 and MABEG=MAEND=mh, NMA=1.
Earlier diagnostics caught incorrect HIGGS=5 (MSSM) and SLHAOUT=1
(SLHA instead of BR tables); these were corrected without modifying
HDECAY's physics source. The exact original run and parsed BR sums are
retained by the standalone workflow.

At native h125 mass 125.25 GeV:

| SM partial width | Official HDECAY [MeV] | HXSWG/HiggsTools SM [MeV] | HDECAY / HXSWG |
|---|---:|---:|---:|
| bb | 2.360536 | 2.38477 | 0.98984 |
| cc | 0.116673 | 0.11841 | 0.9853 |
| tau tau | 0.256189 | 0.25693 | 0.9971 |
| WW* | 0.902220 | 0.89800 | 1.0047 |
| ZZ* | 0.110604 | 0.11048 | 1.0011 |
| gg | 0.336200 | 0.33675 | 0.9984 |
| gamma gamma | 0.010015 | 0.009359 | **1.0700** |
| **Total SM width** | **4.101** | **4.12289** | **0.99469** |

The HDECAY SM reference has not been retuned to all HXSWG input
masses, couplings, and QCD/EW conventions. In particular the gamma-gamma
SM difference (~7%) remains an input/order distinction. The VV* reference
may differ in detail because the HXSWG prediction uses Prophecy4f
for four-fermion final states.

### Zee partial-width comparison using SM × native kappa^2 proxies

**HDECAY does not calculate our Zee model**. To get a transparent
independent cross-check of tree-dominated partial widths only,
rescale each SM HDECAY channel by the archived source-supplied
Zee effective kappa squared and compare it with native SPheno and
HiggsTools effective-kappa partial widths.

| Source-comparison ratio, s00017 NO | bb | cc | WW* | ZZ* | gg | gamma gamma |
|---|---:|---:|---:|---:|---:|---:|
| Native SPheno Zee / HXSWG-kappa Zee | 0.93449 | 0.88777 | 0.96075 | 0.85082 | 0.97101 | 1.02944 |
| Native SPheno Zee / HDECAY-kappa Zee proxy | **0.94408** | **0.90101** | **0.95626** | **0.84985** | **0.97260** | **0.96208** |

**Conclusion:** The \(\sim15\%\) native SPheno \(ZZ^*\) difference persists
across both SM reference approaches. The bb difference remains \(\sim5.6\%\).
The full explanation requires matching SPheno decay order/parameters and
the off-shell four-fermion treatment. The HDECAY gamma-gamma discrepancy
changes sign when the SM reference is changed, so it is not a basis for
claiming a model-specific 3% or 7% uncertainty. No ad hoc correction is
applied to physical SPheno source outputs.

## Working scan policy and remaining scope

For a reasonable LO-oriented phenomenology paper:
- **Primary diagnostic**: official HiggsTools v1.3
  SM HXSWG h125 production normalization × source-supplied Zee loop
  effective factors, and effective-coupling HXSWG h125 branching ratios,
  with calcggH=False and calcHgamgam=False.
- **Explicit validation variations**: old SPheno absolute 13 TeV ggF,
  and SPheno native explicit BRs, each evaluated separately with
  exact provenance; never convert method differences into Gaussian
  theory errors.
- **Keep UNRESOLVED:** full charged-Higgs production, heavy H/A
  source-matched production, off-shell/exotic widths and scalar LFV,
  HiggsSignals dataset age/theory covariance, HiggsBounds cascade
  search coverage, physical-model 95% CL and globally viable points.

**DO NOT claim the currently tested s00017 NO/IO benchmarks pass any
full Zee Higgs 95% CL constraint.** Numerical experimental-tool execution
and publication-quality full-model acceptance are separate milestones.
