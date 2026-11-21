# h125 Higgs production and decay — source-matched validation

**Scientific status:** reproducible production/decay comparisons intended for
exploratory Zee HiggsSignals inputs, **not** a full Zee 95% CL exclusion.
No 2HDMC or ScannerS. SPheno 4.0.5 (SARAH 4.15.4) remains the full Zee
spectrum source. HiggsTools 1.3 (HiggsPredictions, HiggsSignals 3.1) is
the experimental reference; the standalone public **HDECAY 6.62** is used
ONLY as an independent SM partial-width cross-check.

Relevant sources:
- LHCHXSWG production cross sections: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/HiggsXSBR
- LHCHXSWG HDECAY + Prophecy4f branching fractions:
  https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGBRs
- HiggsTools physics and reference model:
  https://higgsbounds.gitlab.io/higgstools/alternative_input.html
- Independent HDECAY source (SM mode only):
  https://gitea.psi.ch/ltpth/HDECAY
  pinned to commit 4884213286e1edbdf3b2bb9fe30095b29e336dff
- HDECAY/Prophecy4f role in SM Higgs decay widths:
  https://www.sciencedirect.com/science/article/pii/S0010465520301454

## Correction: 13 TeV h125 gluon fusion

The source-attested SPheno effective h125 gluon coupling at s00017 NO
is kappa_g=1.03142251, so kappa_g^2=1.063832394.
This must be preserved (including the model's full loop amplitude
according to native source), with calcggH=False and
calcHgamgam=False in HiggsTools 1.3.

| Source / policy | Zee h125 ggF [pb] | SM equivalent [pb] |
|---|---:|---:|
| old native SPheno HiggsLHC13 (SusHi1.5 rescaled) | 48.3793449 | 45.476... |
| official HiggsTools SMHiggsInterp times kappa_g^2 | 51.4948070 | 48.405... |

The independent LHCHXSWG literature gives an approximately 48.6 pb SM
reference for 13 TeV ggF at 125 GeV. The 6.4% discrepancy is primarily
the older SM production normalization; the native kappa_g agrees between
both. This is a scheme/order/source difference **not** a 6.4% Gaussian
theory error. Source-linked exact mass-specific results are computed in
ZeeHiggs/validate_reference.py.

**Adopt** for exploratory h125 signal strengths:
- HiggsTools/HXSWG SM 13 TeV production normalization times native
  source-bound Zee kappa_g^2. Preserve the native heavy H/A approximate
  production inputs unchanged for this one-parameter comparison.
- Keep source-native effective h125 gammagamma and gg coupling inputs:
  calcHgamgam=False and calcggH=False.
- Construct h125 BRs by effective-coupling reweighting of HiggsTools'
  HXSWG SM partial widths. Use SPheno explicit BRs and legacy native
  h125 ggF only for declared sensitivity variations.
- Do NOT apply the same simple scaling blindly to heavy H/A ggF
  (thresholds, bottom-interference, different K factors) or charged H+
  production. These require separate source-matched calculations.

**Crossed four-policy study:** each NO/IO benchmark is evaluated through
the *same* HiggsSignals v1.1 data and same m_h reference with
(1) native SPheno vs HXSWG h125 13 TeV ggF, and
(2) native explicit vs HXSWG-rescaled h125 BRs.
This separates rate normalization from decay-width effects. The
differences are **not** theoretical error distributions or additive
chi2 contributions.

## Independent decay validation

Standalone HDECAY 6.62 is built from the pinned official Fortran source.
The same source is executed in SM mode at m_h=125.250000 (NO)
and m_h=125.250135 GeV (IO), recording its unmodified SM predictions
br.sm1/br.sm2, total width, source hashes and exact input files.

This *does not calculate the Zee model*. Its SM partial widths, rescaled
by the source-bound Zee kappa factors, provide a transparent third
channel-by-channel reference for bb, cc, tau tau, WW*, ZZ*, gg, gamma gamma.
These proxies are only matched-order approximations. In particular
HDECAY's stand-alone WW/ZZ off-shell prediction may differ from the
LHCHXSWG Prophecy4f-based reference, and the SM input masses and
electroweak parameters are not identical across the three codes.
The comparison *must* record these limitations, not fit a correction
factor onto the actual native SPheno widths.

The official SM branching fractions are based on HDECAY with
Prophecy4f for four-fermion decays. This is why we prefer the
HiggsTools/XSWG effective-coupling reference over an arbitrary
SPheno-to-SM total-width rescaling for a 125 GeV Higgs.
A physically complete NLO electroweak Zee calculation (including
two charged scalars) remains outside this validation.

## Execute

The CI workflow compiles official HiggsTools and HDECAY from pinned
sources and runs the diagnostic. From a Linux machine with gfortran,
CMake, a C++ compiler, Git and a compatible HiggsTools Python venv:

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/higgs-production-decay-validation
    bash Higgs/install_higgstools.sh
    python3 -m unittest ZeeHiggs.test_validate_reference \
      ZeeHiggs.test_hdecay_sm_reference -v
    ~/zee-higgstools/venv/bin/python -m ZeeHiggs.validate_reference \
      --hsdataset ~/zee-higgstools/hsdataset \
      --output /tmp/zee-h125-HXSWG-native-compare.json
    bash Higgs/run_hdecay_SM_reference.sh "$HOME/zee-hdecay-sm-reference"

The CI job archives full JSON outputs. It must be executed before any
new numerical claims; untested external build instructions are not
numerical evidence. Do not overwrite original s00017 card files or
historical receipts.

## What the result does NOT establish

- HiggsSignals p-values or a validated 95% CL cutoff;
- complete 2026 measurement coverage (HSDataSet v1.1 is older);
- model-specific theory covariance or uncertain higher-order Zee rates;
- heavy neutral Higgs ggF/bbH QCD accuracy;
- both charged Higgs cross sections and direct-search topology coverage;
- full EW corrections to exotic Higgs decays or the h125 Zgamma amplitude;
- a globally allowed parameter point.

Physical Higgs acceptance remains UNKNOWN pending these audits.
