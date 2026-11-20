# B->Xs gamma: defensible exploratory-screen prescription (2026)

**Implemented: exploratory screening ONLY. Physical exclusion UNKNOWN for
every point.** This method does not provide a complete new-physics
uncertainty likelihood and must not be presented as a 95% CL exclusion.

## Published inputs

Reference: Misiak et al., arXiv:2608.29863v1 (2026), Eqs. (1), (3), (4),
https://arxiv.org/abs/2608.29863.

CP/isospin averaged B->Xs gamma at photon E>1.6 GeV:
- Experimental average: (3.49 +/- 0.19)*1e-4.
- Modern SM NNLO: (3.54 +/- 0.14)*1e-4.
- Generic NP linear approximation, with real LO ordinary Wilson coefficients
  at mu0=160 GeV:
  1e4 * BR_linear = 3.54 - 8.51*DeltaC7(0) - 2.16*DeltaC8(0).

The paper's linear formula neglects quadratic NP terms and NP NLO matching;
its full Type-II 670 GeV 95%-CL bound includes distinct higher-order
contributions and is NOT reproduced here. Do not apply mH1>670 GeV as a
hard exclusion in a Zee model with two mixed charged states.

## Screening algorithm

The archived executed source is:
Flavor/reference/hepfit_restored/s00017/bsgamma_nlo.json.
Raw spectra and cards, rebuilt driver and 37 historical validations are
documented in Flavor/HEPFIT_RESTORATION.md.

For each point, compute:
1. Delta BR(NLO) = BR(HEPfit Zee) - BR(HEPfit SM), **same-scheme**.
2. Hybrid diagnostic H = 3.54 + 1e4*Delta BR(NLO), in 1e-4 units.
   H is not a consistently matched NNLO BSM prediction.
3. Published linear diagnostic L = 3.54 - 8.51*C7LO - 2.16*C8LO.
   Use our source-bound real BMU/NDR LO coefficients at mu0=160 GeV.
4. For orientation ONLY, nominal Gaussian interval
   I = [3.49 +/- 1.96*sqrt(0.19**2+0.14**2)]
     = [3.0274234,3.9525766].
   The errors refer to the modern SM and experimental averages.
   These are not a validated full-model BSM likelihood or the paper's
   reproduced confidence construction.
5. Define method *difference* d=abs(H-L). The engineer's triage guard
   G=[min(H,L)-d,max(H,L)+d]. The multiplier 1 on d is a declared
   non-probabilistic heuristic, **not** a theory 1-sigma or rigorous
   missing-higher-orders bound. Neither method is independent of our
   common scalar matching.
6. SCREEN_PASS only if G is wholly inside I;
   SCREEN_FAIL only if G is wholly outside I on one side;
   SCREEN_BORDERLINE for every overlap/straddle case.
   ALL screen labels are preliminary priority tags, not physical
   passes or exclusions. Never discard a scan point on this basis.
7. Enforce B_PHYSICAL_UNKNOWN, theory_sigma=null,
   experimental_exclusion=null, globally_allowed=null for every output.

The older FlavorKit prediction of 3.634718e-4 has a historical fixed
SM baseline 3.15e-4 and **is not an independent uncertainty estimate**.
It is not included in the screening guard.

## s00017 results (rates in units 1e-4)

| Quantity | NO | IO |
|---|---:|---:|
| Old-scheme HEPfit NLO SM | 3.4914241274 | 3.4914241281 |
| Old-scheme HEPfit Zee NLO | 3.9214198286 | 3.9214193541 |
| Same-scheme NP shift | +0.4299957012 | +0.4299952259 |
| Hybrid H = modern SM + NLO shift | 3.9699957012 | 3.9699952259 |
| Linear published L | 4.0465445806 | approximately 4.046544 |
| d = |H-L|, NOT a sigma | 0.0765488795 | approximately 0.076549 |
| Method guard G | [3.893447,4.123093] | approximately same |
| Nominal exp+SM Gaussian band I | [3.027423,3.952577] | same |
| Smallest downward change to enter I | 0.0174191 | approximately same |
| Change / computed NP shift | 4.05% | 4.05% |
| Hypothetical independent Gaussian NP 1-sigma to reach 1.96 | 0.0653753 | approximately same |
| **Exploratory screen** | **B_SCREEN_BORDERLINE** | **B_SCREEN_BORDERLINE** |
| **Physical B constraint** | **B_PHYSICAL_UNKNOWN** | **B_PHYSICAL_UNKNOWN** |

The hypothetical sigma is a sensitivity calculation, **not an assigned
uncertainty**. It is approximately 15% of our NP shift, but there is no
basis in the current results for assigning this Gaussian distribution,
magnitude or independence. Likewise the 4.05% correction is the signed
change sufficient to cross the illustrative boundary, NOT a confidence
error estimate.

## Requirements for actual physical PASS/FAIL

1. Reproduce an appropriate published same-order one-state rate and
   confidence construction. In particular, our old HEPfit historical
   2HDM table is not a reproducible implementation of the new 2026
   670 GeV 95%-CL one-state limit.
2. Quantify the missing finite strange-quark/primed, internal light-quark,
   NP four-quark and Zee EW contributions in the real top-only NLO
   matching domain and test extrapolation across scan points.
3. Establish a documented Wilson-dependent theory treatment for
   perturbative scales, HQE/kinetic-input correlations, nonlocal/resolved
   photons, photon-cut extrapolation and any NP matching remainder.
   Avoid double counting effects included in HEPfit, and do not reuse
   a generic NNLO SM 4% or historical Type-II 7% nuisance without proof.
4. Only then construct a model-consistent 95% CL test and decide
   the physical gate. Uncertainties must be source-qualified.

## Reproduction

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/bsgamma-conservative-screen
    python3 -m unittest Flavor.test_bsgamma_screen -v
    python3 -m Flavor.bsgamma_screen --output /tmp/zee-s00017-bsgamma-screen.json

Output refuses overwrite and records the exact upstream archived JSON SHA,
Gaussian reference, both rates, model spread, screening guard and UNKNOWN
physics gates. The library function screen_prediction can later process
new valid scan points, but the archived CLI deliberately authenticates
only fixed s00017 NO/IO native cards and source-bound HEPfit runs.
