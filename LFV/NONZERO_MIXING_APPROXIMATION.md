# LFV with nonzero charged-singlet mixing: reproducible approximation audit

Date: 18 November 2026.
**Verdict: `TRUNCATED_LFV_PHOTON_DIAGNOSTIC_READY`; full physical LFV `UNKNOWN`.**

The 2017 general Zee study and Heeck–Thapa 2023 use leading one-loop
and selected two-loop Barr–Zee terms, with assumptions and neglected
electroweak contributions described in those works.
The published *full* general-2HDM two-loop computation
[Altmannshofer et al. (2025), erratum (2026)](https://arxiv.org/html/2410.17313v4)
does not include the charged-singlet/mixed (F) Zee sector.
The 18-point pure-2HDM source control in
[APPROXIMATION_STUDY.md](APPROXIMATION_STUDY.md) found neglected
source amplitudes above 10% of retained amplitudes in all 18 sampled
points for each of three radiative LFV channels.
Thus a **uniform 10% missing-term uncertainty is not justified**.

## Physical model scope of new code

New module [nonzero_mixing_study.py](nonzero_mixing_study.py) evaluates
*actual* restricted Zee tree-level scalar spectra with
(gamma
e 0), (F
e 0), and kinematically open
(H_2^pm	o H_1^pm Z). At each fixed grid node:
1. Invert prescribed (m_h,m_H,m_A,m_{H_1^pm},m_{H_2^pm},alpha,eta,gamma)
   and independently declared (m_{12}^2,lambda_C,K_1,K_2,K_{12}) to
   (lambda_{1..5},mu_4) using audited scalar inversion.
   No parameter is adjusted in response to LFV.
2. Keep the complex (D) texture fixed; vary (F) by a predeclared
   overall multiplier. This **changes neutrino predictions**:
   it is a *sensitivity scan only*, not a neutrino-compatible scan.
3. Compute coherent physical-source-diagnostic
   (C=C_{m scalar,1loop}+C_{m fermion,photon,BZ}
     +C_{m charged,scalar,photon,BZ})
   from existing independently source-audited code. The complete
   potential is supplied at every point; no missing trilinear is set
   silently to zero.
4. Record per-channel complex retained norm, partial BR, cancellation
   metric, per-topology hierarchy, the exact (F^2) one-loop norm,
   and the exact off-diagonal charged-scalar (Z) coupling.
   The (F^2) control checks scale (Fmapsto 2F) increases its
   amplitude by exactly four.
5. **Do not** calculate or bound the complete electroweak (W,Z),
   Goldstone, ghost, pinch/LSZ and mixed-singlet (F)-dependent
   contributions. Their amplitude upper bound stays null and all
   physical LFV gates are `UNKNOWN`.

Source-matched known-photon terms are *not* a complete two-loop
gauge-invariant Zee amplitude, and the retained branching ratio cannot
be called a model-exclusion prediction.

## Predeclared toy diagnostic grid (not physical Zee scan)

| Axis | Values |
| --- | --- |
| tan beta | 2, 8 |
| gamma (rad) | 0.15, 0.50, 0.95 |
| mH2+ (GeV) | 900, 1200 |
| F multiplier | 1, 100 |
| mH1+, h, H, A (GeV) | 610, 125, 640, 720 |
| scalar inversion | m12sq = 180000 GeV^2; K1=.2, K2=.4, K12=.01, lambdaC=.25 |
| lepton Yukawas | fixed synthetic complex D and antisymmetric nonzero F |
| other masses | explicit illustrative non-matched test inputs |

**24 Zee scalar points × three radiative decays = 72 retained LFV
amplitude evaluations.** All sampled gamma and F values are nonzero,
and the charged cascade is open at all nodes. Points are *not*
filtered for neutrino oscillations, Higgs signal strengths, B physics,
vacuum stability, unitarity, perturbativity or experimental LFV.
The h=125 GeV input is a declared scalar-sector control, not
evidence that other Higgs constraints pass.

[GitHub Actions run](https://github.com/aalida0/Zee-Pheno/actions/runs/44400524298)
completed both the actual-mixing retained-amplitude sweep and a second
auxiliary pure-2HDM reference. Complete per-point CSV and JSON files
are saved as workflow artifacts:
[Zee 24-point retained LFV diagnostics](https://github.com/aalida0/Zee-Pheno/actions/runs/44400524298/artifacts/14470279911);
[24-point pure-2HDM doublet reference](https://github.com/aalida0/Zee-Pheno/actions/runs/44400524298/artifacts/14470154934).

## Exact doublet-sector *proxy*, not a Zee gauge approximation

Module [doublet_proxy_study.py](doublet_proxy_study.py) reconstructs
for every Zee point an auxiliary pure-2HDM Lagrangian with **identical
neutral doublet quartics** (lambda_{1ldots5}), mass spectrum and
lepton (D) texture. The charged doublet mass *prior* to singlet
mixing is exactly

[
m_{H_D^+}^2=cos^2gamma,m_{H_1^+}^2
           +sin^2gamma,m_{H_2^+}^2 .
]

Set (gamma=F=mu_4=K_1=K_2=K_{12}=0) **inside this separate
reference only**; reconstruct and verify the original doublet quartics.
The public author source is run using its own renormalization convention,
with the original 2HDM result decomposed into the same retained/omitted
source pieces as the 18-point decoupling control. Nothing from the
auxiliary 2HDM is added to the full Zee decay amplitude.
In particular, this cannot bound the mixed (ZH_1H_2), (F)-dependent,
charged-current, Goldstone or counterterm corrections.

**Auxiliary pure-2HDM omitted/retained amplitude norms**, 24 points
for each channel (the full mixed Zee amplitude was *not* calculated):

| Process | median omitted/retained | range | points above 10% |
| --- | ---: | ---: | ---: |
| mu->e gamma | 2.76090 | 0.35090–8.59622 | 24/24 |
| tau->e gamma | 2.31069 | 0.34961–5.77328 | 24/24 |
| tau->mu gamma | 2.32478 | 0.35011–5.83899 | 24/24 |

The 10% value is merely a stress-test diagnostic threshold. Ratios
may be inflated by interference cancellations in retained amplitudes.
They are **not** estimates of true Zee two-loop truncation errors.
Even the *proxy* points fail any blanket assertion that the 2HDM
non-photonic sector is numerically unimportant for this region.

## Decision for the phenomenology manuscript

**Defensible now:** plots explicitly titled e.g.
"One-loop + photon Barr–Zee LFV sensitivity (incomplete two-loop
approximation)", with the exact retained diagram classes and numerical
phase interference stated. Use the photon amplitude hierarchy,
the model's exact charged-mixing interactions, and both 2HDM source
stress tests as evidence of what has and has not been approximated.

**Not defensible:** plots or text claiming the *full* LFV-allowed
charged-cascade parameter space, an arbitrary 10% theory uncertainty,
or true exclusions based on only retained LFV amplitudes.
Potentially large omitted 2HDM terms **cannot** be dismissed on
alignment alone, and the genuine Zee mixed/F-dependent remainder is
unbounded. Nonzero mixing needed for the cascade rules out simply
adopting a strict pure-2HDM amplitude for the model.

**Required for unconditional Zee LFV exclusion:** complete the missing
gauge-invariant Zee two-loop amplitude, or derive a rigorously
source-qualified per-point bound in a separately justified restricted
region of parameter space. A review-quality paper may instead state
LFV as an intentionally approximate sensitivity overlay and **not
claim full experimental constraint coverage**.

For local checks, after pulling this branch:

~~~bash
python3 -m unittest LFV.test_nonzero_mixing_study LFV.test_doublet_proxy_study -v
python3 -m LFV.nonzero_mixing_study --help
python3 -m LFV.doublet_proxy_study --help
~~~

The input JSON format is the generic `LFV.constraints` point plus
`scalar_inversion` controls supplied separately. The GitHub workflow
generates standalone **diagnostic-only** example inputs using test
fixtures; for physical points the source must be the model-wide
neutrino/sector-constraint pipeline, which is not connected yet.
