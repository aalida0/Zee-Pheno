# Heavy neutral-Higgs LFV — bounded feasibility audit (2026-11-21)

**Status:** `ANALYTIC_ONLY; GLOBAL_PHYSICS_UNKNOWN; NOT A PUBLICATION ACCEPTANCE`.

The reusable scoped widths implementation is `Higgs/neutral_lfv_tree.py` and its 4 focused tests are `Higgs/test_neutral_lfv_tree.py`. The complete numerical pilot (12 tests, benchmark inputs, plots, scripts) is available as a separate downloadable ZIP in the related ChatGPT conversation, and is **not** claimed to be executed by GitHub CI. This branch was **not** merged into `main`.

## Benchmark and calculation scope

- Restricted CP-conserving Zee with Type-II **quark texture** and general complex lepton Yukawas.
- Scalar spectrum from previous interference feasibility input: `m_h=125, m_H=603.60393898, m_A=622.40011587, m_H1+=600, m_H2+=620` GeV; `alpha=0.6088373713, beta=0.5963066963, gamma=pi/4`.
- Sufficient BFB/tree unitarity passed in prior scalar-only audit; global vacuum and full flavor/Higgs/EWPO unknown.
- NO and IO analytic complex `Fzee,D` witnesses from prior interference pilot, `max|F|=1e-3`; complex neutrino mass closure at leading-log. Replacing leading log with exact charged-lepton finite kernel changes singular masses by at most `7.8e-5` relatively; fixed-F inversion re-solves the original target at `~1e-15` relative matrix precision.
- Yukawa convention: `S=c_beta*Ye+s_beta*Ye2=sqrt(2)M_l/v`, `D=-s_beta*Ye+c_beta*Ye2`; `Y_H=(sin(beta-alpha)S+cos(beta-alpha)D)/sqrt(2)`, `Y_A=D/sqrt(2)`. `Gamma(phi->li lj) = mphi*(|Yij|^2+|Yji|^2)/(8*pi)` for the two charge-conjugate modes at negligible lepton masses.
- The partial-width ledger includes tree `tt,bb,ee,mumu,tautau,emu,etau,mutau`, and `H->WW,ZZ,hh`, `A->Zh`. All are **partial**, so BR = Gamma_mutau/sum(included) is an *optimistic upper bound*. Other channels/QCD/EW corrections and experimental acceptance are absent.

## Calculated results

| Yukawa configuration | H->mu tau BR (upper) | A->mu tau BR (upper) |
|---|---:|---:|
| NO, max|F|=1e-3, max|D|=6.34e-4 | 1.549e-7 | 1.013e-7 |
| IO, max|F|=1e-3, max|D|=9.87e-4 | 3.739e-7 | 2.445e-7 |
| NO, max|F|=1e-5, max|D|=0.0634 | 0.001545 | 0.001011 |
| IO, max|F|=1e-5, max|D|=0.0987 | 0.003713 | 0.002434 |

Because `Mnu` is linear in the product `FD`, `F -> r F, D -> D/r` preserves the analytic neutrino target exactly. The complex fixed-F inverse has rank 5 and 4 complex null directions; neutrino data alone therefore **do not predict a definite LFV branching ratio**.

Using the **single-observable 125 GeV LFV proxy** `BR(h125->mu tau)<0.0018` (ATLAS, JHEP 07 (2023) 166, https://arxiv.org/abs/2302.05225) and a 4.07 MeV otherwise-SM width, the saturation points have:

| Ordering | max|F| | max|D| | optimistic H BR | optimistic A BR |
|---|---:|---:|---:|---:|
| NO | 4.638e-6 | 0.13674 | 0.007126 | 0.004676 |
| IO | 7.207e-6 | 0.13694 | 0.007104 | 0.004666 |

These **ARE NOT PHYSICALLY ACCEPTED**. The h125 proxy omits changes to its total width and production signal strength; low-energy LFV including Barr–Zee, B physics, HiggsTools, complete heavy-Higgs channels, global vacuum and the native neutrino comparison remain unresolved. The pilot's dimensionless LO gluon-fusion amplitude rescalings must not be interpreted as absolute LHC cross sections.

## Native SPheno blocker

From `Model/VALIDATION.md`, corrected native input propagation has an attested user-reported PASS, but the same historical card yields analytic neutrino singular masses `[0.0093018294,0.01272100695,0.05095806075]` eV versus native SPheno loop-enabled `[0.012785582,0.0414498832,0.0625655793]` eV. **Native radiative neutrino masses are NOT VALIDATED**.

The already-committed `Model/probe_neutrino_scaling.py` checks the F=0 and F->2F limits on the attested SPheno binary. It must be executed in the user's local environment, as that binary and its generated Majorana self-energy sources are not present in this runtime. Inspect generated `Sigma1LoopFv` assembly, PMNS/UV basis, input/output running scales and charged-lepton rotations; do not relabel a successful executable run as neutrino-physics acceptance.

## Decision

**NO_GO_CURRENT_HA_MU_TAU_PAPER_CENTER.** Tiny original witness BRs; larger possible rates exploit unconstrained lepton-Yukawa freedom and need complete physical LFV and Higgs constraints. No full scan, MG5/SusHi native production run, recast, or physical Zee benchmark acceptance was performed here.

References: Herrero-García et al. 2017 (arXiv:1701.05345), Julio/Primulando/Uttayarat 2022 (arXiv:2201.01960), Heeck/Thapa 2023 (arXiv:2303.13383), ATLAS 2023 (arXiv:2302.05225).
