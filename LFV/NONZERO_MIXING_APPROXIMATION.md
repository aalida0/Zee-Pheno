# Nonzero-mixing Zee LFV: approximation and materiality audit

18 November 2026 | Status: TRUNCATED_LFV_PHOTON_DIAGNOSTIC_READY
Full experimental LFV constraint: UNKNOWN.

This is a **methodological validation**, not a scan of allowed Zee model
parameter space. There is no retuning of points in response to LFV limits.

## Published theory basis

The 2017 general Zee phenomenology study and 2023 Zee LFV analysis
computed one-loop and selected photon Barr-Zee contributions with
specified flavor and scalar-sector assumptions. Neither supplies a
complete mixed-singlet two-loop electroweak amplitude for our Type-II
quark, general charged-Higgs-mixing Zee parameter space.

The corrected full general-2HDM radiative LFV calculation
(Altmannshofer et al., arXiv:2410.17313v4 and July 2026 erratum) is a
valuable *reference*, not an automatic Zee prediction. The earlier
18-point unmodified-author-source comparison found neglected diagrams
above a 10% amplitude diagnostic threshold at all sampled controls:
[APPROXIMATION_STUDY.md](APPROXIMATION_STUDY.md).

## Implemented on actual Zee parameters (gamma != 0, F != 0)

[nonzero_mixing_study.py](nonzero_mixing_study.py) accepts an arbitrary
Zee scalar/Yukawa LFV input, an explicit set of restricted-potential
inversion inputs, and a predeclared four-axis control grid. At each point:

1. Reconstruct Lambda1..5, K1, K2, K12 and Mu4 from the physical
   charged/neutral masses, alpha, beta, gamma, declared soft mass and
   scalar parameters. No potential parameter is changed to satisfy LFV.
2. Preserve all complex D entries and vary antisymmetric complex F by a
   declared multiplier. This is a texture sensitivity study: multiplying
   F generally changes the neutrino masses and is NOT a neutrino-fit scan.
3. Calculate the coherent source-checked neutral/charged scalar
   one-loop and fermion-/charged-scalar-photon Barr-Zee contributions.
4. Record the retained complex-amplitude norm, partial branching
   ratio, interference/cancellation metric, known F-only one-loop
   contribution, the retained photon two-loop/one-loop ratio and exact
   charged mixing Z(H1,H2) current.
5. Preserve a null theoretical bound on omitted electroweak gauge
   amplitudes and UNKNOWN physical LFV for every point. No experimental
   PASS/FAIL is inferred from a retained partial branching ratio.

The actual scalar theory checks are only diagnostic: all quartic
magnitudes below 4*pi and the **necessary** CP-conserving 2HDM
doublet-subpotential bounded-from-below inequalities:
Lambda1>0, Lambda2>0,
Lambda3 > -sqrt(Lambda1*Lambda2),
Lambda3+Lambda4-|Lambda5| > -sqrt(Lambda1*Lambda2).
These are NOT a full Zee BFB, global minimum, unitarity or Higgs-rate gate.

## Fixed revised scalar consistency control (24 points)

The *initial* trial grid used an inconsistent fixed m12sq at high
tan(beta), producing nonperturbatively large quartics. Those initial
numbers are SUPERSEDED and must not be used as physical evidence.
The revised parameters are declared **before** evaluating LFV.

| Control | Revised fixed input |
| --- | --- |
| tan(beta) | 2, 8 |
| gamma [rad] | 0.15, 0.35, 0.55 |
| m(H2+) [GeV] | 900, 1200 |
| F overall multiplier | 1, 100 |
| m(H1+) [GeV] | 650 |
| m(h), m(H), m(A) [GeV] | 125, 640, 720 |
| CP-even alignment | alpha = beta |
| Soft mass | Msoft = 640 GeV, m12sq = Msoft^2 sin(beta) cos(beta) |
| Quartic controls | LambdaC=0.25, K1=0.2, K2=0.4, K12=0.01 |
| D and base F | Fixed, illustrative complex flavor textures |
| Other mass scheme | Explicit but NOT matched physical benchmark scheme |

24 genuine mixed Zee scalar input points x 3 radiative decays = 72
included-topology amplitude evaluations. The on-shell charged cascade
H2+ -> H1+ Z is open at every point. The original 125 GeV h input
and the scalar quartic diagnostics do not certify phenomenological
viability; neutrino, B, full Higgs, vacuum and collider constraints are
NOT evaluated.

## Exact doublet-sector auxiliary pure-2HDM reference (NOT Zee two-loop)

[doublet_proxy_study.py](doublet_proxy_study.py) reconstructs an auxiliary
**different model** from each Zee input. It preserves the exact neutral
doublet quartics and D lepton Yukawas, but removes F, charged-singlet
mixing and singlet portals only in that auxiliary control.
Its unmixed doublet charged mass obeys:

    m(HD+)^2 = cos(gamma)^2 * m(H1+)^2
             + sin(gamma)^2 * m(H2+)^2.

The doublet quartics are independently reconstructed to verify this
identity. The public full 2HDM source is then evaluated using its own
mass/Yukawa renormalization scheme. The additional 2HDM electroweak
contributions remain **outside** the physical Zee amplitude.
These auxiliary values are *neither numerical estimates nor bounds*
for mixed H1/H2 and F-dependent Zee gauge graphs.

Revised 24-point auxiliary pure-2HDM omitted/retained amplitude norms:

| Channel | median | minimum | maximum | >10% |
| --- | ---: | ---: | ---: | ---: |
| mu -> e gamma | 0.50909 | 0.42992 | 0.54478 | 24/24 |
| tau -> e gamma | 0.51984 | 0.43354 | 0.56362 | 24/24 |
| tau -> mu gamma | 0.52022 | 0.43366 | 0.56430 | 24/24 |

These are coherent complex-amplitude *norm ratios*, NOT BR errors.
Even in the aligned auxiliary model, the omitted 2HDM amplitudes
are not generically negligible in this control region. For the
genuine Zee theory the mixed charged-scalar and F-dependent
electroweak loops are still UNKNOWN.

Reproduction is automated in GitHub Actions: run both the unmodified
source and the mixed Zee photon subset. The job uploads three separate
CSV/JSON artifacts: legacy strict-2HDM controls, revised 24-point
nonzero-mixing Zee partial amplitudes, and revised 24-point auxiliary
doublet-reference amplitudes. The public author source is pinned to
Git commit 6a54cebcca3bd5bc13b77c15766fc03fffda5199 and thdm.py
SHA256 9b4af7aa4bf7bf533a7b82c8c2603d118346424e3b7b5f9f06ab8cb41bc4e193.

## Defensible first-paper language

Use: "LFV sensitivity at one loop plus selected fermion-/charged-scalar-
photon Barr-Zee topologies (incomplete electroweak two-loop treatment)."
Show interference and a materiality/approximation warning. Interpret
LFV overlays as *conditional sensitivity*, never as full-model
acceptance or exclusion contours.

Do not: assert a uniform 10% error, merge the auxiliary 2HDM result
into the Zee amplitude, claim experimentally allowed charged-cascade
parameter regions, or hide missing F-dependent gauge/pinch terms.

A complete model-wide LFV exclusion still requires a validated Zee
gauge-invariant two-loop calculation or a defensible per-point bound
for all omitted material amplitudes. No paper can legitimately claim
a full LFV constraint based only on the retained subset.

## Local reproduction

From research/lfv-approximation-study:

    python3 -m unittest LFV.test_nonzero_mixing_study LFV.test_doublet_proxy_study -v
    python3 -m LFV.nonzero_mixing_study --help
    python3 -m LFV.doublet_proxy_study --help

The GitHub Actions workflow creates explicit *synthetic* source inputs
under /tmp for this audit. To apply the method to real fitted points,
first export the source-verified neutrino-fit/SPheno point with the
right lepton and quark mass schemes, then run the per-point LFV
constraint evaluator and preserve its UNKNOWN full-gauge status.
