# Zee scalar BFB, unitarity and EWPO source audit — 21 November 2026

**Summary:** Native s00017 NO/IO passes a **rigorous sufficient quartic BFB
certificate**, the complete published high-energy **tree-level scalar
unitarity** matrices, and an **exploratory correlated PDG 2025 S,T,U
screen**. Full global-vacuum stability and complete EWPO convention
validation are not yet certified. Neither point is globally allowed.

## Sources and conventions

- [Zee-type charged cascade study, EPJC 81 (2021),
doi:10.1140/epjc/s10052-021-09956-2](https://link.springer.com/article/10.1140/epjc/s10052-021-09956-2):
restricted scalar potential Eq.17; 2HDM necessary BFB Eq.41, Zee
necessary conditions Eq.45, high-energy unitarity Eq.57 and full
Appendix E scattering matrices E.22–E.29 plus analytic eigenvalues E.31.
- [PDG 2025, Electroweak Model and Constraints on New Physics,
Table 10.8](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf):
right-hand three-parameter fit with U FLOATING.
- Native SPheno 4.0.5 / SARAH 4.15.4; original SPheno binary
SHA256 2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52.

The code verifies the original input card and SPheno production binary
via the existing saved source-attestation, and checks that native
MINPAR[1–9] matches the nine low-scale ZEEQUARTIC entries at Q=160 GeV.
MINPAR[11–13] similarly match ZEEPOT M12Sq,MCS2,Mu4.

## BFB: proof for a large sufficient region, not all possible cases

The previous conservative condition K1 >= |K12|, K2 >= |K12| is replaced by
a less restrictive, rigorous *sufficient* bound. For the exact potential
in Eq.17,

    V4 = V4_2HDM + LambdaC |chi|^4
         + |chi|² (Phi1,Phi2)† [[K1,-K12],[-K12,K2]] (Phi1,Phi2).

If the doublet quartic is strictly BFB, LambdaC > 0 and
K1 >= 0, K2 >= 0, K1*K2 >= K12², each term is nonnegative for
arbitrary complex doublets and charged singlet. This proves quartic
BFB along all field directions. It is NOT a necessary condition:
failure of portal positive-semidefiniteness alone is **UNKNOWN**,
not a physical exclusion. We separately check the necessary
doublet, pure charged-singlet and singlet-doublet axis conditions.
Eq.41 non-strict equality (flat quartic directions) is never
misreported as a necessary-condition failure.

A quartic BFB PASS **does not imply a global neutral vacuum**:
the Zee cubic Mu4 may produce a deeper charge-breaking stationary point
(2021 paper, Sec.4.1.2).

## Full scalar high-energy tree unitarity

We construct exactly five published matrices (Q,Y):
(2,1) 6x6, (1,1) 6x6, (1,0) 6x6, (0,1) 3x3,
(0,0) 13x13. The complete real symmetric eigenvalue spectra are
obtained with numpy.linalg.eigvalsh. Fifteen simple eigenvalues
are checked separately against the published Eq. E.31; the
remaining roots are incorporated via the full 13x13 matrix.

The actual condition is **max |Lambda_i| < 8 pi**, not the
conventional |each quartic| < 4 pi magnitude prior. A regression
demonstrates a counterexample to replacing scattering eigenvalues
with a 4pi prior.

This is the *full high-energy tree scalar* limit at the declared
renormalization scale. It does not replace finite-energy trilinear
exchange, scalar loop unitarity, Yukawa perturbativity, RGE running
or a charged vacuum search.

## One-loop S,T,U: existing native outputs, correlated PDG2025

Original SPhenoLowEnergy entries are index 1=T, 2=S, 3=U, each
explicitly marked 1-loop BSM. The SPheno low-energy enable flag
57=1, comments, numeric values and block provenance are checked.

We use the right-hand PDG 2025 Table 10.8 3D fit:

    mean (S,T,U)=(0.021, 0.040, 0.008)
    std  (S,T,U)=(0.096, 0.120, 0.092)
    rho(S,T)=+0.91; rho(S,U)=-0.62; rho(T,U)=-0.83
    chi2=(STU-mu)^T Cov^-1 (STU-mu)
    three-dimensional 95% ellipsoid chi2 <= 7.814727903

We do NOT substitute the different U=0 fit or three independent
one-dimensional two-sigma bands. The SPheno separate delta(rho)
entry 39 is recorded only; we do not equate it to alpha*T
without checking its source/scheme.

**Limitation:** The native SPheno model-specific one-loop self-energy
implementation and reference SM subtraction have not yet been
independently matched to the full two-charged-state analytic Zee
calculation. The correlated experimental screen is operational,
but complete physical EWPO is still not independently certified.

## Actually executed NO/IO audit

| Quantity | s00017 NO | s00017 IO |
|---|---:|---:|
| Native scalar input/output binding | PASS | PASS |
| Full restricted quartic BFB sufficient proof | PASS | PASS |
| Portal K determinant | 0.0485689782 | 0.0485689782 |
| Maximum absolute high-energy scattering eigenvalue | 5.942698295 | 5.942698295 |
| High-energy scalar tree unitarity < 8pi (~25.133) | PASS | PASS |
| S (native SPheno BSM) | 0.00296941567 | 0.00296941567 |
| T (native SPheno BSM) | 0.0120201035 | 0.0120201035 |
| U (native SPheno BSM) | 0.00005590245 | 0.00005590245 |
| PDG2025 correlated STU chi2 | 0.49490954 | 0.49490954 |
| PDG2025 conditional 3D screen (threshold 7.81473) | PASS | PASS |
| Global charge-breaking vacuum | UNKNOWN | UNKNOWN |
| Full physical one-loop oblique convention crosscheck | PENDING | PENDING |
| Full model global acceptance | UNKNOWN | UNKNOWN |

The two reference points share a scalar-sector setup, hence their
identical native quartics and BSM oblique values. No *new* SPheno card
was generated for this reference audit.

## Deterministic 30-scalar proposal dry run

The original declared physical_pilot_config.json and native snapshot
were reused as independent scalar sampling priors. These are
**ANALYTICAL PROPOSALS**; NO SPheno re-execution occurred.

| Analytic count | Value |
|---|---:|
| Scalar proposals | 30 |
| Passes stricter declared scalar preselection | 3 |
| Screened by one or more declared scalar priors | 27 |
| High-energy scalar unitarity FAIL | 18 |
| Singlet-doublet axis necessary BFB FAIL | 15 |
| Full quartic BFB sufficient certificate PASS | 12 |
| New native spectra | 0 |
| Neutrino or LFV physical results | 0 |

Counts overlap; never sum the failed counts. The previous pilot
cutflow is archived and not silently overwritten. Reproduce
new cutflows in a NEW output directory before scheduling SPheno.
Passing theory prior screens is not a full Zee experimental pass.

## Reproduction and outputs

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/zee-theory-ewpo-audit
    python3 -m unittest Scan.test_zee_scalar_theory Scan.test_declared_theory -v
    python3 -m Scan.zee_ewpo --output /tmp/zee-native-theory-ewpo.json
    python3 -m Scan.zee_theory_cutflow --output /tmp/zee-theory-dryrun30.json

CI executes all tests on original saved native input files and
publishes the complete source-bound numerical JSON and scalar
dryrun. Output files refuse overwrites. These are not results of
a new full-native SPheno scan or a global fit.

**Next:** cross-check SPheno STU one-loop self energies and the
SM reference convention, audit deeper charge-breaking vacua and
appropriate Yukawa/RGE constraints, then perform native new-card
pilot scans and retain full LFV/B/Higgs exclusions separately.
