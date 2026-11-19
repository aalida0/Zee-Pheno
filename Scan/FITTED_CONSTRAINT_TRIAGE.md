# Fixed NO/IO s00017 — source-bound constraint triage

19 November 2026. Starting point: [native benchmark calibration](NATIVE_BENCHMARK.md)
(commit cb4416b). **This is NOT an allowed-point certificate.** The fitted
points were selected and optimized for native mh and neutrino observables.

The reproducible archived-artifact evaluator is
[\`fitted_benchmark_triage.py\`](fitted_benchmark_triage.py). It re-runs
the validated scalar one-loop \`LFV.pilot\` on the archived, SHA-bound native
cards and reads the actual SPheno 4.0.5 spectrum, FlavorKit, Higgs and decay
blocks. Its direct independent two-body Z width uses **GAUGE at Q=160 GeV,
pole scalar masses**, and CHARGEMIX. This is a tree/running/pole mixed-order
comparison, not a loop-corrected width certification.

## Quantitative baseline

| Quantity | NO | IO | Status |
|---|---:|---:|---|
| mh pole, GeV | 125.250000 | 125.250135 | Higgs mass necessary screen PASS |
| sum neutrino masses, eV | 0.11605997 | 0.13196239 | Cosmology UNKNOWN, model dependent |
| mH1+, GeV | 666.305092 | 666.305622 | Native SPheno |
| mH2+, GeV | 1184.553650 | 1184.553680 | Native SPheno |
| Gamma(H2+), GeV | 13.0165809 | 13.0166628 | SPheno decay table; total-width coverage not certified |
| BR(H2+ -> H1+ Z) | 0.278192521 | 0.278190181 | SPheno decay table |
| Gamma(H2+ -> H1+ Z), GeV | 3.621115455 | 3.621107780 | Native BR × width |
| Independent analytic tree partial width, GeV | 3.621115501 | ~3.6211 | NO verified using printed CHARGEMIX, running g1/g2 and pole masses |
| BR(H2+ -> H1+ h125) | 0.256728625 | 0.256726786 | SPheno decay table |
| BR(H2+ -> H W+) | 0.255540967 | 0.255538721 | SPheno decay table |
| BR(H2+ -> A W+) | 0.203775171 | 0.203773138 | SPheno decay table |
| Gamma(h125), MeV | 3.45213110 | 3.45187339 | SPheno prediction; Higgs-rate gate UNKNOWN |
| BR(B -> Xs gamma) | 3.63471789e-4 | 3.63471789e-4 | SPheno FlavorKit output; physical flavor gate UNKNOWN |
| S | 0.002969416 | 0.002969416 | SPheno; no correlated electroweak fit |
| T | 0.01202010 | 0.01202010 | SPheno; no correlated electroweak fit |
| U | 0.0000559024 | 0.0000559024 | SPheno; no correlated electroweak fit |

The independent tree width uses the physical off-diagonal derivative Z
coupling \(g_{ZH_2H_1}=g_2\sin(2\gamma)/(4c_W)\) and
\(\Gamma=g_{ZH_2H_1}^2\lambda^{3/2}(m_2^2,m_1^2,m_Z^2)/
(16\pi m_2^3m_Z^2)\). At the NO benchmark, the
relative difference to the native BR × total width is 1.27e-8. This is
consistent at SLHA precision, not independent validation of the native
total width or the other branching fractions.

## LFV discrepancy — highest-priority unresolved constraint

The [MEG II result](https://arxiv.org/abs/2504.15711) gives
BR(mu -> e gamma) < 1.5e-13 (90% CL) from 2021–22 data.
The available amplitudes are **not equivalent**:

| Calculation | NO | IO |
|---|---:|---:|
| Independently stored scalar one-loop subset [NATIVE_BENCHMARK](NATIVE_BENCHMARK.md) | 3.6400e-19 | 2.8396e-13 |
| SPheno generated FlavorKitLFV 701 | 2.14943950e-20 | 1.36647538e-14 |
| Ratio (scalar subset / FlavorKit) | 16.93 | 20.78 |
| Subset alone exceeds observed MEG II bound? | No | **Yes** |
| Full Zee LFV experimental gate | **UNKNOWN** | **UNKNOWN** |

**Never mark the IO benchmark excluded on the subset alone.**
The full physical amplitude is coherent. Omitted electroweak,
F-dependent, mixed charged-current and other two-loop terms have no
source-qualified upper bound here. A partial BR can exceed the bound while
the full BR does not; the converse also holds. Equally, the smaller
FlavorKit result is not a validated Zee full-model prediction until its
generated source, diagram/operator coverage, matching, normalizations and
F-dependent contributions are audited. The two numbers **must not** be
interchangeably compared as independent full calculations.

Additional raw SPheno FlavorKit diagnostics, without physical PASS:
- BR(mu -> 3e): 9.63148860e-22 (NO), 1.13448336e-16 (IO).
- CR(mu-e, Al): 2.36665870e-15 (NO), 1.25198714e-11 (IO).
The project's source-matching caveats for conversion and three-body
amplitudes still apply.

**Decisive next LFV work:** trace
\`ZeeTypeII/FlavorKit_*.f90\` and associated generated
\`LowEnergy*.f90\`/source operators in the actual local SPheno installation,
identify exactly what FlavorKitLFV 701 includes, and compare signed complex
CL/CR and chirality conventions to the existing scalar one-loop component
ledger at these two fixed cards. Audit and source-match missing terms rather
than guessing a universal two-loop remainder.

## Other constraints — do not overclaim

- **Higgs:** The native Higgs total widths are numerical inputs, not a
  validated HiggsSignals/HiggsBounds/ATLAS/CMS likelihood. The native
  hbb effective scalar-coupling ratio is approximately 0.892, so current
  channel-dependent signal strengths merit explicit testing. Check mixing,
  production (ggF, VBF, VH, ttH), full partial/total widths and limits
  for H/A/H1+/H2+; a mass-only pass is insufficient.
- **B flavor:** SPheno gives B->Xs gamma = 3.635e-4. Its generated
  FlavorKit scheme, theory uncertainty, operator content and matching
  to the restricted charged-scalar doublet fraction must be independently
  verified before using a significance or exclusion. The previous
  native B_FULL_RATE blocker remains open.
- **EWPO:** S, T, U are numerically small but a proper correlated
  reference SM-subtracted fit and conventions are required. No PASS.
- **Scalar theory:** A sufficient quartic BFB screen and |lambda|<4pi
  do not prove full tree unitarity or global charge-preserving vacuum.
  Both UNKNOWN.
- **Cosmology:** The fitted neutrino sums are 0.1161/0.1320 eV.
  Cosmological limits depend strongly on the assumed cosmological model,
  data combination and statistical prescription. These were not part of
  the neutrino oscillation fit. No universal cosmology PASS/FAIL.
  The lightest neutrino mass is a scan variable and need not remain fixed
  to these numerical control values.

## Release and reproduction

~~~bash
cd ~/Zee-Pheno
python3 -m unittest Scan.test_fitted_benchmark_triage -v
python3 -m Scan.fitted_benchmark_triage \
  --output /tmp/zee-fitted-s00017-triage.json
~~~

\`/tmp/zee-fitted-s00017-triage.json\` stores card/native/spectrum hashes,
neutrino mass sums, the independent cascade width comparison, scalar
one-loop subset versus FlavorKit mu->e gamma, native FlavorKit B->Xs gamma,
S/T/U and partial Higgs widths. Every genuinely pending full-sector
gate remains **UNKNOWN** and \`n_globally_allowed=null\`.

**Verdict:** Cascades have a reproducible nonzero native tree-level
width, and calibrated spectra are useful fixed controls. The full
benchmark viability remains UNKNOWN. Do not generate exclusion contours,
a discovery projection, or a production allowed-region scan from these
two optimized points. Prioritize LFV source coverage and model-compatible
Higgs/direct-search checks, then review B/theory/EWPO/cosmology.
