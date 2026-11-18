# LFV omitted-amplitude materiality decision — 2026-11-18

**Decision:** \`NO_UNIFORM_REMAINDER_BOUND\`. The restricted Zee model's full
radiative LFV prediction is **UNKNOWN**. The validated one-loop + selected
photon Barr–Zee amplitude is a reproducible, useful **scoped approximation**,
but cannot be assigned a general, source-justified fractional error or
used as a model-independent LFV PASS/FAIL throughout a parameter scan.
No parameter point was tuned and no diagrams were invented.

## What the primary sources actually support

| Reference / location | Actual assumption or result | Transfer to our model |
| --- | --- | --- |
| Herrero-García et al., [JHEP04(2017)130](https://arxiv.org/abs/1701.05345), Sec.3.3.2, Eqs.(3.13)–(3.19) | Gives leading Zee radiative formulas and selected Barr–Zee terms; explicitly assumes \`f^(e mu)=0\` for a vanishing charged contribution (Eq.3.15), while Eq.(3.19) retains F; Sec.3.3.1 uses leading tree trileptons and notes loop dipoles | Leading term cross-check only; not arbitrary-F complete 2-loop/mu3e/conversion oracle |
| Heeck–Thapa, [PLB841(2023)137910](https://arxiv.org/abs/2303.13383), Sec.2–3 | Sets couplings of second doublet to quarks to zero, assumes neutral CP-even alignment; explicitly neglects Z, charged-scalar and W-associated Barr–Zee contributions, arguing suppression in its slice | **Not** a bound for generally misaligned/type-II-quark Zee points |
| Altmannshofer et al., [JHEP06(2025)156](https://arxiv.org/html/2410.17313v4), Secs.4–6, especially kinetic Eq.(36) | Gauge-independent **sum** of neutral-kinetic terms, including charged-Higgs-mass-dependent \`f8e(xWh,xWH+)\`; removing individual terms breaks gauge independence | Complete **pure 2HDM limit** reference; no charged singlet F and only one physical charged Higgs |
| [Erratum JHEP07(2026)184](https://doi.org/10.1007/JHEP07(2026)184) | Corrects the 2025 general-2HDM result; source-locked earlier in-house audit identified changed signs/prefactors and kinetic function | Only use erratum-aware code, then match the entire gauge-invariant set to Zee before claiming inclusion |
| Earlier [Zee LFV_COVERAGE](https://github.com/aalida0/zee/blob/research/zee-full-constraints/research/full_constraints/LFV_COVERAGE.md) and [MATCHING_AUDIT](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module4/MATCHING_AUDIT.md) | Complete operator/topology inventory, signed photon subset, W transverse candidate discrepancy, uncancelled F-dependent gauge subgraphs, generated FlavorKit conversion nuclear-index discrepancy | These are **positive source evidence for missing work**, not numerical bounds on its total |

The 2017 source's model assumptions and leading formula are sufficient for
a **declared leading approximation** in a phenomenology paper. They do not
supply a rigorous, general remainder estimate or an arbitrary-parameter-point
physical exclusion. The 2023 alignment assumptions cannot be imposed
silently. The corrected public 2HDM amplitude must be treated as a *whole*;
a W-only or mixing-angle-rescaled subset is not gauge invariant in general.

## Exact missing operator/tensor blocks

1. **Radiative decays:** electroweak kinetic W/Goldstone/ghost/pinch/
   reducible/LSZ contributions; doublet–singlet mixed physical charged
   current; \`F†F\` gauge dressing; Zγ fermion and mixed charged-scalar
   currents; WZ and remaining potential-sector terms. The prior isolated
   \`GWW\` subdivergence and F-neutrino self-energy gauge dependence
   require the missing partners to cancel.
2. **Three-lepton decays:** K1 non-dipole photon, K2 full dipole,
   scalar/vector penguins and all relevant boxes, with complex signs and
   crossed identical-fermion interference. A coherent rate *functional*
   is available upstream but its missing operator coefficients are not.
3. **Muon conversion:** complete electromagnetic/non-photonic quark
   operators and heavy-quark/nuclear matching. The source-locked unmodified
   FlavorKit scalar proton/neutron assignment differs from the direct
   \`Z*g_p + N*g_n\` form, so its output is not automatically a validated
   independent check.
4. **Scheme/remainder:** matched Yukawa/quark mass scale, leading-QED
   evolution, physical pole vs running inputs, finite external invariants.
   Numerical quadrature error alone is NOT a physical theory uncertainty.

A complete two-loop **2HDM** result does not make the mixed charged-singlet
**Zee** correction complete; the latter contains new charged spectra,
off-diagonal Z transitions and F-dependent electroweak topologies.

## Quantitative materiality at the native NO/IO method fixtures

Both fixed cards had an invalid Higgs spectrum (PDG25≈484 GeV), so these are
**mathematical LFV diagnostics**, not accepted physical points. Their native
SPheno lepton masses and input hashes were verified in Step 1. This table
uses the initial **scalar one-loop** LFV pilot BRs, *not* the newer coherent
one-loop+Barr–Zee output, which has not been run on these cards.

Source: MEG II [EPJC85(2025), arXiv:2504.15711](https://arxiv.org/abs/2504.15711),
\`B(mu->e gamma)<1.5e-13\` at 90% CL; normalization from the archived
Module3 scalar pilot \`alpha=0.0072973525692838015\`,
\`m_mu=0.105658372\`, \`m_e=0.00051099893\` GeV,
\`Gamma_mu=2.995983701908041e-19\` GeV.

Define \`K=alpha*m_mu**5/(4*Gamma_mu)*(1-m_e**2/m_mu**2)**3\`,
and \`B=K*||C||²\` with \`||C||=sqrt(|CL|²+|CR|²)\`.
Then \`K=8.01779578456e10 GeV^4\`,
\`C_limit=sqrt(1.5e-13/K)=1.36778593736e-12 GeV^-2\`.

| Quantity | NO | IO |
| --- | ---: | ---: |
| Scalar one-loop BR(mu→eγ) only | 1.23598633497e−17 | 2.42199337497e−17 |
| Retained scalar one-loop \`||C||\` (GeV^-2) | 1.24159324010e−14 | 1.73803685017e−14 |
| Sufficient absolute omitted \`||ΔC||\` bound for passing MEG II (GeV^-2) | **<1.35537000495e−12** | **<1.35040556885e−12** |
| \`C_limit/||C_scalar||\` | 110.16 | 78.70 |
| 2026 source-justified bound on **combined** omitted terms | **NONE** | **NONE** |
| Experimental physical gate | **UNKNOWN** | **UNKNOWN** |

These bounds are **requirements on a future verified amplitude bound**,
not estimates of omitted amplitudes. No assumption that \`ΔC=0\` is made.
Likewise, the scalar one-loop \`CR\` h/H/A amplitudes nearly cancel:
the magnitude of their sum divided by sum of magnitudes is approximately
0.06084 (NO) and 0.09116 (IO), from the user-run component output.
Hence *relative* omitted/retained ratios can blow up near cancellations,
even when the absolute missing contributions are perturbatively ordinary.

The same fixtures have beta−alpha≈0.05267 rad, not exact alignment;
\`|sin(beta-alpha)cos(beta-alpha)|≈0.05257\`. This is a **coupling
suppression indicator** for some neutral-W LFV structures, not a bound
on the complete two-loop Zee gauge aggregate.

## What can be defended for a modest phenomenology paper?

**Allowed scope now:** parameter-space mapping with an explicitly
labeled leading/subset LFV sensitivity overlay, together with separately
certified scalar, neutrino, B-flavor, Higgs and vacuum gates. Show the
retained contribution ledger and emphasize that the LFV overlay is
*approximation dependent*. Never report \`LFV_ALLOWED\` based on this.

**Not allowed:** a plot called the full Zee experimentally allowed region
or any claimed model-wide LFV exclusion/acceptance without a demonstrated
complete leading gauge-invariant contribution (or a reliable per-point
bound on every omitted material amplitude). A scan of arbitrary Yukawa
textures has no universal finite, small **relative** error because
interference can drive the retained amplitude toward zero.

**Minimal decisive next technical task:** independently match the corrected
2HDM *full gauge-invariant* lepton-dipole aggregate in the strict
charged-singlet-decoupling (\`F=0\`, zero charged singlet mixing) limit,
including source-normalized \`CL,CR\` and updated code revision. Only then
derive and test the genuine Zee correction from two physical charged
eigenstates, off-diagonal Z currents and F-dependent W/neutrino diagrams,
including \`R_xi\` cancellation and field renormalization. The existing
code is missing that proof; neither a guessed numerical sigma nor one
transverse W form factor closes this gate.

**Release decision:** \`LFV_LO_SCOPED_READY\`,
\`LFV_FULL_PHYSICAL_UNKNOWN\`, \`SCANS_MAY_ONLY_SHOW_CONDITIONAL_LFV\`.
No additional arbitrary-BSM two-loop amplitudes implemented in this audit.
