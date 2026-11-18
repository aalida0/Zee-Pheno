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

**New strict 2HDM neutral-Yukawa map regression:** The corrected public
2HDM convention in arXiv:2410.17313v4, Eqs.(8),(11),(12), maps onto
the restricted Zee lepton mass basis via
`rho_e^dagger = D/sqrt(2)`,
`q_h=(cos(beta-alpha),-sin(beta-alpha))`,
`q_H=(sin(beta-alpha),cos(beta-alpha))`,
and `q_A2=-i` in the chosen real-CP scalar phase convention.
The corresponding up/down Type-II quark `rho` values are
`(m_u/v)cot(beta)` and `-(m_d/v)tan(beta)`.
`LFV/test_constraints.py` checks this complex 3x3 neutral Yukawa
dictionary and the quark scalar/axial couplings for all three neutral
states. This is **vertex-level** transfer validation, **not** the
missing full gauge-invariant two-loop amplitude.

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


## Strict decoupling reference and genuinely new Zee currents (2026-11-18)

A **separate** public-source reference has been implemented at
[LFV/strict_2hdm_reference.py](strict_2hdm_reference.py).
It dynamically imports an **unmodified** author-supplied
[general-2hdm-pheno/thdm.py](https://gitlab.com/jbrod/general-2hdm-pheno)
without vendoring any author source code. Input conditions are frozen:

\[
F=0,\quad \gamma=0,\quad
\mu_4=K_1=K_2=K_{12}=0.
\]

All four singlet–doublet portal parameters are set to zero **solely for
this mathematical reference limit**. Otherwise a singlet-scalar/photon
Barr–Zee term survives even when the charged mixing vanishes.
This is not tuning a parameter point into physical acceptance.

The CP-conserving author convention is
\(\theta_{12}=\beta-\alpha\), \(\theta_{13}=\theta_{23}=0\),
\(\rho_\ell^\dagger=D/\sqrt 2\);
the CP-odd Zee field differs by an overall sign. Type-II quark matrices
use \(\rho_u=(m_u/v)\cot\beta\),
\(\rho_d=-(m_d/v)\tan\beta\) with the **same source-evolved**
weak-scale masses as the author code. The doublet charged trilinears
are projected into its \(Z_3,Z_7\) notation using the independently
validated native restricted scalar potential. This maps neutral vertices;
it does **not** match unseen higher-order topologies.

The source reference splits raw \(c_L,c_R\) functions into 13 source-local
contributions, verifies their sum reproduces the unmodified total, and
retains the original source-file SHA-256. Published dimensionless
\(c_{L,R}\) are mapped to the frozen GeV\(^{-2}\) operator via
\(C_{L,R}=-c_{L,R}/(8\pi^2v^2)\), using the existing source-audited
convention bridge.

**HW normalization caveat:** the public source's
\`cL_llpg_HW\` / \`cR_llpg_HW\` contain an extra \`/self.MW**2\`
relative to **arXiv:2410.17313v4 Eqs.(65,66)**, which are manifestly
dimensionless. The reference evaluator emits *both* the raw author total
and the separate paper-normalized alternative (multiplying only the HW
component by \(M_W^2\)). It does **not** alter the author source or silently
replace the result. This correction is a source-versus-paper dimensional
check; independent amplitude proof of the actual HW coefficient remains
outside this test. Full 2HDM comparisons must distinguish these versions.

The exact charged Zee current decomposition is implemented in
[LFV/gauge_vertices.py](gauge_vertices.py):
\[
R_\gamma=\begin{pmatrix}c_\gamma&s_\gamma\\-s_\gamma&c_\gamma\end{pmatrix},
\quad d=(c_\gamma,-s_\gamma),\ s=(s_\gamma,c_\gamma),
\]
\[
\frac{G_Z}{e}=\frac{\frac12 dd^T-s_W^2 I}{s_Wc_W},\qquad
\frac{G_\gamma}{e}=I,\qquad
\frac{(G_Z)_{12}}e=-\frac{\sin2\gamma}{4s_Wc_W}.
\]
The native charged Yukawa chiral tensors in this same convention are
\[
L_i=-s_i F,\qquad R_i=-d_i D^\dagger.
\]
(The \(R\) matrix written here has neutrino index first and charged-lepton
index second.) The \(\sum_i d_i s_i=0\) unweighted completeness
identity is exact; it does NOT imply cancellation of different-mass loop
functions. These facts exclude naive one-state angle rescaling, especially
for off-diagonal \(ZH_1H_2\) currents.

**Evidence levels:** these current identities and all zero/covariance
limits are algebraic tree-level checks. Successful evaluation of an
unmodified author source is a *2HDM numerical reference*, not a
gauge-parameter cancellation test for the Zee extension. The additional
mixed charged \(W/Z\), \(F\)-dependent reducible/counterterm and
Goldstone/pinch two-loop aggregates remain
\`ZEE_GAUGE_TWO_LOOP_UNVALIDATED\`. Do not promote \`full_lfv_gate\`
above UNKNOWN.


## Independent strict-2HDM amplitude matching (2026-11-18)

External CI used **unmodified** `thdm.py` at public Git commit
`6a54cebcca3bd5bc13b77c15766fc03fffda5199`, SHA-256
`9b4af7aa4bf7bf533a7b82c8c2603d118346424e3b7b5f9f06ab8cb41bc4e193`.
The mathematical control has F=gamma=Mu4=K1=K2=K12=0 and a **self-consistent**
tree scalar potential reconstructed from the declared physical spectrum.
It is not asserted to pass physical Higgs/flavor/neutrino constraints.

The first source run identified a real **chirality-label bug inherited
from the old Module3/4 subset:** with identical effective operator
`L=e m_i/2 bar(l_j) sigma (CL PL+CR PR) l_i F+h.c.`, the direct
same-chirality mu→eγ one-loop residual was **1.02047**, versus
**0.00260525** with L/R exchanged. The physically correct external-mass
flip follows the source's Eq.(48),(51),(52). The standalone
`LFV/core.py` and `LFV/barr_zee.py` on this branch were corrected to
use physical CL and CR consistently. **Previously saved local
`step2_{NO,IO}_lfv_one_loop.json` chiral coefficient *labels* are
obsolete**, while the scalar-subset BRs remain unchanged under swapping
both chirality components.

After that fix, source-level external run
[GitHub Actions 44387509108](https://github.com/aalida0/Zee-Pheno/actions/runs/44387509108)
reported:

| Physical amplitude subset, mu→eγ | Independent comparison | Maximum relative coefficient residual |
| --- | --- | ---: |
| Scalar one-loop CL/CR, including different subleading daughter-mass terms | Public full one-loop vs our Lavoura approximation | 0.00260525 when normalized to dominant chiral coefficient |
| Charged scalar/photon Barr–Zee CL/CR | Public `HA` vs our independently integrated scalar-QED + seagull | 1.44e−13 |
| Top fermion/photon Barr–Zee CL/CR | Public `th` with Z piece algebraically set to zero vs our photon source kernel | 8.75e−15 |
| Charm fermion/photon Barr–Zee | As above, `ch` | 1.52e−12 |
| Bottom fermion/photon Barr–Zee | As above, `bh` | 8.79e−14 |
| Tau fermion/photon Barr–Zee | As above, `tauh` | 1.86e−14 |

The author group `th,ch,bh,tauh` includes photon+Z contributions.
Setting the internal source object's *external Z vector-current factor*
`vZe=0` **only inside an isolated numerical comparison object** projects
onto its electromagnetic form factor. It is not a physical change to the
model and not a way to compute the missing Z or gauge aggregates. The
original author source file is unmodified.

These strong comparisons establish **source consistency of the retained
photon subsets in the strict 2HDM limit**, not their gauge completions
in the mixed charged-singlet Zee model. All `W/Z/Goldstone/ghost/pinch`,
`F†F`, charged mixed-current, `HW` and full 3l/conversion blocks
remain `UNKNOWN`; no full physical LFV experimental PASS/FAIL has
been obtained. Only one representative complex off-diagonal Yukawa
control has been externally executed so far, not a parameter-grid
source comparison.
