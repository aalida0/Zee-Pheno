# LFV approximation materiality control — 18 November 2026

**Release decision:** `NO_GENERIC_10PCT_APPROXIMATION`. An explicitly
restricted 2017-style **scalar one-loop + photon Barr–Zee** approximation
was checked against the **complete source-defined general-2HDM** calculation.
The omitted sectors are not guaranteed small even in the decoupled 2HDM
theory. This is a *falsification of a universal neglect assumption*, not a
numerical estimate of the missing Zee charged-singlet gauge amplitude.

## Independently executed control

Source: [jbrod/general-2hdm-pheno](https://gitlab.com/jbrod/general-2hdm-pheno)
at Git commit `6a54cebcca3bd5bc13b77c15766fc03fffda5199`,
unmodified `thdm.py` SHA256
`9b4af7aa4bf7bf533a7b82c8c2603d118346424e3b7b5f9f06ab8cb41bc4e193`.
These exact inputs were checked automatically.

[Successful external GitHub Actions run](https://github.com/aalida0/Zee-Pheno/actions/runs/44394862734)
archives all **18 points × 3 radiative channels = 54** complex-amplitude
comparisons as `materiality.csv` and `materiality.json`:
[results artifact](https://github.com/aalida0/Zee-Pheno/actions/runs/44394862734/artifacts/14470050068).

| Control-grid variable | Values |
| --- | --- |
| tan beta | 1.5, 5, 12 |
| beta-alpha (rad) | -0.15, 0, +0.15 |
| heavy CP-even scalar mass (GeV) | 350, 700 |
| h SM-like mass (GeV) | 125 (fixed) |
| mA / mH+ doublet (GeV) | mH+70 / mH+35 |
| physical charged singlet (GeV) | mH+535 |
| singlet sector | F=gamma=Mu4=K1=K2=K12=0, exact 2HDM |
| scalar Lagrangian | reconstructed at each grid node with m12²=mH² sinβ cosβ |
| flavor textures | fixed complex D with three nonzero LFV transition channels; no fitting/tuning |

All 18 grid nodes were **retained irrespective of their LFV output**.
No neutrino fit, LFV acceptance, HiggsSignals, B-physics, vacuum or direct
searches were applied. Hence these are theory/approximation controls,
**not viable Zee scan points**.

The retained author-source groups are
`1loop, th, ch, bh, tauh, HA` (one-loop, neutral-fermion/photon, and
charged-scalar/photon). The omitted groups are
`tHp, cHp, tauHp, hZ, kin, HZ, HW`. Each complex contribution is extracted
from the **same** author source, and the sum must reconstruct its full
complex chiral coefficient. This avoids cross-scheme numerical artifacts
when quantifying *source-level* neglected terms. It is not an independent
calculation of the full author gauge amplitude.

For each channel:
`R = sqrt(|CL_omitted|²+|CR_omitted|²) /
sqrt(|CL_retained|²+|CR_retained|²)`. Both the unmodified public HW
coefficient and the alternate normalization implied by the published
Eqs.(65-66) are preserved. The alternative does not remove the result.

## Numerical outcome, no cuts or adaptive selection

| Channel | median R | min R | max R | nodes R>0.1 |
| --- | ---: | ---: | ---: | ---: |
| mu->e gamma | **0.69194** | 0.15035 | 15.3554 | 18/18 |
| tau->e gamma | **0.71090** | 0.16271 | 54.4686 | 18/18 |
| tau->mu gamma | **0.71277** | 0.16280 | 57.1533 | 18/18 |

**Dominant single omitted term:** in the same 18-point public-source
control grid, the largest individual omitted complex-amplitude group is
`kin` in **10/18** points and `tHp` in **8/18** points,
**for each of the three radiative channels**. The remaining five omitted
source groups are never the largest in this frozen grid. This refers to
individual term norms, **not** the coherent omitted sum; do not add
their absolute values incoherently or infer a Zee gauge-completion
hierarchy from a source label.

The separate published-HW-normalization alternative also has **18/18**
nodes above 10% for all three channels; the maximum alternative R values
are 14.5568, 51.6358, 54.1810 respectively. Large ratios may result
from cancellations of the retained amplitude. The 10% threshold is an
*explicit arbitrary diagnostic*, never an experimental or theoretical
error confidence interval.

**Interpretation:** The 2017-/2023-style neglect cannot be adopted
universally for a general Type-II-quark, potentially misaligned Zee
parameter scan. The actual 2023 paper also sets second-doublet quark
couplings to zero; that condition is **not true** of the Type-II model
used in our paper. This grid does not estimate genuine Zee
mixed H1/H2/Z or F-dependent W/Goldstone/ghost loops.

## Executable code and release scope

`LFV/approximation_study.py` accepts an *explicit decoupled* point JSON,
public-2HDM SM-parameter dictionary, and the pinned unmodified thdm.py.
The deterministic 18-point control grid and mathematical scalar-potential
reconstruction are in one small function. CI uploads both JSON and CSV.
`LFV/test_approximation_study.py` verifies the exact retained/omitted
source partition, coherent complex cancellation, fixed-grid enumeration
and UNKNOWN physical gate.

**Permitted in paper:** approximation-validation control plots and
**clearly labeled** leading LFV sensitivity overlays conditional on
retained terms. The results above **rule out** describing a universal
10% LFV truncation error, or calling the incomplete model scan the
fully Zee-experimentally allowed parameter space.

**Required before a genuine all-constraints claim:** either match the
corrected **full gauge-invariant** 2HDM electroweak result plus
the Zee singlet extension (including charged mixing, off-diagonal
Z, F-dependent loops/counterterms), or demonstrate quantitatively for
the intended *predeclared model subspace* that the omitted amplitude
cannot materially change the experimental decisions.

The current output status remains `STRICT_2HDM_CONTROL_MATERIALITY_ONLY`;
the **full Zee LFV gate is UNKNOWN**. See
[REMAINDER_AUDIT.md](REMAINDER_AUDIT.md).
