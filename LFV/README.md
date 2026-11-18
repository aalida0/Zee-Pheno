# Step 2 — lepton-flavor violation

**Promoted scope:** preliminary scalar one-loop photon dipoles for
mu→e gamma, tau→e gamma, tau→mu gamma using the actual SPheno **running**
lepton Yukawas and charged mixing of validated NO/IO neutrino reference
points. **No physical LFV acceptance, exclusion, or complete prediction.**

The small core.py reuses the earlier source-audited
[aalida0/zee module3/core.py](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module3/core.py),
its [operator conventions](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module3/CONVENTIONS.md)
and Lavoura [hep-ph/0302221](https://arxiv.org/abs/hep-ph/0302221).
The coherent dipole has seven scalar one-loop pieces: h/H/A-lepton, plus
H1/H2 doublet-D and singlet-F charged-scalar diagrams. Signed complex
C_L and C_R are retained in the declared dipole operator. The leading
kernel neglects finite neutrino masses and external invariants in its
loop denominators.

**Not calculated here:** remaining Barr-Zee/two-loop contributions,
mixed F-dependent charged-current/gauge diagrams, mu→3e complete
penguins/boxes/dipoles, mu-e conversion and interference.
Subset BR values cannot be used as an allowed/excluded verdict.
The older [LFV_COVERAGE.md](https://github.com/aalida0/zee/blob/research/zee-full-constraints/research/full_constraints/LFV_COVERAGE.md)
documents these outstanding pieces and their published source matching.

## Execute on NO and IO native reference cards

~~~bash
cd ~/Zee-Pheno
python3 -m unittest LFV.test_pilot -v
WORK="$HOME/zee-lfv64-serial-20261112-094549"
~~~

The native spectrum contains neutral scalar masses (PDG 25,35,36).
The CP-even mixing angle **alpha is extracted automatically** from the same
SPheno `SCALARMIX` block as the neutral pole masses. For the supplied NO
spectrum, `ZH(1,:)=(0.725651649,0.688062268)` gives
`alpha=atan2(ZH12,ZH11)=0.7588153265 rad` in the frozen convention.
The historical example value 1.398 is NOT applicable.
Mixed running couplings and pole masses are an explicitly mixed-order
leading-order scheme diagnostic.

~~~bash
# CP-even alpha is parsed from SCALARMIX; no manual angle is needed.
python3 -m LFV.pilot \
  --native "$WORK/step1_NO_native/result.json" \
  --card "$WORK/step1_NO_proposal.in" \
  --spectrum "$WORK/step1_NO_native/output.spc" \
  --output "$WORK/step2_NO_lfv_one_loop.json"

python3 -m LFV.pilot \
  --native "$WORK/step1_IO_native/result.json" \
  --card "$WORK/step1_IO_proposal.in" \
  --spectrum "$WORK/step1_IO_native/output.spc" \
  --output "$WORK/step2_IO_lfv_one_loop.json"
~~~

You can inspect the native neutral masses/mixing blocks with:

~~~bash
grep -niE '^[[:space:]]*Block[[:space:]]+(MASS|SCALARMIX|ZH|ALPHA)' \
  "$WORK/step1_NO_native/output.spc"
~~~

The tool fails unless the native diagnostic's input SHA-256 matches the
candidate card, and the 2×2 `SCALARMIX` rotation is complete and orthogonal.
It uses native running F, Ye, Ye2, ZEL/ZER and ZP, with
matched internal charged-scalar masses. Neutral masses are taken from
the explicitly supplied native spectrum. Alpha(0) and muon/tau total widths
are taken from the archived Module3 illustrative PDG-2024 normalization;
the tool does NOT apply experimental LFV bounds.

**Next:** audit the signed full complex mu→e gamma amplitude at the same
card, its uncertainty and omitted material contributions. Reuse existing
published and older in-house results first. Do not launch an unnecessary
precision two-loop calculation.

**Release:** SCALAR_ONE_LOOP_SUBSET_ONLY; physical_lfv=UNKNOWN;
global_zee_acceptance=null. Step 2 is started, not complete.



## General LFV engine — any parameter point

**LFV/constraints.py** accepts one explicit restricted-Zee model point in the
charged-lepton mass basis. It does NOT fit/tune Yukawas, sample a benchmark,
or infer global viability.

| Computed topology | Provenance | Gate |
| --- | --- | --- |
| h/H/A scalar one-loop, general D | Earlier Zee Module3, Lavoura (2003) | VALIDATED SUBSET |
| H1/H2 scalar one-loop, both F and D | Earlier Zee Module3 | VALIDATED SUBSET |
| h/H/A + nine charged-fermion photon Barr-Zee | Earlier Zee Module4 matching audit | VALIDATED PHOTON SUBSET |
| h/H/A + physical charged-scalar photon Barr-Zee | Earlier Module4 scalar-QED + seagull audit | CONDITIONAL ON FULL POTENTIAL |
| Neutral h/H/A tree mu→3e, tau→3l (optional) | Earlier Module3 massive crossed spinor rate | VALIDATED TREE SUBSET |
| Remaining complete gauge, LFV three-body and conversion | Explicit missing-amplitude ledger | UNKNOWN |

**Sources:** [2017 general Zee scan](https://arxiv.org/abs/1701.05345),
[2023 Zee study](https://arxiv.org/abs/2303.13383),
[Lavoura 2003](https://arxiv.org/abs/hep-ph/0302221),
[Ilisie 2015](https://arxiv.org/abs/1502.04199),
[2025/2026 general 2HDM + erratum](https://arxiv.org/abs/2410.17313);
earlier in-house [signed fermion Barr-Zee matching](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module4/MATCHING_AUDIT.md)
and [physical scalar-QED matching](https://github.com/aalida0/zee/blob/research/zee-full-constraints/module4/CHARGED_SCALAR_AUDIT.md).
Heeck–Thapa 2023 uses the neutral alignment limit and neglects second-doublet
quark couplings and additional bosonic diagrams. The complete general-2HDM
gauge terms are not automatically transferable to Zee singlet F and mixed
H1/H2. Missing terms have NOT been shown negligible.

### Explicit JSON input contract

Required: point_id; scalar parameters v,beta,alpha,gamma,mh1,mh2,ma,mH1,mH2
(GeV and radians); complex 3x3 F and D arrays with real/imag matrices;
internal_lepton_masses_GeV and pole_lepton_masses_GeV (e,mu,tau);
fermion_loop_masses_GeV (up,charm,top,down,strange,bottom,electron,muon,tau);
quark_yukawa_masses_GeV (independently, six quarks);
alpha_em, total_lepton_widths_GeV (mu,tau), mass_yukawa_scheme (explicit string).
The complex matrices use paired real/imag numeric arrays.

Optional potential: all Lambda1..Lambda5, K1, K2, K12, Mu4 from the
restricted Zee Lagrangian, same spectrum and angle convention. With it,
the physical charged-scalar/photon amplitude is calculated from trilinears;
without it this piece is UNKNOWN, NOT zero. Potential-mass consistency still
requires an independent scalar validation gate.

Neutral scalar LFV partial widths are always calculated for h/H/A to
e-mu, e-tau and mu-tau, with both charge orientations. These are
**partial widths, not branching ratios**; Higgs total widths are separate
required inputs in the Higgs analysis. Optional GF_GeV_minus2 adds the
leading quadratic singlet-F lepton-universality shifts, marked SUBSET_DIAGNOSTIC;
doublet-D/mixed operators are not included.

Optional radiative_limits: a dictionary keyed by mu_to_e_gamma,
tau_to_e_gamma or tau_to_mu_gamma. Each entry requires experimental
BR_90CL (number 0<BR<1) and source (nonempty string). The tool reports
the maximum omitted complex dipole-norm budget that would be sufficient
for compatibility, but **never assumes the omitted amplitude is zero**.
If a separately justified theoretical bound is available, provide
missing_norm_upper_GeV_minus2 (nonnegative) and bound_source. Then the
triangle inequality gives an interval for the full BR and a conditional
PASS_IF_DECLARED_BOUND_VALID, FAIL_IF_DECLARED_BOUND_VALID or UNKNOWN.
These labels are conditional mathematical implications: they are NOT a
certification that the supplied theoretical remainder bound is correct.
Without that bound, the comparison remains UNKNOWN_UNBOUNDED_MISSING.

Optional tree_quadrature_order: integer 8–48, evaluates neutral-tree mu→3e
and tau→3l partial rates with N and 2N nodes (can be slow); missing
dipole/penguin/box interference is NOT assumed zero.

**Quark scheme warning:** the quark masses in the inner propagators and the
Yukawa couplings are separate required inputs. The caller must specify
their actual mass scale, renormalization scheme, and matching procedure.
A label by itself is not a theory-validation certificate.

### Export existing native NO/IO fixtures to generic points

The adapter reconstructs the running F/Ye/Ye2 matrices, lepton rotations and CP-even scalar mixing from the original SPheno diagnostic, card and spectrum. It verifies the card SHA-256 directly; no pre-exported pilot JSON is trusted.
It requires two explicit files of six named masses:
quark_loop.json for the inner propagators and quark_yukawa.json for the
neutral scalar Yukawa couplings. Both have keys:
up,charm,top,down,strange,bottom. The source/scheme of each must be
documented; do NOT fill with arbitrary numbers.

~~~bash
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m LFV.native_point \
  --native "$WORK/step1_NO_native/result.json" \
  --card "$WORK/step1_NO_proposal.in" \
  --spectrum "$WORK/step1_NO_native/output.spc" \
  --quarks "$WORK/quark_loop.json" \
  --quark-yukawas "$WORK/quark_yukawa.json" \
  --scheme "native SPheno Q=160 GeV lepton + documented quark schemes" \
  --point-id NO_diagnostic \
  --output "$WORK/step2_NO_general_input.json"

python3 -m LFV.constraints \
  --input "$WORK/step2_NO_general_input.json" \
  --output "$WORK/step2_NO_general_amplitudes.json"
~~~

If a source-verified restricted potential is available, add the
--potential /path/to/potential.json option. Repeat for IO.
These NO/IO reference cards have neither CP-even scalar near 125 GeV;
the results are method validation only, not allowed collider benchmarks.

Every LFV subprocess returns its complex chirality components,
interference-sensitive included-subset rate and UNKNOWN experimental gate.
The full LFV gate always remains UNKNOWN pending the unvalidated
charged-current/gauge, Z-photon, additional bosonic and three-body/conversion
sectors. Never use partial subset BRs as physical LFV exclusions or passes.

### Tests

~~~bash
python3 -m unittest LFV.test_constraints LFV.test_pilot -v
~~~

Requires NumPy, SciPy, mpmath. Port validation checks source-kernel signs,
finite-quadrature behavior, potential-sourced charged trilinears, zero LFV,
native provenance and fail-closed gates. This is not independent validation
of all higher-order physics or matching-scale variation.


## Native-matched pilot results (2026-11-18)

The user executed the NO/IO LFV pilot against the **card-bound, NuFIT-matched**
SPheno reference points. The outputs remain **SCALAR_ONE_LOOP_SUBSET_ONLY**,
`physical_lfv=UNKNOWN`:

| Radiative decay | NO scalar one-loop subset BR | IO scalar one-loop subset BR |
| --- | ---: | ---: |
| mu -> e gamma | 1.23598633e-17 | 2.42199337e-17 |
| tau -> e gamma | 4.79428599e-19 | 9.56584856e-20 |
| tau -> mu gamma | 6.47172764e-17 | 1.00841541e-16 |

Native angle from NO SCALARMIX: alpha=0.7588153265256222 rad;
IO: 0.758815317317687 rad. Reported native running mixing
gamma≈−2.28665445 rad (the charged one-loop pure-F and pure-D pieces use
its squared singlet/doublet fractions). Original candidate-card SHA-256:
NO `71bb94b67a4db185a6d44fa56c96da751799211f5da23761408de75a5acbbdad`,
IO `26e705a9115cafd0d20833ac5420803f998f48013f97b98649526fb6c35b6e4d`.
Saved LFV JSON: `$WORK/step2_NO_lfv_one_loop.json` and
`$WORK/step2_IO_lfv_one_loop.json` (not committed).

**Critical Higgs-spectrum warning:** the same user-run spectra report
native masses PDG 25 = **484.284433 GeV** (NO),
**484.284454 GeV** (IO), PDG 35 ≈ **767.531 GeV** and
PDG 36 ≈ **718.038 GeV**. The source must be checked for
pole/tree identification and SPheno warnings. If these are the physical
CP-even pole masses, **neither is the observed 125 GeV Higgs**.
Therefore these are *neutrino-fit/LFV-method test cards*, NOT Higgs-viable
benchmark points. No full-scan/global acceptance is justified.
Do not repair by silently substituting historical 125 GeV values
into the native LFV calculation.

These are the coherent seven-piece scalar subset, **not** the full branching
ratios. Missing gauge/two-loop pieces are **not numerically bounded** by a
small scalar subset. The first user-run regression found one outdated test
keyword `alpha_rad`, corrected in the Step 2 branch; the LFV pilot itself
completed for both orderings. A fresh local test rerun remains to be observed.


## Experimental context (not yet a comparison)

- MEG II: B(mu→e gamma) < **1.5×10⁻¹³**, 90% C.L. using the
  2021–2022 sample ([arXiv:2504.15711](https://arxiv.org/abs/2504.15711)).
- Belle II: B(tau→mu gamma) < **9.5×10⁻⁸**, 90% C.L. with
  428 fb⁻¹ ([Belle II public result](https://docs.belle2.org/pub_data/publications/4644/)).

These limits refer to complete physical branching fractions. A number labeled
BR_scalar_one_loop_subset from this module is **not** a prediction of that
quantity: omitted amplitudes and interference may be significant. This
module records these sources only as future comparison targets.
