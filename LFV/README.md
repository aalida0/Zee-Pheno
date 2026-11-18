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
Supply the **physical CP-even scalar mixing angle alpha, radians**, from
the *same* SPheno scalar convention; the historical value 1.398 is NOT a
safe default for this new benchmark. Mixing and pole masses are an
explicitly mixed-order leading-order scheme diagnostic.

~~~bash
# Substitute the real angle for ALPHA_RAD; do not use the literal token.
python3 -m LFV.pilot \
  --native "$WORK/step1_NO_native/result.json" \
  --card "$WORK/step1_NO_proposal.in" \
  --spectrum "$WORK/step1_NO_native/output.spc" \
  --alpha-rad ALPHA_RAD \
  --output "$WORK/step2_NO_lfv_one_loop.json"

python3 -m LFV.pilot \
  --native "$WORK/step1_IO_native/result.json" \
  --card "$WORK/step1_IO_proposal.in" \
  --spectrum "$WORK/step1_IO_native/output.spc" \
  --alpha-rad ALPHA_RAD \
  --output "$WORK/step2_IO_lfv_one_loop.json"
~~~

Find mixing blocks and neutral-scalar masses with:

~~~bash
grep -niE '^[[:space:]]*Block[[:space:]]+(MASS|SCALARMIX|ZH|ALPHA)' \
  "$WORK/step1_NO_native/output.spc"
~~~

The tool fails unless the native diagnostic's input SHA-256 matches the
candidate card. It uses native running F, Ye, Ye2, ZEL/ZER and ZP, with
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


## Experimental context (not yet a comparison)

- MEG II: B(mu→e gamma) < **1.5×10⁻¹³**, 90% C.L. using the
  2021–2022 sample ([arXiv:2504.15711](https://arxiv.org/abs/2504.15711)).
- Belle II: B(tau→mu gamma) < **9.5×10⁻⁸**, 90% C.L. with
  428 fb⁻¹ ([Belle II public result](https://docs.belle2.org/pub_data/publications/4644/)).

These limits refer to complete physical branching fractions. A number labeled
BR_scalar_one_loop_subset from this module is **not** a prediction of that
quantity: omitted amplitudes and interference may be significant. This
module records these sources only as future comparison targets.
