# Step 1 — neutrino oscillations

Minimal **pilot**, not a physical scan. **fit.py** uses the Zee one-loop
finite-mass formula checked in [Model/NEUTRINO_ROOT_CAUSE.md](../Model/NEUTRINO_ROOT_CAUSE.md)
and the saved **native running parameters** to invert an NO or IO NuFIT target.

**Experimental input:** [NuFIT 6.1 (2025), IC24 *with SK atmospheric*](https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf),
lower block, separate NO/IO best fits and one-dimensional 3-sigma intervals.
The convention is Δm²3ℓ = Δm²31 for NO and Δm²32 for IO. Selected best fits are
**targets**, not an official correlated NuFIT likelihood. Majorana phases and
the lightest neutrino mass are user inputs.

## Generate NO and IO proposals

Requires Python 3 and NumPy; native validation needs the prebuilt SPheno.
Run from the repository root:

~~~bash
python3 -m unittest discover -s Model -p 'test_*.py' -v
python3 -m unittest Neutrino.test_fit -v

WORK="$HOME/zee-lfv64-serial-20261112-094549"
CARD="$WORK/model-neutrino-loop-20261118/loop_enabled.in"

python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --reference-card "$CARD" --write-card "$WORK/step1_NO_proposal.in" \
  --output "$WORK/step1_NO_proposal.json"

python3 -m Neutrino.fit --ordering IO --mlight-ev 0.001 \
  --reference-card "$CARD" --write-card "$WORK/step1_IO_proposal.in" \
  --output "$WORK/step1_IO_proposal.json"
~~~

The inverse uses the Zee null-cone condition and a rank-five linear
least-squares solution. This takes the **minimum-norm** member of the
four-complex-parameter inverse family; it does not sample that family.
All outputs remain **ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN** until native
revalidation. Files use exclusive creation: use new filenames when rerunning.

The proposal uses the *reference* running masses and mixing in the saved
[one-card native snapshot](../Model/reference/native_neutrino_snapshot.json).
It changes **only** the four SLHA real/imag blocks for Fzee and Ye2. SPheno
reconstructs Ye, running leptons and scalar inputs during matching; the first
proposal **must not** be assumed to pass oscillation constraints natively.

## Check native masses and PMNS

Run instrumented SPheno in a fresh directory without changing the original
installation:

~~~bash
python3 -m Model.diagnose_neutrino_native \
  --spheno "$WORK/SPheno-4.0.5" \
  --card "$WORK/step1_NO_proposal.in" \
  --output-dir "$WORK/step1_NO_native"

python3 -m Neutrino.fit --ordering NO --mlight-ev 0.0093 \
  --native-result "$WORK/step1_NO_native/result.json" \
  --output "$WORK/step1_NO_native_screen.json"
~~~

Repeat for IO with new output directories. The screen uses the *one-loop*
native Takagi rotation and charged-lepton rotation to extract mixing angles,
δCP and mass splittings, checking published one-dimensional 3-sigma ranges.
**NATIVE_3SIGMA_SCREEN_PASS** is not a correlated χ² or full Zee acceptance.
A failed first proposal requires retuning in the native matching scheme.

**STOP before LFV/scans:** certify both NO and IO under native matching and
assess stability under reasonable matching-scale variations.
