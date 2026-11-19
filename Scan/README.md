# Scan — scalar-first candidates, not an allowed-region claim

**Physics question:** Where can the previously proposed charged cascade
`H2± → H1± Z` survive the restricted Type-II-quark Zee model constraints?
This repository currently **generates candidates and validates available
evidence**. It does not establish a globally viable parameter region.

## 1. Generate scalar and neutrino candidates

The recommended generator is **[physical_proposals.py](physical_proposals.py)**.
It supersedes the old fixed-scalar `proposal_catalog.py`, which was removed.
For an isolated NO/IO target at a *fixed* spectrum, the existing
`python3 -m Neutrino.fit` CLI remains available.

The [versioned example configuration](physical_pilot_config.json) declares
30 independent scalar proposals, three predeclared arms weighted 3:5:2
(sub-TeV, 1000–1350 GeV, broad control), and four neutrino families
per scalar. All distributions and sampling ranges are explicit;
there is **no LFV-driven adjustment, rejection/resampling or optimization**.

At every independent scalar attempt:

1. Draw physical scalar masses, `tanbeta`, charged mixing `gamma`,
   `beta-alpha`, soft-mass ratio and quartics from the declared arm.
   **Tree-level `mh1=125.25 GeV` is a requested input**, not a native
   mass or Higgs-rate PASS. Reject a closed cascade as a declared
   proposal-prior screen, but keep that attempt in the denominator.
2. Reconstruct the restricted scalar potential using
   [Model/model.py](../Model/model.py): `ZEEPOT` codes 1–5 map
   `M112,M222,M12Sq,MCS2,Mu4`; `ZEEQUARTIC` codes 1–9
   map `Lambda1..5,LambdaC,K1,K2,K12`, matching the
   [original SARAH model](https://github.com/aalida0/zee/blob/research/zee-full-constraints/models/SARAH/ZeeTypeII/parameters.m).
   The SHA-attested neutrino reference may omit the two scalar
   blocks; the generator adds complete placeholder rows to the **derived
   output cards only**, then replaces their values with the sampled
   potential. Existing partial scalar blocks are rejected. `MINPAR[3]`
   and all four complex Yukawa input matrices remain mandatory.
3. Record necessary doublet BFB and a deliberately conservative
   **sufficient** Zee *quartic* BFB condition, plus `|quartic|<4pi`
   as a conventional analysis prior. These are **not** a proof of
   tree global-minimum stability, unitarity, EWPO or scalar rate viability.
4. At preselected scalar points, build separate NO/IO targets using the
   published NuFIT 6.1 central oscillation coordinates. Sample
   lightest mass, Majorana phases, the complex Zee null-cone orientation,
   both algebraic roots, several projective charts, a declared `F`
   normalization/phase, and **four independent complex homogeneous
   `D` coefficients**. Use the *current scalar masses and gamma* in
   the one-loop inverse. Check the complex matrix closure each time.
5. Rotate `F,D` into the archived charged-lepton reference gauge basis.
   Construct the proposed `Ye2` using the *reference* running charged
   lepton mass-direction Yukawa. Update the full scalar SLHA blocks,
   `MINPAR[3]`, `FZEEIN/IMFZEEIN`, `YEE2IN/IMYEE2IN`
   through the existing validated card editor.

The starting `MFe/ZEL/ZER/Ye/Ye2` correspond to **one archived
running reference**. They are *proposal coordinates*, not native
matching at new scalar points. SPheno can shift all masses, mixings,
lepton rotations and effective couplings; **every generated card must
be rerun and matched independently**. Native FZEEIN propagation must
also be checked on the actual binary. No candidate is automatically
called neutrino-compatible or LFV-allowed.

Run from the repository root:

~~~bash
cd ~/Zee-Pheno
git fetch origin
git switch research/scalar-neutrino-scan-ponytail
git pull --ff-only

WORK="$HOME/zee-lfv64-serial-20261112-094549"
TEMPLATE="$WORK/model-neutrino-loop-20261118/loop_enabled.in"

python3 -m Scan.physical_proposals \
  --config Scan/physical_pilot_config.json \
  --template "$TEMPLATE" \
  --snapshot Model/reference/native_neutrino_snapshot.json \
  --output-dir "$WORK/scalar_neutrino_pilot_01"
~~~

The template SHA-256 is pinned to the original loop-enabled source
card; that file is never rewritten. Missing `ZEEPOT`/`ZEEQUARTIC`
blocks are synthesized **after** source attestation for derived cards.
Missing Yukawa or tan-beta inputs, partial scalar blocks, and any altered
source hash are refused. All output paths are exclusive: choose a new
output directory for a rerun. Newly added scalar blocks still require
native SPheno input-propagation validation before physics claims.

Read `proposals.json` **before launching SPheno**:
`n_scalar_attempted` is the denominator for scalar sampling;
`n_scalar_preselected` is the number passing the declared cheap
proposal screens; `n_native_cards_written` counts conditional
Yukawa families, **not distinct scalar points**. Every card is
`NATIVE_REVALIDATION_REQUIRED`, and `n_globally_allowed=null`.

## 2. Preflight the native scalar spectrum before instrumenting

The original SPheno executable may print a fatal `CalculateMhh` error
and exit with code zero **without producing any spectrum**. Run every
candidate card through the *original, SHA-pinned executable* first:

~~~bash
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Scan.spectrum_preflight \\
  --proposals "$WORK/scalar_neutrino_pilot_01/proposals.json" \\
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \\
  --output-dir "$WORK/scalar_spectrum_preflight_01"
~~~

Review `scalar_spectrum_preflight_01/preflight.json`. Each case has a
retained `run.log` and optional `spectrum.spc`. `SPECTRUM_PRESENT_UNVALIDATED`
means only that SPheno produced a scalar MASS block with positive neutral
and charged scalar entries. No neutrino, LFV, Higgs or vacuum PASS is implied.
If all four scalar candidates fail the CP-even mass calculation, investigate
native matching/potential inputs before drawing a new scan.

## 3. Validate native outputs on the matching source card


For a *single* generated card, use the existing instrumented SPheno
runner on your machine, once the executable's source/binary attestation
and `FZEEIN` perturbation test have been verified:

~~~bash
CARD=$(find "$WORK/scalar_neutrino_pilot_01/cards" -name '*.in' -print -quit)
test -n "$CARD"
python3 -m Model.diagnose_neutrino_native \
  --spheno "$WORK/SPheno-4.0.5" \
  --card "$CARD" \
  --output-dir "$WORK/scalar_native_test_01"
~~~

**Do not** reuse the original NO/IO native results for newly
generated scalar or neutrino cards. Use one distinct native output per
card and verify its SHA-256. The instrumented SPheno setup is local;
GitHub CI runs *synthetic* numerical/software tests only.

## 4. Report separate physics gates

[Scan/native_pipeline.py](native_pipeline.py) takes an explicit
`schema_version=1` manifest listing actual `card`, `native`,
`spectrum`, `point_id` and `ordering` for each native run.
Relative paths are relative to that manifest. Example:

~~~json
{
  "schema_version": 1,
  "points": [{
    "point_id": "scan_00001_NO",
    "ordering": "NO",
    "card": "/absolute/path/to/generated_card.in",
    "native": "/absolute/path/to/result.json",
    "spectrum": "/absolute/path/to/output.spc"
  }]
}
~~~

~~~bash
python3 -m Scan.native_pipeline \
  --manifest /path/to/actual_native_manifest.json \
  --output-dir /path/to/new_gate_results
~~~

This produces `points.json` and `points.csv` with exact source
hashes and PASS/FAIL/UNKNOWN gates. A wrong Higgs mass or failed
individual NuFIT 3-sigma screen is an explicit necessary-screen FAIL;
a missing input or mismatched card is **NOT_EVALUATED_SOURCE_ERROR**.
A passing set of necessary screens means only
`UNKNOWN_PENDING_GATES`, never `globally_allowed=true`.

Optional quark propagator and Yukawa running masses **must be provided
together with an explicit matching scheme** to enable known photon
Barr–Zee LFV terms. Without them, the pipeline computes the
source-bound scalar one-loop subset only. An optional native
`FZEEIN/Ye2` perturbation report must match the original card SHA
and executed binary before that propagation gate is PASS.

Optional `potential`, `potential_source`, `scalar_inversion` can
record a tree-consistent scalar reconstruction and enable the charged
scalar-photon partial LFV term. It does *not* prove that the input
potential matches the **running native** Lagrangian at a different
renormalization scale. Avoid interpreting a tree-to-pole mismatch as
a physical scalar exclusion without a matched scheme.

**Still unknown:** full Zee electroweak LFV two-loop amplitudes,
`mu→3e`/conversion, correlated oscillation likelihood, full Higgs
rates/searches, B flavor, EWPO, tree-unitarity and global charged
vacuum, and collider signal sensitivity. Approximation-labeled
LFV overlays are permitted; **full Zee exclusion/allowed contours
are not**.

The original source-matched NO and IO controls both pass their six
NuFIT 6.1 one-dimensional screens but fail the necessary h125 mass
screen (native PDG25 ~484 GeV). These are *regressions*, not surviving
charged-cascade benchmarks.

## Tests and code retention

~~~bash
python3 -m unittest Scan.test_physical_proposals Scan.test_native_pipeline -v
python3 -m unittest Neutrino.test_fit -v
~~~

The old `Scan.proposal_catalog` and its duplicate test were removed.
The older LFV/2HDM source comparison scripts are **retained** because
they are executable validation evidence, not disposable demos. The
native instrumentation, source-propagation checks and exact original
reference fixture are also intentionally retained.
