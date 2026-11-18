# Scan — governing physics policy

**Goal:** identify regions of the restricted Zee model where the charged
cascade H₂± → H₁± Z survives independently justified constraints.
A successful neutrino-fit benchmark is **not** an accepted parameter point.

- Separate **candidate generation**, **observable computation**, and
  **experimental/theoretical gates**. Each gate returns PASS, FAIL or UNKNOWN
  with inputs, source/version, approximation, numerical tolerance and hash.
- **Do not retune scalar masses, angles or Yukawa textures to rescue a
  failing card**. A failed physical condition rejects that card. Uncomputed
  mandatory constraints remain UNKNOWN, never silently treated as passes.
- The exact Zee inverse may be used to **parameterize the neutrino-compatible
  manifold** (a physics condition), but the first minimum-norm D solution
  and a single null-cone root are *not* the whole solution space. A future
  scan must state the freedom in F, the four-complex-dimensional homogeneous
  D family, lightest mass and Majorana phases, and sample or restrict it
  transparently. Do not advertise such a restricted slice as a global result.
- Compute LFV with source-matched full complex amplitudes and coherent
  interference. The seven-piece scalar one-loop subset in LFV/ is a
  diagnostic only, not a pass/fail constraint until omitted contributions
  have been justified, bounded or included.
- Enforce the 125-GeV observed-Higgs requirement, masses/mixing, vacuum,
  perturbativity, EWPO, flavor and direct searches for **every** proposed
  physical point; keep theory-scheme uncertainty explicit. In particular
  the currently tested NO/IO cards have native PDG25 ≈484 GeV and PDG35
  ≈768 GeV, so they are not Higgs-viable benchmarks.
- A **diagnostic/exploratory** scan MAY overlay the source-matched LO +
  selected photon Barr-Zee partial LFV sensitivity, flagged
  `LFV_APPROXIMATION_DEPENDENT`; it must not call those points
  LFV-allowed. The [materiality decision](../LFV/REMAINDER_AUDIT.md)
  shows why no uniform source-qualified missing-amplitude bound exists
  for general Zee points. Physical combined LFV gates remain `UNKNOWN`.
- A deterministic 18-point **strict 2HDM** benchmark-independent control
  sweep [LFV/APPROXIMATION_STUDY.md](../LFV/APPROXIMATION_STUDY.md) finds
  omitted electroweak contributions larger than 10% of the retained
  one-loop/photon amplitude in **18/18 points for each of three channels**.
  This is a source-matched falsification of universal suppression,
  not a Zee two-loop result or an LFV experimental exclusion.
- A [24-point nonzero-mixing Zee LFV study](../LFV/NONZERO_MIXING_APPROXIMATION.md)
  now evaluates the **actual retained** one-loop and photon Barr–Zee
  amplitudes at gamma ≠ 0 and F ≠ 0 with an open H2±→H1± Z cascade.
  Its separate exact-doublet **auxiliary 2HDM control** finds omitted
  electroweak terms above 10% in 24/24 corresponding control points
  for each LFV channel. These numbers are **not transferable Zee
  gauge amplitudes**, and none of the toy inputs is neutrino/constraint
  certified. All physical LFV gates remain UNKNOWN.
- Sample/scoring machinery is **not yet implemented**. There is no claim
  of a physical allowed region, exclusion, or cascade sensitivity.


## Run the existing NO/IO native reference checks now

The references can be screened **without any quark mass dictionaries**.
In this mode the pipeline runs the source-matched scalar one-loop LFV
pilot, not the photon Barr-Zee calculation. Do not invent QCD schemes.

~~~bash
cd ~/Zee-Pheno
git fetch origin
git switch research/scan-native-gated-pipeline
git pull --ff-only

WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Scan.prepare_reference_manifest \
  --work "$WORK" \
  --output "$WORK/native_scan_manifest_20261118.json"

python3 -m Scan.native_pipeline \
  --manifest "$WORK/native_scan_manifest_20261118.json" \
  --output-dir "$WORK/native_scan_gates_20261118"
~~~

Both outputs use exclusive file/directory creation; choose fresh
filenames on rerun. The checked source spectra previously reported
PDG25 approximately 484 GeV, so the expected necessary h125 mass
screen is FAIL for both. The six individual neutrino 3sigma results
should remain PASS if the same native snapshots/cards are supplied.
The full Zee LFV gate will remain UNKNOWN.

When separately checked quark running/loop masses are available,
supply **all three** optional manifest fields quark_loop,
quark_yukawa, mass_yukawa_scheme. The pipeline will then evaluate
the additional fermion/photon Barr-Zee pieces. Without them, none
of the scalar, neutrino or one-loop LFV checks is blocked.

## Implemented evidence bridge: neutrino proposals → native SPheno → LFV

This repository now contains a **versioned point workflow**, not a
physical allowed-region scan:

- [proposal_catalog.py](proposal_catalog.py) generates a declared NO/IO
  lightest-mass/Majorana-phase target slice as separate SLHA cards using
  the attested reference-card and Neutrino.fit. **Every card is pending
  independent native SPheno revalidation.** The minimum-norm/single-root
  inverse does not sample the full four-complex-dimensional D freedom.
- [native_pipeline.py](native_pipeline.py) consumes a list of already
  executed native SPheno cards/spectra. It verifies the original card
  SHA-256, recomputes the six individual NuFIT 6.1 3sigma tests through
  the native one-loop Takagi output, then evaluates source-bound full
  complex one-loop and selected photon Barr-Zee LFV amplitudes.
- [native_gates.py](native_gates.py) records independent PASS/FAIL/UNKNOWN
  records. It tests the loose *necessary-only* h125 mass screen
  (PDG25 within 3 GeV of 125.25) and the kinematic on-shell charged
  cascade. All B-flavor, Higgs signal-rate, EWPO, direct-search,
  global-vacuum and *full Zee LFV* gates remain UNKNOWN.
- [theory_gates.py](theory_gates.py) optionally compares a declared
  scalar potential with the mass/angle inverse, plus **necessary**
  CP-conserving 2HDM doublet-subpotential BFB tests and a conventional
  magnitude of quartics below 4pi. These do not certify full Zee
  vacuum stability, native potential matching, or unitarity.

The SPheno charged angle may differ by integer pi from the principal
gamma used in scalar inversion. The bridge treats these as identical
charged-field sign conventions **only for the tree potential identity**.
The native gamma and spectra are never modified. This is not a retune.

### A. Declared analytic candidate slice (not a native pass)

Create a JSON file with schema_version=1, snapshot path, reference_card
path and a targets array. One example target is:

~~~json
{
  "point_id": "NO_mass0093",
  "ordering": "NO",
  "mlight_eV": 0.0093,
  "alpha21_rad": 0.0,
  "alpha31_rad": 0.0
}
~~~

Use Neutrino's pinned loop-enabled source card:

~~~bash
python3 -m Scan.proposal_catalog \
  --catalog /path/to/targets.json \
  --output-dir /path/to/new_proposal_folder
~~~

Each output is clearly marked NATIVE_REVALIDATION_REQUIRED. A native
SPheno run must actually be executed for every candidate before physical
screens may run. New scalar points cannot inherit old neutrino passes.

### B. Per-point native + partial LFV manifest

A manifest has schema_version 1 and a points array. Example:

~~~json
{
  "schema_version": 1,
  "points": [{
    "point_id": "NO_reference",
    "ordering": "NO",
    "card": "/home/amine/zee-lfv64-serial-20261112-094549/step1_NO_proposal.in",
    "native": "/home/amine/zee-lfv64-serial-20261112-094549/step1_NO_native/result.json",
    "spectrum": "/home/amine/zee-lfv64-serial-20261112-094549/step1_NO_native/output.spc",
    "quark_loop": "/path/to/author-scheme-quark-loop-masses.json",
    "quark_yukawa": "/path/to/scale-matched-quark-yukawas.json",
    "mass_yukawa_scheme": "DOCUMENTED_SCALE_AND_RENORMALIZATION_SCHEME"
  }]
}
~~~

The two quark dictionaries must each provide positive masses (GeV)
for up, charm, top, down, strange and bottom. A string describing a
scheme is not by itself a demonstrated running/Yukawa matching. Do
not use arbitrary test values as physical inputs.

~~~bash
python3 -m Scan.native_pipeline \
  --manifest /path/to/native_manifest.json \
  --output-dir /path/to/new_results_directory
~~~

The runner creates a machine-readable points.json and points.csv with
full source SHA-256 and partial LFV rates. It never overwrites an
existing output folder. A tampered source is isolated as an
unevaluated provenance error, **not** an excluded physics point.

Optional fields are native_input_report (a native-executed, same-card
Model.check_native FZEEIN/Ye2 perturbation report), plus the three
joint fields potential, potential_source and scalar_inversion.
The inversion JSON requires m12sq_GeV2, lambda_c, K1, K2, K12.
The potential JSON requires Lambda1..5, K1, K2, K12, Mu4.
Source-declared scalar consistency is **not** an independent
confirmation of the input potential in SPheno. If omitted, the
charged-scalar photon contribution cannot be included.

### Interpretation and current blockers

The known original NO/IO native reference cards pass the individual
NuFIT screens, but have PDG25 around 484 GeV. They must **FAIL**
the restricted model's observed-h125 mass prerequisite. Passing
neutrino targets is not a globally viable scalar/Higgs point.

A model point with all computed necessary screens passing remains
UNKNOWN_PENDING_GATES, not allowed. Full LFV, conversion, flavor,
Higgs signal rates, EWPO, unitarity and the Zee global vacuum are
not yet certified. The old FZEEIN native-input propagation issue is
UNKNOWN until the card-matched perturbation gate is passed.

Run the regression suite:

~~~bash
python3 -m unittest Scan.test_native_pipeline Scan.test_proposal_catalog -v
~~~

All test-native data are synthetic one-loop Takagi control matrices,
**not new native SPheno production**.
