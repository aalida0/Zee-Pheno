# Ponytail code audit and cleanup — 18 November 2026

**Scope:** the entire ~60-file `Zee-Pheno` branch tree, with detailed
inspection of the scan generator, inverse, native input writer, CI and
their callers. Expected load: one researcher generating thousands of
independent SLHA proposal cards, then evaluating separate native runs
on a batch system. No web server or shared mutable production database.

## Must fix (completed in this branch)

1. **Fixed-scalar-only candidates:** `Scan/proposal_catalog.py` used the
   same native reference spectrum for all targets; a passing neutrino
   point could not test the cascade region with a 125 GeV Higgs.
   **Fix:** `Scan/physical_proposals.py` creates independent physical
   scalar proposals and reconstructs Lagrangian input blocks from each.
   Every native result must be recomputed from its own generated card.
2. **Incomplete inverse family:** `Neutrino.fit.inverse` used only the
   smallest-magnitude root, one projective chart, and minimum-norm D.
   **Fix:** expose both roots, alternate charts, complex F phase and all
   four complex homogeneous D directions, preserving the original
   default for single-target backward compatibility.
3. **Unsafe scalar-card generation:** the original row writer only
   supported MINPAR and four Yukawa matrix blocks.
   **Fix:** extend the *same* `Model.check_native.fix_card` path for
   the SARAH `ZEEPOT` and `ZEEQUARTIC` codes. Preflight verifies
   all real/complex entries **before creating any files**; missing
   blocks and changed source-card hashes terminate cleanly.

## Lean changes

- **Deleted:** `Scan/proposal_catalog.py` (fixed-spectrum-only CLI)
  and `Scan/test_proposal_catalog.py` (tests for that obsolete CLI).
  Their purpose is superseded by the scalar-first generator; the
  independent `Neutrino.fit` CLI still provides a single-target
  neutrino proposal.
- **Removed:** corresponding obsolete CI test references; replaced
  repeated `Scan/README.md` instructions with one current workflow.
- **Kept intentionally:** `Scan/prepare_reference_manifest.py`, because
  the user actually runs it for NO/IO native validation.
  `Model.diagnose_neutrino_native.py`, `Model.check_native.py`,
  `Model/check_sarah.wls`, the reference snapshot and its tests:
  removing them would break reproducibility or source attestation.
- **Kept intentionally:** `LFV/strict_2hdm_reference.py`,
  `LFV/approximation_study.py`, `LFV/nonzero_mixing_study.py`,
  `LFV/doublet_proxy_study.py` and their focused tests. These look
  exploratory but reproduce **published source-level approximation
  comparisons** and substantiate that the omitted electroweak amplitudes
  can be material. They are not scan production code, but are not dead
  research evidence.

## Review limits and physics caveats

- Scalar BFB is only necessary doublet + sufficient quartic *screens*.
  The code does **not** certify global vacuum stability, one-loop
  EWPO, perturbative scattering unitarity, Higgs rates, flavor or full
  Zee LFV.
- The scalar + neutrino proposal uses archived native charged-lepton
  running masses/rotations. A newly generated card is analytically
  closed at those inputs only; native SPheno spectrum/matching may move
  the point, and the recorded `FZEEIN` propagation problem must be
  rechecked on the executing binary.
- CI uses a **synthetic SLHA template** with the correct SARAH block
  contract, not the user's original loop-enabled card or actual native
  SPheno. Production execution and full numeric real-card validation
  must take place on the user's machine.
- This branch has **not** reviewed every formula in the LFV two-loop
  implementation or every prior historical note, only the critical
  scan and its source-linked dependencies.

**Verdict:** ship the audited proposal generator for a native pilot,
not a full physics exclusion scan. Follow-up work should focus on
actual SPheno card runs and physical gates, not more speculative
software layers.
