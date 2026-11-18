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
- Sample/scoring machinery is **not yet implemented**. There is no claim
  of a physical allowed region, exclusion, or cascade sensitivity.
