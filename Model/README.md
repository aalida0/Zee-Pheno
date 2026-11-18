# Step 0 — Model definition and SARAH/SPheno

**Frozen scope for the restricted low-energy model:** two Higgs doublets,
charged singlet, real CP-conserving scalar potential with λ₆=λ₇=0 at the
input scale; Type-II quark Yukawa texture; unrestricted complex Ye/Ye2 and
antisymmetric Fzee. The quark texture is **not** protected by a conventional
exact/soft Z₂ in the full model. RGE closure is not claimed.

Minimal reproducible components:
- **model.py**: tree scalar inversion and antisymmetric complex F.
- **check_sarah.wls**: SARAH CheckModel warnings/rotation audit.
- **check_native.py**, **check_basis.py**: native matrix ingestion and tree
  charged-lepton basis checks; Fortran input entries require row comments.
- **audit_neutrino_native.py**: read-only source SLHA helpers and archived
  native loop-mass availability audit.
- **diagnose_neutrino_native.py**: instrument an *isolated* generated
  SPheno copy to expose its actual one-loop Majorana self-energy and Takagi
  rotation. Uses the one-loop formula from Neutrino.fit.
- **reference/native_neutrino_snapshot.json**: original matched native card
  self-energy and parameter values, used by a permanent regression test.

See [VALIDATION.md](VALIDATION.md) for evidence and
[NEUTRINO_ROOT_CAUSE.md](NEUTRINO_ROOT_CAUSE.md) for the detailed matching
investigation. Historical one-use loop probes were removed after a
regression-grade native matrix fixture replaced them.

~~~bash
python3 -m unittest discover -s Model -p 'test_*.py' -v
~~~

**Step 0 complete for this restricted model definition only.**
A native one-loop neutrino-mass **calculator** is validated at the matched
benchmark, while neutrino oscillation *fits* belong to
[Step 1](../Neutrino/README.md). No LFV, Higgs, B-physics, vacuum-globality,
RGE or full-scan acceptance is implied.
