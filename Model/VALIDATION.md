# Model validation record — 2026-11-18

Scope: **restricted, CP-conserving Zee scalar potential at tree level**, Type-II quark Yukawas, general charged-lepton Yukawas and a complex antisymmetric singlet coupling `f=-f^T`. This is **not** a globally accepted Zee benchmark.

## Executed checks

- `python3 -m unittest Model.test_model -v`: **4/4 passed**, Python 3.13.5. Tests cover exact antisymmetry, original scalar spectrum, tadpoles, full charged Goldstone/invariants, and signed charged mixing.
- Independent potential audit (SymPy differentiation of the 10 real scalar degrees of freedom, NumPy diagonalization): benchmark `mh=125, mH=580.7, mA=633.7, mH1+=300, mH2+=700 GeV`, for `gamma=+1.089` and `gamma=-1.089`. Three massless Goldstones; largest |Goldstone squared mass| `6.8e-11 GeV^2`, maximum physical squared-mass relative residual `2.8e-15`, maximum tree-level tadpole `6.3e-10 GeV^3`.
- Independent deterministic stress test of 2,000 unconstrained parameter points (seed `20261118`), covering both charged-mixing signs: maximum relative mass-squared reconstruction residual `1.24e-13`; maximum normalized tadpole residual `6.23e-16`; maximum absolute charged-Goldstone mass-squared residual `1.93e-9 GeV^2`. **These are algebraic points, not phenomenologically accepted scan points.**

The 10x10 symbolic audit and stress test were executed independently of the GitHub production model tests; only the four small standard-library regression tests are retained in this minimal repository.

## Unvalidated / do not promote

- **SARAH/SPheno canonical model:** not established. The older repository's `models/SARAH/ZeeTypeII/ZeeTypeII.m` and the generated-model snapshot at `/home/amine/zee-conditional-cascade-scan/validation/module14/local/models/ZeeTypeII/` differ in their neutrino rotation definitions. Do not silently select either snapshot.
- **Native input propagation:** one corrected native SPheno benchmark and two orthogonal probes now report PASS on the user's machine (details below). The original `.spc` and log artifacts have not been independently inspected here.
- **Radiative neutrino masses, LFV, EDM, Yukawa basis transformations, loop normalization and global vacuum:** not validated by these tree-level tests.
- **BFB/unitarity/EWPO/Higgs/flavor/collider constraints:** no acceptance or exclusion claim from this module.

## Required native-model acceptance gate

On the user's machine with SARAH 4.15.4 and existing compiled SPheno 4.0.5:

1. Reconcile the two historical SARAH model snapshots, preserving neutrino conventions; run `CheckModel[]` and independently inspect `Fzee` antisymmetry and charged scalar vertices.
2. On one corrected four-field SLHA benchmark, verify independent changes to `FZEEIN` and `YEE2IN` affect only their intended output parameters, reconstruct `f+f^T=0`, and match the scalar spectrum and physical Yukawa basis.
3. Validate the radiative neutrino mass and source-matched LFV amplitudes separately before promoting any generated model, code, or global exclusions.

Result: **TREE_SCALAR_VALIDATED; NATIVE_SARAH_SPHENO_PENDING**.

## 2026-11-18 user-reported native SPheno input-gate result

**NATIVE_INPUT_GATE_PASS** from the user's Ubuntu machine using `python3 -m Model.check_native`. This is user-provided stdout, not an independently inspected local `.spc` artifact. No new broad physics acceptance is inferred.

- Input SHA-256: `a969dfae8428c352bfb6f70178bfb971c3357a1dadac83de01be52981e88513e`
- Executable SHA-256: `2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52`
- Baseline `FZEE` relative mismatch: `7.760162818345134e-10`.
- Baseline `YEE2` relative mismatch: `4.69220374237888e-10`.
- `max |F + F^T|`: **0.0**.
- Antisymmetric real `F` probe: input perturbation `2.2068930208234563e-06`; observed `F` change `2.2068930000000017e-06`; observed `Ye2` change `0.0`.
- Imaginary `Ye2(2,2)` probe: input perturbation `1.3071601895116883e-05`; observed `Ye2` change `1.3071602000000004e-05`; observed `F` change `0.0`.
- Three output hashes: baseline `3a7a4fab80c9e1837e72ea4a4e0ab70ea222367d97606b879d9d7115bea1c4d8`; F probe `f8f7eb70b76955caea149d81d33c5daaa6a1541a17b4a6328f0f89bc5a803ad0`; Ye2 probe `2443ccc076740132638fb7f41400c9f9cc3f73c806fcd49dfbc34244073f1533`.

**Decision:** Corrected SLHA/Yukawa native input propagation passes for this card and attested binary. Previous reported F diagonal/symmetry violations are explained by the malformed input rows, not shown to be a SARAH interaction defect. SARAH `CheckModel[]` now completes with warnings on the rotation-bearing snapshot (see below). Physical Yukawa basis `D`, scalar spectrum in this native run, radiative neutrino-mass normalization, and complete LFV remain **PENDING/UNKNOWN**.

## 2026-11-18 SARAH/SPheno validation preparation

### Inspectable evidence

- The older repository's `models/SARAH/ZeeTypeII/parameters.m` declares `Fzee` complex and `Form -> AntiSymmetric`. Its `ZeeTypeII.m` uses doublet Yukawas `Ye, Ye2`, quark `Yd, Yu`, and `Fzee/2 l.l.Splus`. This is **static source evidence only**.
- The November 14 native source audit attributed the anomalous output Yukawas to three-field rows of `FZEEIN/IMFZEEIN/YEE2IN/IMYEE2IN` being read with the generated Fortran's four-field `Read(io,*) i1,i2,wert,read_line`. A separate local gfortran reader check confirms `#` is accepted as the fourth character token.
- `Model/check_native.py` generates a corrected baseline and independent antisymmetric-real-`F` and imaginary-`Ye2` probes, keeps the original benchmark untouched, verifies the executable against the previously recorded forced-build SHA, and writes the diagnostics in a new directory. `Model/test_native.py` tests serialization, antisymmetry and a **mock** executable. Mock success is not physical validation.
- `Model/check_sarah.wls` runs the local SARAH `Start["ZeeTypeII"]` and `CheckModel[]`. Its output **requires human review**; no automatic PASS is inferred from the exit status.

### Execute the native gate on the user's machine

From `~/Zee-Pheno` after pulling the current `main`:

```bash
python3 -m unittest Model.test_model Model.test_native -v

# The two historical model snapshots are not identical; compare first.
diff -u \
  ~/zee/models/SARAH/ZeeTypeII/ZeeTypeII.m \
  ~/zee-conditional-cascade-scan/validation/module14/local/models/ZeeTypeII/ZeeTypeII.m || true

# Supply the real SARAH installation's SARAH.m path.
SARAH_M="$(find "$HOME" -name SARAH.m -type f -print -quit)"
wolframscript -file Model/check_sarah.wls "$SARAH_M" \
  "$HOME/zee-conditional-cascade-scan/validation/module14/local/models" \
  2>&1 | tee /tmp/zee_sarah_checkmodel.log

WORK="$HOME/zee-lfv64-serial-20261112-094549"
CARD="$WORK/EDM3ZERO_SOURCE64_20261114_172527/native-source-64-results/Zee_2264_NO_20261111_newEDM3zero/input.in"
python3 -m Model.check_native \
  --input "$CARD" \
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
  --attestation "$WORK/model-rebuild-attestation.json" \
  --output-dir "$WORK/model-validation-native-20261118"
```

The output directory must not exist beforehand. Inspect `native-input-report.json` and the three `.spc`/`.log` files. If `SARAH.m` cannot be found, provide its explicit path rather than using a substitute version.

### Interpretation and release gate

- `NATIVE_INPUT_GATE_PASS` means the **specific one-card input/output propagation checks** passed; other loop-related correctness remains UNKNOWN.
- `NATIVE_INPUT_GATE_FAIL` or an exception requires investigating SPheno reading, RG/output scale, source linkage, or missing matrices before proceeding. The diagnostic tolerances (`5e-4` relative input matching and `5e-9` absolute antisymmetry) may be too restrictive if output and input are evaluated at different scales. Do **not** loosen them simply to obtain a pass.
- `F_probe_unintended_Ye2_delta` and `Ye2_probe_induced_F_delta` are measurements, not zero-coupling requirements: RGE mixing can be physical.
- Keep **native model validation UNKNOWN** until the SARAH model-snapshot/rotation discrepancy, `CheckModel[]` warnings, source-to-binary provenance and one-card native results are reviewed. Then check the charged vertices, charged-lepton basis `D` and radiative neutrino-mass normalization independently.

**Current release decision: TREE_SCALAR_VALIDATED / NATIVE_SARAH_SPHENO_PENDING / PHYSICAL_LFV_UNKNOWN.**

## 2026-11-18 user-reported SARAH 4.15.4 `CheckModel[]` outcome

Evidence: terminal output provided by the user, covering an actual `Start["ZeeTypeII"]; CheckModel[];` on `/home/amine/zee-conditional-cascade-scan/validation/module14/local/models`. SARAH reports `All Done. ZeeTypeII is ready!` and `All checks finished in 4.366489...` seconds. **Initialization and structural checks completed; this is not a warning-free or loop-physics validation.**

**Snapshot diff (exactly two EWSB changes):**
- Generated model adds `{{vL},{VL,UV}},` to `DEFINITION[EWSB][MatterSector]`.
- Generated model uses `Fv -> {VL,conj[VL]}` instead of `Fv -> {vL,conj[vL]}`.

The rotation-bearing, native-generation snapshot is the **provisional canonical candidate**, because it matches the successful SARAH run and historical native SPheno workflow. Do not copy into this minimal repo or assert physical PMNS matching until the radiative neutrino-mass basis is checked. At tree level the `{{vL},{vL}}` mass matrix is identically zero; `UV` is therefore not a physically determined PMNS rotation.

**CheckModel warning review:**
- `MassMatrix::OnlyZero` for `{{vL},{vL}}`: EXPECTED (Zee radiative neutrinos).
- `PossibleTerms::NonSUSY`: SARAH lists nonincluded gauge-compatible operators, including opposite-doublet quark Yukawas, `lambda6/lambda7`-type terms, singlet-related terms and `{d,conj[Splus],conj[u]}`. **OPEN:** classify each warning as forbidden by an explicitly stated symmetry, identically zero after gauge/statistics/Lorentz contraction, or an intentionally omitted operator. SARAH itself warns this checker is imperfect for nonsupersymmetric models.
- EWSB rotations `vL->VL`, `phi1/phi2->hh`, `sigma1/sigma2->Ah`, three charged scalars `->Hm`, and quark/charged-lepton rotations: the reported quantum-number compatibility checks complete without an error.

**Release gate:** `SARAH_CHECKMODEL_COMPLETED_WARNINGS_OPEN`; `NATIVE_INPUT_GATE_PASS` remains valid. Next verify physical charged vertices, `UV`/charged-lepton basis matching, full scalar spectrum and radiative neutrino loop normalization using the generated source and corrected benchmark. The restricted potential and Type-II quark Yukawas need a clearly justified scope/symmetry statement; no RG closure or full-model certification is claimed.



## 2026-11-18 user-provided generated Fortran vertex/mass audit

Source displayed by the user from the existing model-specific SPheno build:
- `/home/amine/zee-lfv64-serial-20261112-094549/SPheno-4.0.5/ZeeTypeII/Couplings_ZeeTypeII.f90`, lines approximately 3880–3910.
- `/home/amine/zee-lfv64-serial-20261112-094549/SPheno-4.0.5/ZeeTypeII/TreeLevelMasses_ZeeTypeII.f90`, lines approximately 620–645.

**Structural PASS:**
- The left-chiral charged-scalar–lepton–neutrino vertex contains
  `resL = ZP(a,3)/2 * sum_{ij} [conj(UV(nu,j))*conj(ZEL(ell,i))*fzee(i,j) - conj(UV(nu,i))*conj(ZEL(ell,j))*fzee(i,j)]`.
  Renaming dummy indices and using `fzee(i,j)=-fzee(j,i)` gives
  `resL = ZP(a,3)*sum_{ij} conj(UV(nu,j))*conj(ZEL(ell,i))*fzee(i,j)`.
  The two half-terms combine: no residual factor of 1/2.
- The right-chiral charged-scalar vertex includes both `Ye` weighted by `ZP(a,1)` and `Ye2` weighted by `ZP(a,2)`, with `ZER` and `UV` rotations.
- The charged-scalar mass matrix has `M13=Mu4*vu/sqrt(2)` and `M23=-Mu4*vd/sqrt(2)`. With `H_D+=-sin(beta)H1+ + cos(beta)H2+`, the physical doublet–singlet mixing is `M_(HD,S)=-Mu4*v/sqrt(2)`, agreeing with `Model/model.py` and the original Gate 6 scalar convention.
- The printed `RXiWm` terms are gauge-fixing terms in the charged-Goldstone sector, not new scalar couplings.

**Not established by this extraction:** independent literature normalization and charge-conjugation/chirality matching at amplitude level, concrete numerical `UV`/charged-lepton basis mapping, one-loop neutrino mass, complete LFV amplitudes, or a well-defined symmetry protecting all omitted gauge-allowed terms. This record is based on user-provided source excerpts; it does not claim an independent recomputation of the generated vertices.

**Release decision:** `NATIVE_YUKAWA_INPUT_PASS + SARAH_CHECKMODEL_COMPLETED_WARNINGS_OPEN + CHARGED_VERTEX_STRUCTURE_PASS + CHARGED_MIXING_SIGN_PASS`; complete model/LFV physics validation remains **UNKNOWN**.


## 2026-11-18 Neutrino and charged-lepton basis conventions audit

### Source-matched leading one-loop neutrino mass

Primary literature:
- J. Herrero-García et al., JHEP 04 (2017) 130, arXiv:1701.05345v2, Eqs. (2.15)–(2.18), (2.30)–(2.31).
- J. Heeck and A. Thapa, Phys. Lett. B 841 (2023) 137910, arXiv:2303.13383, Eqs. (3)–(4).
- Relevant native generated vertex excerpt: `Couplings_ZeeTypeII.f90`, around lines 3889–3909 (provided by user in this project).
- Previous independent neutrino module: `aalida0/zee` `module2/core.py` and `module2/interface.py`, branch `research/zee-full-constraints`.

Our scalar VEV convention is `v1=v*cos(beta)`, `v2=v*sin(beta)`. In the **charged-lepton mass basis** define

```text
S = cos(beta)*Ye + sin(beta)*Ye2 = sqrt(2)*diag(m_e,m_mu,m_tau)/v
D = -sin(beta)*Ye + cos(beta)*Ye2
Ye2 = sin(beta)*S + cos(beta)*D
```

The 2017 paper writes the singlet lepton term without the explicit one-half that appears in our SARAH `Fzee/2 * l.l.Splus` term, and its physical charged vertex contains `2*f_paper`. The generated SPheno `resL` contains `Fzee` (the antisymmetric terms combine without an extra 1/2). For equal physical vertices, **`Fzee = 2*f_paper`**. The published 2023 Higgs-basis formula is

```text
Mnu = kappa * [f_paper * diag(m_l) * D + D.T * diag(m_l) * f_paper.T]
kappa = sin(2*gamma)/(16*pi^2) * log(mH1^2/mH2^2)
```

and thus the SARAH-normalized form has coefficient `sin(2*gamma)/(32*pi^2)` multiplying the analogous `Fzee` terms. This agrees with the **leading small-charged-lepton-mass limit** of the older `module2/core.py` kernel. Here `f.T` is a transpose (not Hermitian conjugate). Absolute loop/sign normalization relative to native SPheno remains an independent numerical gate.

**Executed locally in the assistant's isolated Python environment:** four standard-library tests in `Model.test_conventions` and `Model.test_basis` pass. The former compares the 2017 generic-basis formula, the 2023 Higgs-basis formula, and the SARAH factor-of-two convention for nonzero complex, off-diagonal entries. This is an **analytic/algebraic identity check**, not a calculation of the generated SPheno neutrino self energy. The latter tests a synthetic charged-lepton benchmark and deliberately corrupted matrix. The exact combined GitHub test suite still requires the user's local execution.

The tree-level `vL -> VL` SARAH `UV` rotation belongs to a **zero** neutrino mass matrix and cannot be identified with the physical PMNS matrix. The physical PMNS matrix is determined only after the loop-generated Majorana matrix is Takagi diagonalized in the properly fixed charged-lepton basis.

### One-card native charged-lepton-basis gate (user-reported PASS)

`Model/check_basis.py` reads the **corrected existing** SPheno `baseline.spc` from the prior input gate, and checks that `S=c_beta*Ye+s_beta*Ye2` is real, diagonal, ordered and positive; it reconstructs the complete complex `D` and checks `Fzee+Fzee.T=0`. It **does not** compute PMNS, neutrino masses, LFV or an experimental acceptance. This is a tree-level Yukawa matching check. It does not equate running Yukawa eigenvalues to pole lepton masses.

```bash
cd ~/Zee-Pheno
git pull --ff-only origin main
python3 -m unittest Model.test_model Model.test_native Model.test_conventions Model.test_basis -v

WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.check_basis --spectrum "$WORK/model-validation-native-20261118/baseline.spc"
```

If the diagnostic fails, **do not force diagonalization or loosen tolerances**: check the output renormalization scale and SARAH charged-lepton rotations `ZEL/ZER` first. In a non-diagonal gauge basis, transform the Yukawas before applying the mass-basis test.

### Restricted model / symmetry warning register

The 2026-11-18 SARAH `CheckModel[]` reported gauge-invariant terms not written in the input potential. The checked implementation is a **restricted effective Lagrangian**, not yet a demonstrated radiatively closed complete Zee model.

| Warning class | Disposition | Confidence |
| --- | --- | --- |
| `H1.H1 S-` and `H2.H2 S-` | Vanish identically for commuting identical SU(2) scalar doublets under antisymmetric contraction | Algebraic, not added |
| Opposite-doublet up/down quark Yukawas | Omitted by the chosen Type-II quark texture; no complete exact symmetry assignment has yet been demonstrated with the general lepton sector | Scope restriction |
| `|H1|^2 H1^dag H2`, `|H2|^2 H1^dag H2` and conjugates | `lambda6`, `lambda7`-type allowed quartics absent at the specified input scale; **RGE closure UNKNOWN** with both lepton Yukawas and soft terms | Scope/STOP for high-scale claim |
| Quark coupling to charged singlet (schematically `bar(u_R) d_R S+`) | Listed by SARAH's gauge-only possible-term checker; independent Lorentz/color/Fermi-statistics and lepton-number verification is required before treating this as a genuine renormalizable interaction | Scope/STOP for fully general model |
| Other checker-listed operators | Check individually for redundancies and symmetry before asserting forbidden terms | UNKNOWN |

The absence of `lambda6,lambda7` and other allowed operators must **not** be advertised as protected by an exact `Z2` without explicitly checking all lepton and scalar couplings.

**Current release:** tree-level scalar and charged vertices structurally checked, one-card native `Fzee/Ye2` input PASS; leading published neutrino-formula **convention bridge PASS (algebra only)**; native charged-lepton matching **user-reported PASS**; native one-loop Majorana mass matching **UNKNOWN**; symmetry/RGE closure **UNKNOWN**.


## 2026-11-18 Native mass-basis and conventional Z2 consistency audit

### User-provided native charged-lepton result

Result of `Model.check_basis` applied to the corrected `baseline.spc`:
- **`TREE_LEPTON_BASIS_PASS`**, `beta=0.8114830062749087` rad.
- `S=c_beta Ye+s_beta Ye2` diagonal, real and ordered: `[2.9350282036806365e-06, 6.068707417274492e-04, 1.0205763472088182e-02]`.
- Relative off-diagonal `2.2850369020644617e-10`; diagonal imaginary relative `4.0188679963578914e-11`.
- `max|Fzee+Fzee.T|=0.0`; the physical Higgs-basis `D=-s_beta Ye+c_beta Ye2` is a complex, non-diagonal matrix as intended.
- The script explicitly reports `pmns_from_tree_UV=NOT_ESTABLISHED` and `radiative_neutrino_mass=NOT_EVALUATED`. This is a **user-reported successful one-card test**, not independent review of the underlying `.spc`.

### Why our current model is not an exact or softly broken conventional Type-II Z2

Using the ordinary Type-II assignment
`H1: +, H2: -, Q_L: +, d_R: +, u_R: -`, choose `L_L=e_R=S+: +` so that `Ye L H1 e` and `F LL S+` are allowed. Then:

| Operator | Parity under conventional quark Type-II Z2 | Interpretation |
| --- | --- | --- |
| `Q H1 d`, `Q H2 u`, `L H1 e`, `LL S+` | even | Allowed |
| `L H2 e` (`Ye2`) | **odd** | **Dimension-4 hard breaking** for general nonzero `Ye2` |
| `|S|^2 H1^dag H2` (`k12`) | **odd** | **Dimension-4 hard breaking** |
| `(H1^dag H1)(H1^dag H2)` (`lambda6`), `(H2^dag H2)(H1^dag H2)` (`lambda7`) | **odd** | Set to zero at chosen scale, but not symmetry-protected |
| `m12sq H1^dag H2`, `Mu4 S- H1.H2` | **odd** | Dimension-2/3 **soft breaking** |

With generic `Ye` and `Ye2`, no conventional family-universal `Z2` with `H1` and `H2` opposite parities can simultaneously allow both lepton matrices. Reassigning `S+` parity changes whether the `LLS` or trilinear term is forbidden but does not remove the `Ye2` contradiction.

**Crucial precedent:** the 2021 charged cascade paper [arXiv:2106.08332, EPJC 81 (2021) 1148](https://link.springer.com/article/10.1140/epjc/s10052-021-09956-2), §3.1 Eq. (17), includes `k12` in precisely this restricted scalar potential. Its §8 explicitly acknowledges that imposing the quark `Z2`-motivated texture while dropping the symmetry in the leptons is *not strictly symmetry-consistent* (their wording); they leave the combined quark+lepton flavor fit as future work. The 2017 general Zee scan [arXiv:1701.05345](https://arxiv.org/abs/1701.05345) also requires both lepton Yukawa matrices to reproduce oscillation data. Oredsson/Rathsman [arXiv:1810.02588](https://arxiv.org/abs/1810.02588) show that hard `Z2` breaking in 2HDM Yukawa/scalar sectors can radiatively spread and alter FCNCs. Their quantitative conclusions cannot be transferred to the charged-singlet extension without checking it.

**Scoped model decision:** retain the 2021-inspired **low-scale restricted Zee ansatz** (quark Type-II *texture*, general `Ye/Ye2`, `k12` as specified, `lambda6=lambda7=0` at the declared reference scale), not an exact/softly broken globally `Z2`-invariant theory. Do not claim `Z2` enforces the absence of wrong-Higgs quark couplings or `lambda6/7` under RG. For an experimentally allowed scan, specify the input scale, estimate radiatively generated terms or associated model-dependence, and test sensitivity of the main charged cascade predictions to those omissions. The approximate texture is **not automatically a full renormalizable high-scale theory**.

**Release status:** `NATIVE_TREE_LEPTON_BASIS_PASS` (user-reported); `CONVENTIONAL_SOFT_Z2_INCOMPATIBLE_WITH_CURRENT_LAGRANGIAN` (analytic); `RGE_CLOSURE_UNKNOWN`; `NEUTRINO_1LOOP_NATIVE_UNKNOWN`.


## 2026-11-18 Neutrino native self-energy execution gate (NOT YET PASSED)

**Necessary distinction:** `CheckModel[]` found the neutrino tree-level mass matrix identically zero. The user-validated SPheno charged-lepton basis does **not** prove that the generated program evaluates the one-loop Majorana neutrino mass. In SARAH/SPheno, input switch `SPhenoInput[55]` selects loop mass corrections: `0` disables them, `1` enables them (see Florian Staub, *Exploring New Models in All Detail with SARAH*, 2015). Its treatment must be checked for this specific generated model; simply finding a `SelfEnergyFv` symbol is not enough.

The minimal read-only diagnostic `Model/audit_neutrino_native.py` inspects the corrected benchmark and the already-generated Fortran; it makes **no mass prediction or PASS claim**. Its two targeted `Model.test_neutrino_native` tests passed on synthetic fixtures. No native one-loop neutrino run has been executed in this environment.

```bash
cd ~/Zee-Pheno && git pull --ff-only origin main
python3 -m unittest Model.test_neutrino_native -v

WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.audit_neutrino_native \
  --model "$WORK/SPheno-4.0.5/ZeeTypeII" \
  --card "$WORK/model-validation-native-20261118/baseline.in" \
  --spectrum "$WORK/model-validation-native-20261118/baseline.spc"
```

Interpretation: `LOOP_MASSES_DISABLED` blocks a numerical loop-mass comparison from this original spectrum. `SOURCE_OR_OUTPUT_INSPECTION_REQUIRED` also remains **UNKNOWN**, even if `SelfEnergyFv` or nonzero PDG 12/14/16 masses appear. The decisive next step is to inspect the full generated neutrino self-energy, whether it is actually called in the SPheno mass calculation, and then compare its properly scaled complex symmetric Majorana matrix against the independently established full-finite Zee expression on a single fixed card.

A native SPheno loop-mass comparison is **not** automatically mandatory if the generator omits it: the already-tested published analytic Zee formula can be the validated neutrino backend once externally benchmarked, with SPheno used for a separate scalar/lepton/LFV calculation. Do not rebuild SPheno solely to manufacture a neutrino mass if the physical calculation is already independently validated.

**Release decision:** `NATIVE_MNU_1LOOP_UNKNOWN`; Step 0 remains open only for the requested native loop-mass equivalence gate. All earlier tree-level and Yukawa-input gates retain their scoped PASS statuses.


## 2026-11-18 Native loop-mass gate: disabled baseline identified, one-card probe available

**User-provided output of `Model.audit_neutrino_native`: `LOOP_MASSES_DISABLED`** on the previously validated baseline spectrum.
- Original benchmark input SHA-256 `f2efe759aa3196978d466739fc788b1dd7fb01d71cbc6e883c8a87830da454bf`.
- Original spectrum SHA-256 `3a7a4fab80c9e1837e72ea4a4e0ab70ea222367d97606b879d9d7115bea1c4d8`.
- Both occurrences of `SPhenoInput[55]` were `0`; three output neutrino pole masses (PDGs 12,14,16) were `-0.0 GeV`. These zeros cannot test radiative neutrino-mass generation.
- The generated `LoopMasses_ZeeTypeII.f90` contains calls of `Sigma1LoopFv` at lines 606, 666, 757, with the routine defined at lines 854–1024; these are **source-symbol evidence**, not yet a confirmed physical Majorana mass.
- The previously validated compiled executable SHA-256 is `2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52`.

Official SARAH/SPheno flag documentation: Florian Staub, *Exploring New Models in All Detail with SARAH* (2015), `SPhenoInput[55]=0` disables loop-corrected masses; `1` enables them. See https://onlinelibrary.wiley.com/doi/10.1155/2015/840780.

**Controlled follow-up (not executed here on the user's SPheno binary):** `Model/probe_neutrino_loop.py` verifies both SHA-256 pins, changes **exactly two** `SPhenoInput[55]` records from `0` to `1`, writes a new card in a previously nonexistent directory, runs the identical executable, and collects the `MASS` entries and logs. The three synthetic tests in `Model/test_loop_probe.py` pass. It never overwrites the original benchmark or output and never returns a physics PASS.

```bash
cd ~/Zee-Pheno && git pull --ff-only origin main
python3 -m unittest Model.test_loop_probe -v
WORK="$HOME/zee-lfv64-serial-20261112-094549"
python3 -m Model.probe_neutrino_loop \
  --card "$WORK/model-validation-native-20261118/baseline.in" \
  --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
  --output-dir "$WORK/model-neutrino-loop-20261118"
cat "$WORK/model-neutrino-loop-20261118/result.json"
grep -nEi 'error|warn|problem|fail' "$WORK/model-neutrino-loop-20261118/loop_enabled.log" | head -n 35
# Inspect the actual generated neutrino Majorana self energy when diagnosing the output:
sed -n '854,1024p' "$WORK/SPheno-4.0.5/ZeeTypeII/LoopMasses_ZeeTypeII.f90"
```

**Interpretation:** `NEUTRINO_POLE_MASS_NONZERO_UNVALIDATED` means only that the loop-enabled executable returned nonzero pole-mass entries. The full complex symmetric `Mnu`, phase, basis rotation, finite loop integrals, sign and normalization still require physical matching to the published Zee expression and an independent calculation. `NEUTRINO_MASS_ZERO_OR_MISSING` requires examining the self-energy/call site and whether loop masses were actually evaluated. Numerically small eV masses may be suppressed or rounded in generic spectrum-output handling. Any SPINFO warnings or unexpected scalar-spectra shifts must also be investigated.

**Status remains:** `STEP0_NATIVE_MAJORANA_MASS_UNKNOWN`; baseline loop mass setting was the observed blocker.
