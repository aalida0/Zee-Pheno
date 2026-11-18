# Native Zee neutrino discrepancy — resolved, 2026-11-18

**Cause: different renormalized Yukawas and internal running lepton masses.**
The original analytic comparison used the loop-disabled output; the native
loop-enabled calculation recomputes `Ye` by matching the measured lepton inputs.
Consequently its physical combination `D=-sb*Ye+cb*Ye2` differs. A cancellation in
`D(3,3)` amplifies a small change in `Ye(3,3)` into a large flavor-dependent mass
change. Neither an incorrect Zee loop coefficient nor a broken Takagi solver was
found. No production Fortran, model assumption, executable, or reference card was
changed. The correction is to the comparison's parameter matching.

## Reproduced results

All five supplied executable/card/spectrum SHA-256 values were verified before
and after the investigation. Provenance, full complex input arrays, generated
vertex arrays, native self-energies, returned rotation, and numerical results are
retained in [reference/native_neutrino_snapshot.json](reference/native_neutrino_snapshot.json).

| Calculation | m1 (eV) | m2 (eV) | m3 (eV) | Delta m21^2 (eV^2) | Delta m31^2 (eV^2) |
|---|---:|---:|---:|---:|---:|
| Original analytic, unmatched baseline | .009301829404 | .012721006947 | .050958060748 | 7.52999875e-5 | 2.51019992e-3 |
| Original native SLHA | .012785582000 | .041449883200 | .062565579300 | 1.55462171e-3 | 3.75098061e-3 |
| Analytic, actual native running inputs | .012785581989 | .041449883225 | .062565579316 | 1.554621713e-3 | 3.750980608e-3 |
| Native, before SLHA rounding | .012785581989 | .041449883225 | .062565579316 | 1.554621713e-3 | 3.750980608e-3 |

Matched analytic/native maximum matrix residual divided by the largest native
entry: **1.721e-14**; largest relative singular-mass residual **<2e-14**.
SLHA rounding permits only about 1e-9 relative comparisons. Regression thresholds
are 1e-11 for the full matrix and 1e-8 for rounded masses; these measure numerical
agreement, not an estimate of uncalculated physical corrections.

| Parameter | Original comparison | Actual native self-energy |
|---|---:|---:|
| v (GeV) | 246.219650039 | 246.709282247 |
| me (GeV) | .00051099893 | .000503190433 |
| mmu (GeV) | .105658372 | .104043812 |
| mtau (GeV) | 1.77686 | 1.74970486 |
| charged scalar masses (GeV) | 674.274427, 1165.97482 | 673.95831715, 1166.61113182 |
| sin(2 gamma) | .990465251816 | .990343945182 |
| Re Ye(3,3) | .00692697912 | .00667144383 |
| Re D(3,3) | .000136367865 | .000321709868 |
| Im D(3,3) | -.000146371656 | -.000146371649 |

`F` and `Ye2` retain their input values to SLHA serialization accuracy, respectively
7.76e-10 and 4.69e-10 relative. Beta is unchanged. The loop-output lepton MASS block
is explicitly reset to pole masses; it conceals the running masses used inside
the loop. The loop-output charged scalar poles, 696.792942 and 1272.12599 GeV,
are also **not** those internal masses.

Substituting only native scalar masses/mixing into the old analytic expression
gives (.00931792, .01274301, .05104620) eV; substituting only running lepton masses
gives (.00915968, .01252662, .05017936) eV. Substituting only native `D` gives
(.01296157, .04202045, .06342679) eV. These diagnostic substitutions isolate the
large effect of `D`; they are not alternative consistent renormalization schemes.

## Source trace and amplitude check

Line references below are to the original installation at
`/home/amine/zee-lfv64-serial-20261112-094549/SPheno-4.0.5/`, not the instrumented copy.
Its SARAH headers identify 4.15.4; SPheno is 4.0.5, compiler gfortran 11.4.0.
Source and linked-library hashes accompany the snapshot.

- `ZeeTypeII/SPhenoZeeTypeII.f90:176-204`: automatic matching `19=-2` chooses
  `MatchingOrder=-1` without loop masses and `2` with loop masses. Thus changing
  `55` also changes matching. The pole-input path constructs SM Yukawas from
  pole lepton masses at lines 250-268; the matched path calls `CalculateSpectrum`.
- `ZeeTypeII/Boundaries_ZeeTypeII.f90:1013-1033,1182-1204`: running lepton
  masses, self-energy matching via `Yukawas3`, and reconstruction of `Ye`.
- `ZeeTypeII/Model_Data_ZeeTypeII.f90:2099`:
  `Ye=-(vu/vd)*Ye2+(vSM/vd)*YeSM`. At fixed `Ye2`, modifying `YeSM` and the vev
  modifies `D`. This exactly accounts for the measured parameter changes.
- `ZeeTypeII/LoopMasses_ZeeTypeII.f90:174-179,295-298`: `TreeMasses`, then
  `CouplingsForLoopMasses`, then `OneLoopFv`. Internal scalar masses and `ZP`
  are tree values at the matched running parameters. Scalar pole calculations
  occur afterward; lines 375-384 reset external lepton masses to pole values.
- `ZeeTypeII/Couplings_ZeeTypeII.f90:11068-11083,11243-11258`: the two
  antisymmetric half-couplings combine into the **full** `F` vertex. With
  identity charged-lepton rotations,
  `L_ai=-F_ai*ZP(k,3)`, `R_ai=-conj(Ye_ia*ZP(k,1)+Ye2_ia*ZP(k,2))`.
  The opposite charge orientation supplies the transpose contraction.
- `ZeeTypeII/LoopMasses_ZeeTypeII.f90:884-903,966-986,1019-1022`:
  each scalar orientation has `2*m_l*B0` times `1/2`, and the total has
  `1/(16*pi^2)`. Their sum gives `A+A.T`, never `A+A.dagger`.
- The W diagrams have one zero chiral coupling and no Majorana mass insertion;
  the Z diagram has an internal zero tree neutrino mass. Wave terms multiply
  the zero tree matrix. Therefore these diagrams generate no extra mass.
- `ZeeTypeII/LoopMasses_ZeeTypeII.f90:610-613,670-673,761-764`:
  `M=-[conj(SigSL)+SigSR]/2` for zero tree mass. Numerically
  `conj(SigSL)=SigSR`, and `M=M.T` exactly within printed precision.

Let `q_k=ZP(k,3)` and `h_k=-sb*ZP(k,1)+cb*ZP(k,2)`.
The observed convention has `2*q_2*h_2=sin(2 gamma)` (Fortran index 2 is the
lighter physical scalar). Orthogonality cancels the common B0 constant/scale;
`B0(0,ml^2,m1^2)-B0(0,ml^2,m2^2)=-I_l`. The minus sign from mass reconstruction
then yields precisely

```
M = sin(2 gamma)/(32 pi^2) * [F diag(ml I_l) D + transpose] .
```

This normalization is derived from both generated charge orientations, without
fitting a coefficient. Expanding the literature term without the explicit 1/2
gives two identical antisymmetric terms; hence `F_SARAH=2*f_literature` in the
matched field convention. The result agrees with Heeck–Thapa
[arXiv:2303.13383v2, Eqs. (3)-(4)](https://arxiv.org/pdf/2303.13383), where the
coefficient multiplying `f` is `sin(2 phi)/(16 pi^2)` times the mass logarithm.
The Higgs rotation and charged-scalar convention agree with
[arXiv:1701.05345v2, Eqs. (2.1), (2.7)-(2.8), (2.15)-(2.17), (2.30)-(2.31)](https://arxiv.org/pdf/1701.05345)
after accounting for its Yukawa daggers and scalar ordering. The finite kernel
was independently integrated at 45 digits for all three actual lepton masses;
agreement with `finite_kernel` is better than 1e-13 absolute. Its massless limit
is the published logarithm. Conversion is 1 GeV = 1e9 eV.

## Bases, diagonalization, and controls

Generated `UFv` external legs are in the unrotated neutrino gauge basis. For
`e_mass=ZEL e_gauge`, `ec_mass=ZER ec_gauge`, the exact contraction is
`F ZEL.dagger diag(ml I_l) conj(ZER) D + transpose`.
For comparison with the shared analytic solver we use
`F_mass=conj(ZEL) F ZEL.dagger`,
`D_mass=conj(ZER) D ZEL.dagger`, and
`M_mass=conj(ZEL) M_gauge ZEL.dagger`.
The charged-lepton mass-diagonalization residual is 1.27e-16. Dropping the tiny
rotations would leave a 1.8e-7 relative matrix error, so they were retained.

`OneLoopFv:605-653,665-724,756-815` first evaluates at p^2=0, then iterates each
pole using tree internal masses. Every evaluated matrix in the spectrum call
is identical in binary64; maximum p^2 is 3.915e-21 GeV^2. Complex matrices are
diagonalized through `M.dagger M`, with phases fixed so
`conj(UV_1L) M UV_1L.dagger=diag(m)`. The observed Takagi residual is 2.45e-15.
Positive eigenvalues never enter the negative-small-eigenvalue repair branch.
SLHA states 12,14,16 are ascending masses. Tree `UV=i*identity` is not PMNS.

A Fourier unitary substituted for tree `UV` before generation of the loop
couplings leaves the complete native Majorana matrix **bitwise unchanged**.
A separate random complex congruence check preserves singular values to 1e-13.
No PMNS oscillation-angle or Majorana-phase fit is inferred from these tests.

| Check | Measured result | Status |
|---|---|---|
| Original hashes before/after | all five unchanged | PASS |
| Full matched complex matrix | relative residual 1.721e-14 | PASS |
| Independent sum of dumped generated vertices | relative residual 2.22e-15 | PASS |
| Symmetry / Takagi | 0 / 2.45e-15 | PASS |
| F=0, rerun native | exactly zero matrix and masses | PASS |
| F doubled, rerun native | masses (.025571163995,.082899766487,.125131158649) eV; largest scaling residual <3.2e-9 | PASS |
| Doubled full matrix versus twice baseline | normalized residual 1.25e-9 | PASS |
| Arbitrary tree UV control | full matrix unchanged | PASS |
| Scalar/lepton output preservation | every diagnostic MASS entry equals original loop spectrum | PASS |
| Existing Model regression suite | 25 tests passed | PASS |
| Intended oscillation splittings in native matched scheme | not reproduced | FAIL |
| PMNS angle / Majorana-phase phenomenology | not evaluated | UNKNOWN |
| Complete Step 0 beyond this calculator issue | previously open gates not re-certified | UNKNOWN |

The invariants also distinguish the original and native matrices: respectively
`Tr(M.dagger M)=.002845072003/.005796015641 eV^2` and
`abs(det M)=6.029797846e-6/3.315730950e-5 eV^3`.

There are **two** native `OneLoopFv` invocations in this run: a later Higgs-coupling
calculation uses modified local parameters and gives (.01284507,.04226890,.06307493)
eV. It does not supply the written neutrino MASS block. The diagnostic selects
the snapshot matching that block; using the last call would introduce a new
spurious discrepancy. An exploratory card changing only matching to explicit
`19=-1` terminated with a tadpole NaN; it is not evidence for or against the
resolution, and that alternative card is not validated. SPheno's exit code was
zero despite termination, so complete snapshots and masses are required.

## Reproduction and decision

```
python3 -m Model.diagnose_neutrino_native \
  --spheno /home/amine/zee-lfv64-serial-20261112-094549/SPheno-4.0.5 \
  --card /home/amine/zee-lfv64-serial-20261112-094549/model-neutrino-loop-20261118/loop_enabled.in \
  --output-dir /tmp/zee-neutrino-fresh
# Repeat in another fresh directory with --rotate-tree-neutrinos for the UV control.
python3 -m unittest discover -s Model -p 'test_*.py' -v
```

The diagnostic copies two sources, adds output only, and links against the
existing libraries in a fresh directory. No SARAH regeneration or production
recompilation is necessary. Investigation artifacts are in `/tmp/zee-neutrino-*`;
the essential reference snapshot is committed, not binaries or build trees.

**The calculator discrepancy is closed.** Native SPheno is consistent with the
analytic Zee amplitude under its actual matched inputs. The original analytic
values describe the loop-disabled parameter definition; their agreement with
target splittings does not validate that fit in the native matched definition.
Next, choose and freeze the renormalized input scheme and perform neutrino
phenomenology in that scheme. No refit was performed here. Full Step 0 remains
UNKNOWN pending its other existing gates; this result alone does not authorize
a physically accepted parameter scan or establish higher-order precision in
this cancellation-sensitive benchmark.
