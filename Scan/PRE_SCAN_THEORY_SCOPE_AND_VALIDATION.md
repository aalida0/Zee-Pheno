# Before scanning: completed LOW-SCALE Zee theory engine and exact scope

**Date: 21 November 2026. Branch: research/pre-scan-theory-completion.**

This report freezes the **low-scale tree-scalar + one-loop oblique**
theory definition. It is a physics scope, not a statement that *all*
possible higher-order/UV/vacuum observables in the full Zee model are
known. The branch contains a single per-native-spectrum gate usable
for any NO/IO card (no S00017 hardcoding).

**Do NOT run a full physics scan yet if the paper requires certified
one-loop effective-potential global stability, a high-energy cutoff,
or non-oblique EW corrections.** Such a theory scope must be stated
before scanning; it cannot be established from one native SLHA.

## Physics inputs and published calculations

1. **BFB**: exact restricted Zee Eq.17 all-complex-field sufficient
   quartic positivity proof using the positive-semidefinite
   charged-singlet portal matrix, plus separate necessary fail tests.
   No unknown portal configurations are artificially excluded.

2. **Complete scalar high-energy tree perturbative unitarity**: all
   2021 EPJC Appendix E.22–E.29 scattering matrices, including full
   13x13 sector and 15 independent analytic eigenvalue checks.
   Criterion is max |eigenvalue| < 8*pi, not |lambda| < 4*pi.
   The latter remains a separately named conventional prior.

3. **Local-scale running Yukawa perturbativity**: read **native**
   complex Yu,Yd,Ye,YEE2,FZEE, compute largest singular value for
   each (basis independent). Require max singular value < sqrt(4*pi).
   FZEE antisymmetry and source SPheno running Q are checked.
   At each point additionally test native scalar quartics and tree
   unitarity at the same Q. Native SPheno flag 38=2 states which beta
   functions the source generates, but is NOT proof of 10 TeV
   perturbativity.

4. **Independent full one-loop scalar S,T,U**: Grimus et al.,
   Nucl.Phys.B801 (2008) 81–96, Eqs. 28, 30, 31, Appendices B/C.
   Uses two charged eigenstates plus Goldstone and all three
   physical neutral scalars; full source-matched mixing matrices,
   explicit SM subtraction and stable one-dimensional numerical
   Feynman-parameter quadratures for the G and Ghat kernels.
   Matched to source SPheno in its **running gauge** convention.
   On-shell weak mixing angle predictions are included only as
   sensitivity variations; they are not erroneously required to
   reproduce a different source convention.
   Compare the independently calculated STU with the actual
   SPhenoLowEnergy values at numerical tolerance
   |delta S|<0.0005, |delta T|<0.0005, |delta U|<0.0001.
   The full correlated PDG2025 Table 10.8 U-FLOATING 3x3
   covariance is then evaluated; chi2_3 95% threshold 7.81473.
   This is a one-loop **oblique scalar** screen, not a claim of
   precision-complete non-oblique electroweak observables.

5. **Rigorous SUFFICIENT global EW tree-vacuum certificate**:
   Exact lower bound by completing the charged singlet's quadratic
   square and applying the known global-minimum theorem for a
   softly-broken Z2 2HDM (Barroso et al. EPJC 73 (2013) 2537).
   Define d=mu4^2/M_C^2 for M_C^2>0. If lambdaC>=0 and the
   2x2 scalar portal K is PSD,

     V_Zee(Phi1,Phi2,chi) >= V_aux_2HDM(Phi1,Phi2),
     lambda3_aux=lambda3-d, lambda4_aux=lambda4+d.

   Equality holds at our neutral EW vacuum with chi=0 and
   det(Phi1,Phi2)=0. If the AUXILIARY 2HDM potential is strictly
   BFB, has a normal EW local minimum (positive h/H/A/H+ tree
   masses squared) and satisfies the published global normal-minimum
   discriminant

     D=m12² (m11²-k²*m22²) (tanbeta-k)>0,
     k=(lambda1/lambda2)^(1/4),

   then the Zee EW vacuum is a **rigorous GLOBAL minimum of the
   DEFINED TREE potential over all complex field configurations**.
   It needs no random minimizer or merely missing deeper witness.
   This is a sufficient criterion, not necessary: failure leaves
   global-tree status UNKNOWN unless an actual deeper witness
   is demonstrated.

   The theorem is validated by an exact gauge-orbit parametrization
   with arbitrary CP phase and charged directions, an independent
   1800-complex-field lower-bound numerical challenge, and an
   intentionally constructed charge-breaking deeper witness.

   **CRITICAL SPheno matching caveat**:
   Original SPheno ZEEPOT m11/m22 are loop-tadpole-corrected and do
   **not** satisfy the tree stationarity relations when combined
   directly with the printed quartics/VEVs. NO residual ~0.139848,
   IO ~0.139849. Consequently the original SPheno running
   loop-tadpole-potential global minimum stays UNKNOWN.
   The rigorous proof applies to the explicitly *tree-tadpole
   reconstructed* reference potential, the appropriate defined
   tree-level screening scheme. On 21 Nov 2026, the original source
   executable with input option 510=1 was independently run for both
   s00017 NO and IO and reported source-matched tree-stationary ZEEPOT
   values meeting this sufficient theorem (user-reported run summary;
   raw artifacts retained locally). This is not a certification of the
   DEFAULT loop-tadpole ZEEPOT values. Never turn the tree result into
   a 1-loop effective-potential or metastability certification.

## Actual source-matched NO/IO numerical evidence

| Item | s00017 NO | s00017 IO |
|---|---:|---:|
| All-direction quartic BFB sufficient | PASS | PASS |
| Complete high-energy tree scalar unitarity | PASS | PASS |
| Max absolute scattering eigenvalue | 5.942698295 | 5.942698295 |
| Max running native Yukawa singular value | 0.9623703849 | 0.9623703859 |
| Local Yukawa < sqrt(4*pi) ~ 3.5449 | PASS | PASS |
| SPheno STU S/T/U | 0.00296942 / 0.01202010 / 0.00005590 | same |
| Independent Grimus S/T/U | 0.00303747 / 0.01217905 / 0.00002765 | ~same |
| STU numerical consistency at source gauge scheme | PASS | PASS |
| Correlated PDG2025 3D chi2 | 0.49490954 | 0.49490954 |
| Correlated 3D 95% experimental screen | PASS | PASS |
| Rigorous **defined TREE** EW global-minimum certificate | PASS | PASS |
| Native SPheno loop-tadpole global minimum certification | UNKNOWN | UNKNOWN |
| RG flow to 10 TeV | NOT REQUESTED IN LOW-SCALE SCOPE | same |
| Full all-sector model viability | UNKNOWN | UNKNOWN |

The **HiggsSignals reference screen independently disfavors both
s00017 points**, so theoretical PASS does NOT imply their experimental
survival. Their flavor/LFV sectors are also not fully validated.

## What code to use on every future native scan point

The reusable source-bound evaluator is
[theory_gate_native.py](theory_gate_native.py). Input is the SAME
schema_version=1 manifest used by Scan/native_pipeline.py, with
point_id, ordering, original card, native result JSON and full
SPheno spectrum. Source card SHA256 and original SPheno binary
SHA256 are validated. The nine running quartics, all native mixing
matrices, Yukawas, physical masses, SPheno output STU, full tree
EW potential and charged singlet quantities are parsed from
**that specific point**.

    cd ~/Zee-Pheno
    git fetch origin
    git switch research/pre-scan-theory-completion

    python3 -m unittest Scan.test_zee_oblique_STU_independent \
      Scan.test_zee_vacuum_global_sufficient \
      Scan.test_theory_gate_native \
      Scan.test_verify_tree_tadpoles \
      Scan.test_prescan_completion -v

    python3 -m Scan.theory_gate_native \
      --manifest /absolute/path/to/actual_native_manifest.json \
      --output /absolute/new/path/zee-native-theory.json

Do not supply --cutoff-Gev 10000 unless you are **explicitly requiring
UV validity up to 10 TeV**, with an attested continuous model-specific
RG trajectory. The default criterion is at the native SPheno
renormalization Q (160 GeV for original s00017).

The pass name is deliberately
DECLARED_LOW_SCALE_THEORY_SCOPE_PASS, not
ALL_PHYSICAL_CONSTRAINTS_PASS. It is safe as one theory-stage
preselection for a pheno parameter scan with the same scope.

## One-time SPheno tree-tadpole source reproduction (already scripted)

The original source SPhenoInput options are:
- 510=1: request tree-level tadpole solutions;
- 530=1: request additional tree/loop information for Vevacious.

Do NOT assume either flag guarantees a particular SLHA block name or
overwrites ZEEPOT in this custom generated model. The original cards set
both to zero. The verifier inspects source output, records available
TREE/LOOP blocks and requires explicit neutral tree stationarity before
any matching certificate. Missing source-matched tree output remains
UNKNOWN, not an exception or an automatic PASS.
[verify_tree_tadpoles.py](verify_tree_tadpoles.py) creates isolated
temporary variants for both NO and IO, **changing only one option
at a time**, verifies the original SPheno executable SHA256, runs
the real binary, archives logs/spectrum hashes, and checks
the tree-level stationarity/global certificate on the *actual*
source-output potential blocks.

**Original SPheno execution (21 Nov 2026; user-reported summary):**

| Case | Source-matched stationary blocks | Result |
|---|---|---|
| s00017_NO_tree510 | ZEEPOT | SOURCE_MATCHED_GLOBAL_TREE_CERTIFICATE |
| s00017_IO_tree510 | ZEEPOT | SOURCE_MATCHED_GLOBAL_TREE_CERTIFICATE |
| s00017_NO_tree530 | none | UNKNOWN_NO_MATCHED_NATIVE_TREE_TADPOLES |
| s00017_IO_tree530 | none | UNKNOWN_NO_MATCHED_NATIVE_TREE_TADPOLES |

Thus **510=1 is the demonstrated source-matched TREE-output mode** in this
specific custom SPheno executable for both reference controls. The
530=1 result is *not* a global-vacuum FAIL; it simply lacks matched tree
parameters in the inspected blocks. The outputs and logs remain on
the user's machine in `pre_scan_tree_tadpole_2026_v2`; they have not
been independently downloaded to this repository or rerun in CI.
The original default (510=0,530=0) spectra remain loop-tadpole
outputs. Future native scan points must still have their own
tree stationarity check, and their spectrum/observable differences
between tree and default tadpole modes must not be silently ignored.

    WORK="$HOME/zee-lfv64-serial-20261112-094549"
    python3 -m Scan.verify_tree_tadpoles \
      --exe "$WORK/SPheno-4.0.5/bin/SPhenoZeeTypeII" \
      --output-dir "$WORK/pre_scan_tree_tadpole_2026"

This uses ONLY original s00017 cards as independent scheme controls,
not a parameter-space scan. An option producing a valid SPheno spectrum
but no explicitly matched tree values is reported UNKNOWN with hashes
and retained raw logs, not claimed as a physics PASS. Inspect
`native-tree-audit.json` and each output `spectrum.slha` before
using any tree-potential scheme for *new* cards. The verifier cannot
execute the private local SPheno binary in GitHub CI.

## Other limitations do not disappear

- 1-loop effective-potential global stability and tunnelling lifetimes
  remain out of the declared TREE theory scope; no vacua are proven
  globally stable at one loop.
- High-cutoff model RG running is **not** inferred from local couplings
  and is outside the declared low-energy phenomenology theory scope.
- The simple Yukawa spectral criterion is conventional
  perturbativity, not a complete fermionic tree unitarity calculation
  or high-scale UV completion.
- Non-oblique EW corrections and a precision global likelihood
  (including updated experimental reference data) remain out of scope.
- FULL LFV, flavor and Higgs experimental/SM-covariance model
  constraints remain separate; the final physics scan should not
  claim all accepted points physical until those gates are resolved.

**No new parameter-space scan was started here.** Native original
s00017 test points and synthetic perturbations were used only to
validate the theoretical machinery and its source-binding contract.
