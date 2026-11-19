# Restored HEPfit backend and executed s00017 NLO B→Xsγ

19 November 2026; branch `research/spheno-hepfit-bsgamma`.

**Execution PASS. Physical exclusion UNKNOWN.** A persistent, source-qualified
HEPfit executable is working at `/home/amine/zee-hepfit-build/driver`. All
37 historical HEPfit runs reproduce every recorded observable exactly. The
actual two-state s00017 calculations have completed, with raw input/log
receipts. No frozen historical calculation or native benchmark was modified.

## Numerical results, Eγ > 1.6 GeV

| Prediction | BR NLO | NP shift from same-scheme SM |
|---|---:|---:|
| Historical HEPfit SM control, mh=125.09 GeV | 3.491423243740102e-4 | — |
| s00017 NO same-scheme SM, mh=125.25 GeV | 3.491424127389079e-4 | 0 |
| **s00017 NO Zee** | **3.921419828565858e-4** | **4.299957011767789e-5** |
| s00017 IO same-scheme SM, mh=125.250135 GeV | 3.491424128134351e-4 | 0 |
| **s00017 IO Zee** | **3.921419354070257e-4** | **4.299952259359058e-5** |
| Original SPheno FlavorKit, NO | 3.634717890000000e-4 | Not a same-scheme shift |
| Original SPheno FlavorKit, IO | 3.634717890000000e-4 | Not a same-scheme shift |

The relative NP shifts are +12.3157682% (NO) and +12.3157546% (IO).
HEPfit exceeds the archived FlavorKit numbers by 2.86701939e-5 and
2.86701464e-5, approximately 7.8879%. This difference is **not an uncertainty
estimate**. FlavorKit uses its preserved `NNLO_SM=3.15` normalization in
units of 1e-4 and its polynomial response, whereas the HEPfit calculation
uses its own source-qualified NLO/EM/HQE prescription. No intercept was
fitted, no SM normalization was replaced, and no exclusion was inferred.

[Actual s00017 result](reference/hepfit_restored/s00017/bsgamma_nlo.json)
contains exact source hashes and separate SM/zero-NP/Zee execution receipts.

## Environment repair and source qualification

Initial inspection found no HEPfit source/build directory or executable to
reuse, including the missing `/tmp/zee-hepfit/driver`. The machine had 22 GiB
RAM, about 8.7 GiB available, 12 logical CPUs and 242 GiB free disk. The build
used two C++ jobs, serial SuperIso make, and occupies approximately 678 MiB.

Reused installations: GCC 11.4.0, CMake 3.22.1, ROOT 6.32.12 and Boost 1.74.
GSL headers are 2.7.1; the original Herwig GSL runtime is 2.7. The latter's
actual bytes, ROOT config and MathCore/MathMore bytes match the historical
build receipt. ROOT still required the old `libtbb.so.2` and `libgsl.so.25`
SONAMEs, so the SHA-pinned Ubuntu libtbb2 package was extracted locally and
the existing Herwig GSL library copied into a local compatibility directory.
No global installation or system-library replacement occurred.

The new [restore script](restore_hepfit.sh) follows the old REPRODUCE.md,
reuses verified archives on reruns, and builds only the required library
clients. The four downloaded archives were all SHA256-verified. Before
compilation, all 346 locked build/rate/configuration files matched the
source lock after applying exactly the required constructor patch.

| Component | Source identity |
|---|---|
| HEPfit | `322bd73a2dadd94ce94a6d71cb79d616e569120a` |
| gslpp | `f5b059e058a47f6797823a3846fd2cf9d17ff12c` |
| RGESolver | `d1eb75559e03ae74a9ed0532760560757ef3beec` |
| Initialization patch SHA256 | `fbf9c7daa170da46506785e6074bf2648299f110867788abd7a6b36096aa44ec` |
| Independent BMU matcher SHA256 | `0c2dfb32a7138b02fa4a52c163c63f942f313bb5d697545e138ef1915c200c4e` |
| Historical executable SHA256 | `2e622645491118041361607faf040120fec59103cf0dfd85dbb6925de5fee0a3` |
| **Restored executable SHA256** | **`a9ce45a05c77d12b8cf2259ea7757e4d56a2c91c144dbea46d67130ac00a0786`** |

The exact unchanged initialization fix zero-initializes CWbsgArrayLO/NLO/NNLO
and initializes four cache keys to −1. Every fresh non-SM driver run passes
the placement-construction witness with previously nonzero storage: unused
LO and NNLO C1 slots return exactly zero. No matching kernel or physics
formula was patched. The historical unpatched falsifier remains preserved.

The older `~/zee` checkout was actually on
`research/zee-lfv64-closure-20261113`, with unrelated user edits. Its HEPfit
adapter and BMU source files have no differences from the specified
`research/zee-full-constraints` branch. Source hashes qualify their use; the
checkout was not switched, cleaned or modified. Exact input identities and
resource observations are in [INPUT_PROVENANCE.json](reference/hepfit_restored/INPUT_PROVENANCE.json).

## Historical numerical validation

The unchanged original `execute.py` and original source-qualified scalar
cards were used, including its independent SuperIso comparison clients.
Outputs were written into the new persistent `historical-replay/` directory,
not over any historical receipt.

| Historical central control | Required value | Reproduced value |
|---|---:|---:|
| D7 | 4.17955384241e-4 | 4.1795538424101365e-4 |
| F13 | 3.83711422745e-4 | 3.837114227449281e-4 |
| SM | 3.491423243740102e-4 | 3.491423243740102e-4 |

All 37 run input dictionaries, modes and complete output vectors match the
historical receipt, with maximum absolute numerical difference **zero**.
The original source/receipt verifier passes on the new execution. No expected
value or original test was changed. The original 23 receipt/uncertainty tests
also pass; those tests read the preserved historical data, while the separate
fresh replay establishes the rebuilt executable's numerical equivalence.

The historical suite includes SM/zero-NP equality, both one-state Type-II
limits, equal masses, charged-state exchange, native Type-II matching,
ordinary-to-effective conversion, the NLO alpha_s factor, interpolation of
the published table, signed C4/C7/C8 responses, an independent held-out
bilinear interference test, and 120/200 GeV scale controls. The old HEPfit vs
SuperIso prescription mismatch is reproduced and remains a limitation; it
was not removed by a normalization adjustment.

[Fresh full execution receipt](reference/hepfit_restored/historical-replay/HEPFIT_EXECUTION.json)
and all referenced raw logs/cards are committed. Build warnings are retained
in [build.log.gz](reference/hepfit_restored/build.log.gz).

## Provenance enforcement for the rebuilt executable

The original executable SHA check was retained as the legacy path. A changed
binary without an attestation is rejected. The added
[validator](validate_hepfit.py) does all of the following before issuing or
accepting a rebuilt-runtime attestation:

1. Verify the source lock, patch, BMU matcher, all locked HEPfit source files,
   actual driver source, and actual executable SHA256.
2. Verify every build-receipt artifact, including static libraries and
   reused ROOT/GSL/TBB artifacts. Bind the CMake cache to the supplied source.
3. Apply the unchanged historical `provenance.py` execution verifier to both
   historical and fresh receipts; authenticate every raw stdout/input file.
4. Reparse fresh raw stdout and cards, rather than trusting JSON claims.
5. Require all 37 exact input/mode/output-key sets and compare every number
   against the historical run: 2e-14 absolute for operator coefficients,
   2e-13 for other outputs, using the original operator/implementation bounds.
   Observed differences are zero. These are implementation checks, not
   theory errors; the original suite's other tighter algebraic checks also
   reproduce because their output vectors are identical.
6. Recompute this verification when the bridge consumes the attestation;
   a self-declared PASS or stale receipt hash is insufficient.

The fixed new-point parser still requires μw=160 GeV. The validator explicitly
supplies each historical run's matching scale when checking the original
120/200 GeV controls. Its first development attempt rejected these controls
because the shared parser assumed 160 GeV; this parser issue was corrected
without changing any input, expected value or tolerance. No backend numerical
run failed. Regression falsifiers cover changed binaries, missing runs,
modified inputs/coefficients/rates, NaNs and tampered receipts.

[RESTORED_BACKEND.json](reference/hepfit_restored/RESTORED_BACKEND.json),
[BUILD_RECEIPT.json](reference/hepfit_restored/BUILD_RECEIPT.json),
[SOURCE_LOCK.json](reference/hepfit_restored/SOURCE_LOCK.json) and
[RUNTIME_LIBRARIES.json](reference/hepfit_restored/RUNTIME_LIBRARIES.json)
retain build, source and shared-library identities. The executable stays in
persistent local storage; no linked binary or third-party source tree is
committed. GPL notices and complete corresponding-source obligations still
apply if someone later distributes linked binaries.

## New native inputs and matching conventions

The bridge reads the frozen native `input.in`, `spectrum.slha` and
`native.json` for each fitted s00017 benchmark. Card hashes and the original
SPheno production-binary identity are checked. Both charged **pole** masses
and the physical rows of the native **CHARGEMIX** are used.

| Native input | NO | IO |
|---|---:|---:|
| mH1± pole, GeV | 666.305092 | 666.305622 |
| mH2± pole, GeV | 1184.553650 | 1184.553680 |
| tanβ from vu/vd | 4.6927909112 | 4.6927909112 |
| H1± squared doublet norm | .9776069374014071 | .9776069374014071 |
| H2± squared doublet norm | .02239306202065159 | .02239306202065159 |

The weights are sums of squares of the two doublet columns in each native
physical charged row, without a fitted angle or silent renormalization.
Their sum differs from one by only 5.8e-10 due to printed SLHA precision.
βH in HMIX[10] is a scalar angle; the Yukawa tanβ is read from the running
vevs. Their small loop-level difference is retained as documented by the
existing bridge, not treated as a new matching convention.

The **same HEPfit prescription** supplies the top mass and SM normalization.
The original 173.2 GeV top input is explicitly interpreted as a pole mass
through `FlagMtPole=true`, giving mt_MSbar(160)=164.29247164143644 GeV.
The SPheno running top mass is not substituted into this prescription.
Other SM/hadronic inputs remain those of the validated HEPfit cards; only
mh is aligned between each Zee point and its corresponding SM control.

Scales: μw=μt=160 GeV, μb=2 GeV, μc=1.99999999 GeV, with Eγ>1.6 GeV.
The measured alpha_s(160)/(4π) is **0.008689679872399127**.
For each operator, the ordinary matching coefficients are summed coherently:

```text
Delta C_i = w1 Delta C_i^II(mH1) + w2 Delta C_i^II(mH2).
```

The source-qualified BMU evaluator retains real, top-only, color-singlet
charged-Higgs matching through NLO QCD in the ordinary NDR basis. It does
not average two branching ratios. Both charged-state contributions enter
before evolution, interference and evaluation of the finite-cut rate.

| Ordinary NP matching coefficient at μw | NO | IO |
|---|---:|---:|
| ΔC7^(0) | −.0464087353537987 | −.046408688203482575 |
| ΔC8^(0) | −.051669556843575014 | −.05166951018447357 |
| ΔC4^(1), raw | .001669812282783969 | .001669810736247945 |
| ΔC7^(1), raw | 1.0000746811312704 | 1.0000740166958761 |
| ΔC8^(1), raw | .5588120202020941 | .5588117685455243 |

Raw NLO entries have **no** alpha_s/(4π) included. In the executed C++ client,
with a=alpha_s(μw)/(4π), the nonzero order-tagged effective NP entries are

```text
LO:  C7eff = C7^(0), C8eff = C8^(0)
NLO: C4 = a C4^(1)
     C7eff = a [C7^(1) − (4/9) C4^(1)]
     C8eff = a [C8^(1) − (1/6) C4^(1)].
```

These actual high-scale entries are independently checked against driver
stdout, including the nonzero C4 contribution and exactly one NLO alpha_s
factor. The SM matching vector is added once. General C3–C6 conversion and
RG mixing are retained by the original client. NP NNLO matching is zero;
`BR_stock_NNLO` is intentionally ignored for the Zee prediction.

## Scientific status

| Claim / gate | Status | Evidence and limit |
|---|---|---|
| Pinned source and initialization fix restored | PASS | Four archives, 346 files, patch and runtime zero-slot witness |
| Historical SM/D7/F13 and full implementation suite | PASS | 37 fresh numerical vectors identical; original 23 tests pass |
| New native pole masses/doublet projections | PASS | Archived SPheno artifacts and existing bridge guards |
| SM/zero-NP equality for each new point | PASS | Actual separately executed controls |
| Two-state matching/effective conversion/NLO normalization | PASS | Full original algebraic suite and actual high-scale vectors |
| Actual finite-cut central NO/IO prediction | PASS | Six fresh HEPfit processes with source-bound logs |
| Same-order matched independent physical rate validation | UNKNOWN | Historical SuperIso prescription mismatch persists |
| Complete model-specific uncertainty | UNKNOWN | Not calibrated by this restoration |
| Physical B→Xsγ acceptance/exclusion | UNKNOWN | theory_sigma=null; experimental_exclusion=null |
| Global Zee viability | UNKNOWN | Outstanding flavor/LFV/Higgs/vacuum and other gates remain |

The retained calculation includes HEPfit's native partial EW/HQE and
semileptonic normalization prescription. It is not a uniformly truncated full
Zee QCD+QED/EW calculation. General primed/light-quark matching,
model-specific electroweak corrections, resolved/nonlocal photon effects and
correlated NP uncertainties remain incompletely covered. No uncertainty was
fabricated, no experimental limit was converted into a verdict, and no scan
or new benchmark fit was performed.

## Reproduction

The already restored backend can be used directly with the README command.
To rebuild on this machine, use the following with fresh receipt/output names:

```bash
set -e
cd /home/amine/Zee-Pheno
HEPFIT_BUILD="$HOME/zee-hepfit-build"
HEPFIT_SOURCE="$HEPFIT_BUILD/HEPfit-322bd73a2dadd94ce94a6d71cb79d616e569120a"
HEPFIT_TASK="$HOME/zee/research/full_constraints/bsgamma_hepfit"
bash Flavor/restore_hepfit.sh "$HEPFIT_BUILD" "$HOME/zee" > /tmp/zee-hepfit-rebuild.log 2>&1
mkdir "$HEPFIT_BUILD/historical-replay-02"
PYTHONDONTWRITEBYTECODE=1 python3 "$HEPFIT_TASK/execute.py" \
  --source "$HEPFIT_SOURCE" --binary "$HEPFIT_BUILD/driver" \
  --superiso "$HEPFIT_BUILD/superiso_probe" \
  --superiso-inputs "$HEPFIT_BUILD/superiso_inputs" \
  --compat-lib-dir "$HEPFIT_BUILD/compat/usr/lib/x86_64-linux-gnu" \
  --root-lib-dir "$HOME/root/lib" \
  --output-dir "$HEPFIT_BUILD/historical-replay-02"
python3 -m Flavor.validate_hepfit \
  --upstream "$HOME/zee" --driver "$HEPFIT_BUILD/driver" \
  --source "$HEPFIT_SOURCE" --build-receipt "$HEPFIT_BUILD/BUILD_RECEIPT.json" \
  --execution-receipt "$HEPFIT_BUILD/historical-replay-02/HEPFIT_EXECUTION.json" \
  --output "$HEPFIT_BUILD/RESTORED_BACKEND_02.json"
python3 -m Flavor.bsgamma \
  --upstream "$HOME/zee" --driver "$HEPFIT_BUILD/driver" \
  --hepfit-source "$HEPFIT_SOURCE" \
  --compat-lib-dir "$HEPFIT_BUILD/compat/usr/lib/x86_64-linux-gnu" \
  --root-lib-dir "$HOME/root/lib" \
  --attestation "$HEPFIT_BUILD/RESTORED_BACKEND_02.json" \
  --output-dir /tmp/zee-bsgamma-s00017-nlo-02
python3 -m unittest Flavor.test_bsgamma Flavor.test_validate_hepfit -v
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s "$HEPFIT_TASK" -p 'test_*.py' -v
```

The new attestation/output writers refuse existing paths. The historical
runner itself permits overwriting, so the explicit `mkdir` under `set -e`
guards its fresh output path. Do not point it at the frozen historical data.
The restore script reuses already hash-verified archives and CMake targets;
it does not repeat downloads or force complete recompilation.

Validation completed: **100 Zee-Pheno regression tests PASS**, plus **23 old
HEPfit receipt/uncertainty tests PASS**. CI checks archived actual output and
synthetic falsifiers, not a live ROOT/HEPfit build. The next physics task is
uncertainty and remaining model coverage, not assigning an exclusion from
these central numbers.
