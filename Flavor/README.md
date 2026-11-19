# Flavor — physical constraint evaluation

## Source-bound SPheno → HEPfit B→Xsγ bridge

The executable module [`bsgamma.py`](bsgamma.py) reuses the validated
**independent BMU NLO-QCD matching** and patched **HEPfit bridge client** from
`aalida0/zee` (`research/zee-full-constraints`). No new NLO kernels are
invented or copied, and the old FlavorKit rate polynomial is unchanged.

This adapter is for the two **fixed fitted native controls** in
`Scan/reference/native_benchmark/s00017_{NO,IO}_benchmark`, not yet an
arbitrary-card production scanner. It:

- verifies each spectrum/card/native-snapshot association and reconstructs
  charged **pole** masses, native `tanβ` (`MINPAR[10]`), the input soft mass,
  and squared doublet projections of both charged physical eigenstates from
  `CHARGEMIX`. The spectrum's `HMIX[10]` is **β in radians**, not tanβ;
  these are independently cross-checked against the snapshot's `vd,vu`;
- locks the existing BMU matching source SHA-256
  `0c2dfb32a7138b02fa4a52c163c63f942f313bb5d697545e138ef1915c200c4e`,
  the previous HEPfit execution receipt, its patched driver source and exact
  executable SHA-256. HEPfit source commit is
  `322bd73a2dadd94ce94a6d71cb79d616e569120a`;
- executes, for each ordering, a **same-Higgs-mass HEPfit SM control** and
  **zero-NP bridge**. Both must agree. It obtains the **running top mass at
  μ=160 GeV from this same HEPfit executable**; no SPheno running top mass
  is silently mixed into the source-validated HEPfit normalization;
- computes LO ordinary `C7,C8` and unweighted NLO ordinary `C4,C7,C8`
  for **each charged state**, weights them **coherently at the Wilson level**
  by its squared native doublet projection, and uses the HEPfit client to
  perform the BMU→effective conversion and NLO evolution exactly once;
- independently verifies the actual high-scale order-tagged coefficient
  values reported by the executed HEPfit client, extracts **`BR_NLO`**,
  and archives source hashes, raw input/log files, the same-scheme SM
  baseline and the original FlavorKit 200 value **separately**.

The `BR_stock_NNLO` output is a known mixed-order diagnostic and is **never**
used as a Zee NNLO rate. The HEPfit central rate retains its native HQE/partial
EW prescription; it is not a uniformly truncated, full Zee calculation.
The old validation covers real top-only NLO NP matching, with omitted u/c,
primed/light-quark, Zee EW and Wilson-dependent nonlocal/resolved-photon
uncertainties. `theory_sigma=null`, `experimental_exclusion=null`, and
`globally_allowed=null` are **mandatory** for these runs.

### Run on the local machine

The separate historical `~/zee` checkout and compiled patched driver must
exist; GitHub CI executes native-artifact and synthetic contracts, not HEPfit.
Build instructions for a missing driver are in
[`research/full_constraints/bsgamma_hepfit/REPRODUCE.md`](https://github.com/aalida0/zee/blob/research/zee-full-constraints/research/full_constraints/bsgamma_hepfit/REPRODUCE.md).

~~~bash
cd ~/Zee-Pheno
git fetch origin
git switch research/spheno-hepfit-bsgamma
python3 -m unittest Flavor.test_bsgamma -v

# Source-bound, built HEPfit client from the previous project stage:
test -x /tmp/zee-hepfit/driver
sha256sum /tmp/zee-hepfit/driver
# Expected SHA256: 2e622645491118041361607faf040120fec59103cf0dfd85dbb6925de5fee0a3

python3 -m Flavor.bsgamma \
  --upstream ~/zee \
  --driver /tmp/zee-hepfit/driver \
  --hepfit-source /tmp/zee-hepfit/HEPfit-322bd73a2dadd94ce94a6d71cb79d616e569120a \
  --compat-lib-dir /tmp/zee-hepfit/compat/usr/lib/x86_64-linux-gnu \
  --root-lib-dir ~/root/lib \
  --output-dir /tmp/zee-bsgamma-s00017-nlo-01
~~~

Read `/tmp/zee-bsgamma-s00017-nlo-01/bsgamma_nlo.json`. The runner refuses
an already existing output directory. If the earlier HEPfit scratch was
removed, **restore the pinned build first**; a fabricated central value or
an unrelated HEPfit executable is not an acceptable substitute.

Source for the archived FlavorKit polynomial:
`SPheno-4.0.5/ZeeTypeII/Observables/FlavorKit_Observables_ZeeTypeII.f90:1073–1109`.
Its fixed `NNLO_SM=3.15` (×10⁻⁴), Lunghi–Matias response and Wilson
normalization are not the HEPfit NLO calculation. The old output entry
`FlavorKitQFV[200]` is retained only as a provenance-labeled diagnostic.

## Other flavor constraints

Full meson and B physics likelihoods, NP-dependent theoretical uncertainties,
Higgs and LFV constraints remain **UNKNOWN**, even when the above central
calculation executes. See the historical
[coverage ledger](https://github.com/aalida0/zee/blob/research/zee-full-constraints/research/full_constraints/BSGAMMA_COVERAGE.md).
