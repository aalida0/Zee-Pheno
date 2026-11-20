# Executed HiggsTools pilot — native s00017 NO/IO

**Date:** 20 November 2026. **Status:** OFFICIAL EXPERIMENTAL TOOL EXECUTED;
FULL ZEE-MODEL HIGGS CONSTRAINT UNKNOWN.

The source-qualified native SPheno spectra for s00017 NO and IO were
passed through the unmodified upstream HiggsTools v1.3 (HiggsBounds 6.1 /
HiggsSignals 3.1). No 2HDMC or ScannerS was used.

Evidence:
- Official executed GitHub Actions workflow:
  https://github.com/aalida0/Zee-Pheno/actions/runs/44732975159
- Full machine-readable result artifact:
  https://github.com/aalida0/Zee-Pheno/actions/runs/44732975159/artifacts/14627701307
- Audited SPheno input module: ../ZeeHiggs/spheno_audit.py
- Official HiggsTools API adapter: ../ZeeHiggs/higgstools_pilot.py
- Input limitations and local reproduction: README.md

## HiggsSignals results (same HS dataset and same h125 mass)

| Quantity | NO | IO |
|---|---:|---:|
| SM reference chi2, evaluated at the point's own mh | 151.55556002183712 | 151.55556065245997 |
| Zee chi2, effective-coupling h125 BR prescription | 163.01017406056917 | 166.7511144599829 |
| Zee - SM delta chi2 (effective-coupling BR) | **+11.454614038732046** | **+15.195553807522913** |
| Zee chi2, explicit native SPheno branching ratios | 178.13293606131302 | 181.55994766632205 |
| Zee - SM delta chi2 (explicit SPheno BR) | **+26.577376039475894** | **+30.004387013862072** |

Both pairs of calculations were performed with the same official
HiggsSignals v1.1 dataset and a SM-like reference Higgs at the respective
native mh=125.25 GeV. The reference calculation uses the official
HiggsTools SM-like effective couplings, not a hand-chosen chi2 baseline.

Delta chi2 is a comparison diagnostic, **not** by itself a justified
frequentist 95% confidence exclusion. The implied degree-of-freedom and
likelihood-ratio interpretation, theory errors and data coverage have not
been established. In particular the marked difference between explicit
SPheno BRs and effective-coupling-recomputed h125 BRs is material.
Do not assign a p-value or declare either benchmark excluded.

## HiggsBounds results for supplied, incomplete model predictions

The official HiggsBounds dataset returned allowed=true for both native
controls. Its selected limits only concern neutral Higgs bosons:

| PDG | Selected topology | Citation | NO observed/expected ratios | IO observed/expected ratios |
|---|---|---|---:|---:|
| 25 (h125) | LHC8 associated HW/HZ -> bb | ATLAS arXiv:1409.6212 | 0.7298757462 / 1.1735638705 | 0.7299803928 / 1.1737360121 |
| 35 (H) | LHC13 H -> (tautau)(bb) via intermediate scalars | ATLAS arXiv:2209.10910 | 0.2094020918 / 0.2692624444 | 0.2088824150 / 0.2685950875 |
| 36 (A) | LHC13 H -> Z+(scalar -> bb) | ATLAS arXiv:2207.00230 | 0.0554104584 / 0.0496001369 | 0.0553146098 / 0.0495146979 |

No selected limit for either charged state (PDG 37 or 9000037)
was recorded. This does not imply that the official dataset has no
charged-Higgs analyses at all: it means none was selected for the
currently supplied input set. Charged-Higgs production cross sections
and dedicated H2+ -> H1+ Z experimental topology/acceptance coverage
are missing. HiggsBounds allowed=true must not be promoted to a
full Zee-model PASS.

## Source and dataset provenance

- Spectrum: SARAH 4.15.4 / SPheno 4.0.5, both original native
  NO/IO cards and their authenticated source hashes, physical h,H,A,
  H1+,H2+ masses, widths and branching ratios.
- HiggsTools official upstream v1.3.
- HiggsBounds data: tag v1.7, commit
  3539f9adf3724f8752701208bf059d1658824cd2,
  310 JSON files, manifest SHA256
  a2ca0729f3c2d1d5bc2769177978f683e96a81c7e9d81cf93c77781cb070adf5.
- HiggsSignals data: tag v1.1, commit
  3c6c2538d9ba0549ea0559cf9b38c18eb80c3af1,
  22 JSON files, manifest SHA256
  d84ad943f9fb05f50f11f3861eab721ad04c6c02866191db733fedce65476d68.

The HiggsSignals v1.1 measurements date from an older dataset release;
the 2026 HiggsTools executable alone does NOT imply complete 2026
experimental measurement coverage. Neutral 13 TeV ggF rates are SPheno
approximate rescaled SusHi-1.5.0 estimates; full production theory is
not yet source-validated. SPheno labels its neutral Z-gamma loop
couplings as uncalculated. The actual imported HiggsTools dictionary
exposed no Zgamma keys to remove; their absence is still missing input,
not a validated zero.

## Physics gates and priority

- Historical SPheno source/decay audit: PASS.
- Official HiggsTools runtime and observable output: PASS.
- Neutral HiggsBounds numerical allowed for *provided input set*: TRUE.
- Physical complete Zee-model HiggsBounds acceptance: UNKNOWN.
- HiggsSignals validated p-value/95% exclusion: UNKNOWN.
- Global Zee benchmark acceptance: UNKNOWN.

Immediate priority: understand the h125 BR/total-width treatment and
identify dominant HiggsSignals chi2 components, since that ambiguity
changes chi2 by about 15 for NO. Next validate neutral production
schemes, supply both charged-Higgs production cross sections, audit
charged analysis/topology coverage, and test updated HiggsSignals data.
The aim is a defensible Higgs-sector screen for the phenomenology scan,
not an unjustified hard veto.
