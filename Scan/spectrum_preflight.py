"""Run a source-bound SPheno spectrum preflight, not a physical scan gate.

SPheno can terminate early with exit code zero. Retain all logs and spectrum
files; only complete, positive scalar MASS entries count as spectrum-present.
"""
import argparse
import json
import subprocess
from pathlib import Path

from Model.audit_neutrino_native import parse_entries, sha
from Model.diagnose_neutrino_native import BINARY_SHA

SCALAR_PDG = (25, 35, 36, 37, 9000037)


def preflight(proposals_path, executable, output_dir):
    if sha(executable) != BINARY_SHA:
        raise ValueError("Preflight requires the attested original SPheno executable")
    catalog = json.loads(proposals_path.read_text())
    if (catalog.get("schema_version") != 2 or
            catalog.get("status") != "SCALAR_NEUTRINO_PROPOSALS_NATIVE_UNKNOWN"):
        raise ValueError("Unrecognized candidate catalog")
    cases = []
    seen = set()
    for scalar in catalog["outputs"]:
        for family in scalar["families"]:
            if family.get("status") != "NATIVE_REVALIDATION_REQUIRED":
                continue
            name = family["card"]
            if name in seen or Path(name).name != name or not name.endswith(".in"):
                raise ValueError("Duplicate or unsafe scan card name")
            seen.add(name)
            card = proposals_path.parent / "cards" / name
            if sha(card) != family["card_sha256"]:
                raise ValueError("Input card changed since proposal generation")
            cases.append((scalar["scalar_id"], family["ordering"], card, family["card_sha256"]))
    if not cases:
        raise ValueError("No pending native candidate cards")
    output_dir.mkdir(parents=True, exist_ok=False)
    rows = []
    for scalar_id, ordering, card, card_sha in cases:
        case = output_dir / card.stem
        case.mkdir()
        spectrum = case / "spectrum.spc"
        try:
            run = subprocess.run(
                [str(executable), str(card.resolve()), str(spectrum)],
                cwd=case, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=300, check=False)
            log = run.stdout
            code = run.returncode
        except subprocess.TimeoutExpired as exc:
            log = (exc.stdout or b"").decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            log += "\nTIMEOUT_AFTER_300_SECONDS\n"
            code = None
        (case / "run.log").write_text(log)
        masses = {}
        if spectrum.is_file() and spectrum.stat().st_size:
            try:
                masses = dict(parse_entries(spectrum.read_text(), "MASS"))
            except (ValueError, UnicodeError):
                masses = {}
        complete = all(k in masses and masses[k] > 0 for k in SCALAR_PDG)
        status = ("SPECTRUM_PRESENT_UNVALIDATED"
                  if code == 0 and complete and "TerminateProgram" not in log
                  else "NATIVE_SPECTRUM_FAILED_OR_INCOMPLETE")
        rows.append(dict(scalar_id=scalar_id, ordering=ordering,
                         card=card.name, card_sha256=card_sha,
                         status=status, exit_code=code,
                         scalar_pole_masses_GeV={str(k): masses.get(k) for k in SCALAR_PDG},
                         log=str(case / "run.log"),
                         spectrum=str(spectrum) if spectrum.is_file() else None))
        print(f"{card.name}: {status}", flush=True)
    report = dict(status="NATIVE_SCALAR_SPECTRUM_PREFLIGHT_ONLY",
                  executable_sha256=sha(executable),
                  catalog_sha256=sha(proposals_path),
                  n_cases=len(rows),
                  n_spectrum_present=sum(x["status"]=="SPECTRUM_PRESENT_UNVALIDATED" for x in rows),
                  n_failed=sum(x["status"]=="NATIVE_SPECTRUM_FAILED_OR_INCOMPLETE" for x in rows),
                  source_note="Spectrum presence is NOT neutrino/LFV/Higgs/vacuum validation.",
                  cases=rows)
    (output_dir / "preflight.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proposals", type=Path, required=True)
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    a = ap.parse_args()
    report = preflight(a.proposals.resolve(), a.exe.resolve(), a.output_dir.resolve())
    print(json.dumps({k:v for k,v in report.items() if k!="cases"}, indent=2))


if __name__ == "__main__":
    main()
