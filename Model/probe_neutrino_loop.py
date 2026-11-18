"""Run one attested SPheno card with loop masses on; no physics PASS claim."""
import argparse
import json
import re
import subprocess
from pathlib import Path
from Model.audit_neutrino_native import parse_entries, sha

BASELINE_SHA = "f2efe759aa3196978d466739fc788b1dd7fb01d71cbc6e883c8a87830da454bf"
BINARY_SHA = "2c1a73c797b2f21b3266f3558dae3196a6b2514a5fd3d2e05fc83b0ce918af52"


def enable_loop_masses(source):
    current = None
    changed = 0
    output = []
    for line in source.splitlines(keepends=True):
        header = re.match(r"^[ \t]*BLOCK[ \t]+(\w+)\b", line, re.I)
        if header:
            current = header[1].upper()
        elif re.match(r"^[ \t]*DECAY\b", line, re.I):
            current = None
        if current == "SPHENOINPUT" and not header:
            body = line.rstrip("\r\n")
            suffix = line[len(body):]
            match = re.match(r"^([ \t]*55[ \t]+)([^ \t#\r\n]+)(.*)$", body)
            if match:
                if float(match[2].replace("D", "E").replace("d", "E")) != 0.:
                    raise ValueError("SPhenoInput[55] is not zero: refusing to alter card")
                line = match[1] + "1" + match[3] + suffix
                changed += 1
        output.append(line)
    result = "".join(output)
    if changed != 2 or [v for k, v in parse_entries(result, "SPHENOINPUT") if k == 55] != [1., 1.]:
        raise ValueError("Expected exactly two disabled SPhenoInput[55] entries")
    return result


def probe(card, exe, directory):
    if sha(card) != BASELINE_SHA or sha(exe) != BINARY_SHA:
        raise ValueError("Source card or executable SHA-256 mismatch: refusing to run")
    changed = enable_loop_masses(card.read_text())
    directory.mkdir(parents=True, exist_ok=False)
    new_card = directory / "loop_enabled.in"
    spectrum = directory / "loop_enabled.spc"
    log = directory / "loop_enabled.log"
    new_card.write_text(changed)
    outcome = subprocess.run([str(exe.resolve()), str(new_card.resolve()), str(spectrum.resolve())],
                             cwd=directory.resolve(), stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, timeout=300)
    log.write_text(outcome.stdout)
    masses = {str(k): v for k, v in parse_entries(spectrum.read_text(), "MASS")
              if k in (12, 14, 16)} if spectrum.exists() else {}
    result = {"status": "NATIVE_LOOP_RUN_FAILED" if outcome.returncode or not spectrum.exists() else
              "NEUTRINO_POLE_MASS_NONZERO_UNVALIDATED" if any(abs(v) > 0 for v in masses.values()) else
              "NEUTRINO_MASS_ZERO_OR_MISSING",
              "native_one_loop_neutrino_mass_validated": False,
              "baseline_input_sha256": BASELINE_SHA, "binary_sha256": BINARY_SHA,
              "loop_input_sha256": sha(new_card),
              "loop_output_sha256": sha(spectrum) if spectrum.exists() else None,
              "exit_code": outcome.returncode, "neutrino_mass_entries_GeV": masses,
              "model_note": "Pole mass output is not a substitute for a matched complex Majorana self-energy."}
    (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--card", required=True, type=Path)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(probe(args.card, args.exe, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
