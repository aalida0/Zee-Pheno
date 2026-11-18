"""Read-only audit of native Zee neutrino loop-mass availability, not validation."""
import argparse
import hashlib
import json
import re
from pathlib import Path

NEUTRINO_PDG = (12, 14, 16)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_entries(text, block):
    rows = []
    active = False
    for line in text.splitlines():
        hit = re.match(r"\s*BLOCK\s+(\w+)\b", line, re.I)
        if hit:
            active = hit[1].upper() == block.upper()
            continue
        if re.match(r"\s*DECAY\b", line, re.I):
            active = False
        if not active:
            continue
        fields = line.split("#", 1)[0].split()
        if len(fields) < 2:
            continue
        try:
            key = int(fields[0])
            val = float(fields[1].replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        rows.append((key, val))
    return rows


def audit(model, card, spc):
    for path in (model, card, spc):
        if not path.exists():
            raise FileNotFoundError(path)
    flags = [v for k, v in parse_entries(card.read_text(), "SPHENOINPUT") if k == 55]
    if any(v not in (0., 1.) for v in flags):
        loop_state = "UNRECOGNIZED_FLAG55"
    elif not flags:
        loop_state = "FLAG55_NOT_EXPLICIT"
    elif len(set(flags)) > 1:
        loop_state = "CONFLICTING_FLAG55"
    else:
        loop_state = "ENABLED" if flags[0] == 1 else "DISABLED"

    spectra = {k: v for k, v in parse_entries(spc.read_text(), "MASS")
               if k in NEUTRINO_PDG}
    sources = sorted(p for p in model.rglob("*.f90")
                     if p.name.startswith(("LoopMasses", "SPheno", "TreeLevelMasses")))
    pattern = re.compile(
        r"\b\w*(?:SelfEnergy\w*Fv|Sigma\w*Fv|Mass\w*Fv|Fv\w*Mass|Fv\w*SelfEnergy)\w*\b",
        re.I)
    matches = []
    for path in sources:
        for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if pattern.search(line):
                matches.append(f"{path.relative_to(model)}:{n}: {line.strip()[:180]}")
    status = ("LOOP_MASSES_DISABLED" if loop_state == "DISABLED" else
              "SOURCE_OR_OUTPUT_INSPECTION_REQUIRED")
    return {
        "status": status,
        "native_one_loop_neutrino_mass_validated": False,
        "card_sha256": sha(card), "spectrum_sha256": sha(spc),
        "spheno_input_55_occurrences": flags,
        "loop_mass_switch": loop_state,
        "neutrino_mass_entries_GeV": {str(k): v for k, v in spectra.items()},
        "loop_mass_source_files": [str(p.relative_to(model)) for p in sources],
        "source_lines_to_inspect": matches[:35],
        "next": "Check the actual Fv Majorana self-energy and call site against the analytic Zee mass, not merely the source symbol."
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--card", type=Path, required=True)
    p.add_argument("--spectrum", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(audit(args.model, args.card, args.spectrum), indent=2))


if __name__ == "__main__":
    main()
