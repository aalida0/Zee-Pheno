"""Two one-card Zee controls: F=0 and F*2, existing loop-enabled SPheno."""
import argparse
import json
import re
import subprocess
from pathlib import Path

from Model.audit_neutrino_native import parse_entries, sha
from Model.check_native import blocks, mat, antisym
from Model.probe_neutrino_loop import BINARY_SHA

LOOP_CARD_SHA = "e83f315928a3906bc14e48aec39b2bd842abf66949f53b6cd78824d4c45a8e07"


def rescale_f(source, scale):
    if scale not in (0, 2):
        raise ValueError("Only the zero and factor-two controls are supported")
    current = None
    changed = {name: 0 for name in ("FZEEIN", "IMFZEEIN")}
    out = []
    for line in source.splitlines(keepends=True):
        hdr = re.match(r"\s*BLOCK\s+(\w+)\b", line, re.I)
        if hdr:
            current = hdr[1].upper()
        elif re.match(r"\s*DECAY\b", line, re.I):
            current = None
        if current in changed and not hdr and line.strip() and not line.lstrip().startswith("#"):
            match = re.fullmatch(
                r"([ \t]*\d+[ \t]+\d+[ \t]+)([-+.\dEeDd]+)([ \t]+#.*)?(\r?\n)?",
                line)
            if match is None:
                raise ValueError("Unexpected native matrix row: " + repr(line))
            number = float(match[2].replace("D", "E").replace("d", "e"))
            line = match[1] + f"{scale*number:.17e}" + (match[3] or " # probe entry") + (match[4] or "")
            changed[current] += 1
        out.append(line)
    if changed != {"FZEEIN": 9, "IMFZEEIN": 9}:
        raise ValueError("Expected exactly 9 full complex Zee matrix entries per block")
    result = "".join(out)
    baseline, modified = blocks(source), blocks(result)
    for name in ("FZEEIN", "IMFZEEIN"):
        for key, value in baseline[name].items():
            if abs(modified[name][key]-scale*value)>1e-15:
                raise ArithmeticError("Rescaled Zee matrix does not match source")
    if antisym(mat(modified, "FZEEIN")) > 1e-12:
        raise ValueError("Rescaled Fzee is not antisymmetric")
    return result


def probe(input_card, exe, baseline_spectrum, output_dir):
    if sha(input_card) != LOOP_CARD_SHA or sha(exe) != BINARY_SHA:
        raise ValueError("Only the attested loop-enabled card and executable are supported")
    source = input_card.read_text()
    baseline_masses = dict(parse_entries(baseline_spectrum.read_text(), "MASS"))
    if any(baseline_masses.get(pdg, 0.) <= 0 for pdg in (12, 14, 16)):
        raise ValueError("Missing original loop-enabled neutrino pole masses")
    output_dir.mkdir(parents=True, exist_ok=False)
    masses = {}
    for label, factor in (("zero_F", 0), ("double_F", 2)):
        card = output_dir / (label + ".in")
        spectrum = output_dir / (label + ".spc")
        logfile = output_dir / (label + ".log")
        card.write_text(rescale_f(source, factor))
        result = subprocess.run(
            [str(exe.resolve()), str(card.resolve()), str(spectrum.resolve())],
            cwd=output_dir.resolve(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, timeout=300)
        logfile.write_text(result.stdout)
        if result.returncode or not spectrum.exists():
            raise RuntimeError(f"SPheno {label} returned {result.returncode}; inspect {logfile}")
        parsed = dict(parse_entries(spectrum.read_text(), "MASS"))
        if not all(pdg in parsed for pdg in (12, 14, 16)):
            raise ValueError(f"Missing neutrino pole masses for {label}")
        masses[label] = sorted(parsed[pdg]*1e9 for pdg in (12, 14, 16))
    baseline = sorted(baseline_masses[pdg]*1e9 for pdg in (12, 14, 16))
    info = {
        "status": "F_SCALING_CONTROL_EXECUTED_PHYSICS_UNVALIDATED",
        "loop_card_sha256": LOOP_CARD_SHA,
        "exe_sha256": BINARY_SHA,
        "baseline_eV": baseline,
        "zero_F_eV": masses["zero_F"],
        "double_F_eV": masses["double_F"],
        "zero_F_over_baseline_max": max(abs(x) for x in masses["zero_F"])/max(baseline),
        "double_F_over_2baseline": [x/(2*y) for x,y in zip(masses["double_F"], baseline)],
        "note": "Tests lepton-number limit and linearity; not a validation of Majorana matrix or PMNS."
    }
    (output_dir / "result.json").write_text(json.dumps(info, indent=2)+"\n")
    return info


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(probe(a.input, a.exe, a.baseline, a.output_dir), indent=2))


if __name__ == "__main__":
    main()
