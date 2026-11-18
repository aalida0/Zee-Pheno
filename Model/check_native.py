"""One-card ZeeTypeII SPheno input-propagation audit (not an LFV validation).

Runs the existing executable with corrected matrix row comments; never edits
an input card or existing output. Only standard-library Python is required.
"""
import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

MATRICES = ("FZEEIN", "IMFZEEIN", "YEE2IN", "IMYEE2IN")
TRACKED = set(MATRICES) | {"MINPAR", "FZEE", "IMFZEE", "YEE2", "IMYEE2", "YE", "IMYE"}


def blocks(text):
    """Parse only the blocks needed for the input/output matrix gate."""
    result, current = {}, None
    for line in text.splitlines():
        hdr = re.match(r"\s*BLOCK\s+(\w+)\b", line, re.I)
        if hdr:
            name = hdr[1].upper()
            current = name if name in TRACKED else None
            if current:
                if current in result:
                    raise ValueError(f"Repeated block {current}")
                result[current] = {}
            continue
        if re.match(r"\s*DECAY\b", line, re.I):
            current = None
            continue
        if current not in result:
            continue
        words = line.split("#", 1)[0].split()
        if not words:
            continue
        try:
            key = tuple(map(int, words[:-1]))
            val = float(words[-1].replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        if not key or not math.isfinite(val):
            raise ValueError(f"Invalid numeric line in {current}: {line}")
        if key in result[current] and current != "SPHENOINPUT":
            raise ValueError(f"Duplicate {current}{key}")
        result[current][key] = val
    return result


def mat(data, stem):
    indices = {(i, j) for i in range(1, 4) for j in range(1, 4)}
    for name in (stem, "IM" + stem):
        if set(data.get(name, {})) != indices:
            raise ValueError(f"Missing/incomplete complex 3x3 block {name}")
    return {(i, j): complex(data[stem][(i, j)], data["IM" + stem][(i, j)])
            for i, j in sorted(indices)}


def antisym(f):
    return max(abs(f[i, j] + f[j, i]) for i in range(1, 4) for j in range(1, 4))


def distance(a, b):
    return max(abs(a[k] - b[k]) for k in a)


def fix_card(text, update=None):
    """Add Fortran's required final string token; preserve all other rows."""
    update = update or {}
    current, seen = None, set()
    out = []
    for line in text.splitlines(keepends=True):
        header = re.match(r"\s*BLOCK\s+(\w+)\b", line, re.I)
        if header:
            current = header[1].upper()
        elif re.match(r"\s*DECAY\b", line, re.I):
            current = None
        if current not in (*MATRICES, "MINPAR") or header:
            out.append(line)
            continue
        raw, sep, comment = line.rstrip("\r\n").partition("#")
        fields = raw.split()
        n = 2 if current == "MINPAR" else 3
        if len(fields) != n:
            out.append(line)
            continue
        key = tuple(map(int, fields[:-1]))
        if (current, key) in update:
            fields[-1] = f"{update[(current, key)]:.17e}"
            seen.add((current, key))
        tail = (" #" + comment) if sep else " # generated entry"
        out.append("  " + " ".join(fields) + tail + ("\n" if line.endswith("\n") else ""))
    if set(update) != seen:
        raise ValueError(f"Missing probe inputs {set(update) - seen}")
    result = "".join(out)
    parsed = blocks(result)
    for name in MATRICES:
        if len(parsed.get(name, {})) != 9:
            raise ValueError(f"Native input {name} is not complete")
    if antisym(mat(parsed, "FZEEIN")) > 1e-10:
        raise ValueError("Input FZEEIN is not antisymmetric")
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate(inp, exe, output):
    spc = output.with_suffix(".spc")
    run = subprocess.run([str(exe), str(inp), str(spc)], text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         timeout=180)
    output.with_suffix(".log").write_text(run.stdout)
    if run.returncode or not spc.is_file():
        raise RuntimeError(f"SPheno failed for {inp}: exit={run.returncode}; see {output}.log")
    data = blocks(spc.read_text())
    return mat(data, "FZEE"), mat(data, "YEE2"), digest(spc)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--exe", type=Path, help="existing, already-built SPhenoZeeTypeII executable")
    ap.add_argument("--attestation", type=Path, help="receipt containing binary_sha256")
    args = ap.parse_args()
    source = args.input.read_text()
    original = blocks(source)
    f = mat(original, "FZEEIN")
    y = mat(original, "YEE2IN")
    if antisym(f) > 1e-10:
        raise ValueError("Source benchmark itself is not antisymmetric")
    if args.attestation and not args.exe:
        ap.error("--attestation requires --exe")
    if args.exe:
        if not args.exe.is_file():
            raise FileNotFoundError(args.exe)
        if not args.attestation:
            ap.error("--attestation is required for an attested native result")
        receipt = json.loads(args.attestation.read_text())
        if (receipt.get("binary_sha256") != digest(args.exe) or
                receipt.get("force_complete_rebuild") is not True or
                receipt.get("build_exit_code") != 0):
            raise ValueError("SPheno binary does not match a successful forced-build attestation")
    df = max(1e-6, abs(f[1, 2].real)*0.1)
    dy = max(1e-6, abs(y[2, 2].imag)*0.1)
    variants = {
        "baseline": {},
        "f_probe": {("FZEEIN", (1, 2)): f[1, 2].real + df,
                    ("FZEEIN", (2, 1)): f[2, 1].real - df},
        "ye2_probe": {("IMYEE2IN", (2, 2)): y[2, 2].imag + dy},
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for name, changes in variants.items():
        card = args.output_dir / (name + ".in")
        card.write_text(fix_card(source, changes))
    report = {"status": "NATIVE_PENDING" if not args.exe else "NATIVE_FAIL",
              "input_sha256": digest(args.input), "delta_f": df, "delta_imye2": dy,
              "binary_sha256": digest(args.exe) if args.exe else None,
              "native_executed": bool(args.exe)}
    if args.exe:
        results = {}
        for name in variants:
            inp = args.output_dir / (name + ".in")
            results[name] = evaluate(inp, args.exe, args.output_dir / name)
        fb, yb, _ = results["baseline"]
        ff, yf, _ = results["f_probe"]
        fy, yy, _ = results["ye2_probe"]
        fscale = max(1e-12, *(abs(z) for z in f.values()))
        yscale = max(1e-12, *(abs(z) for z in y.values()))
        checks = {
            "baseline_F_relative": distance(fb, f) / fscale,
            "baseline_Ye2_relative": distance(yb, y) / yscale,
            "F_antisymmetric_max_abs": antisym(fb),
            "F_probe_output_delta": distance(ff, fb),
            "F_probe_unintended_Ye2_delta": distance(yf, yb),
            "Ye2_probe_output_delta": distance(yy, yb),
            "Ye2_probe_induced_F_delta": distance(fy, fb),
        }
        # Off-target responses are diagnostics: genuine RG mixing is possible.
        passed = (checks["baseline_F_relative"] <= 5e-4 and
                  checks["baseline_Ye2_relative"] <= 5e-4 and
                  checks["F_antisymmetric_max_abs"] <= 5e-9 and
                  checks["F_probe_output_delta"] >= 0.5*df and
                  checks["Ye2_probe_output_delta"] >= 0.5*dy)
        report.update(status="NATIVE_INPUT_GATE_PASS" if passed else "NATIVE_INPUT_GATE_FAIL",
                      checks=checks,
                      output_sha256={k: v[2] for k, v in results.items()})
    (args.output_dir / "native-input-report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] in ("NATIVE_PENDING", "NATIVE_INPUT_GATE_PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
