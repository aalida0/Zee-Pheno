"""One-time SOURCE-MATCHED SPheno 510/530 tree-tadpole audit.

The native loop-corrected SPheno spectrum has ZEEPOT m11/m22 which
do NOT satisfy tree neutral tadpoles. SARAH/SPheno option 530=1 writes
tree-level tadpole values INSTEAD of loop-corrected values, while
510=1 writes both source tree values for Vevacious.

This runner changes only these flags in frozen input cards, invokes
the original SHA-pinned SPhenoZeeTypeII executable, retains raw logs,
spectrum sha hashes and rechecks exact potential parameters. Never
interpret a missing TREEZEEPOT or changed physics spectrum as a PASS.
Requires the authenticated original SPheno executable on the machine.

The real native loop-effective global vacuum remains UNKNOWN; this
validates the TREE SCHEME used in the analytic sufficient certificate.
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

from ZeeHiggs.spheno_audit import NATIVE,read_slha,sha256
from Scan.zee_chargebreaking_vacuum import source_inputs_from_blocks
from Scan.zee_vacuum_global_sufficient import from_native_point
from Scan.theory_gate_native import ORIGINAL_BINARY_SHA256


def patch_card(text,flag):
    if flag not in (510,530):
        raise ValueError("Only separate tree-tadpole output flags permitted")
    lines=text.splitlines(keepends=True)
    start=next((i for i,l in enumerate(lines) if
                re.match(r"^\s*Block\s+SPhenoInput\b",l,re.I)),None)
    if start is None:
        raise ValueError("Original SPhenoInput block absent")
    stop=next((i for i in range(start+1,len(lines)) if
               re.match(r"^\s*Block\s+",lines[i],re.I)),len(lines))
    hits=[]
    for i in range(start+1,stop):
        stripped=lines[i].split("#",1)[0].split()
        if len(stripped)>=2 and stripped[0].isdigit() and int(stripped[0]) in (510,530):
            hits.append(i)
    if len(hits)!=2:
        raise ValueError("Expected both original tree input flags")
    numbers={int(lines[i].split()[0]) for i in hits}
    if numbers!={510,530}:
        raise ValueError("Ambiguous option numbering")
    for i in hits:
        key=int(lines[i].split()[0])
        old_comment=lines[i].split("#",1)
        newline="\n" if lines[i].endswith("\n") else ""
        lines[i]=f" {key} {int(key==flag)}  # source-bound tree tadpole audit, only {flag} enabled{newline}"
    return "".join(lines)


def read_reference(original_order):
    if original_order not in ("NO","IO"):
        raise ValueError("Invalid ordering")
    root=NATIVE/f"s00017_{original_order}_benchmark"
    inp=root/"input.in"
    native=json.loads((root/"native.json").read_text())
    if sha256(inp)!=native["input_sha256"] or native["production_binary_sha256"]!=ORIGINAL_BINARY_SHA256:
        raise ValueError("Original input source not authenticated")
    return inp,root


def run_one(exe,order,flag,output_dir):
    card,ref=read_reference(order)
    label=f"s00017_{order}_tree{flag}"
    d=output_dir/label
    d.mkdir()
    generated=d/"input.in"
    generated.write_text(patch_card(card.read_text(),flag))
    out=d/"spectrum.slha"
    command=[str(exe.resolve()),str(generated.resolve()),str(out.resolve())]
    try:
        r=subprocess.run(command,cwd=d,text=True,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                         timeout=300,check=False)
    except subprocess.TimeoutExpired as e:
        (d/"run.log").write_text(str(e)+"\nTIMEOUT_AFTER_300_SECONDS\n")
        raise ValueError("Original SPheno tree-tadpole audit timed out")
    (d/"run.log").write_text(r.stdout)
    if (r.returncode!=0 or not out.is_file() or out.stat().st_size<=1000 or
        "TerminateProgram" in r.stdout):
        raise ValueError(f"Source-matched SPheno tree tadpole {flag} failed; inspect {d}/run.log")
    blocks,_=read_slha(out.read_text())
    p=source_inputs_from_blocks(blocks)
    answer=from_native_point(p)
    # Flag 530 replaces loop tadpole parameters with tree ones;
    # flag 510 should be treated as a separate origin test, not as a
    # guarantee that unchanged ZEEPOT has become a tree output.
    expected=answer["matched_native_potential"]
    if flag==530 and not expected.get("native_potential_certified"):
        raise ValueError(
            "SPheno 530=1 did not produce stationarity-matched tree values; "
            "tree global vacuum is NOT certified for this generated output")
    if flag==510 and not any(k.startswith("TREE") for k in blocks):
        raise ValueError("SPheno 510=1 generated no explicit TREE blocks")
    return {
       "point_id":label,
       "tadpole_option":flag,
       "original_card_sha256":sha256(card),
       "generated_card_sha256":sha256(generated),
       "original_source_binary_sha256":ORIGINAL_BINARY_SHA256,
       "executed_spectrum_sha256":sha256(out),
       "raw_log_sha256":sha256(d/"run.log"),
       "reported_tree_blocks":sorted(k for k in blocks if k.startswith("TREE")),
       "matched_native_tadpole_certificate":expected,
       "tree_consistent_reference_certificate":answer["defined_tree_tadpole_potential"],
       "on_this_distinct_output_tree_vacuum_certified":bool(
            expected.get("native_potential_certified") and flag==530),
       "one_loop_effective_vacuum_global":"UNKNOWN",
    }


def run(exe,output_dir):
    if sha256(exe)!=ORIGINAL_BINARY_SHA256:
        raise ValueError("Wrong SPheno executable: original production SHA mismatch")
    output_dir.mkdir(parents=True,exist_ok=False)
    result={"status":"ACTUAL_NATIVE_SPHENO_TREE_TADPOLE_SOURCE_CHECK",
            "binary_sha256":ORIGINAL_BINARY_SHA256,
            "cases":[],"full_one_loop_vacuum":"UNKNOWN",
            "globally_allowed":None}
    for ordering in ("NO","IO"):
        for flag in (510,530):
            result["cases"].append(run_one(exe,ordering,flag,output_dir))
    (output_dir/"native-tree-audit.json").write_text(
        json.dumps(result,indent=2,allow_nan=False)+"\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    r=run(a.exe.resolve(),a.output_dir.resolve())
    print("SOURCE_ATTESTED_TREE_TADPOLE_VALIDATIONS",
          len(r["cases"]),"one_loop_vacuum=UNKNOWN")


if __name__=="__main__":main()
