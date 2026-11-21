"""Source-bound SPheno tree-tadpole output diagnostic (options 510/530).

510 requests tree-level tadpole solutions; 530 requests Vevacious
tree/loop output blocks. Neither flag alone guarantees a particular
block name in a custom SARAH-generated model. Only numerically matched
tree-stationary potential parameters can support a tree certificate.

Use the original SHA-pinned executable and archived NO/IO cards.
Unmatched or absent source tree values remain UNKNOWN, never PASS.
The one-loop effective-potential global vacuum is not certified here.
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
            if stripped[1] not in ("0","1","0.","1."):
                raise ValueError("Original tree option has unexpected value")
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



def evaluate_tree_output(blocks):
    """Check actual source values, never infer them from an option number."""
    native=from_native_point(source_inputs_from_blocks(blocks))
    candidates={"ZEEPOT":native["matched_native_potential"]}
    parse_errors={}
    if "TREEZEEPOT" in blocks:
        try:
            tree_blocks={**blocks,"ZEEPOT":blocks["TREEZEEPOT"]}
            candidates["TREEZEEPOT"]=from_native_point(
                source_inputs_from_blocks(tree_blocks))["matched_native_potential"]
        except (ValueError,KeyError,IndexError,TypeError) as exc:
            parse_errors["TREEZEEPOT"]=str(exc)
    stationarity=[name for name,result in candidates.items()
                  if result.get("source_matching")=="MATCHED_SOURCE_TREE_TADPOLES"]
    certified=[name for name in stationarity
               if candidates[name].get("native_potential_certified") is True]
    return {
        "status":("SOURCE_MATCHED_GLOBAL_TREE_CERTIFICATE" if certified else
                  "SOURCE_MATCHED_TREE_STATIONARITY_ONLY" if stationarity else
                  "UNKNOWN_NO_MATCHED_NATIVE_TREE_TADPOLES"),
        "reported_tree_blocks":sorted(k for k in blocks if k.startswith("TREE")),
        "reported_loop_blocks":sorted(k for k in blocks if k.startswith("LOOP")),
        "source_potential_stationarity":candidates,
        "tree_block_parse_errors":parse_errors,
        "matched_stationary_source_blocks":stationarity,
        "certified_global_tree_source_blocks":certified,
        "matched_native_tadpole_certificate":native["matched_native_potential"],
        "tree_consistent_reference_certificate":native["defined_tree_tadpole_potential"],
        "on_this_distinct_output_tree_vacuum_certified":bool(certified),
        "original_loop_corrected_global_vacuum":"UNKNOWN",
        "one_loop_effective_vacuum_global":"UNKNOWN",
    }


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
    return {
       "point_id":label,
       "tadpole_option":flag,
       "original_card_sha256":sha256(card),
       "generated_card_sha256":sha256(generated),
       "original_source_binary_sha256":ORIGINAL_BINARY_SHA256,
       "executed_spectrum_sha256":sha256(out),
       "raw_log_sha256":sha256(d/"run.log"),
       **evaluate_tree_output(blocks),
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
