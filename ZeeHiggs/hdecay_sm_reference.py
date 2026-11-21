"""Parse independent source-pinned HDECAY SM br.sm1/br.sm2.

HDECAY is used ONLY in its original SM mode as an independent
partial-width benchmark.  No 2HDM/ScannerS/2HDMC or fake Zee widths.
The HDECAY WW/ZZ 4f treatment may differ from the HXSWG
Prophecy4f reference even when inputs and masses match.
"""
import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

COMMIT="4884213286e1edbdf3b2bb9fe30095b29e336dff"
C1=("bb","tautau","mumu","ss","cc","tt")
C2=("gg","gamgam","Zgam","WW","ZZ")
MASS={"NO":125.250000,"IO":125.250135}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_table(path, columns, mass):
    rows=[]
    for raw in Path(path).read_text().splitlines():
        values=[]
        for token in raw.split():
            try:
                v=float(token.replace("D","E").replace("d","e"))
                if not math.isfinite(v):
                    raise ValueError("Nonfinite HDECAY result")
                values.append(v)
            except ValueError:
                values=[]
                break
        if len(values)==len(columns)+1 and abs(values[0]-mass)<.01:
            rows.append(values)
    if len(rows)!=1:
        raise ValueError(f"Expected one HDECAY SM row at mass {mass}: {rows}")
    vals=rows[0]
    if any(x<0 for x in vals[1:]):
        raise ValueError("Negative HDECAY rate")
    return {"printed_mass_GeV":vals[0],
            "values":dict(zip(columns,vals[1:]))}


def parse_outputs(directory,mass):
    d=Path(directory)
    a=parse_table(d/"br.sm1",C1,mass)
    b=parse_table(d/"br.sm2",C2+("total_width_GeV",),mass)
    if abs(a["printed_mass_GeV"]-b["printed_mass_GeV"])>.00001:
        raise ValueError("HDECAY two files have different masses")
    br={**a["values"], **{k:v for k,v in b["values"].items() if k!="total_width_GeV"}}
    width=b["values"]["total_width_GeV"]
    if width<=0 or not(.995<=sum(br.values())<=1.005):
        raise ValueError(f"HDECAY total BR sum or width invalid: {width}, {sum(br.values())}")
    return {"printed_mass_GeV":a["printed_mass_GeV"],
            "width_GeV":width,"BR":br,
            "partial_width_GeV":{k:v*width for k,v in br.items()},
            "BR_sum":sum(br.values()),
            "source_files_sha256":{
                str(f):sha(d/f) for f in ("hdecay.in","br.sm1","br.sm2","run.stdout.log")
            }}


def evaluate(source,runs):
    head=subprocess.check_output(
        ["git","-C",str(source),"rev-parse","HEAD"],text=True).strip()
    if head!=COMMIT:
        raise ValueError("HDECAY source version mismatch: "+head)
    output={"status":"INDEPENDENT_HDECAY_SM_REFERENCE_ONLY",
            "reference_code":"HDECAY, source-pinned 2025-11-17 snapshot, internal 6.62",
            "HDECAY_git_commit":head,
            "HDECAY_source_sha256":sha(Path(source)/"zzz/hdecay.f"),
            "HDECAY_makefile_sha256":sha(Path(source)/"zzz/makefile"),
            "physics_scope":"SM ONLY; no Zee particles, no charged loops, no BSM electroweak corrections",
            "results":{},
            "physical_higgs":"UNKNOWN","globally_allowed":None}
    for order,dirpath in runs.items():
        if order not in MASS:
            raise ValueError("Unsupported HDECAY point")
        row=parse_outputs(dirpath,MASS[order])
        row["native_Zee_mass_GeV"]=MASS[order]
        row["order"]=order
        row["widths_are_Zee_predictions"]=False
        output["results"][order]=row
    if set(output["results"])!=set(MASS):
        raise ValueError("Both NO and IO HDECAY SM controls required")
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--NO",type=Path,required=True)
    ap.add_argument("--IO",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=evaluate(a.source,{"NO":a.NO,"IO":a.IO})
    with a.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False); f.write("\n")
    for o,r in result["results"].items():
        print("HDECAY_SM",o,"mH",r["printed_mass_GeV"],"Gamma",r["width_GeV"],
              "bb",r["partial_width_GeV"]["bb"],"WW",r["partial_width_GeV"]["WW"],
              "ZZ",r["partial_width_GeV"]["ZZ"],"physical=UNKNOWN")


if __name__=="__main__":
    main()
