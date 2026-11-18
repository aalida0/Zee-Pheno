"""Generate a declared *restricted neutrino-inverse slice* of candidate cards.

Uses the archived reference running point and the existing minimum-norm,
single-null-root Zee inverse. This is not the full solution family and NOT
a native neutrino fit. Every proposed card must be run through SPheno and
the source-bound Scan.native_pipeline before any scan gates can be applied.
"""
import argparse
import json
import re
from pathlib import Path

from Model.audit_neutrino_native import sha
from Model.check_native import fix_card
from Neutrino.fit import (fit,complex_matrix,REFERENCE_CARD_SHA)

FIELDS={"point_id","ordering","mlight_eV","alpha21_rad","alpha31_rad"}


def build(catalog_path,output_dir):
    d=json.loads(catalog_path.read_text())
    if not isinstance(d,dict) or set(d)!={"schema_version","snapshot","reference_card","targets"}:
        raise ValueError("Unexpected neutrino target catalog keys")
    if d["schema_version"]!=1 or not isinstance(d["targets"],list) or not d["targets"]:
        raise ValueError("Invalid neutrino catalog schema or empty targets")
    if len(d["targets"])>1000:
        raise ValueError("Too many targets")
    template=(catalog_path.parent/d["reference_card"]).resolve()
    snapshot=(catalog_path.parent/d["snapshot"]).resolve()
    if sha(template)!=REFERENCE_CARD_SHA:
        raise ValueError("Unrecognized attested loop-enabled reference card")
    ids=set()
    generated=[]
    for t in d["targets"]:
        if not isinstance(t,dict) or set(t)!=FIELDS:
            raise ValueError("Each target must specify ordering, mass and phases")
        point_id=t["point_id"]
        if not isinstance(point_id,str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,80}",point_id):
            raise ValueError("point_id must be an ASCII filename-safe identifier")
        if point_id in ids:
            raise ValueError("Duplicate target ID")
        ids.add(point_id)
        if t["ordering"] not in ("NO","IO"):
            raise ValueError("Unknown neutrino ordering")
        report=fit(snapshot,t["ordering"],t["mlight_eV"],
                   t["alpha21_rad"],t["alpha31_rad"])
        if report["status"]!="ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN":
            raise AssertionError("Analytic result masquerades as native")
        names={"FZEEIN":complex_matrix(report["F_gauge"]),
               "YEE2IN":complex_matrix(report["Ye2_gauge"])}
        changes={(("IM" if imag else "")+name,(i+1,j+1)):
                 float(matrix[i,j].imag if imag else matrix[i,j].real)
                 for name,matrix in names.items()
                 for imag in (False,True) for i in range(3) for j in range(3)}
        card=fix_card(template.read_text(),changes)
        generated.append((point_id,card,report))
    output_dir.mkdir(parents=True,exist_ok=False)
    queue=[]
    for point_id,card,report in generated:
        filename=f"{point_id}.in"
        (output_dir/filename).write_text(card)
        report_path=f"{point_id}.proposal.json"
        report["proposed_card_sha256"]=sha(output_dir/filename)
        (output_dir/report_path).write_text(json.dumps(report,indent=2)+"\n")
        queue.append({
            "point_id":point_id,"ordering":report["ordering"],
            "card":filename,"report":report_path,
            "card_sha256":report["proposed_card_sha256"],
            "analytic_forward_residual":report["forward_residual"],
            "status":"NATIVE_REVALIDATION_REQUIRED",
            "solution_family":"ONE_ROOT_MINIMUM_NORM_D_SLICE_ONLY",
            "scalar_and_LFV_gates":"NOT_EVALUATED",
        })
    out={"schema_version":1,"status":"ANALYTIC_PROPOSALS_NOT_NATIVE_POINTS",
         "template_sha256":sha(template),
         "snapshot_sha256":sha(snapshot),
         "n_proposals":len(queue),
         "neutrino_inverse_scope":"MINIMUM_NORM_ONE_NULL_CONE_ROOT_NOT_GLOBAL",
         "all_proposals":queue}
    (output_dir/"proposal_catalog.json").write_text(
        json.dumps(out,indent=2)+"\n")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    summary=build(a.catalog,a.output_dir)
    print(json.dumps({k:v for k,v in summary.items() if k!="all_proposals"},indent=2))


if __name__=="__main__":
    main()
