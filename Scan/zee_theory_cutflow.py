"""Deterministic scalar-only preselection cutflow on declared pilot priors.

NO SPheno event, spectrum, neutrino fit or physical globally allowed
point is manufactured. This cheap dry run audits precisely the new
necessary BFB / full-Zee high-energy tree unitarity decisions before
an expensive source-attested native scan.
"""
import argparse
import json
import math
import random
from pathlib import Path

from Model.audit_neutrino_native import sha
from Scan.physical_proposals import scalar_point, validate


ROOT=Path(__file__).resolve().parents[1]


def run(config_path,snapshot_path):
    cfg=validate(json.loads(config_path.read_text()))
    snapjson=json.loads(snapshot_path.read_text())
    snap=(snapjson["snapshots"][snapjson["spectrum_snapshot_index"]]
          if "snapshots" in snapjson else snapjson["snapshot"])
    v=math.hypot(float(snap["vd"]),float(snap["vu"]))
    arms=[a for a in cfg["arms"] for _ in range(a["weight"])]
    rows=[]
    for i in range(cfg["scalar_count"]):
        arm=arms[i%len(arms)]
        rng=random.Random(f"{cfg['seed']}:{i}")
        sc,lag,gates=scalar_point(arm,rng,v)
        rows.append({
            "scalar_id":f"s{i:05d}",
            "arm":arm["name"],
            "status":gates["status"],
            "tree_theory":gates,
            "native_SPheno_evaluated":False,
        })
    names=(
        "ANALYTIC_SCALAR_PRESELECTION",
        "SCREENED_AT_DECLARED_SCALAR_PRIORS",
        "FAIL_PREDECLARED_KINEMATIC_OR_MASS_PRIOR",
        "INVALID_SCALAR_INPUT",
    )
    dist={name:sum(r["status"]==name for r in rows) for name in names}
    if sum(dist.values())!=len(rows):
        raise ValueError("Unaccounted scalar proposal status")
    return {
        "status":"SCALAR_PRIOR_THEORY_DRYRUN_NO_NATIVE_VALIDATION",
        "config_sha256":sha(config_path),
        "native_snapshot_sha256":sha(snapshot_path),
        "n_scalar_proposals":len(rows),
        "counts":dist,
        "full_Zee_unitarity_failed_analytic":sum(
            r["tree_theory"].get("full_tree_unitarity")=="FAIL" for r in rows),
        "BFB_axis_necessary_failed_analytic":sum(
            r["tree_theory"].get("singlet_doublet_axis_BFB_necessary")=="FAIL"
            for r in rows),
        "BFB_sufficient_certified_analytic":sum(
            r["tree_theory"].get("BFB_exact_sufficient_certificate")=="CERTIFIED_SUFFICIENT_PASS"
            for r in rows),
        "spheno_run_count":0,
        "neutrino_validated_count":0,
        "EWPO_native_validated_count":0,
        "globally_allowed":None,
        "rows":rows,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config",type=Path,default=ROOT/"Scan/physical_pilot_config.json")
    ap.add_argument("--snapshot",type=Path,
                    default=ROOT/"Model/reference/native_neutrino_snapshot.json")
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    result=run(args.config,args.snapshot)
    with args.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k!="rows"},indent=2))


if __name__=="__main__":
    main()
