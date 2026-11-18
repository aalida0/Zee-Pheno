"""Convert a native Zee/SPheno LFV pilot into a point-independent LFV input.

No benchmark tuning or automatic quark mass scheme selection. Quark input
masses and their renormalization scheme must be supplied explicitly.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from LFV.barr_zee import QUARKS
from Model.audit_neutrino_native import sha
from Neutrino.fit import complex_matrix, encode


def from_native(pilot_path,native_path,quark_path,yukawa_path,scheme,point_id,potential_path=None):
    pilot=json.loads(pilot_path.read_text())
    native=json.loads(native_path.read_text())
    if (pilot["status"]!="SCALAR_ONE_LOOP_SUBSET_ONLY" or
        pilot["native_diagnostic_sha256"]!=sha(native_path) or
        pilot["native_input_sha256"]!=native["input_sha256"]):
        raise ValueError("Pilot/native source hash or status mismatch")
    snapshot=native["snapshots"][native["spectrum_snapshot_index"]]
    vd,vu=float(snapshot["vd"]),float(snapshot["vu"])
    v=math.hypot(vd,vu)
    cb,sb=vd/v,vu/v
    el,er=complex_matrix(snapshot["ZEL"]),complex_matrix(snapshot["ZER"])
    f=el.conj()@complex_matrix(snapshot["fzee"])@el.conj().T
    d=-sb*complex_matrix(snapshot["Ye"])+cb*complex_matrix(snapshot["Ye2"])
    d=er.conj()@d@el.conj().T
    ml=np.asarray(pilot["internal_running_lepton_masses_GeV"],float)
    pole=np.asarray(pilot["external_pole_lepton_masses_GeV"],float)
    if not np.allclose(ml,snapshot["MFe"],rtol=1e-12,atol=1e-15):
        raise ValueError("Pilot lepton masses differ from native snapshot")
    q=json.loads(quark_path.read_text())
    if set(q)!=set(QUARKS) or any(not math.isfinite(x) or x<=0 for x in q.values()):
        raise ValueError("Explicit six positive quark loop masses required")
    qy=json.loads(yukawa_path.read_text())
    if (set(qy)!=set(QUARKS) or
            any(not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0 for x in qy.values())):
        raise ValueError("Explicit six positive quark Yukawa masses required")
    if not isinstance(scheme,str) or not scheme.strip():
        raise ValueError("Mass/renormalization scheme label is mandatory")
    p=pilot["scalar_inputs"]
    data={"point_id":point_id, "scalar":p,
          "F":encode(f),"D":encode(d),
          "internal_lepton_masses_GeV":ml.tolist(),
          "pole_lepton_masses_GeV":pole.tolist(),
          "fermion_loop_masses_GeV":dict(q,electron=float(ml[0]),muon=float(ml[1]),tau=float(ml[2])),
          "quark_yukawa_masses_GeV":qy,
          "alpha_em":pilot["alpha_em"],
          "total_lepton_widths_GeV":pilot["total_lepton_widths_GeV"],
          "mass_yukawa_scheme":scheme}
    if potential_path is not None:
        data["potential"]=json.loads(potential_path.read_text())
    return data


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("pilot","native","quarks","quark-yukawas","scheme","point-id","output"):
        ap.add_argument("--"+key,type=Path if key in ("pilot","native","quarks","quark-yukawas","output") else str,
                        required=True)
    ap.add_argument("--potential",type=Path)
    a=ap.parse_args()
    data=from_native(a.pilot,a.native,a.quarks,a.quark_yukawas,a.scheme,a.point_id,a.potential)
    with a.output.open("x") as file:
        json.dump(data,file,indent=2,allow_nan=False)
        file.write("\n")
    print(f"Prepared {data['point_id']} from hash-matched native running Yukawas")


if __name__=="__main__":
    main()
