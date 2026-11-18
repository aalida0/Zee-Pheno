"""Build a general LFV point from an actual native SPheno card and spectrum.

No fitting or assumption about physical viability. Quark propagator masses
and Yukawa masses must be supplied with an explicit scheme.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from LFV.barr_zee import QUARKS
from LFV.pilot import evaluate as native_pilot
from Neutrino.fit import complex_matrix, encode


def from_native(native_path,card_path,spectrum_path,quark_path,yukawa_path,
                scheme,point_id,potential_path=None):
    # The existing pilot verifies card SHA, spectrum blocks and mixing matrix.
    pilot=native_pilot(native_path,card_path,spectrum_path)
    native=json.loads(native_path.read_text())
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
        raise ValueError("Native lepton running masses changed")
    q=json.loads(quark_path.read_text())
    qy=json.loads(yukawa_path.read_text())
    for label,values in (("loop",q),("Yukawa",qy)):
        if (not isinstance(values,dict) or set(values)!=set(QUARKS) or
            any(not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0
                for x in values.values())):
            raise ValueError(f"Explicit six positive quark {label} masses required")
    if not isinstance(scheme,str) or not scheme.strip():
        raise ValueError("Explicit mass and renormalization scheme label required")
    point={"point_id":point_id,"scalar":pilot["scalar_inputs"],
           "F":encode(f),"D":encode(d),
           "internal_lepton_masses_GeV":ml.tolist(),
           "pole_lepton_masses_GeV":pole.tolist(),
           "fermion_loop_masses_GeV":dict(q,electron=float(ml[0]),
                                         muon=float(ml[1]),tau=float(ml[2])),
           "quark_yukawa_masses_GeV":qy,
           "alpha_em":pilot["alpha_em"],
           "total_lepton_widths_GeV":pilot["total_lepton_widths_GeV"],
           "mass_yukawa_scheme":scheme}
    if potential_path:
        point["potential"]=json.loads(potential_path.read_text())
    return point


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("native","card","spectrum","quarks","quark-yukawas","scheme",
                 "point-id","output"):
        ap.add_argument("--"+name,required=True,
                        type=Path if name in ("native","card","spectrum","quarks",
                                              "quark-yukawas","output") else str)
    ap.add_argument("--potential",type=Path)
    args=ap.parse_args()
    point=from_native(args.native,args.card,args.spectrum,args.quarks,
                      args.quark_yukawas,args.scheme,args.point_id,args.potential)
    with args.output.open("x") as out:
        json.dump(point,out,indent=2,allow_nan=False)
        out.write("\n")
    print(f"Prepared {point['point_id']} from native SPheno card and spectrum")


if __name__=="__main__":
    main()
