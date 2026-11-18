"""Focused general-LFV tests: conventions, arbitrary points, fail-closed gates."""
import json
import tempfile
from pathlib import Path
import math
import unittest

import numpy as np

from LFV.barr_zee import (photon_kernels, fermion_photon,
                          charged_scalar_photon, physical_trilinears)
from LFV.constraints import evaluate, input_parameters, amplitude_budget
from LFV.native_point import from_native
from Model.audit_neutrino_native import sha
from LFV.core import combine, dipole_ledger


def cmat(a):
    a=np.asarray(a,complex)
    return {"real":a.real.tolist(),"imag":a.imag.tolist()}


def sample():
    # Illustrative algebraic point, NOT a scalar/neutrino viable benchmark.
    return {
        "point_id":"unit-test-parameter-point",
        "scalar":dict(v=246.2,beta=1.1,alpha=.7,gamma=.8,
                      mh1=125.,mh2=640.,ma=720.,mH1=610.,mH2=1030.),
        "F":cmat(np.array([[0,2e-5,3e-5],[-2e-5,0,4e-5],
                          [-3e-5,-4e-5,0]],complex)),
        "D":cmat(np.array([[.002,1e-4+3e-5j,0],
                           [1e-5j,.004,0],[0,0,.01]],complex)),
        "internal_lepton_masses_GeV":[.000503,.1041,1.75],
        "pole_lepton_masses_GeV":[.00051099893,.105658372,1.77686],
        "fermion_loop_masses_GeV":dict(
            up=.0022,charm=1.3,top=172.,down=.0047,strange=.095,bottom=4.2,
            electron=.000503,muon=.1041,tau=1.75),
        "quark_yukawa_masses_GeV":dict(up=.0022,charm=1.3,top=172.,down=.0047,strange=.095,bottom=4.2),
        "alpha_em":.0072973525692838015,
        "total_lepton_widths_GeV":{"mu":2.995983701908041e-19,
                                   "tau":2.267350867723045e-12},
        "mass_yukawa_scheme":"explicit nominal example; no matched-scale certification",
    }


class GeneralLFVTests(unittest.TestCase):
    def test_signed_photon_kernels_and_source_limits(self):
        (f,g,f2),errors=photon_kernels(1.)
        self.assertTrue(0<f<g and f2<0)
        self.assertTrue(all(e<1e-8 for e in errors))
        for z in (1e-9,.1,10.):
            a,e=photon_kernels(z)
            self.assertGreater(a[0],0)
            self.assertGreater(a[1],0)
            self.assertLess(a[2],0)
            self.assertTrue(all(x>=0 for x in e))
        a,_=fermion_photon(np.array([.002,.003],complex),.2,172.,640.,
                          .105658372,.0073,3,2/3)
        self.assertLess(a[0].real,0)
        self.assertLess(a[1].real,0)
        b,_=charged_scalar_photon(np.array([.002,.003]),120.,630.,640.,
                                  .105658372,.0073)
        self.assertLess(b[0].real,0)
        self.assertLess(b[1].real,0)

    def test_tiny_z_scalar_kernel_uses_stable_identity(self):
        # Tiny-z direct F2 quadrature triggered SciPy IntegrationWarning
        # despite an excellent stable f,g evaluation.
        from scipy.integrate import quad
        from LFV.barr_zee import _ratio
        for z in (1e-24, 6.18e-13, 9.34e-12, 1e-9, 0.1, 1., 10.):
            (f,g,f2), errors=photon_kernels(z)
            self.assertTrue(math.isfinite(f2) and f2<0)
            self.assertLess(errors[2],max(2e-10,abs(f2)*2e-9))
            self.assertAlmostEqual(f2, (f-g)/(2*z),places=10)
        # Independent direct integral at moderate mass ratios.
        for z in (0.1, 1., 10.):
            def integrand(x):
                u=x*(1-x)
                return -u*_ratio(u,z)/(2*z)
            direct=quad(integrand,0.,.5,epsabs=1e-12,epsrel=1e-11)[0]*2
            self.assertAlmostEqual(photon_kernels(z)[0][2],direct,places=10)

    def test_arbitrary_point_is_partial_not_physical_exclusion(self):
        p=sample()
        out=evaluate(p)
        self.assertEqual(out["full_lfv_gate"],"UNKNOWN")
        self.assertIsNone(out["global_zee_acceptance"])
        self.assertEqual(out["charged_scalar_photon"],"UNKNOWN_POTENTIAL_NOT_SUPPLIED")
        for channel in out["processes"].values():
            self.assertEqual(len(channel["included_components"]),34)
            self.assertIsNone(channel["physical_branching_ratio"])
            self.assertEqual(channel["experimental_gate"],"UNKNOWN")
            self.assertTrue(math.isfinite(channel["branching_ratio_included_topologies_only"]))
            self.assertGreaterEqual(channel["branching_ratio_included_topologies_only"],0.)
            self.assertLessEqual(channel["cancellation_ratio"],1+1e-12)
            parts=channel["included_components"].values()
            for key in ("C_L","C_R"):
                summed=sum(complex(row[key]["real"],row[key]["imag"]) for row in parts)
                ret=complex(channel[key+"_GeV_minus2"]["real"],
                            channel[key+"_GeV_minus2"]["imag"])
                self.assertLess(abs(summed-ret),1e-25)

    def test_neutral_higgs_partial_widths_and_singlet_universality(self):
        p=sample()
        p["GF_GeV_minus2"]=1.1663787e-5
        out=evaluate(p)
        self.assertEqual(set(out["neutral_higgs_lfv_partial_widths"]),{"h","H","A"})
        self.assertGreater(out["neutral_higgs_lfv_partial_widths"]["h"]["e_mu"]
                           ["partial_width_two_charges_GeV"],0)
        self.assertIsNone(out["neutral_higgs_lfv_partial_widths"]["h"]["e_mu"]
                          ["physical_branching_ratio"])
        uni=out["universality_F_subset"]
        self.assertEqual(uni["physical_status"],"SUBSET_DIAGNOSTIC")
        p2=sample()
        f=np.array(p2["F"]["real"])+1j*np.array(p2["F"]["imag"])
        p2["F"]=cmat(2*f)
        p2["GF_GeV_minus2"]=p["GF_GeV_minus2"]
        uni2=evaluate(p2)["universality_F_subset"]
        self.assertAlmostEqual(uni2["rate_shifts_before_GF_subtraction"]["mu_e"]/
                               uni["rate_shifts_before_GF_subtraction"]["mu_e"],4.)
        p3=sample()
        self.assertIsNone(evaluate(p3)["universality_F_subset"])

    def test_potential_adds_six_photon_scalar_rows_without_loss(self):
        p=sample()
        p["potential"]=dict(Lambda1=.25,Lambda2=.45,Lambda3=.1,
                            Lambda4=-.1,Lambda5=-.05,K1=.5,K2=.4,
                            K12=.01,Mu4=110.)
        gs=physical_trilinears(p["scalar"],p["potential"])
        self.assertEqual(set(gs),{"h","H"})
        self.assertEqual(np.shape(gs["h"]),(2,))
        out=evaluate(p)
        self.assertEqual(out["charged_scalar_photon"],
                         "INCLUDED_FROM_RESTRICTED_POTENTIAL")
        self.assertEqual(len(out["processes"]["mu_to_e_gamma"]["included_components"]),40)
        self.assertEqual(out["full_lfv_gate"],"UNKNOWN")
        pieces=out["processes"]["mu_to_e_gamma"]["included_components"]
        for x in ("H1","H2"):
            row=pieces[f"2loop_A_{x}_charged_gamma"]
            self.assertEqual(row["C_L"]["real"],0.)
            self.assertEqual(row["C_R"]["real"],0.)

    def test_identical_zero_lfv_point_has_zero_radiative_rates(self):
        p=sample()
        p["F"]=cmat(np.zeros((3,3)))
        p["D"]=cmat(np.diag([.002,.004,.01]))
        out=evaluate(p)
        for channel in out["processes"].values():
            self.assertEqual(channel["branching_ratio_included_topologies_only"],0.)
            self.assertIsNone(channel["cancellation_ratio"])
            self.assertEqual(channel["experimental_gate"],"UNKNOWN")

    def test_native_adapter_is_hash_bound_and_has_no_hidden_fit(self):
        source=json.loads((Path(__file__).resolve().parents[1]/
                           "Model/reference/native_neutrino_snapshot.json").read_text())
        snap=source["snapshot"]
        ex=sample()
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            card=root/"card.in"
            card.write_text("native input fixture\n")
            native_path=root/"native.json"
            native_path.write_text(json.dumps({
                "input_sha256":sha(card),"snapshots":[snap],
                "spectrum_snapshot_index":0}))
            spectrum=root/"spectrum.spc"
            spectrum.write_text(
                "BLOCK MASS\n 11 0.00051099893\n 13 0.105658372\n"
                " 15 1.77686\n 25 125.\n 35 600.\n 36 640.\n"
                "BLOCK SCALARMIX Q= 1.600E+02\n"
                " 1 1 7.25651649E-01\n 1 2 6.88062268E-01\n"
                " 2 1 -6.88062268E-01\n 2 2 7.25651649E-01\n")
            quarks=root/"quarks.json"
            quarks.write_text(json.dumps({k:ex["fermion_loop_masses_GeV"][k]
                                         for k in ("up","charm","top","down","strange","bottom")}))
            point=from_native(native_path,card,spectrum,quarks,quarks,
                              "explicit test scheme","fixture")
            self.assertEqual(point["point_id"],"fixture")
            self.assertEqual(point["mass_yukawa_scheme"],"explicit test scheme")
            self.assertEqual(evaluate(point)["full_lfv_gate"],"UNKNOWN")
            spectrum.write_text("BLOCK MASS\n 25 125.\n")
            with self.assertRaises(ValueError):
                from_native(native_path,card,spectrum,quarks,quarks,
                            "explicit test scheme","fixture")
            card.write_text("tampered native input\n")
            with self.assertRaises(ValueError):
                from_native(native_path,card,spectrum,quarks,quarks,
                            "explicit test scheme","fixture")

    def test_optional_neutral_tree_three_lepton_gate(self):
        p=sample()
        p["tree_quadrature_order"]=8
        rows=evaluate(p)["three_lepton_tree"]
        self.assertEqual(set(rows),{
            "mu_to_3e","tau_to_3e","tau_to_3mu",
            "tau_to_mu_ee","tau_to_e_mumu"})
        for data in rows.values():
            self.assertGreaterEqual(data["BR_neutral_tree_only"],0.)
            self.assertTrue(math.isfinite(data["BR_neutral_tree_only"]))
            self.assertIsNone(data["physical_branching_ratio"])
            self.assertEqual(data["experimental_gate"],"UNKNOWN")
        p["tree_quadrature_order"]=0
        with self.assertRaises(ValueError):
            evaluate(p)

    def test_strict_2hdm_neutral_yukawa_dictionary(self):
        # Altmannshofer et al. arXiv:2410.17313v4, Eqs.(8),(11),(12).
        # Our paper conventions: rho_e^dagger = D/sqrt(2),
        # q_h=(cos(beta-alpha),-sin(beta-alpha));
        # q_H=(sin(beta-alpha),cos(beta-alpha)); q_A2=-i.
        # Exact neutral-sector map; NOT a Zee charged-gauge completion.
        from LFV.core import vertices
        from LFV.barr_zee import quark_coupling
        p=sample()
        scal=p["scalar"]
        ml=np.asarray(p["internal_lepton_masses_GeV"])
        D=np.array(p["D"]["real"])+1j*np.array(p["D"]["imag"])
        neutral=vertices(D,scal,ml)
        delta=scal["beta"]-scal["alpha"]
        q={"h":(math.cos(delta),-math.sin(delta)),
           "H":(math.sin(delta),math.cos(delta)),
           "A":(0.,-1j)}
        for scalar,(C,_) in neutral.items():
            q1,q2=q[scalar]
            expected=q1*np.diag(ml/scal["v"])+q2*D/np.sqrt(2)
            np.testing.assert_allclose(C,expected,rtol=1e-13,atol=1e-16)
        for species,mass in (("top",172.),("bottom",4.2)):
            rho=(mass/scal["v"]/math.tan(scal["beta"])
                 if species=="top" else -mass/scal["v"]*math.tan(scal["beta"]))
            for name,(q1,q2) in q.items():
                expected=(q1*mass/scal["v"]+
                          rho*(q2.conjugate() if species=="top" else q2))
                self.assertAlmostEqual(abs(quark_coupling(species,mass,name,scal)
                                            -expected),0.,places=13)

    def test_MEGII_budget_for_native_NO_IO_scalar_controls(self):
        # Source: MEG II arXiv:2504.15711 (1.5e-13, 90% CL).
        # These are *required upper bounds* on omitted amplitudes,
        # not estimates that have been established for a Zee point.
        alpha=0.0072973525692838015
        mu=0.105658372
        electron=0.00051099893
        gamma=2.995983701908041e-19
        K=alpha*mu**5*(1-electron**2/mu**2)**3/(4*gamma)
        limit=1.5e-13
        self.assertAlmostEqual(K,8.017795784561368e10,delta=2.)
        for label,br,nref,remaining in (
            ("NO",1.2359863349731284e-17,1.2415932400992839e-14,1.3553700049540977e-12),
            ("IO",2.421993374968115e-17,1.7380368501672127e-14,1.3504055688534183e-12)):
            with self.subTest(ordering=label):
                norm=math.sqrt(br/K)
                self.assertAlmostEqual(norm,nref,delta=5e-25)
                budget=amplitude_budget(norm,mu,electron,gamma,alpha,
                         {"BR_90CL":limit,
                          "source":"MEG II 2025 arXiv:2504.15711"})
                self.assertAlmostEqual(
                    budget["required_max_missing_norm_for_sufficient_pass_GeV_minus2"],
                    remaining,delta=5e-23)
                self.assertEqual(budget["status"],"UNKNOWN_UNBOUNDED_MISSING")

    def test_conditional_missing_amplitude_envelope(self):
        # factor = alpha * mi^5/(4 Gamma) = 1 exactly.
        limit={"BR_90CL":.25,"source":"synthetic experimental limit"}
        for retained,omitted,expected in ((.1,.2,"PASS_IF_DECLARED_BOUND_VALID"),
                                           (1.,.2,"FAIL_IF_DECLARED_BOUND_VALID"),
                                           (.8,.4,"UNKNOWN_EVEN_WITH_DECLARED_BOUND")):
            x=dict(limit,missing_norm_upper_GeV_minus2=omitted,
                   bound_source="synthetic rigorous norm bound")
            status=amplitude_budget(retained,1.,0.,.025,.1,x)
            self.assertEqual(status["status"],expected)
        no_bound=amplitude_budget(.1,1.,0.,.025,.1,limit)
        self.assertEqual(no_bound["status"],"UNKNOWN_UNBOUNDED_MISSING")
        oversized=amplitude_budget(1.,1.,0.,.025,.1,limit)
        self.assertIsNone(
            oversized["required_max_missing_norm_for_sufficient_pass_GeV_minus2"])
        self.assertEqual(oversized["status"],"UNKNOWN_UNBOUNDED_MISSING")
        p=sample()
        p["radiative_limits"]={"mu_to_e_gamma":dict(limit)}
        out=evaluate(p)
        self.assertEqual(out["full_lfv_gate"],"UNKNOWN")
        self.assertEqual(out["processes"]["mu_to_e_gamma"]["experimental_gate"],"UNKNOWN")
        self.assertIsNotNone(out["processes"]["mu_to_e_gamma"]
                             ["conditional_missing_amplitude_budget"])
        p["radiative_limits"]["mu_to_e_gamma"]["missing_norm_upper_GeV_minus2"]=-1.
        p["radiative_limits"]["mu_to_e_gamma"]["bound_source"]="invalid"
        with self.assertRaises(ValueError):
            evaluate(p)

    def test_rejects_bad_input_not_fabricated_rates(self):
        p=sample()
        p["F"]["real"][0][1]=1.
        with self.assertRaisesRegex(ValueError,"antisymmetric"):
            evaluate(p)
        p=sample();p["scalar"]["mh1"]=0.
        with self.assertRaises(ValueError):
            evaluate(p)
        p=sample();p["fermion_loop_masses_GeV"].pop("top")
        with self.assertRaises(ValueError):
            evaluate(p)
        p=sample();p["mass_yukawa_scheme"]=""
        with self.assertRaises(ValueError):
            evaluate(p)
        p=sample();p["potential"]={"Lambda1":.1}
        with self.assertRaises(ValueError):
            evaluate(p)


if __name__=="__main__":
    unittest.main()
