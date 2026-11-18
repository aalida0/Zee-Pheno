"""No-tuning LFV omitted-term study regressions (strict 2HDM only)."""
import copy
import math
import unittest

import numpy as np

from LFV.approximation_study import (
    DEFAULT_GRID, RETAINED_SOURCE_BLOCKS, OMITTED_2HDM_SOURCE_BLOCKS,
    SOURCE_BLOCKS, assess_reference_channel, grid_points
)
from LFV.test_strict_2hdm_reference import decoupled_sample
from Model.model import scalar_parameters


def z(x):
    return {"real":x.real,"imag":x.imag}


def mock_reference(kept=1.+0j, omitted=.3+0j, hw=.001+0j):
    parts={}
    for name in SOURCE_BLOCKS:
        x=kept if name=="1loop" else omitted if name=="kin" else hw if name=="HW" else 0j
        parts[name]={"C_L":z(x),"C_R":z(0j)}
    full=kept+omitted+hw
    return {"reference_dipoles":{"mu_to_e_gamma":{
        "C_L_GeV_minus2":z(full),"C_R_GeV_minus2":z(0j),
        "raw_public_source_components":parts,
    }}}


class MaterialityTests(unittest.TestCase):
    def test_component_partition_and_ratio_exact(self):
        self.assertEqual(RETAINED_SOURCE_BLOCKS & OMITTED_2HDM_SOURCE_BLOCKS,
                         frozenset())
        self.assertEqual(RETAINED_SOURCE_BLOCKS | OMITTED_2HDM_SOURCE_BLOCKS,
                         SOURCE_BLOCKS)
        row=assess_reference_channel(mock_reference(), "mu_to_e_gamma",80.)
        self.assertAlmostEqual(row["author_approx_norm_GeV_minus2"],1.)
        self.assertAlmostEqual(row["author_omitted_norm_GeV_minus2"],.301)
        self.assertAlmostEqual(row["author_omitted_over_approx"],.301)
        self.assertAlmostEqual(row["author_full_BR_over_approx_BR"],1.301**2)
        self.assertAlmostEqual(row["paper_HW_variant_omitted_over_approx"],6.7)
        self.assertEqual(row["physical_Zee_LFV_gate"],"UNKNOWN")
        self.assertEqual(row["genuine_Zee_F_gauge_remainder"],"UNBOUNDED")

    def test_cancelled_retained_amplitude_never_inferred_accurate(self):
        report=mock_reference(kept=0j,omitted=.1+0j,hw=0j)
        row=assess_reference_channel(report, "mu_to_e_gamma",80.)
        self.assertTrue(row["approx_zero"])
        self.assertTrue(row["omitted_over_approx_above_10pct"])
        self.assertEqual(row["physical_Zee_LFV_gate"],"UNKNOWN")

    def test_default_grid_no_adaptive_lfv_selection(self):
        template=decoupled_sample()
        rows=list(grid_points(template))
        self.assertEqual(len(rows),18)
        self.assertEqual(len(set(p["point_id"] for p,_ in rows)),18)
        self.assertEqual(len(set(v["tan_beta"] for _,v in rows)),3)
        self.assertEqual(len(set(v["beta_minus_alpha_rad"] for _,v in rows)),3)
        self.assertEqual(len(set(v["heavy_neutral_mass_GeV"] for _,v in rows)),2)
        self.assertEqual(template["scalar"]["mh2"],640.)
        for point,meta in rows:
            p=point["scalar"]
            self.assertEqual(point["potential"]["K1"],0.)
            self.assertEqual(point["potential"]["K2"],0.)
            self.assertEqual(point["potential"]["Mu4"],0.)
            self.assertAlmostEqual(math.tan(p["beta"]),meta["tan_beta"])
            self.assertAlmostEqual(p["beta"]-p["alpha"],
                                   meta["beta_minus_alpha_rad"],places=14)
            self.assertEqual(p["mh1"],125.)
            self.assertGreater(p["mH2"],p["mH1"])

    def test_malformed_grid_fails_not_selectively_accepted(self):
        g=copy.deepcopy(DEFAULT_GRID)
        g["tan_beta"]=[2.,float("nan")]
        with self.assertRaises(ValueError):
            list(grid_points(decoupled_sample(),g))
        g=copy.deepcopy(DEFAULT_GRID)
        g["heavy_neutral_mass_GeV"]=[]
        with self.assertRaises(ValueError):
            list(grid_points(decoupled_sample(),g))


if __name__=="__main__":
    unittest.main()
