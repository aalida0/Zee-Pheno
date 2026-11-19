"""Archived-card regression; never relabel partial LFV as physical acceptance."""
import math
import tempfile
import unittest
from pathlib import Path

from Scan.fitted_benchmark_triage import evaluate, decay_block, mixing_entries


class FixedBenchmarkTriageTests(unittest.TestCase):
    def test_two_archived_benchmarks(self):
        for order, reference_br, flavor in (
            ("NO", 0.278192521, 3.6400e-19),
            ("IO", 0.278190181, 2.8396e-13),
        ):
            with self.subTest(order=order):
                result = evaluate(order)
                self.assertIsNone(result["globally_allowed"])
                self.assertEqual(result["lfv"]["mu_to_e_gamma"]["full_Zee_lfv_gate"],
                                 "UNKNOWN")
                self.assertAlmostEqual(
                    result["cascade"]["branching_fraction_spheno"],
                    reference_br, places=8)
                self.assertLess(result["cascade"]["tree_width_relative_residual"], 1e-5)
                self.assertAlmostEqual(
                    result["lfv"]["mu_to_e_gamma"]["scalar_one_loop_subset_BR"] /
                    flavor, 1.0, delta=0.001)
                self.assertAlmostEqual(result["flavor"]["spheno_BR_B_to_Xs_gamma"],
                                       3.63471789e-4, places=10)
                self.assertGreater(result["neutrino"]["sum_masses_eV"], 0.1)

    def test_malformed_decay_and_mixing_are_refused(self):
        with self.assertRaises(ValueError):
            decay_block("DECAY 9000037 1.0\n 0.5 3 37 23\n", 9000037)
        with self.assertRaises(ValueError):
            mixing_entries("Block CHARGEMIX\n 1 1 1.0\n", "CHARGEMIX")


if __name__ == "__main__":
    unittest.main()
