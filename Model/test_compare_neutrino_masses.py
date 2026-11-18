"""Check singular-value comparison and deliberately distorted native result."""
import math
import tempfile
import unittest
from pathlib import Path

from Model.compare_neutrino_masses import compare, finite_kernel
from Model.test_basis import output


class NeutrinoComparisonTests(unittest.TestCase):
    def test_kernel_limits(self):
        self.assertEqual(finite_kernel(300., 300., 0.105), 0.)
        self.assertAlmostEqual(finite_kernel(300., 700., 0.),
                               math.log((300./700.)**2), places=14)
        # Independent numerical quadrature of the full finite muon integral.
        self.assertAlmostEqual(finite_kernel(300., 700., .105658),
                               -1.6945941491221044, places=12)

    def test_synthetic_native_agreement_and_bad_result(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            baseline, input_card, loop = (root/x for x in ("tree.spc","tree.in","loop.spc"))
            baseline.write_text(output() + "Block MASS\n"
                                " 11 0.000511\n 13 0.105658\n 15 1.77686\n"
                                " 37 300.\n 9000037 700.\n")
            tanb = math.tan(1.24)
            mu4 = -.8*(300.**2-700.**2)/(math.sqrt(2)*246.21965)
            input_card.write_text(f"Block MINPAR\n 10 {tanb:.17e}\n 13 {mu4:.17e}\n")
            loop.write_text("Block MASS\n 12 0.\n 14 0.\n 16 0.\n")
            reference = compare(input_card, baseline, loop)
            self.assertEqual(reference["charged_mass_relative_shifts"], {})
            self.assertFalse(reference["matched_native_running_parameter_scheme"])
            values = reference["analytic_mass_singular_values_eV"]
            self.assertGreater(values[-1], 0.)
            loop.write_text("Block MASS\n" + "".join(
                f" {i} {v/1e9:.17e}\n" for i,v in zip((12,14,16),values)))
            matched = compare(input_card, baseline, loop)
            self.assertEqual(matched["status"], "SINGULAR_VALUES_AGREE_ONLY")
            self.assertFalse(matched["one_loop_majorana_matrix_validated"])
            loop.write_text("Block MASS\n" + "".join(
                f" {i} {2*v/1e9:.17e}\n" for i,v in zip((12,14,16),values)))
            self.assertEqual(compare(input_card, baseline, loop)["status"],
                             "NATIVE_ANALYTIC_MASS_MISMATCH")


if __name__ == "__main__":
    unittest.main()
