"""Fail closed: scaling must touch Fzee only and preserve complex antisymmetry."""
import unittest
from Model.probe_neutrino_scaling import rescale_f
from Model.check_native import blocks, mat, antisym


def sample():
    lines = ["Block MINPAR\n 10 2 # beta\n", "Block SPHENOINPUT\n 55 1\n"]
    for block in ("FZEEIN", "IMFZEEIN", "YEE2IN", "IMYEE2IN"):
        lines.append(f"Block {block}\n")
        for i in range(1, 4):
            for j in range(1, 4):
                value = (1e-4 if (i,j)==(1,2) else
                         -1e-4 if (i,j)==(2,1) else 0.)
                lines.append(f" {i} {j} {value:.17e} # keep me\n")
    return "".join(lines)


class ScalingTests(unittest.TestCase):
    def test_controls_preserve_input_others_and_antisymmetry(self):
        original = sample()
        for factor in (0,2):
            updated = rescale_f(original, factor)
            original_data, new_data = blocks(original), blocks(updated)
            for name in ("MINPAR", "YEE2IN", "IMYEE2IN"):
                self.assertEqual(original_data[name], new_data[name])
            self.assertEqual(antisym(mat(new_data,"FZEEIN")), 0.)
            for name in ("FZEEIN","IMFZEEIN"):
                for key,value in original_data[name].items():
                    self.assertAlmostEqual(new_data[name][key], factor*value)

    def test_missing_matrix_records_are_rejected(self):
        with self.assertRaises(ValueError):
            rescale_f(sample().replace("Block IMFZEEIN", "Block LOST"), 0)
        with self.assertRaises(ValueError):
            rescale_f(sample(), 1)


if __name__ == "__main__":
    unittest.main()
