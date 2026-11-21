"""SPheno source-locking 510/530 tree tadpole reporter without running fake binary."""
import unittest

from ZeeHiggs.spheno_audit import NATIVE
from Scan.verify_tree_tadpoles import patch_card,read_reference,evaluate_tree_output
from ZeeHiggs.spheno_audit import read_slha


class NativeTreeTadpoleBridgeTests(unittest.TestCase):
    def test_exact_native_original_sha_and_flags(self):
        for order in ("NO","IO"):
            orig,path=read_reference(order)
            card=orig.read_text()
            a=patch_card(card,510)
            b=patch_card(card,530)
            self.assertNotEqual(a,b)
            self.assertNotEqual(card,a)
            self.assertNotEqual(card,b)
            self.assertIn("510 1  # source-bound",a)
            self.assertIn("530 0  # source-bound",a)
            self.assertIn("510 0  # source-bound",b)
            self.assertIn("530 1  # source-bound",b)
            # Preserve untouched physics model input cards beyond
            # the two documented SPhenoInput control flags.
            original_input=card.split("Block FZEEIN",1)[1]
            self.assertEqual(a.split("Block FZEEIN",1)[1],original_input)
            self.assertEqual(b.split("Block FZEEIN",1)[1],original_input)

    def test_no_tree_block_reports_unknown_not_error(self):
        for order in ("NO","IO"):
            _,directory=read_reference(order)
            blocks,_=read_slha((directory/"spectrum.slha").read_text())
            result=evaluate_tree_output(blocks)
            self.assertEqual(result["status"],"UNKNOWN_NO_MATCHED_NATIVE_TREE_TADPOLES")
            self.assertFalse(result["on_this_distinct_output_tree_vacuum_certified"])
            self.assertEqual(result["matched_stationary_source_blocks"],[])

    def test_tree_block_requires_stationarity_and_global_certificate(self):
        from Scan.zee_chargebreaking_vacuum import (
            source_inputs_from_blocks,normal_tree_mass_squared)
        from Scan.zee_ewpo import block_1d
        _,directory=read_reference("NO")
        blocks,_=read_slha((directory/"spectrum.slha").read_text())
        p=source_inputs_from_blocks(blocks)
        masses=normal_tree_mass_squared(p["quartics"],p["vevs"],p["m12sq"])
        tree_rows=[]
        for row,raw in blocks["ZEEPOT"]:
            if row[0] in ("1","2"):
                value=masses[int(row[0])-1]
                row=[row[0],f"{value:.15E}"]
                raw=" ".join(row)
            tree_rows.append((row,raw))
        result=evaluate_tree_output({**blocks,"TREEZEEPOT":tree_rows})
        self.assertIn("TREEZEEPOT",result["matched_stationary_source_blocks"])
        self.assertTrue(result["on_this_distinct_output_tree_vacuum_certified"])
        malformed=evaluate_tree_output({**blocks,"TREEZEEPOT":[(["garbage"],"garbage")]})
        self.assertFalse(malformed["on_this_distinct_output_tree_vacuum_certified"])
        self.assertIn("TREEZEEPOT",malformed["tree_block_parse_errors"])

    def test_reject_unexpected_source_or_flag(self):
        inp,_=read_reference("NO")
        with self.assertRaises(ValueError):
            patch_card(inp.read_text(),530.1)
        with self.assertRaises(ValueError):
            patch_card("Block MINPAR\n 1 2\n",510)
        with self.assertRaises(ValueError):
            patch_card(inp.read_text().replace("510 0","510 X"),510)
        with self.assertRaises(ValueError):
            read_reference("neither")


if __name__=="__main__":unittest.main()
