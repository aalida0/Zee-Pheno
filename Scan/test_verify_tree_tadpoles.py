"""SPheno source-locking 510/530 tree tadpole reporter without running fake binary."""
import unittest

from ZeeHiggs.spheno_audit import NATIVE
from Scan.verify_tree_tadpoles import patch_card,read_reference


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
