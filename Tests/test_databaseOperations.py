from unittest import TestCase
from VEnCode_FANTOM5 import Classes_2

promoters = Classes_2.Promoters("hg19.cage_peak_phase1and2combined_tpm.osc.txt",)
prom_data = promoters

class TestDatabaseOperations(TestCase):
    def test_at_least_one_node_calculator_2(self):
        self.assertEqual(Classes_2.DatabaseOperations.at_least_one_node_calculator())
        self.fail()
