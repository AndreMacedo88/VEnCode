"""
test_extensions.py: File containing a set of unittest.TestCase runnable tests for the classes and functions
in the internals_extensions.py file.
"""

import os
import sys
import unittest
from pathlib import Path

file_dir = str(Path(__file__).resolve().parents[2])
sys.path.append(file_dir)

from VEnCode import internals_extensions as iext
from VEnCode import common_variables as cv


class GettingVencodesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        cls.inputs = cv.expression_data1
        cls.celltype_analyse = "celltypetarget"
        cls.replicate_suffix = "_donor"
        cls.algorithm = "heuristic"
        cls.k = 4
        cls.thresholds = (0, 0.5, 0)  # inact, act, and sparseness, respectively
        cls.files_path = ""


class GetVencodesHeuristicTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencodes_obj = iext.GetVencodes(inputs=self.inputs,
                                             files_path="test",
                                             cell_type=self.celltype_analyse,
                                             algorithm=self.algorithm,
                                             n_regulatory_elements=self.k,
                                             number_vencodes=4,
                                             thresholds=self.thresholds, n_samples=10000,
                                             merge={"replicate_suffix": self.replicate_suffix})
        self.vencodes = self.vencodes_obj.coordinates

    def test_first_vencode(self):
        expected = ['chr12:109554241..109554255,+', 'chr12:109568972..109568988,+', 'chr12:109569139..109569148,+',
                    'chr5:42175384..42175396,-']
        self.assertCountEqual(expected, self.vencodes[0])

    def test_second_vencode(self):
        expected = ['chr12:109554241..109554255,+', 'chr12:109568972..109568988,+', 'chr12:109569139..109569148,+',
                    'chr7:112614228..112614232,+']
        self.assertCountEqual(expected, self.vencodes[1])

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencodes_obj.data:
            vencode_data.drop(self.vencodes_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencodes_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)


class GetVencodesSamplingTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencodes_obj = iext.GetVencodes(inputs=self.inputs,
                                             files_path="test",
                                             cell_type=self.celltype_analyse,
                                             algorithm="sampling",
                                             n_regulatory_elements=self.k,
                                             number_vencodes=4,
                                             thresholds=self.thresholds, n_samples=1000,
                                             merge={"replicate_suffix": self.replicate_suffix})
        self.vencodes = self.vencodes_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencodes_obj.data:
            vencode_data.drop(self.vencodes_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencodes_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)


class GetVencodesUnmergedTest(GettingVencodesTest):

    def builder(self, cell_type, replicates):
        """
        Builds the general skeleton for the next set of tests.

        Parameters
        ----------
        cell_type : str
            The celltype to get VEnCodes for.
        replicates : bool
            True or False whether the data contains replicates for the target celltype.
        """
        self.vencodes_obj = iext.GetVencodes(inputs=self.inputs,
                                             files_path="test",
                                             cell_type=cell_type,
                                             algorithm=self.algorithm,
                                             n_regulatory_elements=self.k,
                                             number_vencodes=4,
                                             thresholds=self.thresholds, n_samples=10000,
                                             merge=None, replicates=replicates)
        self.vencodes = self.vencodes_obj.coordinates
        for vencode_data in self.vencodes_obj.data:
            vencode_data.drop(self.vencodes_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencodes_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_merged(self):
        self.builder("celltypetarget", True)

    def test_merged_no_replicates(self):
        self.builder("celltypetarget_donor1", False)


class GetVencodesAddCelltypeTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        add_celltype_file = cv.expression_data2
        self.add_celltype_ctp = "celltypetarget2"
        add_celltype_ = [add_celltype_file, self.add_celltype_ctp, {"files_path": "test"}]
        self.vencodes_obj = iext.GetVencodes(inputs=self.inputs,
                                             files_path="test",
                                             cell_type=self.add_celltype_ctp,
                                             algorithm=self.algorithm,
                                             n_regulatory_elements=self.k,
                                             number_vencodes=4,
                                             thresholds=self.thresholds, n_samples=10000,
                                             merge={"replicate_suffix": self.replicate_suffix},
                                             add_celltype=add_celltype_)
        self.vencodes = self.vencodes_obj.coordinates

    def test_correct_addition(self):
        for vencode_data in self.vencodes_obj.data:
            with self.subTest(i=vencode_data.index.values.tolist()):
                self.assertIn(self.add_celltype_ctp, vencode_data.columns)

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencodes_obj.data:
            vencode_data.drop(self.vencodes_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencodes_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)


if __name__ == "__main__":
    unittest.main(verbosity=2)
