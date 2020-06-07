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
from VEnCode.utils import dir_and_file_handling as dh
from VEnCode.utils import validation_utils as val
from VEnCode import outside_data


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
        cls.thresholds = (0.5, 0, 0)  # act, inact, and sparseness, respectively
        cls.files_path = "test"


class GetVencodesHeuristicTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencode_obj = iext.GetVencodes(inputs=self.inputs,
                                            files_path=self.files_path,
                                            cell_type=self.celltype_analyse,
                                            algorithm=self.algorithm,
                                            n_regulatory_elements=self.k,
                                            number_vencodes=4,
                                            thresholds=self.thresholds, n_samples=10000,
                                            merge={"replicate_suffix": self.replicate_suffix})
        self.vencodes = self.vencode_obj.coordinates

    def test_first_vencode(self):
        expected = ['chr12:109554241..109554255,+', 'chr12:109568972..109568988,+', 'chr12:109569139..109569148,+',
                    'chr5:42175384..42175396,-']
        self.assertCountEqual(expected, self.vencodes[0])

    def test_second_vencode(self):
        expected = ['chr12:109554241..109554255,+', 'chr12:109568972..109568988,+', 'chr12:109569139..109569148,+',
                    'chr7:112614228..112614232,+']
        self.assertCountEqual(expected, self.vencodes[1])

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))

    def test_export_ven_and_e(self):
        self.vencode_obj.export("vencodes", "e-values", verbose=False)
        folder_path = self.vencode_obj._parent_path
        paths = [os.path.join(folder_path, "celltypetarget_vencode.csv"),
                 os.path.join(folder_path, "celltypetarget_vencode-1.csv"),
                 os.path.join(folder_path, "celltypetarget_vencode-2.csv"),
                 os.path.join(folder_path, "celltypetarget_vencode-3.csv"),
                 os.path.join(folder_path, "celltypetarget_evalues.csv")]
        for i in paths:
            with self.subTest(i=i):
                self.assertTrue(os.path.exists(i))
                dh.remove_file(i)


class GetVencodesSamplingTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencode_obj = iext.GetVencodes(inputs=self.inputs,
                                            files_path=self.files_path,
                                            cell_type=self.celltype_analyse,
                                            algorithm="sampling",
                                            n_regulatory_elements=self.k,
                                            number_vencodes=4,
                                            thresholds=self.thresholds, n_samples=1000,
                                            merge={"replicate_suffix": self.replicate_suffix})
        self.vencodes = self.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
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
        self.vencode_obj = iext.GetVencodes(inputs=self.inputs,
                                            files_path=self.files_path,
                                            cell_type=cell_type,
                                            algorithm=self.algorithm,
                                            n_regulatory_elements=self.k,
                                            number_vencodes=4,
                                            thresholds=self.thresholds, n_samples=10000,
                                            merge=None, replicates=replicates)
        self.vencodes = self.vencode_obj.coordinates
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
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
        self.vencode_obj = iext.GetVencodes(inputs=self.inputs,
                                            files_path=self.files_path,
                                            cell_type=self.add_celltype_ctp,
                                            algorithm=self.algorithm,
                                            n_regulatory_elements=self.k,
                                            number_vencodes=4,
                                            thresholds=self.thresholds, n_samples=10000,
                                            merge={"replicate_suffix": self.replicate_suffix},
                                            add_celltype=add_celltype_)
        self.vencodes = self.vencode_obj.coordinates

    def test_correct_addition(self):
        for vencode_data in self.vencode_obj.expression_data:
            with self.subTest(i=vencode_data.index.values.tolist()):
                self.assertIn(self.add_celltype_ctp, vencode_data.columns)

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)


class GetVencodesValidateTest(GettingVencodesTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.validate_with = outside_data.Bed(source="get_val_ven_test.bed", files_path=self.files_path,
                                              folder="Validation_files")
        self.vencode_obj = iext.GetVencodes(validate_with=(self.validate_with, (":", "\..", ",")),
                                            inputs=self.inputs,
                                            files_path=self.files_path,
                                            cell_type=self.celltype_analyse,
                                            algorithm=self.algorithm,
                                            n_regulatory_elements=self.k,
                                            number_vencodes=4,
                                            thresholds=self.thresholds, n_samples=10000,
                                            merge={"replicate_suffix": self.replicate_suffix})
        self.vencodes = self.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))

    def test_validation(self):
        vencode_data = self.vencode_obj.expression_data[0]
        indexes = vencode_data.index
        _chr, temp = indexes[0].split(":")
        start, end = temp.split("..")
        end, strand = end.split(",")
        coordinates_set = range(int(start) - 200, int(end) + 201)
        condition = False
        for index, row in self.validate_with.data.iterrows():
            condition = (row.Start in coordinates_set) or (row.End in coordinates_set)
            if condition:
                break
        self.assertTrue(condition)


class GettingVencodesFantomTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        cls.celltype_analyse = "Adipocyte - Breast"
        cls.data_type = "promoters"
        cls.sample_type = "primary cells"
        cls.algorithm = "heuristic"
        cls.k = 4
        cls.thresholds = (0.5, 0, 0)  # act, inact, and sparseness, respectively
        cls.parsed = True
        cls.files_path = "test"


class GetVencodesFantomHeuristicTest(GettingVencodesFantomTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencode_obj = iext.GetVencodesFantom(files_path=self.files_path,
                                                  cell_type=self.celltype_analyse,
                                                  algorithm=self.algorithm,
                                                  n_regulatory_elements=self.k,
                                                  number_vencodes=4,
                                                  parsed=self.parsed,
                                                  thresholds=self.thresholds, n_samples=10000,
                                                  data_type=self.data_type, sample_type=self.sample_type)
        self.vencodes = self.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))


class GetVencodesFantomSamplingTest(GettingVencodesFantomTest):
    def setUp(self):
        """ Sets-up variables to be used in the tests. """
        self.vencode_obj = iext.GetVencodesFantom(files_path=self.files_path,
                                                  cell_type=self.celltype_analyse,
                                                  algorithm="sampling",
                                                  n_regulatory_elements=self.k,
                                                  number_vencodes=4,
                                                  parsed=self.parsed,
                                                  thresholds=self.thresholds, n_samples=10000,
                                                  data_type=self.data_type, sample_type=self.sample_type)
        self.vencodes = self.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))


class GetVencodesFantomAddCelltypeTest(GettingVencodesFantomTest):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        super().setUpClass()
        add_celltype_file = cv.test_promoter_file_name
        add_celltype_ = [add_celltype_file, "hIPS", {"sample_types": "time courses",
                                                     "files_path": cls.files_path}]
        cls.vencode_obj = iext.GetVencodesFantom(inputs=cv.test_promoter_file_name,
                                                 files_path=cls.files_path,
                                                 cell_type="hIPS",
                                                 algorithm=cls.algorithm,
                                                 n_regulatory_elements=cls.k,
                                                 number_vencodes=4,
                                                 parsed=False,
                                                 thresholds=cls.thresholds, n_samples=10000,
                                                 data_type=cls.data_type, sample_type=cls.sample_type,
                                                 add_celltype=add_celltype_, merge={"exclude_target": True})
        cls.vencodes = cls.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))


class GetVencodesFantomNotMergedAdd(GettingVencodesFantomTest):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        super().setUpClass()
        add_celltype_file = cv.test_promoter_file_name
        add_celltype_ = [add_celltype_file, "hIPS", {"sample_types": "time courses",
                                                     "files_path": cls.files_path}]
        cls.vencode_obj = iext.GetVencodesFantom(inputs=cv.test_promoter_file_name,
                                                 files_path=cls.files_path,
                                                 cell_type="hIPS",
                                                 algorithm=cls.algorithm,
                                                 n_regulatory_elements=cls.k,
                                                 number_vencodes=4,
                                                 parsed=False,
                                                 thresholds=cls.thresholds, n_samples=10000,
                                                 data_type=cls.data_type, sample_type=cls.sample_type,
                                                 add_celltype=add_celltype_)
        cls.vencodes = cls.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_if_in_target(self):
        for column in self.vencode_obj.target_replicates_data.columns:
            with self.subTest(column=column):
                self.assertRegex(column, "hIPS.*")

    def test_if_added(self):
        self.assertIn("hIPS, biol_rep2", self.vencode_obj.data.columns)


class GetVencodesFantomValidatedTest(GettingVencodesFantomTest):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        super().setUpClass()
        add_celltype_file = cv.test_enhancer_file_name
        add_celltype_ = [add_celltype_file, "hIPS", {"sample_types": "time courses",
                                                     "files_path": cls.files_path}]

        cls.validate_with = val.get_data_to_validate("hIPS", optional=False, files_path=cls.files_path)
        cls.vencode_obj = iext.GetVencodesFantom(validate_with=cls.validate_with,
                                                 inputs=cv.test_enhancer_file_name,
                                                 files_path=cls.files_path,
                                                 cell_type="hIPS",
                                                 algorithm=cls.algorithm,
                                                 n_regulatory_elements=cls.k,
                                                 number_vencodes=2,
                                                 parsed=False,
                                                 thresholds=cls.thresholds, n_samples=10000,
                                                 data_type="enhancers", sample_type=cls.sample_type,
                                                 add_celltype=add_celltype_,
                                                 merge={"exclude_target": True})
        cls.vencodes = cls.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.expression_data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))

    def test_validation(self):
        vencode_data = self.vencode_obj.expression_data[0]
        indexes = vencode_data.index
        _chr, temp = indexes[0].split(":")
        start, end = temp.split("-")
        coordinates_set = range(int(start) - 200, int(end) + 201)
        condition = False
        for index, row in self.validate_with.data.iterrows():
            condition = (row.Start in coordinates_set) or (row.End in coordinates_set)
            if condition:
                break
        self.assertTrue(condition)


class GetVencodesFantomExternalHeuristicTest(GettingVencodesFantomTest):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        super().setUpClass()
        cls.validate_with = val.get_data_to_validate("A549_singlecell", file_name="A549_singlecell.csv",
                                                     files_path=cls.files_path, positions=(0, 0, 0, 1),
                                                     splits=(":", "-"))
        cls.vencode_obj = iext.GetVencodeFantomExternalData(validate_with=cls.validate_with,
                                                            inputs=cv.test_enhancer_file_name,
                                                            files_path=cls.files_path,
                                                            cell_type="A549_singlecell",
                                                            algorithm=cls.algorithm,
                                                            n_regulatory_elements=cls.k,
                                                            number_vencodes=2,
                                                            parsed=False,
                                                            thresholds=cls.thresholds, n_samples=10000,
                                                            data_type="enhancers", sample_type="cell lines",
                                                            merge={"exclude_target": True})
        cls.vencodes = cls.vencode_obj.coordinates

    def test_if_correct_vencodes(self):
        for vencode_data in self.vencode_obj.data:
            vencode_data.drop(self.vencode_obj.target_replicates, axis=1, inplace=True)
            with self.subTest(i=vencode_data.index.values.tolist()):
                condition = self.vencode_obj._assess_vencode_one_zero_boolean(vencode_data)
                self.assertTrue(condition)

    def test_e_values_created(self):
        maximum, minimum = 100, 0
        for i in self.vencode_obj.e_values.values():
            with self.subTest(i=i):
                self.assertTrue(maximum >= i >= minimum, msg="{} is not between {} and {}".format(i, minimum, maximum))

    def test_validation(self):
        vencode_data = self.vencode_obj.data[1]
        indexes = vencode_data.index
        _chr, temp = indexes[0].split(":")
        start, end = temp.split("-")
        coordinates_set = range(int(start) - 200, int(end) + 201)
        condition = False
        for index, row in self.validate_with.data.iterrows():
            condition = (row.Start in coordinates_set) or (row.End in coordinates_set)
            if condition:
                break
        self.assertTrue(condition)


class CheckElementExpression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Sets-up class variables to be used in the tests.
        """
        celltype_analyse = "Adipocyte - breast"
        data_type = "promoters"
        sample_type = "primary cells"
        parsed = False
        files_path = "test"
        cls.element_list = ('chr10:100027943..100027958,-', 'chr10:100174900..100174956,-',
                            'chr10:100204220..100204230,-', 'chr10:100206642..100206717,-')
        expression_obj = iext.CheckElementExpression(inputs=cv.test_promoter_file_name,
                                                     element_list=cls.element_list,
                                                     cell_type=celltype_analyse,
                                                     data_type=data_type, sample_type=sample_type,
                                                     parsed=parsed, files_path=files_path)
        cls.expression = expression_obj.export_expression_data(method="return")

    def test_elements_in_results(self):
        self.assertCountEqual(self.element_list, self.expression.index)

    def test_celltype_in_results(self):
        columns = ["Adipocyte - breast, donor1", "Adipocyte - breast, donor2"]
        self.assertCountEqual(columns, self.expression.columns)

    def test_values(self):
        values = self.expression.isna()
        cond = all(values)
        self.assertTrue(cond)


if __name__ == "__main__":
    unittest.main(verbosity=2)
