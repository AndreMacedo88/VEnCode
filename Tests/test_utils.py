import unittest
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
from Utils import Util
from Utils import pandas_utils as pdutil


class FilterByExpressionAndPercentileTest(unittest.TestCase):
    def setUp(self):
        """Produce some filtered data sets."""
        data_pre = np.fromfunction(lambda x, y: np.sin(x) * np.tan(y) * x * 10, (400, 154))
        for array in data_pre:
            array[array < 0] = 0
        self.data = pd.DataFrame(data_pre)
        self.data.columns = self.data.columns.astype(str)
        print(self.data)

    def test_one(self):
        prediction = (199, 154)
        result = Util.df_filter_by_expression_and_percentile(self.data, "90", 1, 2, threshold=90)
        print(result)
        self.assertEqual(prediction, result.shape, msg="Result not expected")


class AssessVencodeOneZeroBooleanTest(unittest.TestCase):
    def setUp(self):
        """ Define some samples to assess if VEnCode. """
        data_pre = {"names": ["one", "two", "three", "four"],
                    "col1": [0, 1, 1, 0], "col2": [0, 1, 0, 1], "col3": [0, 0, 1, 1],
                    "col4": [0, 1, 0, 0]}
        self.test_pandas_df = pd.DataFrame(data=data_pre).set_index("names").astype("float")

    def test_threshold_zero_success(self):
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0)
        self.assertTrue(test, msg="Result not expected")

    def test_threshold_zero_fail(self):
        pdutil.multi_set_data_frame(self.test_pandas_df, (("one", "col1"), ("four", "col1")), 1)
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0)
        self.assertFalse(test, msg="Result not expected")

    def test_threshold_zero_one_success(self):
        pdutil.multi_set_data_frame(self.test_pandas_df, (("one", "col1"), ("four", "col1")), 0.1)
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0.1)
        self.assertTrue(test, msg="Result not expected")

    def test_threshold_zero_one_fail(self):
        pdutil.multi_set_data_frame(self.test_pandas_df, (("one", "col1"), ("four", "col1")), 1)
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0.1)
        self.assertFalse(test, msg="Result not expected")

    def test_threshold_zero_zero_one_success(self):
        pdutil.multi_set_data_frame(self.test_pandas_df, (("one", "col1"), ("four", "col1")), 0.01)
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0.01)
        self.assertTrue(test, msg="Result not expected")

    def test_threshold_zero_zero_one_fail(self):
        pdutil.multi_set_data_frame(self.test_pandas_df, (("one", "col1"), ("four", "col1")), 0.1)
        test = Util.assess_vencode_one_zero_boolean(self.test_pandas_df, threshold=0.01)
        self.assertFalse(test, msg="Result not expected")

class TestCFNL:
    """ Tests function combinations_from_nested_lists"""
    def __init__(self, lst):
        self.list = lst

    def general_testing(self):
        """
        runs normal function
        """
        for i in Util.combinations_from_nested_lists(self.list):
            print(i)

class TestMakingWritingFiles:
    def test_dicts(self):
        results_dict = {"aa": [1,2], "bb": [2,4]}
        results_directory = Util.check_if_and_makefile(r"/example/testing dicts",
                                                       path_type="parent")
        # Util.check_if_and_makedir(results_directory)
        Util.write_dict_to_csv(results_directory, results_dict, deprecated=False)

if __name__ == "__main__":
    unittest.main()

    """
    test = TestFiles_in_folder_to_csv()
    test.test_file_creation()
    
    lst1 = ["aa", "bb", "cc"]
    lst2 = ["aa", ["bb", "cc"]]
    lst3 = ["aa", ["bb", ["cc"]]]
    lst4 = ["aa", ["bb", ["cc", "dd"]]]
    lst5 = ["aa", ["bb", ["cc", ["dd"]]]]
    lst6 = ["aa", ["bb", ["cc", ["dd", "ee"]]]]
    test = TestCFNL(lst6)
    test.general_testing()
    """
