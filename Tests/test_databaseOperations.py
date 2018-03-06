from unittest import TestCase, main
from .context import Classes
from Common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list


class EvalueNormalizerTest(TestCase):
    def setUp(self):
        """Define a few normalized e-values."""
        self.excess = []
        self.near_maximum = []
        self.zero = []
        self.fifty_percent = None

        def small_database(n=80, k=4):
            """tests a small sized db"""
            self.excess.append(Classes.Promoters.e_value_normalizer(1000, n, k))
            self.near_maximum.append(Classes.Promoters.e_value_normalizer(98, n, k))
            self.zero.append(Classes.Promoters.e_value_normalizer(0, n, k))

        def standard_database(n=154, k=4):
            """tests a standard sized db"""
            self.excess.append(Classes.Promoters.e_value_normalizer(1000, n, k))
            self.near_maximum.append(Classes.Promoters.e_value_normalizer(160, n, k))
            self.zero.append(Classes.Promoters.e_value_normalizer(0, n, k))

        def large_database(n=800, k=8):
            """tests a large sized db"""
            self.excess.append(Classes.Promoters.e_value_normalizer(3000, n, k))
            self.near_maximum.append(Classes.Promoters.e_value_normalizer(2620, n, k))
            self.zero.append(Classes.Promoters.e_value_normalizer(0, n, k))

        def database_cut_half():
            """tests a database when it's half full row- or column-wise"""
            cut_rows = Classes.Promoters.e_value_normalizer(25.3866, 200, 4)
            cut_cols = Classes.Promoters.e_value_normalizer(117.1372, 200, 4)
            self.fifty_percent = (cut_rows, cut_cols)

        small_database(), standard_database(), large_database(), database_cut_half()

    def test_more_than_hundred_percent(self):
        """When raw e-values reach more than the maximum predicted by the model, should return 100."""
        for e_value in self.excess:
            self.assertEqual(e_value, 100)

    def test_around_hundred_percent(self):
        """When raw e-values are near or at the maximum predicted by the model, should return 100 +/-10."""
        for e_value in self.near_maximum:
            self.assertGreater(e_value, 90)
            self.assertLessEqual(e_value, 100)

    def test_e_value_zero(self):
        """ When raw e-values are zero, should normalize to zero"""
        for e_value in self.zero:
            self.assertEqual(e_value, 0)

    def test_fifty_percent(self):
        """ When database is about half full (row- or column-wise) should normalize to a predictable number"""
        for e_value, expected in zip(self.fifty_percent, (12.8831, 59.4446)):
            self.assertAlmostEqual(e_value, expected, places=3)

if __name__ == '__main__':  # run tests if called from command-line
    main()
