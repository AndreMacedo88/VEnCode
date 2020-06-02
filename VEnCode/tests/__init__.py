"""
Module containing all the tests for the VEnCode module.
"""

import unittest


def run_all_tests(verbosity=2):
    """
    Runs all tests in the test module in one go.

    Parameters
    ----------
    verbosity : int
        The amount of verbosity output by the tests.
    """
    from VEnCode.tests import test_internals

    suite = unittest.TestLoader().loadTestsFromModule(test_internals)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
