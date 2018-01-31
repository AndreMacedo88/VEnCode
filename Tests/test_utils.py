from Utils import Util


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


if __name__ == "__main__":
    """
    test = TestFiles_in_folder_to_csv()
    test.test_file_creation()
    """
    lst1 = ["aa", "bb", "cc"]
    lst2 = ["aa", ["bb", "cc"]]
    lst3 = ["aa", ["bb", ["cc"]]]
    lst4 = ["aa", ["bb", ["cc", "dd"]]]
    lst5 = ["aa", ["bb", ["cc", ["dd"]]]]
    lst6 = ["aa", ["bb", ["cc", ["dd", "ee"]]]]
    test = TestCFNL(lst6)
    test.general_testing()
