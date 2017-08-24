import os
from unittest import TestCase
import Utils


class TestFiles_in_folder_to_csv(TestCase):

    def setUp(self):
        self.folder = "/Figure 3-b2/cell lines/"

    def test_file_creation(self):
        file_name = "all_cell_lines.csv"
        Utils.files_in_folder_to_csv(self.folder,file_name)

if __name__ == "__main__":
    test = TestFiles_in_folder_to_csv()
    test.test_file_creation()