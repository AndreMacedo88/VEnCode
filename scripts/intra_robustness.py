#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""intra_robustness.py: Functions for generating intra robustness data """
import os

import utils.directory_handlers as directory_handlers
import utils.writing_files as writing_files
import classes
from common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list


class Setup:
    """
    Sets the variables and other
    """

    def __init__(self):
        self.cell_list = complete_primary_cell_list
        self.vens_to_take = 20
        self.combinations_number = 4
        self.threshold = 90


if __name__ == "__main__":
    Setup = Setup()
    initialize_promoters = classes.Promoters(promoter_file_name, Setup.cell_list,
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types="primary cells", second_parser=None,
                                             conservative=True, skip_raw_data=True)
    """ All cell types 
    # use: cell_list = complete_primary_cell_list
    initialize_promoters.codes_to_csv("codes_all_cells.csv", "list", "/Figure 2/Test codes/")
    initialize_promoters.celltypes_to_csv("celltypes_all.csv", "list", "/Figure 2/Test codes/")
    """

    """ 3 Donors
    # use: cell_list = three_donors_cell_list
    initialize_promoters.codes_to_csv("codes_3_donors.csv", "list", "/Figure 2/Test codes/")
    initialize_promoters.celltypes_to_csv("celltypes_3_donors.csv", "list", "/Figure 2/Test codes/")
    """

    """ 4 Donors
    # use: cell_list = four_donors_cell_list
    initialize_promoters.codes_to_csv("codes_4_donors.csv", "list", "/Figure 2/Test codes/")
    initialize_promoters.celltypes_to_csv("celltypes_4_donors.csv", "list", "/Figure 2/Test codes/")
    """

    results = initialize_promoters.intra_individual_robustness(Setup.combinations_number, Setup.vens_to_take,
                                                               threshold=Setup.threshold)

    results_directory = directory_handlers.check_if_and_makefile(
        os.path.join("Figure 2", "VEnCode E-values {} samples {} VEnCodes".format(len(Setup.cell_list),
                                                                                  Setup.vens_to_take)),
        path_type="parent2")
    writing_files.write_dict_to_csv(results_directory, results, deprecated=False)
