#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""intra_robustness.py: Functions for generating intra robustness data """

import Classes
from Common_variables import file_name, complete_primary_cell_list, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# region "Setup Variables"

cell_list = complete_primary_cell_list
vens_to_take = 20
combinations_number = 4
threshold = 90

# endregion "Global variables"

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name, cell_list,
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types="primary cells",
                                             second_parser=None)
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

    initialize_promoters.intra_individual_robustness(combinations_number, vens_to_take, threshold=threshold)

