#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""inter_robustness.py: Functions for generating inter robustness data """

import Classes
from Common_variables import file_name, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list, \
    three_donors_cell_list, four_donors_cell_list

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name,
                                             three_donors_cell_list,
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types="primary cells",
                                             second_parser=None)
    """ get the percentage of VEnCodes taken for 1 donor that work for all donors: """
    initialize_promoters.ven_diagram_interception(2000, 5, 3, combinations_number=4, threshold=90)

    """ get the percentage of VEnCodes taken for 1, 2, 3, etc donors that work for all and com """
    # initialize_promoters.inter_donor_percentage_difference(2000, 3, 4, combinations_number=4,
    #                                                        threshold=90)

