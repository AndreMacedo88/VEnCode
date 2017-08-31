#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Best VEnCode.py: Functions for generating the best VEnCode """

import Classes
from Common_variables import complete_primary_cell_list, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

file_name = "hg19.cage_peak_phase1and2combined_tpm.osc.txt"

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name, complete_primary_cell_list,
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types="primary cells",
                                             second_parser=None, nrows=20)
    initialize_promoters.best_vencode_generator()
