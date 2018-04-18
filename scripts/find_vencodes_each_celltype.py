#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
find_vencodes_each_celltype.py: Script to finding which cell types have at least one VEnCode
"""
import sys
import os

import utils.directory_handlers as directory_handlers
import utils.writing_files as writing_files

sys.path.append(os.path.abspath(os.path.join('..', '')))

import classes
from common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# Promoters
"""
initialize_promoters = classes.Promoters(promoter_file_name, complete_primary_cell_list,
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=True, log_level="info")

results = initialize_promoters.find_vencodes_each_celltype(stop=5, combinations_number=range(1, 8),
                                                           method="heuristic",
                                                           n_samples=10000, threshold_inactivity=0,
                                                           threshold_activity=0.5)
"""

# Enhancers

initialize_enhancers = classes.Promoters(enhancer_file_name,
                                         complete_primary_cell_list,
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=True, log_level="info", enhancers=enhancer_names_db,
                                         skiprows=None, skip_raw_data=True)

results = initialize_enhancers.find_vencodes_each_celltype(stop=5, combinations_number=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                           method="sampling", n_samples=10000,
                                                           threshold_inactivity=0, threshold_activity=0.5)

results_directory = directory_handlers.check_if_and_makefile(os.path.join("VEnCode Search", "All CellTps VEnC search"),
                                                             path_type="parent2")
writing_files.write_dict_to_csv(results_directory, results, deprecated=False)
