#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
At least one VEnCode.py: Functions for finding which cell types have at least one VEnCode
"""
import sys
import os

sys.path.append(os.path.abspath(os.path.join('..', '')))

import Classes
from Common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# Promoters

initialize_promoters = Classes.Promoters(promoter_file_name, complete_primary_cell_list,
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=True, log_level="info")

initialize_promoters.find_vencodes_each_celltype(stop=5, combinations_number=[1, 2, 3, 4], method="sampling",
                                                 n_samples=10000, threshold_inactivity=0.1)

# Enhancers
"""
initialize_enhancers = Classes.Promoters(enhancer_file_name,
                                         ["CD34+ Progenitors"],
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=False, log_level="info", enhancers=enhancer_names_db,
                                         skiprows=None)

initialize_enhancers.find_vencodes_each_celltype(stop=5, combinations_number=[1, 2, 3, 4, 5], method="heuristic",
                                                 n_samples=10000, expression=0.5)
"""
