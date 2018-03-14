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

initialize_promoters = Classes.Promoters(promoter_file_name, complete_primary_cell_list,
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=False, log_level="info")

initialize_promoters.find_vencodes_each_celltype(stop=5)
