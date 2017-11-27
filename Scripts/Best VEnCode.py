#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Best VEnCode.py: Functions for generating the best VEnCode """

import Classes
from Common_variables import file_name, complete_primary_cell_list, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name, complete_primary_cell_list,
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types="primary cells", second_parser=None,
                                             nrows=2000, log_level="debug")
    initialize_promoters.best_vencode_generator("Hepatocyte")

    # some tests:
    """
    initialize_promoters.test_vencode_data(rows=('chr10:114581583..114581600,-', 'chr10:101841398..101841415,-'),
                                           file_name="test_hepatocyte_vencode.csv")
    
    initialize_promoters.test_vencode_data(('chr10:114581583..114581600,-', 'chr10:101841398..101841415,-'),
                                           (initialize_promoters.codes["Hepatocyte"]))
    """
