#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
sampling_vs_heuristic.py: Script to compare VEnCodes between the sampling and heuristic algorithms.
"""
import matplotlib.pyplot as plt

import utils.writing_files as wfs
import utils.directory_handlers as dhs
import utils.input_handlers as ihs
import classes
from common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

rows_number = ihs.input_integers("Number of rows from the file to open: ")
vencodes_number = ihs.input_integers("Number of VEnCodes to get: ")
algorithm = ihs.input_string("Algorithm(s) to use: ")

initialize_promoters = classes.Promoters(promoter_file_name, complete_primary_cell_list,
                                         celltype_exclude=complete_primary_exclude_list,
                                         not_include=complete_primary_non_include_list,
                                         partial_exclude=complete_primary_jit_exclude_list,
                                         sample_types="primary cells", second_parser=None,
                                         conservative=True, log_level="info", nrows=rows_number)

vencodes = initialize_promoters.vencode_generator("Hepatocyte", algorithm=algorithm, combinations_number=4,
                                                  threshold_activity=1, threshold_inactivity=0,
                                                  number_vencodes=vencodes_number)
data = []
labels = []
for key, value in vencodes.items():
    results_directory = dhs.check_if_and_makefile(r"/sampling_vs_heuristic/vencodes {} method". format(key),
                                                  path_type="parent")
    wfs.write_dict_to_csv(results_directory, value, deprecated=False)
    data.append(list(value.values()))
    labels.append(key)
plt.boxplot(data)
"""
xtickNames = plt.setp(ax1, xticklabels=np.repeat(randomDists, 2))
plt.setp(xtickNames, rotation=45, fontsize=8)
"""
plt.show()
