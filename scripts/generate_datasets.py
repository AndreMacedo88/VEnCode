#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
generate_datasets.py: Script to generate files with parsed data sets.
"""
import os
import pandas as pd
from tqdm import tqdm

import classes
from utils import directory_handlers as dhs
from utils import util
from common_variables import promoter_file_name, enhancer_file_name, enhancer_names_db, complete_primary_cell_list, \
    complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# Promoters:
"""
init = classes.Promoters(promoter_file_name, complete_primary_cell_list,
                         celltype_exclude=complete_primary_exclude_list,
                         not_include=complete_primary_non_include_list,
                         partial_exclude=complete_primary_jit_exclude_list,
                         sample_types="primary cells", second_parser=None,
                         conservative=True, log_level="info", nrows=None)
"""
# Enhancers:

init = classes.Promoters(enhancer_file_name,
                         complete_primary_cell_list,
                         celltype_exclude=complete_primary_exclude_list,
                         not_include=complete_primary_non_include_list,
                         partial_exclude=complete_primary_jit_exclude_list,
                         sample_types="primary cells", second_parser=None,
                         conservative=True, log_level="info", enhancers=enhancer_names_db,
                         skiprows=None, nrows=None)

# Prepare directory:
folder = os.sep.join(["Files", "Dbs"])
parent_path = os.sep.join([os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.pardir)),
                           folder])
dhs.check_if_and_makedir(parent_path)

data_copy = init.data.copy()
init.data = init.merge_donors_into_celltypes()
for celltype in tqdm(complete_primary_cell_list):
    # file name:
    filename = "{}_tpm_enhancers".format(celltype)
    file_directory = dhs.check_if_and_makefile(filename, path=parent_path, file_type=".csv")

    data_celltype = init.data[celltype]
    init.data.drop(celltype, axis=1, inplace=True)
    init.data = pd.concat([init.data, data_copy[init.codes[celltype]]],
                          axis=1)
    data_set = util.df_filter_by_expression(init.data, init.codes[celltype], 0.0001)
    data_set.to_csv(file_directory, sep=";")

    init.data.drop(data_copy[init.codes[celltype]], axis=1, inplace=True)
    init.data = pd.concat([init.data, data_celltype], axis=1)
