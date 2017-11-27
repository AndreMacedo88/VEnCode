#!/usr/bin/env python

"""Tests classes for VEnCode"""

from .context import Classes
from Common_variables import file_name, complete_primary_cell_list, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# region "test enhancer class methods"
"""
x = Classes.Enhancers("human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",
                      "Human.sample_name2library_id.txt", complete_primary_cell_list,
                      celltype_exclude=complete_primary_exclude_list,
                      not_include=complete_primary_non_include_list,
                      partial_exclude=complete_primary_jit_exclude_list)
for key, value in x.partial_exclude_codes.items():
    print("Enhancers", key, value, len(value), sep="\n", end="\n\n")
print(x.celltype, x.raw_data.shape, x.sample_types, x.file, x.parent_path, x.data.shape, x.codes, sep="\n", end="\n\n")
"""
# endregion "test enhancer class methods"

# region "test promoter class methods"
"""
y = Classes.Promoters("hg19.cage_peak_phase1and2combined_tpm.osc.txt", complete_primary_cell_list,
                      celltype_exclude=complete_primary_exclude_list, not_include=complete_primary_non_include_list,
                      partial_exclude=complete_primary_jit_exclude_list, nrows=1000)
for key, value in y.partial_exclude_codes.items():
    print("Promoters", key, value, len(value), sep="\n", end="\n\n")

# print(y.codes, len([item for items in y.codes.values() for item in items]), sep="\n")
# Defs.write_dict_to_csv("test.csv", y.codes, "/Tests/", path="parent")
# print(y.celltype, y.raw_data.shape, y.sample_types, y.file, y.parent_path, y.data.shape, y.codes, sep="\n")
"""

# endregion "test promoter class methods"
