#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import Classes
from Common_variables import file_name, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name,
                                             "acute myeloid leukemia",
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types=["primary cells", "cell lines"],
                                             second_parser="primary cells")

    # get vencodes:
    initialize_promoters.get_vencodes(n=1, write_file=True)
