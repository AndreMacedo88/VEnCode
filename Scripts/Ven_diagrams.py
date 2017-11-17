#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import Classes
from Common_variables import file_name, complete_cancer_cell_type, complete_primary_exclude_list, \
    complete_primary_non_include_list, complete_primary_jit_exclude_list

# region Variables
case_studies = ["small cell lung carcinoma", "testicular germ cell embryonal carcinoma",
                "chronic myelo leukemia", "acute myeloid leukemia (FAB M2)", "teratocarcinoma"]
# endregion Variables

if __name__ == "__main__":
    initialize_promoters = Classes.Promoters(file_name,
                                             case_studies[3],
                                             celltype_exclude=complete_primary_exclude_list,
                                             not_include=complete_primary_non_include_list,
                                             partial_exclude=complete_primary_jit_exclude_list,
                                             sample_types=["primary cells", "cell lines"],
                                             second_parser="primary cells")
    # get files for VEn diagram:
    initialize_promoters.ven_diagrams(50000, 4, threshold=50)

    # test if biologically (statistically) relevant:
    # initialize_promoters.statistics_ven_diagram(5000, 5, 22, threshold=50)
