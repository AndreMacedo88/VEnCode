#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
check_element_expression.py: file used to check expression levels of certain promoter and enhancer regions
in specific cell types of the FANTOM5 data.
"""

import os
import sys

file_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(file_dir)

from VEnCode import internals_extensions as iext

element_list = ("chr9:99976780-99977107", "chr1:213090277-213090696", "chr6:121940262-121940519",
                "chr7:157103071-157103446", "chr7:131290113-13129045", "chr7:51395107-51395355",
                "chr15:65101103-65101281", "chr12:31523933-31524084", "chr5:93697867-93698390",
                "chrX:40016662-40016910", "chrX:40028554-40028949", "chr6:15022090-15022421")

cell_type = "epitheloid carcinoma cell line: HelaS3 ENCODE"
cell_type2 = "embryonic kidney cell line: HEK293/SLAM untreated"
path_out = "D:/Utilizador HDD/OneDrive - Nova Medical School Faculdade de Ciências Médicas da UNL/1-Research/3-Vencode/Fantom5/Expression Levels/"

path_out = path_out + cell_type2.replace(":", "-").replace("/", "-") + ".csv"

if __name__ == "__main__":
    expression = iext.CheckElementExpression(element_list=element_list,
                                             cell_type=cell_type2,
                                             data_type="enhancers",
                                             parsed=False)
    expression.export_expression_data(path=path_out)
