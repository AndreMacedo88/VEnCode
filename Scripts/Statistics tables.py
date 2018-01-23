#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Statistics tables.py: Functions for generating tables to use in statistic analyses."""
import pandas as pd
import numpy as np

from Utils import Util
import Classes

if __name__ == "__main__":
    celltype_number = [20, 80, 100, 154, 200, 250, 350, 450, 550, 650, 750]
    promoter_number = range(3, 6)  # number of promoters ranging from x to y
    e_values = pd.DataFrame(index=promoter_number, columns=celltype_number)
    for i in celltype_number:
        for z in promoter_number:
            data = pd.DataFrame(np.zeros(shape=(z, i)), dtype=np.int8)
            e_value = Classes.Promoters.e_value_calculator(data, reps=100)
            e_values.loc[z, i] = e_value
    file_name = Util.file_directory_handler("Table for e_value statistics_quick.csv", "/Files/", path="parent")
    with open(file_name, 'w') as f:
        e_values.to_csv(f, sep=";")
