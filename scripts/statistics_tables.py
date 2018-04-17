#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""statistics_tables.py: Functions for generating tables to use in statistic analyses."""
import pandas as pd
import numpy as np

import utils.directory_handlers
from utils import util
import classes

if __name__ == "__main__":
    celltype_number = [20, 80, 100, 154, 200, 250, 350, 450, 550, 650, 800, 1000]
    promoter_number = range(1, 11)  # number of promoters ranging from x to y
    e_values = pd.DataFrame(index=promoter_number, columns=celltype_number)
    for i in celltype_number:
        print("Starting number of cell types: {}".format(i))
        for z in promoter_number:
            print("Starting number of promoters: {}".format(z))
            data = pd.DataFrame(np.zeros(shape=(z, i)), dtype=np.int8)
            e_value = classes.Promoters.e_value_calculator(data, reps=1000)
            e_values.loc[z, i] = e_value
    file_name = utils.directory_handlers.file_directory_handler("Table for e_value statistics_final.csv", "/Files/", path_type="parent")
    with open(file_name, 'w') as f:
        e_values.to_csv(f, sep=";")
