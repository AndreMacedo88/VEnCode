#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Classes.py: Classes module for the VEnCode project """

# import csv
import os
import re
import logging

import pandas as pd
import numpy as np
from scipy.special import comb

import Defs


class DatabaseOperations:
    """ A class with classes for each type of database and common methods to all """

    def __init__(self, file, celltype, celltype_exclude=None, not_include=None, sample_types="primary cells",
                 skiprows=None, second_parser=None, nrows=None):
        self.file = file
        self.celltype = celltype
        self.celltype_exclude = celltype_exclude
        self.not_include = not_include
        self.sample_types = sample_types
        self.second_parser = second_parser
        self.parent_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) + "/Files/"
        # First, import only the necessary data
        self.raw_data = pd.read_csv(self.parent_path + self.file, sep="\t", index_col=0, skiprows=skiprows, nrows=nrows)

    @staticmethod
    def sample_category_selector(sample_types_file, types, path="parent", get="index"):
        """
        Returns a list of cell types to keep/drop from a file containing the list of cell types and a 'Sample category'
        column which determines which cell types to retrieve.
        """
        if not isinstance(types, list): types = [types]
        if path == "parent":
            parent_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
            database = pd.read_csv(parent_path + "/Files/" + sample_types_file, sep=",", index_col=0)
        elif path == "normal":
            database = pd.read_csv("./Files/" + sample_types_file, sep=",", index_col=0)
        else:
            raise Exception("path type is not recognized")
        try:
            possible_types = database["Sample category"].drop_duplicates().values.tolist()
        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)
            raise
        assert all(
            sample in possible_types for sample in types), "Sample type is not valid.\nValid sample types: {}".format(
            possible_types)
        celltypes = []
        for sample in types:
            selected = database[database["Sample category"] == sample]
            if get == "index":
                [celltypes.append(value) for value in selected.index.values]
            elif get == "name":
                [celltypes.append(value) for value in selected["Name"].tolist()]
            else:
                pass
        return celltypes

    @staticmethod
    def dataframe_regex_searcher(string, database):
        """ Returns a list containing only the columns of a dataframe which contain the string somewhere
        in its label """
        regular = ".*" + string.replace(" ", ".*").replace("+", "%2b") + ".*"
        idx = database.columns.str.contains(regular, flags=re.IGNORECASE, regex=True, na=False)
        regex_filtered_df = database.loc[:, idx]
        regex_filtered = regex_filtered_df.columns.values.tolist()
        return regex_filtered

    @staticmethod
    def test_codes(codes, celltype, codes_type="list"):
        """ Tests if any codes were generated """
        if codes_type == "list":
            if not codes:
                raise Exception("No codes for {}!".format(celltype))
        elif codes_type == "dict":
            if bool([a for a in codes.values() if a == []]):
                print([item for item, value in codes.items() if not value])
                raise Exception("Some celltypes might not have had codes generated!")
        elif codes_type == "ndarray":
            if codes.size == 0:
                raise Exception("No codes for {}!".format(celltype))
        else:
            raise Exception("Wrong codes type to test for the generation of codes!")

    @staticmethod
    def logging(specific_path):
        parent_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        logging.basicConfig(filename=parent_path + specific_path, level=logging.DEBUG,
                            format="{asctime} - {levelname} - {message}", filemode='w', style="{")
        logger = logging.getLogger(__name__)
        return logger

    @staticmethod
    def sorted_ven_robustness_test(data, codes, celltype, combinations_number, samples_to_take, reps, file_type,
                                   threshold=90, expression=1, multi_plot=False,
                                   include_problems=False, folder="/Figure 1/"):
        if not isinstance(combinations_number, list):
            combinations_number_list = range(1, combinations_number + 1)
        else:
            combinations_number_list = combinations_number
        filter_2 = Defs.fantom_filters(data, codes, expression, 2, threshold)
        if include_problems:
            k_ven_percent, problems = Defs.fantom_sampling(codes, celltype, filter_2, combinations_number_list,
                                                           samples_to_take,
                                                           reps, include_problems=include_problems)
        else:
            k_ven_percent = Defs.fantom_sampling(codes, celltype, filter_2, combinations_number_list,
                                                 samples_to_take,
                                                 reps, include_problems=include_problems)
        # Saving to csv or plotting:
        folder = folder
        if not multi_plot:  # multi_plot is there in case this function is used to generate other plots after.
            file_name = u"/VEnC - {1:s} - {2:s} - exp bigger or = {4:d} - {3:d}x {0:d} samples of k".format(
                samples_to_take, file_type, celltype, reps, expression)
            title = "Probability of VEnCode from sample of size k \n {0:s} expression >= {1:d}".format(
                celltype, expression)
            Defs.write_dict_to_csv(file_name + ".csv", k_ven_percent, folder)
            fig, path = Defs.errorbar_plot(k_ven_percent, folder, file_name, label=celltype, title=title)
            fig.savefig(path)
        if include_problems:
            logging.info("{}: {}".format(celltype, problems))
            new_file_name = u"/Probs for {} - {}x {} samples of k-{}".format(celltype, reps, samples_to_take,
                                                                             combinations_number)
            Defs.write_dict_to_csv(new_file_name + ".csv", problems, folder, path="parent")
        return k_ven_percent

    def not_include_code_getter(self, not_include, df):
        if isinstance(not_include, list):
            not_codes = []
            for item in not_include:
                not_codes_item = self.dataframe_regex_searcher(item, df)
                not_codes.append(not_codes_item)
            not_codes = [item for sublist in not_codes for item in sublist]
        else:
            not_codes = self.dataframe_regex_searcher(not_include, df)
        return not_codes

    def non_combinatory_loop(self, class_instance, combinations_number, expression, threshold, multi=False,
                             mode="count"):
        # Start the loop through all cell types:
        calculated_vencodes = {}
        for celltype in self.celltype:
            print("", "Starting {}".format(celltype), sep="-> ")
            logging.info("Starting %s", celltype)
            codes = class_instance.codes.get(celltype)
            filter_2 = Defs.fantom_filters(class_instance.data, codes, expression, 2, threshold)
            filter_2 = filter_2.applymap(lambda x: 0 if x == 0 else 1)  # change expression to 1 and 0 for quickness
            nodes = Defs.node_calculator(filter_2.drop(codes, axis=1))  # the time consumer!
            if multi:
                for number in range(2, combinations_number + 1):
                    try:
                        calculated_vencodes[number]
                    except KeyError:
                        calculated_vencodes[number] = []
                    ven_combinations = Defs.number_of_combination_from_nodes(nodes, len(filter_2.index), number)
                    logging.info("Number of VEnCodes found: %s for k=%s", ven_combinations, number)
                    if mode == "count":
                        calculated_vencodes[number].append(ven_combinations)
                    elif mode == "percentage":
                        total_comb = comb(len(filter_2.index), combinations_number, exact=False)
                        calculated_vencodes[number] = ven_combinations / total_comb * 100
            else:
                ven_combinations = Defs.number_of_combination_from_nodes(nodes, len(filter_2.index),
                                                                         combinations_number)
                logging.info("Number of VEnCodes found: %s for k=%s", ven_combinations, combinations_number)
                if mode == "count":
                    calculated_vencodes[celltype] = ven_combinations
                elif mode == "percentage":
                    total_comb = comb(len(filter_2.index), combinations_number, exact=False)
                    calculated_vencodes[celltype] = ven_combinations / total_comb * 100
        return calculated_vencodes

    def at_least_one_vencode(self, class_instance, combinations_number, expression, threshold):
        # Start the loop through all cell types:
        vencodes = {}
        for celltype in self.celltype:
            print("", "Starting {}".format(celltype), sep="-> ")
            logging.info("Starting %s", celltype)
            vencodes[celltype] = 0
            codes = class_instance.codes.get(celltype)
            try:
                if class_instance.partial_exclude_codes is not None:
                    if celltype in class_instance.partial_exclude_codes.keys():
                        partial_exclude_codes = class_instance.partial_exclude_codes.get(celltype)
                        data_frame = class_instance.data.drop(partial_exclude_codes, axis=1)
                    else:
                        data_frame = class_instance.data
            except AttributeError:
                data_frame = class_instance.data
                # else:
                # data_frame = class_instance.data
            filter_2 = Defs.fantom_filters(data_frame, codes, expression, 2, threshold)
            filter_2 = filter_2.applymap(lambda x: 0 if x == 0 else 1)  # change expression to 1 and 0 for quickness
            filter_2["sum"] = filter_2.sum(axis=1)
            filter_2.sort_values(["sum"], inplace=True)
            to_get_nodes = filter_2.drop(codes + ["sum"], axis=1)
            promoters = to_get_nodes.index.values
            for promoter in promoters:
                breaks = {}
                for item in range(1, combinations_number):
                    breaks["breaker_" + str(item)] = 0
                nodes = self.at_least_one_node_calculator_2(to_get_nodes, promoter,
                                                            combinations_number=combinations_number, breaks=breaks)
                if nodes > 0:
                    vencodes[celltype] = 1
                    break
                else:
                    to_get_nodes.drop(promoter, axis=0, inplace=True)
            print(nodes)
        return vencodes

    def at_least_one_node_calculator(self, data_frame, promoter, combinations_number=4, counter=1, breaker_1=0,
                                     breaker_2=0,
                                     breaker_3=0,
                                     breaker_4=0):
        # print("counter: " str(counter), end="\n")
        cols = data_frame.loc[promoter] != 0
        cols = data_frame.columns[cols]
        new_df = data_frame[cols].drop(promoter, axis=0)
        nodes = (new_df == 0).all(axis=1)
        node_count = np.sum(nodes)
        if node_count > 0:
            return node_count
        else:
            new_df["sum"] = new_df.sum(axis=1)
            new_df.sort_values(["sum"], inplace=True)
            new_df = new_df.drop(["sum"], axis=1)
            counter = counter
            # region early quit variables, remove in final
            breakr_1 = breaker_1
            breakr_2 = breaker_2
            breakr_3 = breaker_3
            breakr_4 = breaker_4
            # endregion of early quit variables
            while counter < combinations_number:
                counter += 1
                for prom in new_df.index.values:
                    # region "early quit if loop is taking too long"
                    if counter == combinations_number - 3:
                        breakr_1 += 1
                        if breakr_1 == 5:
                            node_count = 0
                            return node_count
                    if counter == combinations_number - 2:
                        breakr_1 += 1
                        if breakr_1 == 5:
                            node_count = 0
                            return node_count
                    elif counter == combinations_number - 1:
                        breakr_2 += 1
                        if breakr_2 == 5:
                            node_count = 0
                            return node_count
                    elif counter == combinations_number:
                        breakr_3 += 1
                        if breakr_3 == 5:
                            node_count = 0
                            return node_count
                    # endregion "early quit if loop is taking too long"
                    node_count = self.at_least_one_node_calculator(new_df, prom,
                                                                   combinations_number=combinations_number,
                                                                   counter=counter, breaker_1=breakr_1,
                                                                   breaker_2=breakr_2, breaker_3=breakr_3,
                                                                   breaker_4=breakr_4)
                    try:
                        if node_count > 1:
                            return node_count
                    except TypeError:
                        node_count = 0
                        continue
            return node_count

    def at_least_one_node_calculator_2(self, data_frame, promoter, combinations_number=4, counter=1, breaks=None):
        # print("counter: " str(counter), end="\n")
        cols = data_frame.loc[promoter] != 0
        cols = data_frame.columns[cols]
        new_df = data_frame[cols].drop(promoter, axis=0)
        nodes = (new_df == 0).all(axis=1)
        node_count = np.sum(nodes)
        if node_count > 0:
            return node_count
        else:
            new_df["sum"] = new_df.sum(axis=1)
            new_df.sort_values(["sum"], inplace=True)
            new_df = new_df.drop(["sum"], axis=1)
            counter = counter
            counter_thresholds = [i for i in range(2, (combinations_number + 1))]
            while counter < combinations_number:
                counter += 1
                for prom in new_df.index.values:
                    # region "early quit if loop is taking too long"
                    if counter in counter_thresholds:
                        breaker_index = str(counter_thresholds.index(counter) + 1)
                        breaks["breaker_" + breaker_index] += 1
                        if breaks["breaker_" + breaker_index] == 3:
                            node_count = 0
                            breaks["breaker_" + breaker_index] = 0
                            return node_count
                    # endregion "early quit if loop is taking too long"
                    node_count = self.at_least_one_node_calculator_2(new_df, prom,
                                                                     combinations_number=combinations_number,
                                                                     counter=counter, breaks=breaks)
                    try:
                        if node_count > 1:
                            return node_count
                    except TypeError:
                        node_count = 0
                        continue
            return node_count

    def best_vencode_generator(self):
        pass


class Enhancers(DatabaseOperations):
    """ A class describing the methods for the enhancers database """

    # unique to each call of Enhancers
    def __init__(self, file, names_db, celltype, celltype_exclude=None, not_include=None, partial_exclude=None,
                 sample_types="primary cells", nrows=None):
        super().__init__(file, celltype, celltype_exclude=celltype_exclude, not_include=not_include,
                         sample_types=sample_types, nrows=nrows)
        self.names_db = pd.read_csv(self.parent_path + names_db, sep="\t", index_col=1, header=None,
                                    names=["celltypes"])
        self.data, self.filtered_names_db = self.first_parser()
        self.codes = self.code_selector(self.celltype, custom=self.filtered_names_db,
                                        remove="fraction", not_include=self.not_include, to_dict=True)
        if partial_exclude:
            self.partial_exclude_codes = {}
            for key, value in partial_exclude.items():
                self.partial_exclude_codes[key] = self.code_selector(value[0], not_include=value[1])

    def first_parser(self):
        """remove universal RNA pools and select the desired sample_type"""
        samples_universal = self.code_selector("universal")
        data_dropped = self.raw_data.drop(samples_universal, axis=1, inplace=False)
        to_keep = self.sample_category_selector("sample types - FANTOM5.csv", self.sample_types, get="name")
        to_keep_codes = self.code_selector(to_keep, match="normal", remove="fraction")
        data = data_dropped[to_keep_codes]
        # Exclude some specific, on-demand, cell-types from the data straight away:
        if self.celltype_exclude is not None:
            codes_exclude = self.code_selector(self.celltype_exclude)
            data.drop(codes_exclude, axis=1, inplace=True)
            # Exclude and select new found codes from the names_db dataframe, to use in self.codes
            new_names_db = self.names_db.loc[to_keep_codes].drop(codes_exclude, axis=0)
        else:
            new_names_db = self.names_db.loc[to_keep_codes]
        return data, new_names_db

    def code_selector(self, celltype, custom=None, match="regex", remove=None, not_include=None, to_dict=False):
        """ Selects codes from a database containing enhancer style indexes. Input is normal celltype nomenclature """
        if custom is not None:
            lines_and_codes = custom
        else:
            lines_and_codes = self.names_db
        if isinstance(celltype, list):
            codes = []
            if to_dict:
                code_dict = {}
            for item in celltype:
                if match == "regex":
                    regular = ".*" + item.replace(" ", ".*").replace("+", "\+") + ".*"
                    idx = lines_and_codes.celltypes.str.contains(regular, flags=re.IGNORECASE, regex=True, na=False)
                elif match == "normal":
                    idx = lines_and_codes.celltypes.str.contains(item, regex=False, na=False)
                elif match == "exact":
                    idx = [True if item == cell_type else False for cell_type in lines_and_codes.celltypes]
                    pass
                else:
                    raise Exception("Match type is not recognized")
                codes_df = lines_and_codes[idx]
                codes.append(codes_df.index.values)
                if to_dict:
                    code_dict[item] = [code for sublist in codes for code in sublist]
                    codes = []
            if not to_dict:
                codes = [code for sublist in codes for code in sublist]
                codes = list(set(codes))
            if remove is not None:
                to_remove_bool = lines_and_codes.celltypes.str.contains(remove, regex=False, na=False)
                to_remove = lines_and_codes[to_remove_bool]
                to_remove = to_remove.index.values
                # remove all codes containing -remove-
                if to_dict:
                    codes = {}
                    for item, values in code_dict.items():
                        codes[item] = [code for code in values if code not in to_remove]
                else:
                    codes = [code for code in codes if code not in to_remove]
            if not_include is not None:
                if isinstance(not_include, dict):
                    for key, values in not_include.items():
                        if key not in codes.keys():
                            continue
                        values_codes = self.code_selector(values, custom=self.filtered_names_db)
                        not_include_codes = self.not_include_code_getter(values_codes, self.data)
                        codes[key] = list(set(codes[key]) - set(not_include_codes))
                else:
                    values_codes = self.code_selector(not_include, custom=self.filtered_names_db)
                    not_include_codes = self.not_include_code_getter(values_codes, self.data)
                    codes = list(set(codes) - set(not_include_codes))
        else:
            regular = ".*" + celltype.replace(" ", ".*").replace("+", "\+") + ".*"
            idx = lines_and_codes.celltypes.str.contains(regular, flags=re.IGNORECASE, regex=True, na=False)
            codes_df = lines_and_codes[idx]
            codes = codes_df.index.values.tolist()
            if not_include is not None:
                if isinstance(not_include, dict):
                    for key, values in not_include.items():
                        if key not in codes.keys():
                            continue
                        values_codes = self.code_selector(values, custom=self.filtered_names_db)
                        not_include_codes = self.not_include_code_getter(values_codes, self.data)
                        codes[key] = list(set(codes[key]) - set(not_include_codes))
                else:
                    values_codes = self.code_selector(not_include, custom=self.filtered_names_db)
                    not_include_codes = self.not_include_code_getter(values_codes, self.data)
                    codes = list(set(codes) - set(not_include_codes))
        self.test_codes(codes, celltype, codes_type=type(codes).__name__)
        return codes


class Promoters(DatabaseOperations):
    """ A class describing the methods for the promoters database """

    def __init__(self, file, celltype, celltype_exclude=None, not_include=None, partial_exclude=None,
                 sample_types="primary cells", second_parser=None, nrows=None):
        super().__init__(file, celltype, celltype_exclude=celltype_exclude, not_include=not_include,
                         sample_types=sample_types, skiprows=1831, second_parser=second_parser, nrows=nrows)
        self.data = self.first_parser()
        self.codes = self.code_selector(self.data, self.celltype, not_include=self.not_include,
                                        to_dict=True)
        if self.second_parser is not None:
            temp_codes = [x for a in self.codes.values() for x in a]
            code_df = self.data.filter(items=temp_codes)
            self.sample_types = self.second_parser
            self.data = self.first_parser()
            self.data = self.data.join(code_df)

        if partial_exclude:
            self.partial_exclude_codes = {}
            for key, value in partial_exclude.items():
                self.partial_exclude_codes[key] = self.code_selector(self.data, value[0], not_include=value[1])

    def first_parser(self):
        data_1 = self.raw_data.drop(self.raw_data.index[[0, 1]])
        universal_rna = self.code_selector(data_1, "universal", not_include=None)
        data_1.drop(universal_rna, axis=1, inplace=True)
        to_keep = self.sample_category_selector("sample types - FANTOM5.csv", self.sample_types)
        data = pd.DataFrame(index=data_1.index.values)
        for sample in to_keep:
            data_temp = data_1.filter(regex=sample)
            data = data.join(data_temp)
        # Exclude some specific, on-demand, cell-types from the data straight away:
        if self.celltype_exclude is not None:
            codes_exclude = self.code_selector(data, self.celltype_exclude)
            data.drop(codes_exclude, axis=1, inplace=True)
        return data

    def code_selector(self, db, celltype, not_include=None, to_dict=False):
        """ Selects codes from database using """
        if isinstance(celltype, list):
            codes = []
            if to_dict:
                code_dict = {}
            for item in celltype:
                codes_df = Defs.dataframe_regex_searcher(item, db)
                codes.append(codes_df.columns.values)
                if to_dict:
                    code_dict[item] = [code for sublist in codes for code in sublist]
                    codes = []
            if not to_dict:
                codes = [item for sublist in codes for item in sublist]
        else:
            codes_df = Defs.dataframe_regex_searcher(celltype, db)
            codes = list(codes_df.columns.values)

        if not_include is not None:
            if isinstance(not_include, dict):
                not_codes = []
                codes = code_dict
                for key, values in not_include.items():
                    if key not in code_dict.keys():
                        continue
                    codes_df = db[code_dict.get(key)]
                    not_codes = self.not_include_code_getter(values, codes_df)
                    # not_codes_df = Defs.dataframe_regex_searcher(values, codes_df)
                    # not_codes.append(not_codes_df.columns.values.tolist())
                    # not_codes = [item for sublist in not_codes for item in sublist]
                    codes[key] = list(set(code_dict[key]) - set(not_codes))
            else:
                not_codes_df = Defs.dataframe_regex_searcher(not_include, codes_df)
                not_codes = not_codes_df.columns.values
                if not codes:
                    codes = [item for sublist in code_dict.values() for item in sublist]
                codes = list(set(codes) - set(not_codes))
        self.test_codes(codes, celltype)
        return codes

# TODO: with the changes in __init__ to the BaseClass, some of these static methods may now be converted to self.xx!
