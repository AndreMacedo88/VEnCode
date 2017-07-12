#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Classes.py: Classes module for the VEnCode project """

# import csv
import os
import re
import logging

import numpy as np
import pandas as pd
from scipy.special import comb
import itertools as iter

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
        self.raw_data = pd.read_csv(self.parent_path + self.file, sep="\t", index_col=0, skiprows=skiprows, nrows=nrows,
                                    engine="python")

    @staticmethod
    def sample_category_selector(sample_types_file, types, path="parent", get="index"):
        """
        Returns a list of cell types to keep/drop from a file containing the list of cell types and a 'Sample category'
        column which determines which cell types to retrieve.
        """
        if not isinstance(types, list): types = [types]
        if path == "parent":
            parent_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
            database = pd.read_csv(parent_path + "/Files/" + sample_types_file, sep=",", index_col=0, engine="python")
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
        for celltype in self.celltype:  # start by cycling through the celltypes to find a VEnCode
            print("", "Starting {}".format(celltype), sep="-> ")
            logging.info("Starting %s", celltype)
            vencodes[celltype] = 0  # init the dict to append a possible VEnCode
            codes = class_instance.codes.get(celltype)  # get codes for this celltype
            try:  # this next section can remove some celltypes from the data just in time before the analysis
                if class_instance.partial_exclude_codes is not None:
                    if celltype in class_instance.partial_exclude_codes.keys():
                        partial_exclude_codes = class_instance.partial_exclude_codes.get(celltype)
                        data_frame = class_instance.data.drop(partial_exclude_codes, axis=1)
                    else:
                        data_frame = class_instance.data
            except AttributeError:
                data_frame = class_instance.data
            filter_2 = Defs.fantom_filters(data_frame, codes, expression, 2, threshold)  # apply efficiency filters
            filter_2 = filter_2.applymap(lambda x: 0 if x == 0 else 1)  # change expression to 1s and 0s for quickness
            filter_2["sum"] = filter_2.sum(axis=1)  # create a extra column with the sum of 1s for each row (promoter)
            filter_2.sort_values(["sum"], inplace=True)  # sort promoters based on the previous sum. Descending order
            to_get_nodes = filter_2.drop(codes + ["sum"], axis=1)  # remove the celltype to analyse and the sum column
            promoters = to_get_nodes.index.values  # get a list (not really a list) of all the promoters, to cycle
            breaks = {}  # this next section creates a dictionary to update with how many times each node is cycled
            for item in range(1, combinations_number):
                breaks["breaker_" + str(item)] = 0
            for promoter in promoters:  # cycle the promoters, starting with these of lower sum of 1s. Founder node.
                """
                Finally, start cycling through nodes - promoters - used one at a time to combine with preceding and 
                subsequent nodes with the purpose of finding a VEnCode - a combination of promoters where 
                for every celltype there's at least one promoter of that combination that does not exhibit any activity
                """
                vencodes_from_nodes = self.at_least_one_node_calculator(to_get_nodes, promoter,
                                                                        combinations_number=combinations_number,
                                                                        breaks=breaks)
                if vencodes_from_nodes > 0:
                    vencodes[celltype] = 1  # if at least one VEnCode was found with the node strategy, return 1
                    break  # and stop the cycling of the founder nodes
                else:
                    to_get_nodes.drop(promoter, axis=0, inplace=True)  # else, drop this founder node and go to next
            print(vencodes_from_nodes)
        return vencodes

    def at_least_one_node_calculator(self, data_frame, promoter, combinations_number=4, counter=1, breaks=None):
        """

        :param data_frame: Data frame containing cage-seq expression profile for several celltypes. Dataframe object
        :param promoter: Previous promoter name(founder node if first time calling this function). str type
        :param combinations_number: Number of combinations of promoters to find VEnCodes. int type
        :param counter: Counter is equal to the depth of the current node. int type
        :param breaks: Dictionary containing keys for the different levels of breaks (one per each combination number)
        and values corresponding to how many times each combination already cycled. dict type
        :return: If the algorithm found a definite VEnCode or not.
        """
        cols = data_frame.loc[promoter] != 0  # create a mask where True marks the celltypes in which the previous
        # node is still expressed
        cols = data_frame.columns[cols]  # apply that mask, selecting the columns that are True
        new_df = data_frame[cols].drop(promoter, axis=0)  # apply the selection and take the prom out of the dataframe
        nodes = (new_df == 0).all(axis=1)  # Check if any VEnCode - if any other promoter have 0 expression in all cells
        node_count = np.sum(nodes)  # if any True (VEnCode) the value of that True becomes 1 and sum gives num VEnCodes
        if node_count > 0:
            return node_count  # found at least one VEnCode so it can return a successful answer
        else:  # if in previous node could not get a definite VEnCode, re-start search with next node
            new_df["sum"] = new_df.sum(axis=1)  # create a extra column with the sum for each row (promoter)
            new_df.sort_values(["sum"], inplace=True)  # sort promoters based on the previous sum. Descending order
            new_df = new_df.drop(["sum"], axis=1)  # now remove the sum column
            counter = counter  # counter is defined with previous counter for recursive use of this function
            counter_thresholds = [i for i in range(2, (combinations_number + 1))]  # set maximum number for counter
            # loop the next area until number of nodes in combination exceeds num of desired proms in comb for VEnCode
            while counter < combinations_number:
                counter += 1  # updates the counter as it will enter the next node depth
                for prom in new_df.index.values:  # next promoter
                    # region "early quit if loop is taking too long"
                    if counter in counter_thresholds:
                        breaker_index = str(counter_thresholds.index(counter) + 1)
                        breaks["breaker_" + breaker_index] += 1
                        if breaks["breaker_" + breaker_index] == 3:  # here, we only test 3 promoters per node level
                            node_count = 0
                            breaks["breaker_" + breaker_index] = 0
                            return node_count
                    # endregion "early quit if loop is taking too long"
                    node_count = self.at_least_one_node_calculator(new_df, prom,
                                                                   combinations_number=combinations_number,
                                                                   counter=counter, breaks=breaks)
                    try:
                        if node_count > 1:
                            return node_count
                    except TypeError:
                        node_count = 0
                        continue
            return node_count


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
            code_dict = {celltype: codes}
        if not_include is not None:
            if isinstance(not_include, dict):
                codes = code_dict
                for key, values in not_include.items():
                    if key not in code_dict.keys():
                        continue
                    codes_df = db[code_dict.get(key)]
                    not_codes = self.not_include_code_getter(values, codes_df)
                    codes[key] = list(set(code_dict[key]) - set(not_codes))
            else:
                not_codes_df = Defs.dataframe_regex_searcher(not_include, codes_df)
                not_codes = not_codes_df.columns.values
                if not codes:
                    codes = [item for sublist in code_dict.values() for item in sublist]
                codes = list(set(codes) - set(not_codes))
        self.test_codes(codes, celltype)
        return codes

    # For figures:

    def ven_diagrams(self, vens_to_take, combinations_number=4, expression=1, threshold=90):
        if isinstance(self.codes, dict):
            codes = [j for i in list(self.codes.values()) for j in i]
        else:
            codes = self.codes
        for code in codes:  # subsequent analysis will be done for each sample (code) separately
            print("donor:", code, sep="\n", end="\n\n")
            codes_2 = [x for x in codes if x != code]
            donors_data = self.data[codes_2]
            data_2 = self.data.drop(codes_2, axis=1)
            print("donors to exclude:", *codes_2, sep="\n", end="\n\n")
            filter_2 = Defs.fantom_filters(data_2, code, expression, 2, threshold)
            ven_diagram = {"None": []}  # create the final dict to later append VEnCode counts:
            for r in reversed(range(1, (len(codes_2) + 1))):  # create dict keys up front because of append
                for z in iter.combinations(codes_2, r):
                    z = list(z)
                    string_z = "  ".join(z)
                    ven_diagram[string_z] = []
            n = 0
            while n < vens_to_take:
                sample = filter_2.sample(n=combinations_number)  # taking a sample of n RE to assess if VEnCode
                sample_dropped = sample.drop(code, axis=1).values  # dropping code celltype from the sample
                assess_if_vencode = np.any(sample_dropped == 0, axis=0)  # assess if at least 1 zero per column
                if all(assess_if_vencode):  # is RE sample a VEnCode for code?
                    n += 1
                    donors_data_sample = donors_data.loc[sample.index.values]  # locating the n RE in the other codes db
                    no_ven = True  # will change to false to prevent double attribution of VEnCodes to the other codes
                    counter = 0
                    none_counter = 0
                    # start iterator through all possible other codes combinations:
                    for i in reversed(range(1, (len(codes_2) + 1))):
                        for y in iter.combinations(codes_2, i):  # take a combination of i number of other codes
                            y = list(y)
                            string_y = "  ".join(y)
                            to_assess = donors_data_sample[y]
                            assess_if_not_vencode_donors = np.any(to_assess.values == 0, axis=0)
                            # is RE sample not a VEnCode for this comb other codes?
                            try:
                                if assess_if_not_vencode_donors:
                                    pass
                                else:
                                    counter += 1
                                    no_ven = False
                                    break
                            except:
                                if any(assess_if_not_vencode_donors):
                                    pass
                                else:
                                    counter += 1
                                    no_ven = False
                                    break
                        ven_diagram[string_y].append(counter)
                        if not no_ven:
                            break
                    if counter == 0:
                        none_counter += 1
                        ven_diagram["None"].append(none_counter)
                else:
                    pass
            for key in ven_diagram:
                ven_diagram[key] = sum(ven_diagram[key])

            folder = "/Figure 3-b2/cell lines/"
            if not os.path.exists(folder):
                os.makedirs(folder)
            position = codes.index(code) + 1
            file_name = u"/Donor{} - {}".format(position, self.celltype)
            Defs.write_one_value_dict_to_csv(file_name + ".csv", ven_diagram, folder)

    def statistics_ven_diagram(self, vens_to_test, sampling_number, number_donors, combinations_number=4,
                               expression=1, threshold=90):
        """ To use to generate statistics based on random sampling the data for functions ven_diagram and, especially,
        ven_diagram_interception """
        codes_1 = Defs.possible_dict_to_list(self.codes)
        number = 0
        ven_diagram_1 = {}
        while number < sampling_number:
            number += 1
            codes = np.random.choice(codes_1, number_donors, replace=False)
            ven_diagram = {}
            for code in codes:
                print("Sample:", code, sep="\n", end="\n\n")
                codes_2 = [x for x in codes if x != code]
                donors_data = self.data[codes_2]
                data_2 = self.data.drop(codes_2, axis=1)
                print("Samples to exclude:", *codes_2, sep="\n", end="\n\n")
                filter_2 = Defs.fantom_filters(data_2, code, expression, 2, threshold)
                ven_diagram[code] = []
                n = 0
                stop_counter = 0
                while n < vens_to_test:
                    stop_counter += 1
                    sample = filter_2.sample(n=combinations_number)
                    sample_dropped = sample.drop(code, axis=1).values
                    assess_if_vencode = np.any(sample_dropped == 0, axis=0)
                    if all(assess_if_vencode):
                        n += 1
                        donors_data_sample = donors_data.loc[sample.index.values]
                        counter = 0
                        to_assess = donors_data_sample[codes_2]
                        assess_if_not_vencode_donors = np.any(to_assess.values == 0, axis=0)
                        try:
                            if assess_if_not_vencode_donors:
                                pass
                            else:
                                counter += 1
                        except:
                            if any(assess_if_not_vencode_donors):
                                pass
                            else:
                                counter += 1
                        ven_diagram[code].append(counter)
                    if stop_counter > 500000:  # in case there is not enough VEnCodes to check them all
                        break
                ven_diagram[code] = sum(ven_diagram[code]) / stop_counter
                print(n)
            codes_string = ''.join(codes.tolist())
            ven_diagram_1[codes_string] = np.mean(list(ven_diagram.values()))
        print(ven_diagram_1)
        folder = "/Figure 3-b2/cell lines/"
        if not os.path.exists(folder):
            os.makedirs(folder)
        if isinstance(self.celltype, str):
            file_name = u"/stats for {} celltype".format(self.celltype)
        else:
            file_name = u"/stats for {} celltypes".format(len(self.codes))
        Defs.write_one_value_dict_to_csv(file_name + ".csv", ven_diagram_1, folder)

    def ven_diagram_interception(self, vens_to_test, sampling_number, number_donors, combinations_number=4,
                                 expression=1, threshold=90):
        for code_1 in self.codes:
            codes_1 = self.codes[code_1]
            number = 0
            ven_diagram_1 = {}
            while number < sampling_number:
                number += 1
                codes = np.random.choice(codes_1, number_donors, replace=False)
                ven_diagram = {}
                for code in codes:
                    print("Sample:", code, sep="\n", end="\n\n")
                    codes_2 = [x for x in codes if x != code]
                    donors_data = self.data[codes_2]
                    data_2 = self.data.drop(codes_2, axis=1)
                    print("Samples to exclude:", *codes_2, sep="\n", end="\n\n")
                    filter_2 = Defs.fantom_filters(data_2, code, expression, 2, threshold)
                    ven_diagram[code] = []
                    n = 0
                    stop_counter = 0
                    while n < vens_to_test:
                        stop_counter += 1
                        sample = filter_2.sample(n=combinations_number)
                        sample_dropped = sample.drop(code, axis=1).values
                        assess_if_vencode = np.any(sample_dropped == 0, axis=0)
                        if all(assess_if_vencode):
                            n += 1
                            donors_data_sample = donors_data.loc[sample.index.values]
                            counter = 0
                            to_assess = donors_data_sample[codes_2]
                            assess_if_not_vencode_donors = np.any(to_assess.values == 0, axis=0)
                            try:
                                if assess_if_not_vencode_donors:
                                    pass
                                else:
                                    counter += 1
                            except:
                                if any(assess_if_not_vencode_donors):
                                    pass
                                else:
                                    counter += 1
                            ven_diagram[code].append(counter)
                        if stop_counter > 500000:  # in case there is not enough VEnCodes to check them all
                            break
                    ven_diagram[code] = sum(ven_diagram[code]) / stop_counter
                    print(n)
                codes_string = ''.join(codes.tolist())
                ven_diagram_1[codes_string] = np.mean(list(ven_diagram.values()))
            print(ven_diagram_1)
            folder = "/Figure 3-b2/cell lines/"
            if not os.path.exists(folder):
                os.makedirs(folder)
            file_name = u"/interception of {}".format(code_1)
            Defs.write_one_value_dict_to_csv(file_name + ".csv", ven_diagram_1, folder)

    def inter_donor_percentage_difference(self, vens_to_test, sampling_number, number_donors, combinations_number=4,
                                          expression=1, threshold=90):
        for code_1 in self.codes:
            codes_1 = self.codes[code_1]
            number = 0
            ven_diagram_1 = {}
            while number < sampling_number:
                number += 1
                codes = np.random.choice(codes_1, number_donors, replace=False)
                ven_diagram = {}
                for code in codes:
                    print("Sample:", code, sep="\n", end="\n\n")
                    codes_2 = [x for x in codes if x != code]
                    donors_data = self.data[codes_2]
                    data_2 = self.data.drop(codes_2, axis=1)
                    print("Samples to exclude:", *codes_2, sep="\n", end="\n\n")
                    filter_2 = Defs.fantom_filters(data_2, code, expression, 2, threshold)
                    ven_diagram[code] = []
                    n = 0
                    stop_counter = 0
                    while n < vens_to_test:
                        stop_counter += 1
                        sample = filter_2.sample(n=combinations_number)
                        sample_dropped = sample.drop(code, axis=1).values
                        assess_if_vencode = np.any(sample_dropped == 0, axis=0)
                        if all(assess_if_vencode):
                            n += 1
                            donors_data_sample = donors_data.loc[sample.index.values]
                            counter = 0
                            to_assess = donors_data_sample[codes_2]
                            assess_if_not_vencode_donors = np.any(to_assess.values == 0, axis=0)
                            try:
                                if assess_if_not_vencode_donors:
                                    pass
                                else:
                                    counter += 1
                            except:
                                if any(assess_if_not_vencode_donors):
                                    pass
                                else:
                                    counter += 1
                            ven_diagram[code].append(counter)
                        if stop_counter > 500000:  # in case there is not enough VEnCodes to check them all
                            break
                    ven_diagram[code] = sum(ven_diagram[code]) / stop_counter
                    print(n)
                codes_string = ''.join(codes.tolist())
                ven_diagram_1[codes_string] = np.mean(list(ven_diagram.values()))
            print(ven_diagram_1)
            folder = "/Figure 3-b2/cell lines/"
            if not os.path.exists(folder):
                os.makedirs(folder)
            file_name = u"/interception of {}".format(code_1)
            Defs.write_one_value_dict_to_csv(file_name + ".csv", ven_diagram_1, folder)

    def get_vencodes(self, combinations_number=4, p=None, n=None, write_file=False):
        if isinstance(self.codes, dict):
            codes = [j for i in list(self.codes.values()) for j in i]
        else:
            codes = self.codes
        if len(codes) == 1:
            code = codes[0]
            others_df = None
            filter_1 = Defs.fantom_filter_1(self.data, codes, 1)
            threshold = 100
            while threshold > 0:
                with_percentage_of_zeros, column_name = Defs.fantom_percentile_calculator(filter_1, codes,
                                                                                          threshold)
                filter_2 = Defs.fantom_filter_2(with_percentage_of_zeros, column_name)
                if p is not None:
                    sorted_1 = filter_2.nlargest(p, codes)
                else:
                    sorted_1 = filter_2.sort(codes.tolist(), ascending=False)
                print("starting %s -> threshold = %s" % (code, threshold))
                success = Defs.reform_vencode_n_combinations_of_k(threshold, sorted_1.drop(column_name, axis=1),
                                                                  code, self.celltype, "Promoters", combinations_number,
                                                                  n, others_df, write_file)
                if not success:
                    threshold -= 10
                if success:
                    threshold = 0
        else:
            code = codes[0]
            codes_2 = np.array(codes).tolist()
            codes_2.remove(code)
            filter_1 = Defs.fantom_filter_1(self.data, codes, 1)
            threshold = 100
            while threshold > 0:
                with_percentage_of_zeros, column_name = Defs.fantom_percentile_calculator(filter_1, codes,
                                                                                          threshold)
                filter_2 = Defs.fantom_filter_2(with_percentage_of_zeros, column_name)
                if p is not None:
                    sorted_1 = filter_2.nlargest(p, codes)
                else:
                    sorted_1 = filter_2.sort(codes, ascending=False)
                codes_2_df = sorted_1[codes_2]
                cropped = sorted_1.drop(codes_2, axis=1)
                print("Starting %s -> threshold = %s" % (self.celltype, threshold))
                success = Defs.reform_vencode_n_combinations_of_k(threshold, cropped.drop(column_name, axis=1),
                                                                  code, self.celltype, "Promoters", combinations_number,
                                                                  n, codes_2_df, write_file)
                if not success:
                    threshold -= 10
                if success:
                    threshold = 0

    def best_vencode_generator(self):
        print(self.codes)
        pass

    def intra_individual_robustness(self, combinations_number, vens_to_take, reps=1, threshold=90, expression=1,
                                    get_vencodes=False):
        final = {}
        final_vencodes = {}
        base_threshold = threshold
        # Starting loop through all cell types:
        for cell in self.codes:
            threshold = base_threshold
            codes = self.codes[cell]
            print("Cell types to get VEnCodes:", *codes, sep="\n", end="\n\n")
            filter_1 = Defs.fantom_filters(self.data, codes, expression, 1)
            false_negatives = []
            counter = 0
            while len(false_negatives) < vens_to_take:
                filter_2, threshold = Defs.fantom_filter_threshold(filter_1, codes, threshold, combinations_number)
                if get_vencodes:
                    false_negatives, vencodes = Defs.fantom_sampling_monte_carlo(codes, filter_2, combinations_number,
                                                                                 vens_to_take, reps,
                                                                                 vencodes=threshold, stop_at=250000)
                    final_vencodes.update(vencodes)
                else:
                    false_negatives = Defs.fantom_sampling_monte_carlo(codes, filter_2, combinations_number,
                                                                       vens_to_take,
                                                                       reps, vencodes=get_vencodes, stop_at=250000)
                threshold -= 5
                counter += 1
                if threshold < 50 or counter == 3:
                    if len(false_negatives) == vens_to_take:
                        break
                    for i in range(len(false_negatives), vens_to_take):
                        false_negatives.append("-")
                    break
            final[cell] = false_negatives

        folder = "/Figure 2/"
        if get_vencodes:
            file_name = u"/VEnCodes {} samples {} VEnCodes.csv".format(len(self.codes), vens_to_take)
            Defs.write_dict_to_csv(file_name, final_vencodes, folder, path="parent")
        file_name = u"/VEnCode E-values {} samples {} VEnCodes.csv".format(len(self.codes), vens_to_take)
        Defs.write_dict_to_csv(file_name, final, folder, path="parent")
        return

    # Tests:

    def test_code_size(self):
        for cell in self.codes:
            codes = self.codes[cell]
            print(cell, len(codes), sep=": ")

    def test_code_names(self, size=None):
        for cell in self.codes:
            codes = self.codes[cell]
            if size is not None:
                if len(codes) == size:
                    print(cell, codes, sep=": ")
                else:
                    pass
            else:
                print(cell, codes, sep=": ")

    def codes_to_csv(self, file_name, type, folder_name):
        if type == "dict":
            Defs.write_dict_to_csv(file_name, self.codes, folder_name, path="parent")
        if type == "list":
            codes_list = Defs.possible_dict_to_list(self.codes)
            Defs.write_list_to_csv(file_name, codes_list, folder_name, path="parent")

    def celltypes_to_csv(self, file_name, type, folder_name):
        cell_list = list(self.codes.keys())
        Defs.write_list_to_csv(file_name, cell_list, folder_name, path="parent")

# TODO: with the changes in __init__ to the BaseClass, some of these static methods may now be converted to self.xx!
