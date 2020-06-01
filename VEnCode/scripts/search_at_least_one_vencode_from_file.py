"""
search_at_least_one_vencode_from_file.py: Script that finds out if there is at least one VEnCode for every
cell type in a data set.
"""
import pandas as pd
import argparse

from tqdm import tqdm

from VEnCode import internals
import VEnCode.utils.dir_and_file_handling as dfh


def main(data_name, sep, out, algorithm, number_of_re, thresholds):
    data = pd.read_csv(data_name, sep=sep, index_col=0, engine="python")
    results = dict()
    for ctp in tqdm(data.columns, desc="Searching for VEnCodes for each cell type"):
        data_tpm = internals.DataTpm(data)
        data_tpm.load_data()
        data_tpm.make_data_celltype_specific(ctp, replicates=False)
        data_tpm.filter_by_target_celltype_activity(threshold=thresholds[0])
        data_tpm.define_non_target_celltypes_inactivity(threshold=thresholds[1])
        data_tpm.filter_by_reg_element_sparseness(threshold=thresholds[2])
        vencode = 0
        vencodes = internals.Vencodes(data_tpm, algorithm=algorithm, number_of_re=number_of_re)
        vencodes.next(1)
        if vencodes.vencodes:
            vencode = 1
        results[ctp] = vencode
    dfh.write_dict_to_csv(out, results, deprecated=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Searches for at least one VEnCode in the supplied data.
    And returns a list of celltypes indicating if they have at least one VEnCode or not.
    Example on how to use in command-line:
    >python search_at_least_one_vencode_from_file.py data.txt -o results.txt -a heuristic -k 4
    """)
    parser.add_argument("file", help="file with the data to search for VEnCodes")
    parser.add_argument("--sep", "-s", type=str, help="separator of the file columns", default=";")
    parser.add_argument("--out", "-o", type=str, help="path for the results, with name and extension of the file",
                        default="results.csv")
    parser.add_argument("--alg", "-a", type=str, help="Algorithm to use in the search", default="heuristic")
    parser.add_argument("--res", "-k", type=int, help="Number of REs (k) that comprise a VEnCode", default="4")
    parser.add_argument('--thr', '-t', nargs='+', type=int,
                        help="Thresholds to apply to the data. You should supply 3 values, first the target cell type "
                             "activity threshold, then the non-target cell types inactivity threshold. Finally the "
                             "sparseness threshold, that filters the data retaining only the sparsest REs.",
                        default=[1, 0, 0])

    args = parser.parse_args()
    main(args.file, sep=args.sep, out=args.out, algorithm=args.alg, number_of_re=args.res, thresholds=args.thr)
