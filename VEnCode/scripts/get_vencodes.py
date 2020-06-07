#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" get_vencodes.py: Script used to generate VEnCodes from user supplied data. """

import sys
import argparse
from pathlib import Path

file_dir = str(Path(__file__).resolve().parents[2])
sys.path.append(file_dir)

from VEnCode import internals_extensions as iext


def main(inputs, cell_type, number_vencodes=1, algorithm="heuristic", number_of_re=4, thresholds=(0.5, 0, 0),
         n_samples=10000, stop=5, files_path=None, out="VEnCodesFantom"):
    def export(data_, path):
        """
        Export the E values and VEnCode data to a file.
        """
        data_.export("vencodes", "e-values", path=path)

    if out is not None:
        pass

    if thresholds is None:
        non_target_celltypes_inactivity = 0
        target_celltype_activity = 0.5
        if algorithm == "heuristic":
            reg_element_sparseness = 0
        elif algorithm == "sampling":
            reg_element_sparseness = 90  # For some, you probably have to reduce sparseness to 0.
        else:
            raise AttributeError("Algorithm - {} - currently not supported".format(algorithm))
        thresholds = (target_celltype_activity, non_target_celltypes_inactivity, reg_element_sparseness)

    data_obj = iext.GetVencodes(inputs=inputs,
                                cell_type=cell_type,
                                algorithm=algorithm,
                                n_regulatory_elements=number_of_re,
                                number_vencodes=number_vencodes,
                                thresholds=thresholds, n_samples=n_samples, stop=stop, files_path=files_path)

    export(data_obj, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Searches for VEnCodes in the supplied data.
        Exports the VEnCodes and their E values to files, which path you can set with -out.
        Example on how to use in command-line:
        >python get_vencodes.py Hepatocyte --algorithm heuristic --number_vencodes 2
        """)
    parser.add_argument("file", help="File with the data to search for VEnCodes.")
    parser.add_argument("cell_type", help="sCelltype to search for VEnCodes.")

    parser.add_argument("--number_vencodes", "-n", type=int, help="Number of VEnCodes to retrieve.", default=1)
    parser.add_argument("--algorithm", "-a", type=str, help="Algorithm to use in the search.", default="heuristic")
    parser.add_argument("--number_of_re", "-k", type=int, help="Number of REs (k) that comprise a VEnCode.", default=4)
    parser.add_argument('--thresholds', '-t', nargs='+', type=int,
                        help="Thresholds to apply to the data. You should supply 3 values, first the target cell type "
                             "activity threshold, then the non-target cell types inactivity threshold. Finally the "
                             "sparseness threshold, that filters the data retaining only the sparsest REs.",
                        default=None)
    parser.add_argument("--out", "-o", type=str, help="Path or partial path to a folder to store the result files."
                                                      "Creates a new folder if non-existent.", default="VEnCodesFantom")
    parser.add_argument("--n_samples", "-c", type=int, help="Optional for sampling method. Number of samples to take "
                                                            "before giving up finding a VEnCode.", default=10000)
    parser.add_argument("--stop", "-w", type=int, help="Optional for heuristic method. Depth of nodes to go before "
                                                       "giving up finding a VEnCode.", default=5)
    parser.add_argument("--files_path", "-f", type=str, help="Folder where the FANTOM5 files are located."
                                                             "Used if full path is not supplied in argument 'file'.",
                        default=None)

    args = parser.parse_args()
    main(args.file, args.cell_type, number_vencodes=args.number_vencodes, algorithm=args.algorithm,
         number_of_re=args.number_of_re, thresholds=args.thresholds, n_samples=args.n_samples, stop=args.stop,
         files_path=args.files_path, out=args.out)
