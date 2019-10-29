import os

from VEnCode.common_variables import promoter_file_name, enhancer_file_name
from VEnCode.utils import dir_and_file_handling as dfh
from scripts import num_validated_elements, get_validated_vencodes, vencodes_percent_active


class SetUp:
    """set up some variables: """
    cell_type = "B lymphoblastoid cell line: GM12878 ENCODE"
    type = "cell lines"
    data_type = "enhancers"
    algorithm = "sampling"
    k = 4
    number_vencodes_to_get = 200

    data_set = "both"  # default: None

    # Paths to export data:
    name_vencodes = "B lymph-both 200 val vencodes"

    name_cell_type = "B lymphoblastoid-both"
    name_heatmap = "heatmap 200 val vens - {}_inferno".format(name_cell_type)

    path_out_ven = os.path.join(
        "D:/Utilizador HDD/OneDrive - Nova Medical School Faculdade de Ciências Médicas da UNL/1-Research/"
        "3-Vencode/Fantom5/VEnCodes", name_vencodes)

    name_evalues = dfh.str_replace_multi(cell_type, {":": "-", "*": "-", "?": "-", "<": "-", ">": "-"})
    path_e_values = os.path.join(path_out_ven, name_evalues + "_evalues.csv")

    path_vencodes = path_out_ven

    path_out_heat = os.path.join(
        "D:/Utilizador HDD/OneDrive - Nova Medical School Faculdade de Ciências Médicas da UNL/1-Research/"
        "3-Vencode/Fantom5/Validations/2- Figures - py", name_heatmap)

    # Next ones you may not need to change:
    non_target_celltypes_inactivity = 0
    if data_type == "enhancers":
        file_name = enhancer_file_name
        target_celltype_activity = 0.1
    elif data_type == "promoters":
        file_name = promoter_file_name
        target_celltype_activity = 0.5
    else:
        raise AttributeError("data_type - {} - currently not supported".format(data_type))

    if algorithm == "heuristic":
        reg_element_sparseness = 0
    elif algorithm == "sampling":
        reg_element_sparseness = 90  # For some, you probably have to reduce sparseness to 0.
    else:
        raise AttributeError("Algorithm - {} - currently not supported".format(algorithm))


setup = SetUp()

elem = num_validated_elements.ValidatedElements(setup)
elem.export()

ven = get_validated_vencodes.ValidatedVEnCodes(setup)
ven.export()

heat = vencodes_percent_active.ValidatedVEnCodesHeatmap(setup, to_drop=len(ven.data.vencodes._celltype_donors)-1)
heat.plot()
