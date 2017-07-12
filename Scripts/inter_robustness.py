#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""inter_robustness.py: Functions for generating inter robustness data """

import Classes_2

# region "Global Variables"
complete_primary_non_include_list = {"Adipocyte - breast": "pre", "Adipocyte - omental": "pre",
                                     "Adipocyte - perirenal": "pre", "Adipocyte - subcutaneous": "pre",
                                     "Endothelial Cells - Vein": "Umbilical",
                                     "Renal Epithelial Cells": "Cortical",
                                     "Skeletal Muscle Cells": ["satellite", "differentiated"]}
complete_primary_exclude_list = ["mesenchymal precursor cell - ovarian", "Osteoblast - differentiated",
                                 "Peripheral Blood Mononuclear Cells", "Whole blood"]
complete_primary_jit_exclude_list = {"CD14+ CD16- Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD14+ CD16+ Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD14- CD16+ Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD4+CD25+CD45RA- memory regulatory T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25+CD45RA+ naive regulatory T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25-CD45RA- memory conventional T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25-CD45RA+ naive conventional T cells": ("CD4+ T Cells", "CD25")}

three_donors_cell_list = ["Adipocyte - omental", "Adipocyte - subcutaneous",
                          "Amniotic Epithelial Cells", "amniotic membrane cells", "Basophils", "Cardiac Myocyte",
                          "CD14+CD16+ Monocytes",
                          "chorionic membrane cells", "Ciliary Epithelial Cells",
                          "Endothelial Cells - Artery", "Endothelial Cells - Lymphatic",
                          "Endothelial Cells - Microvascular", "Endothelial Cells - Umbilical vein",
                          "Endothelial Cells - Vein", "Eosinophils", "Esophageal Epithelial Cells",
                          "Fibroblast - Aortic Adventitial",
                          "Fibroblast - Lung", "Fibroblast - Lymphatic",
                          "Fibroblast - Villous Mesenchymal", "Gingival epithelial cells",
                          "granulocyte macrophage progenitor",
                          "Hepatic Sinusoidal Endothelial Cells", "Hepatic Stellate Cells (lipocyte)", "Hepatocyte",
                          "Keratinocyte - epidermal", "Keratocytes", "Lens Epithelial Cells",
                          "Macrophage - monocyte derived", "Mammary Epithelial Cell",
                          "Meningeal Cells", "mesenchymal precursor cell - adipose",
                          "mesenchymal precursor cell - bone marrow", "Mesenchymal stem cells - hepatic",
                          "Mesothelial Cells", "Natural Killer Cells", "Neurons", "Nucleus Pulposus Cell", "Osteoblast",
                          "Pericytes", "Placental Epithelial Cells", "Preadipocyte - omental",
                          "Preadipocyte - subcutaneous", "Preadipocyte - visceral", "Prostate Epithelial Cells",
                          "Prostate Stromal Cells", "Renal Epithelial Cells", "Renal Proximal Tubular Epithelial Cell",
                          "Schwann Cells", "Sebocyte",
                          "Skeletal muscle cells differentiated into Myotubes - multinucleated",
                          "Skeletal Muscle Satellite Cells", "Smooth Muscle Cells - Carotid",
                          "Smooth Muscle Cells - Colonic", "Smooth Muscle Cells - Coronary Artery",
                          "Smooth Muscle Cells - Prostate",
                          "Smooth Muscle Cells - Tracheal",
                          "Smooth Muscle Cells - Umbilical Vein", "Synoviocyte",
                          "tenocyte", "Trabecular Meshwork Cells", "Tracheal Epithelial Cells"]

four_donors_cell_list = ["Smooth Muscle Cells - Umbilical Artery", "Retinal Pigment Epithelial Cells",
                         "Smooth Muscle Cells - Aortic", "Mesenchymal Stem Cells - umbilical",
                         "Endothelial Cells - Aortic", "Mesenchymal Stem Cells - adipose", "Urothelial Cells",
                         "mature adipocyte", "Mesenchymal Stem Cells - bone marrow", "Olfactory epithelial cells"]

partial_list = ["Adipocyte - omental", "Adipocyte - subcutaneous", "Alveolar Epithelial Cells",
                "Amniotic Epithelial Cells", "amniotic membrane cells", "Astrocyte - cerebellum",
                "Astrocyte - cerebral cortex", "Basophils", "Cardiac Myocyte", "CD14+CD16- Monocytes",
                "CD14+CD16+ Monocytes", "CD14-CD16+ Monocytes",
                "chorionic membrane cells", "Ciliary Epithelial Cells", "Corneal Epithelial Cells",
                "Dendritic Cells - plasmacytoid",
                "Endothelial Cells - Artery", "Endothelial Cells - Lymphatic",
                "Endothelial Cells - Microvascular", "Endothelial Cells - Umbilical vein",
                "Endothelial Cells - Vein", "Eosinophils", "Esophageal Epithelial Cells",
                "Fibroblast - Aortic Adventitial", "Fibroblast - Choroid Plexus",
                "Fibroblast - Lung", "Fibroblast - Lymphatic", "Fibroblast - Mammary",
                "Fibroblast - Villous Mesenchymal", "Gingival epithelial cells",
                "granulocyte macrophage progenitor", "Hair Follicle Dermal Papilla Cells",
                "Hepatic Sinusoidal Endothelial Cells", "Hepatic Stellate Cells (lipocyte)", "Hepatocyte",
                "Keratinocyte - epidermal", "Keratocytes", "Lens Epithelial Cells",
                "Macrophage - monocyte derived", "Mallassez-derived cells", "Mammary Epithelial Cell",
                "Meningeal Cells", "mesenchymal precursor cell - adipose",
                "mesenchymal precursor cell - bone marrow", "Mesenchymal stem cells - hepatic",
                "Mesothelial Cells", "migratory langerhans cells",
                "Myoblast", "Natural Killer Cells", "Neurons", "Nucleus Pulposus Cell", "Osteoblast",
                "Pericytes", "Placental Epithelial Cells", "Preadipocyte - omental",
                "Preadipocyte - subcutaneous", "Preadipocyte - visceral", "Prostate Epithelial Cells",
                "Prostate Stromal Cells", "Renal Epithelial Cells", "Renal Mesangial Cells",
                "Renal Proximal Tubular Epithelial Cell",
                "salivary acinar cells", "Schwann Cells", "Sebocyte",
                "Skeletal muscle cells differentiated into Myotubes - multinucleated",
                "Skeletal Muscle Satellite Cells", "Small Airway Epithelial Cells",
                "Smooth Muscle Cells - Brachiocephalic",
                "Smooth Muscle Cells - Brain Vascular", "Smooth Muscle Cells - Carotid",
                "Smooth Muscle Cells - Colonic", "Smooth Muscle Cells - Coronary Artery",
                "Smooth Muscle Cells - Internal Thoracic Artery", "Smooth Muscle Cells - Prostate",
                "Smooth Muscle Cells - Pulmonary Artery", "Smooth Muscle Cells - Subclavian Artery",
                "Smooth Muscle Cells - Tracheal",
                "Smooth Muscle Cells - Umbilical Vein", "Synoviocyte",
                "tenocyte", "Trabecular Meshwork Cells", "Tracheal Epithelial Cells"]
# endregion "Global variables"

if __name__ == "__main__":
    initialize_promoters = Classes_2.Promoters("hg19.cage_peak_phase1and2combined_tpm.osc.txt",
                                               three_donors_cell_list,
                                               celltype_exclude=complete_primary_exclude_list,
                                               not_include=complete_primary_non_include_list,
                                               partial_exclude=complete_primary_jit_exclude_list,
                                               sample_types="primary cells",
                                               second_parser=None)
    # get the percentage of VEnCodes taken for 1 donor that work for all donors:
    initialize_promoters.ven_diagram_interception(2000, 5, 3, combinations_number=4, threshold=90)

    # get the percentage of VEnCodes taken for 1, 2, 3, etc donors that work for all and com

# TODO: run inter for cancer cell lines: 3 donors, 4 donors, etc..
