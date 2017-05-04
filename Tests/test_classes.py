#!/usr/bin/env python

"""Tests classes for VEnCode"""

from .context import Classes

# region "Global Variables"
complete_primary_cell_list = ["Adipocyte - breast", "Adipocyte - omental", "Adipocyte - perirenal",
                              "Adipocyte - subcutaneous", "Alveolar Epithelial Cells", "Amniotic Epithelial Cells",
                              "amniotic membrane cells", "Anulus Pulposus Cell", "Astrocyte - cerebellum",
                              "Astrocyte - cerebral cortex", "Basophils", "Bronchial Epithelial Cell",
                              "Cardiac Myocyte", "CD133+ stem cells - adult bone marrow derived",
                              "CD133+ stem cells - cord blood derived",
                              "CD14+ monocyte derived endothelial progenitor cells", "CD14+ Monocytes",
                              "CD14+ CD16- Monocytes", "CD14+ CD16+ Monocytes", "CD14- CD16+ Monocytes",
                              "CD19+ B Cells", "CD34 cells differentiated to erythrocyte lineage", "CD34+ Progenitors",
                              "CD34+ stem cells - adult bone marrow derived", "CD4+ T Cells",
                              "CD4+CD25+CD45RA- memory regulatory T cells", "CD4+CD25+CD45RA+ naive regulatory T cells",
                              "CD4+CD25-CD45RA- memory conventional T cells",
                              "CD4+CD25-CD45RA+ naive conventional T cells", "CD8+ T Cells", "Chondrocyte",
                              "chorionic membrane cells", "Ciliary Epithelial Cells", "common myeloid progenitor CMP",
                              "Corneal Epithelial Cells", "Dendritic Cells - monocyte immature derived",
                              "Dendritic Cells - plasmacytoid", "Endothelial Cells - Aortic",
                              "Endothelial Cells - Artery", "Endothelial Cells - Lymphatic",
                              "Endothelial Cells - Microvascular", "Endothelial Cells - Thoracic",
                              "Endothelial Cells - Umbilical vein", "Endothelial Cells - Vein", "Eosinophils",
                              "Esophageal Epithelial Cells", "Fibroblast - Aortic Adventitial", "Fibroblast - Cardiac",
                              "Fibroblast - Choroid Plexus", "Fibroblast - Conjunctival", "Fibroblast - Dermal",
                              "Fibroblast - Gingival", "Fibroblast - Lung", "Fibroblast - Lymphatic",
                              "Fibroblast - Mammary", "Fibroblast - Periodontal Ligament",
                              "Fibroblast - Pulmonary Artery", "Fibroblast - skin", "Fibroblast - Villous Mesenchymal",
                              "gamma delta positive T cells", "Gingival epithelial cells",
                              "granulocyte macrophage progenitor", "Hair Follicle Dermal Papilla Cells",
                              "Hair Follicle Outer Root Sheath Cells", "Hepatic Sinusoidal Endothelial Cells",
                              "Hepatic Stellate Cells (lipocyte)", "Hepatocyte", "immature langerhans cells",
                              "Intestinal epithelial cells (polarized)", "Iris Pigment Epithelial Cells",
                              "Keratinocyte - epidermal", "Keratinocyte - oral", "Keratocytes", "Lens Epithelial Cells",
                              "Macrophage - monocyte derived", "Mallassez-derived cells", "Mammary Epithelial Cell",
                              "Mast cell", "mature adipocyte", "Melanocyte", "Meningeal Cells",
                              "mesenchymal precursor cell - adipose", "mesenchymal precursor cell - bone marrow",
                              "mesenchymal precursor cell - cardiac", "Mesenchymal stem cells - adipose",
                              "Mesenchymal Stem Cells - amniotic membrane", "Mesenchymal Stem Cells - bone marrow",
                              "Mesenchymal stem cells - hepatic", "Mesenchymal stem cells - umbilical",
                              "Mesenchymal Stem Cells - Vertebral", "Mesenchymal Stem Cells - Wharton Jelly",
                              "Mesothelial Cells", "migratory langerhans cells",
                              "Multipotent Cord Blood Unrestricted Somatic Stem Cells", "Myoblast",
                              "nasal epithelial cells", "Natural Killer Cells", "Neural stem cells", "Neurons",
                              "Neutrophil", "Nucleus Pulposus Cell", "Olfactory epithelial cells",
                              "Oligodendrocyte - precursors", "Osteoblast", "Pancreatic stromal cells", "Pericytes",
                              "Perineurial Cells", "Peripheral Blood Mononuclear Cells", "Placental Epithelial Cells",
                              "Preadipocyte - breast", "Preadipocyte - omental", "Preadipocyte - perirenal",
                              "Preadipocyte - subcutaneous", "Preadipocyte - visceral", "promyelocytes",
                              "Prostate Epithelial Cells", "Prostate Stromal Cells", "Renal Cortical Epithelial Cells",
                              "Renal Epithelial Cells", "Renal Glomerular Endothelial Cells", "Renal Mesangial Cells",
                              "Renal Proximal Tubular Epithelial Cell", "Retinal Pigment Epithelial Cells",
                              "salivary acinar cells", "Schwann Cells", "Sebocyte", "Sertoli Cells",
                              "Skeletal Muscle Cells", "Skeletal muscle cells differentiated into Myotubes",
                              "Skeletal Muscle Satellite Cells", "Small Airway Epithelial Cells",
                              "Smooth muscle cells - airway", "Smooth Muscle Cells - Aortic",
                              "Smooth Muscle Cells - Bladder", "Smooth Muscle Cells - Brachiocephalic",
                              "Smooth Muscle Cells - Brain Vascular", "Smooth Muscle Cells - Bronchial",
                              "Smooth Muscle Cells - Carotid", "Smooth Muscle Cells - Colonic",
                              "Smooth Muscle Cells - Coronary Artery", "Smooth Muscle Cells - Esophageal",
                              "Smooth Muscle Cells - Internal Thoracic Artery", "Smooth Muscle Cells - Intestinal",
                              "Smooth Muscle Cells - Prostate", "Smooth Muscle Cells - Pulmonary Artery",
                              "Smooth Muscle Cells - Subclavian Artery", "Smooth Muscle Cells - Tracheal",
                              "Smooth Muscle Cells - Umbilical artery", "Smooth Muscle Cells - Umbilical Vein",
                              "Smooth Muscle Cells - Uterine", "Synoviocyte", "tenocyte", "Trabecular Meshwork Cells",
                              "Tracheal Epithelial Cells", "Urothelial cells", "Whole blood"]
complete_primary_non_include_list = {"Adipocyte - breast": "pre", "Adipocyte - omental": "pre",
                                     "Adipocyte - perirenal": "pre", "Adipocyte - subcutaneous": "pre",
                                     "Endothelial Cells - Vein": "Umbilical",
                                     "Renal Epithelial Cells": "Cortical",
                                     "Skeletal Muscle Cells": ["satellite", "differentiated"]}
complete_primary_exclude_list = ["mesenchymal precursor cell - ovarian", "Osteoblast - differentiated"]
complete_primary_jit_exclude_list = {"CD14+ CD16- Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD14+ CD16+ Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD14- CD16+ Monocytes": ("CD14+ Monocytes", "CD16"),
                                     "CD4+CD25+CD45RA- memory regulatory T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25+CD45RA+ naive regulatory T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25-CD45RA- memory conventional T cells": ("CD4+ T Cells", "CD25"),
                                     "CD4+CD25-CD45RA+ naive conventional T cells": ("CD4+ T Cells", "CD25")}
# endregion "Global Variables"

# region "test enhancer class methods"
x = Classes.Enhancers("human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",
                      "Human.sample_name2library_id.txt", complete_primary_cell_list,
                      celltype_exclude=complete_primary_exclude_list,
                      not_include=complete_primary_non_include_list,
                      partial_exclude=complete_primary_jit_exclude_list)
for key, value in x.partial_exclude_codes.items():
    print("Enhancers", key, value, len(value), sep="\n", end="\n\n")

#
# print(x.celltype, x.raw_data.shape, x.sample_types, x.file, x.parent_path, x.data.shape, x.codes, sep="\n", end="\n\n")

# endregion "test enhancer class methods"

# region "test promoter class methods"
y = Classes.Promoters("hg19.cage_peak_phase1and2combined_tpm.osc.txt", complete_primary_cell_list,
                      celltype_exclude=complete_primary_exclude_list, not_include=complete_primary_non_include_list,
                      partial_exclude=complete_primary_jit_exclude_list, nrows=1000)
for key, value in y.partial_exclude_codes.items():
    print("Promoters", key, value, len(value), sep="\n", end="\n\n")

# print(y.codes, len([item for items in y.codes.values() for item in items]), sep="\n")
# Defs.write_dict_to_csv("test.csv", y.codes, "/Tests/", path="parent")
# print(y.celltype, y.raw_data.shape, y.sample_types, y.file, y.parent_path, y.data.shape, y.codes, sep="\n")

# endregion "test promoter class methods"
