#This code CLEANS two data sets so that we can look at LUAD and LUSC data and TP53 gene expression for our project 
import pandas as pd
from pathlib import Path

#Path to the gene expression file 
file_path = Path(r"C:\Users\logan\Documents\GitHub\desktop-tutorial\New folder\CancerModule\tp53 data\GSE62944_subsample_log2TPM.csv")
expr = pd.read_csv(file_path, index_col=0)

#Filter all tp53 related gene expressions (genes starting with TP53)
#UPDATED: keep ALL genes starting with TP53, ALL BCL, and ALL CASP
tp53_all = expr.index[expr.index.str.startswith("TP53")].tolist()
bcl_all  = expr.index[expr.index.str.startswith("BCL")].tolist()
casp_all = expr.index[expr.index.str.startswith("CASP")].tolist()

genes_to_keep = tp53_all + bcl_all + casp_all

expr_tp53 = expr.loc[expr.index.intersection(genes_to_keep)]

#Fix TCGA ID formatting so it matches expression column format
expr_tp53.columns = expr_tp53.columns.str[:12]

#Path to the metadata file
file_path_meta = Path(r"C:\Users\logan\Documents\GitHub\desktop-tutorial\New folder\CancerModule\meta data\GSE62944_metadata.csv")
meta = pd.read_csv(file_path_meta)

#Filter for BOTH LUSC and LUAD (instead of only LUSC)
meta_filtered = meta[meta["cancer_type"].isin(["LUSC"])]

#Remove any patients with missing barcodes (NaN)
meta_filtered = meta_filtered.dropna(subset=["bcr_patient_barcode"])

#Match filtered patient IDs to expression columns
filtered_barcodes = meta_filtered["bcr_patient_barcode"]
expr_tp53_filtered = expr_tp53[filtered_barcodes]

#Transpose expression so patients are rows instead of one row
expr_tp53_filtered = expr_tp53_filtered.T

#Combine metadata and TP53-family expression into one file
merged = meta_filtered.merge(expr_tp53_filtered, left_on="bcr_patient_barcode", right_index=True)

#Keep only patient id, cancer type, tumor stage, and TP53 family genes
cols_to_keep = ["bcr_patient_barcode", "cancer_type", "ajcc_pathologic_tumor_stage"] + list(expr_tp53_filtered.columns)
merged = merged[cols_to_keep]

#Put clean combined data in another file
output = file_path_meta.parent / "LUSC_TP53_BCL_CASP_combined.csv"
merged.to_csv(output, index=False)
