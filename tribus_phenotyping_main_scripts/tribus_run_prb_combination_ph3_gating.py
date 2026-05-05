import os
import pandas as pd
import numpy as np
import tribus
from visualization import heatmap_for_median_expression, marker_expression, umap_vis, z_score, cell_type_distribution

# -----------------------------
# Paths
# -----------------------------
data_path = "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Data/csv"
output_path = "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_PanCK_0_aSMA_1_pRb_combined_PH3_gated_600"
logic_table_path = "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/logic_tables"
logic_table_name = "logic_table_proliferation_11"

os.makedirs(output_path, exist_ok=True)

# -----------------------------
# Sample list
# -----------------------------
samples = ["S032_iOvaL_b1", "S091_iOme"]

# -----------------------------
# Marker columns
# -----------------------------
cols = ["DAPI_1", "CyclinE_1", "p27_1", "PanCK_1",
        "DAPI_2", "PCNA_2", "CyclinD1_2", "p21_2", "PH3_2",
        "DAPI_3", "Ki-67_3", "Geminin_3", "CyclinB1_failed", "Vimentin_3",
        "DAPI_4", "P-RPA32_4", "Jun_4", "HLA-DPB1_4", "pRb_4",
        "DAPI_5", "LaminB1_5", "PAX8_5", "H2AX_5", "bCatenin_failed",
        "DAPI_6", "N-Cadherin_6", "aSMA_6", "CyclinB1_6", "Iba1_6",
        "DAPI_7", "E-Cadherin_7", "CD3D_7", "Zic2_7", "bCatenin_7"]

# -----------------------------
# Load logic table ONCE
# -----------------------------
logic_xl = pd.ExcelFile(f"{logic_table_path}/{logic_table_name}.xlsx")
logic = pd.read_excel(logic_xl, logic_xl.sheet_names, index_col=0)

depth = 2

# -----------------------------
# Process each sample
# -----------------------------
for sample_name in samples:
    print(f"\nProcessing {sample_name}...")

    file_path = os.path.join(data_path, f"{sample_name}.csv")
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Load CSV
    df = pd.read_csv(file_path, low_memory=False)
    print(f"Original shape: {df.shape}")

    # Remove outliers
    Q = df[cols].quantile(0.999)
    df = df[~((df[cols] > Q)).any(axis=1)]
    print(f"After outlier removal: {df.shape}")

    # -----------------------------
    # Run TRIBUS
    # -----------------------------
    labels, scores = tribus.run_tribus(
        np.arcsinh(df[cols] / 5.0),
        logic,
        depth=depth,
        normalization=z_score,
        tuning=0,
        sigma=1,
        learning_rate=1,
        clustering_threshold=100,
        undefined_threshold=0.0005,
        other_threshold=0.4,
        random_state=42
    )

    # Join labels to the original dataframe
    result_data = df.join(labels)

    # -----------------------------
    # Combine pRb positive cells into Tumor
    # -----------------------------
    mask_pRb = result_data["final_label"] == "Tumor_pRb"
    result_data.loc[mask_pRb, "final_label"] = "Tumor"

    # -----------------------------
    # Gate PH3 positive tumor cells
    # -----------------------------
    if "PH3_2" in result_data.columns:
        ph3_threshold = 600  # adjust if needed
        mask_ph3 = (result_data["PH3_2"] > ph3_threshold) & (result_data["final_label"] == "Tumor")
        result_data.loc[mask_ph3, "final_label"] = "Tumor_PH3"
        print(f"{sample_name}: {mask_ph3.sum()} cells gated as Tumor_PH3")

    # -----------------------------
    # Save results

    labels_new = result_data[[
        "CellID",
        "X_centroid",
        "Y_centroid",
        "final_label"
    ]]

    labels_new.to_csv(
        os.path.join(
            output_path,
            f"{sample_name}_{logic_table_name}_tribus_annotation_pRb_combined_PH3_gated.csv"
        ),
        index=False
    )

    # Full dataset with markers + updated labels
    result_data.to_csv(
        os.path.join(
            output_path,
            f"{sample_name}_{logic_table_name}_raw_tribus_annotated_pRb_combined_PH3_gated.csv"
        ),
        index=False
    )

    print(f"{sample_name} saved with updated labels.")