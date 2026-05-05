# Standard library
import os
import warnings
warnings.filterwarnings('ignore')

# Scientific computing
import pandas as pd


# -----------------------------
# Define ALL input/output pairs
# -----------------------------
path_pairs = [

    (
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_Vimentin_0_pRb_combined",
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_Vimentin_0_pRb_combined_PH3_gated_600"
    ),

    (
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_PanCK_0_pRb_combined",
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_PanCK_0_pRb_combined_PH3_gated_600"
    ),

    (
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_PanCK_0_aSMA_1_pRb_combined",
        "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_PanCK_0_aSMA_1_pRb_combined_PH3_gated_600"
    )
]


# -----------------------------
# Loop through each folder pair
# -----------------------------
for input_folder, output_folder in path_pairs:

    print(f"\nProcessing folder: {input_folder}")

    os.makedirs(output_folder, exist_ok=True)

    all_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

    raw_files = [f for f in all_files if "raw_tribus_annotated" in f]
    annotation_files = [f for f in all_files if "tribus_annotation" in f]

    # Map annotation files by prefix
    annotation_dict = {}
    for ann_file in annotation_files:
        prefix = ann_file.split("tribus_annotation")[0]
        annotation_dict[prefix] = ann_file

    # Process each raw file
    for raw_file in raw_files:

        prefix = raw_file.split("raw_tribus_annotated")[0]

        raw_path = os.path.join(input_folder, raw_file)
        raw_df = pd.read_csv(raw_path)

        if "PH3_2" not in raw_df.columns:
            print(f"Skipping {raw_file} (PH3_2 not found)")
            continue

        # PH3 gating mask
        ph3_mask = (
            (raw_df["PH3_2"] > 600) &
            (raw_df["final_label"].str.contains("Tumor", case=False, na=False))
        )

        n_changed = ph3_mask.sum()

        # Update RAW
        raw_df.loc[ph3_mask, "final_label"] = "Tumor_PH3"

        raw_out = os.path.join(
            output_folder,
            os.path.splitext(raw_file)[0] + "_PH3_gated.csv"
        )
        raw_df.to_csv(raw_out, index=False)

        # Update annotation if matching file exists
        if prefix in annotation_dict:

            ann_file = annotation_dict[prefix]
            ann_path = os.path.join(input_folder, ann_file)
            ann_df = pd.read_csv(ann_path)

            if len(ann_df) != len(raw_df):
                print(f"Row mismatch in {prefix} — skipping annotation update")
            else:
                ann_df.loc[ph3_mask, "final_label"] = "Tumor_PH3"

                ann_out = os.path.join(
                    output_folder,
                    os.path.splitext(ann_file)[0] + "_PH3_gated.csv"
                )
                ann_df.to_csv(ann_out, index=False)

        if n_changed > 0:
            print(f"{prefix}: {n_changed} cells gated as Tumor_PH3")

print("\nAll folders processed successfully.")