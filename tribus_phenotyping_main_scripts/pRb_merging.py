# Input and output folders
input_folder = "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_Vimentin_0"
output_folder = "/home/ad/P-drive/h30492/farkkilab2/4_CellCycle/tcycif/4_10_Expansion/Analysis/01-phenotyping/tribus/labels_Vimentin_0_pRb_combined"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop through CSV files
for filename in os.listdir(input_folder):

    if filename.endswith(".csv"):

        input_path = os.path.join(input_folder, filename)

        # Read CSV
        df = pd.read_csv(input_path)

        # Modify labels
        df.loc[df["final_label"] == "Tumor_pRb", "final_label"] = "Tumor"

        # Create new filename with suffix
        base_name = os.path.splitext(filename)[0]
        new_filename = f"{base_name}_pRb_combined.csv"

        output_path = os.path.join(output_folder, new_filename)

        # Save modified file
        df.to_csv(output_path, index=False)

print("All files processed successfully.")