import os

def process_dat_files(root_directory, save_directory):
    print("Starting processing...")
    print(f"Root Directory: {root_directory}")
    print(f"Save Directory: {save_directory}")
    print("-" * 50)

    total_files = 0
    total_folders = 0
    processed_files = 0
    processed_folders = 0

    # First, count total files and folders
    for folder, subfolders, files in os.walk(root_directory):
        if os.path.basename(folder) == "grain" and "output" in os.path.dirname(folder):
            total_folders += 1
            total_files += len([file for file in files if file.endswith(".dat") and file.startswith("2d")])

    print(f"Total folders to process: {total_folders}")
    print("-" * 50)

    # Traverse all subdirectories in the root directory
    for folder, subfolders, files in os.walk(root_directory):
        if os.path.basename(folder) == "grain" and "output" in os.path.dirname(folder):
            print(f" ({processed_folders + 1}/{total_folders}) {folder} Processing folder")
            processed_folders += 1

            # Get relative path and determine save location
            relative_path = os.path.relpath(folder, root_directory)  # Get relative path
            save_path = os.path.join(save_directory, relative_path)  # Destination folder
            os.makedirs(save_path, exist_ok=True)

            # Process all .dat files in the current folder
            for file in files:
                if file.endswith(".dat") and file.startswith("2d"):
                    input_file = os.path.join(folder, file)
                    output_file = os.path.join(save_path, file)
                    processed_files += 1

                    print(f"Processing file:({processed_files}/{total_files}) {input_file} ")

                    try:
                        # Process the file
                        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
                            for line in infile:
                                columns = line.strip().split()

                                try:
                                    value_to_check = abs(float(columns[-4]))  # Take the absolute value
                                except (ValueError, IndexError) as e:
                                    continue

                                # Add a new column based on the condition
                                new_column = 1 if 0.45 <= value_to_check <= 0.55 else 0
                                columns.append(str(new_column))

                                # Write the modified line to the output file
                                outfile.write("\t".join(columns) + "\n")

                    except Exception as e:
                        print(f"Error processing file {input_file}: {e}")
                        continue

    print("Processing completed.")

# Replace the paths with your directory structure
root_directory = r"C:\sem7\btp\data\run"  # Path containing folders like 20_100, 24_100, etc.
save_directory = "./processed"  # Path to save processed files

process_dat_files(root_directory, save_directory)
