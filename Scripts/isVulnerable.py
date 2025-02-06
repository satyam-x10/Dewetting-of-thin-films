import os

def process_dat_files(root_directory, save_directory):
    print("Starting processing...")
    print(f"Root Directory: {root_directory}")
    print(f"Save Directory: {save_directory}")
    print("-" * 50)

    # Traverse all subdirectories in the root directory
    for folder, subfolders, files in os.walk(root_directory):
        # Debug: Display the current folder being checked
        print(f"Checking folder: {folder}")

        # Look for "output/grain" in the folder path
        if os.path.basename(folder) == "grain" and "output" in os.path.dirname(folder):
            print(f"Found target folder: {folder}")

            # Get relative path and determine save location
            relative_path = os.path.relpath(folder, root_directory)  # Get relative path
            save_path = os.path.join(save_directory, relative_path)  # Destination folder
            print(f"Save path for processed files: {save_path}")
            
            # Ensure the destination folder exists
            os.makedirs(save_path, exist_ok=True)

            # Process all .dat files in the current folder
            for file in files:
                if file.endswith(".dat"):
                    input_file = os.path.join(folder, file)
                    output_file = os.path.join(save_path, file)
                    print(f"Processing file: {input_file}")

                    try:
                        # Process the file
                        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
                            for line in infile:
                                columns = line.strip().split()
                                
                                # Debug: Display the parsed columns
                                
                                # Identify the absolute value of the 4th-to-last column
                                try:
                                    value_to_check = abs(float(columns[-4]))  # Take the absolute value
                                except (ValueError, IndexError) as e:
                                    # print(f"Error parsing value in line: {line.strip()}")
                                    # print(f"Skipping line due to error: {e}")
                                    continue
                                
                                # Add a new column based on the condition
                                new_column = 1 if 0.45 <= value_to_check <= 0.55 else 0
                                columns.append(str(new_column))
                                
                                # Write the modified line to the output file
                                outfile.write("\t".join(columns) + "\n")
                        
                        print(f"Processed file saved: {output_file}")
                    except Exception as e:
                        print(f"Error processing file {input_file}: {e}")
                        continue

    print("Processing completed.")


root_directory = "G:\BTP\data"  # Path containing folders like 20_100, 24_100, etc.
save_directory = "./processed"  # Path to save processed files

process_dat_files(root_directory, save_directory)
