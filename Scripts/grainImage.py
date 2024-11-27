import os
import matplotlib.pyplot as plt
import numpy as np

def create_grain_images_with_boundaries(input_directory, output_directory):
    print("Generating grain boundary images with red regions...")
    print(f"Input Directory: {input_directory}")
    print(f"Output Directory: {output_directory}")
    print("-" * 50)

    for folder, subfolders, files in os.walk(input_directory):
        for file in files:
            if file.startswith("2d_grain") and file.endswith(".dat"):
                input_file = os.path.join(folder, file)
                print(f"Processing file: {input_file}")

                # Determine save path for the image
                relative_path = os.path.relpath(folder, input_directory)
                save_path = os.path.join(output_directory, relative_path)
                os.makedirs(save_path, exist_ok=True)
                output_image = os.path.join(save_path, file.replace(".dat", ".png"))

                # Read the .dat file
                try:
                    data = []
                    with open(input_file, "r") as infile:
                        for line in infile:
                            columns = line.strip().split()
                            try:
                                x = int(columns[0])  # X-coordinate
                                y = int(columns[1])  # Y-coordinate
                                grain_id = int(float(columns[-4]))  # Grain ID (ensure integer)
                                value = int(columns[-1])  # True (1) or False (0)
                                data.append((x, y, grain_id, value))
                            except ValueError:
                                print(f"Skipping invalid row: {line.strip()}")
                                continue

                    # Get grid dimensions
                    max_x = max(row[0] for row in data) + 1
                    max_y = max(row[1] for row in data) + 1

                    # Create grid for grain IDs and visualization
                    grain_grid = np.full((max_x, max_y), -1, dtype=int)  # For grain IDs
                    image_grid = np.zeros((max_x, max_y, 3), dtype=np.uint8)  # For visualization

                    for x, y, grain_id, value in data:
                        grain_grid[x, y] = grain_id  # Populate grain grid
                        if value == 1:
                            # Highlight a 5x5 area for red pixels (indicating specific condition)
                            image_grid[x , y ] = [255, 0, 0]  # Red
                        else:
                            image_grid[x, y] = [255, 255, 0]  # Yellow for other points

                    # Identify and mark grain boundaries
                    for x in range(max_x):
                        for y in range(max_y):
                            if grain_grid[x, y] == -1:
                                continue
                            # Check neighbors for boundary
                            neighbors = [
                                (x + dx, y + dy)
                                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]
                                if 0 <= x + dx < max_x and 0 <= y + dy < max_y
                            ]
                            for nx, ny in neighbors:
                                if grain_grid[nx, ny] != grain_grid[x, y]:
                                    image_grid[x, y] = [0, 0, 0]  # Black boundary

                    # Plot and save the image
                    plt.figure(figsize=(10, 10))
                    plt.imshow(image_grid, origin="upper")
                    plt.axis("off")

                    # Save the image
                    plt.savefig(output_image, bbox_inches="tight", pad_inches=0)
                    plt.close()
                    print(f"Image saved: {output_image}")

                except Exception as e:
                    print(f"Error processing file {input_file}: {e}")
                    continue

    print("Grain boundary image generation with red regions completed.")

# Replace the paths with your directory structure
input_directory = "./processed"  # Directory containing the .dat files
output_directory = "./grain_images_with_red_regions"  # Directory to save generated images

create_grain_images_with_boundaries(input_directory, output_directory)
