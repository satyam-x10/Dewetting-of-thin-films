import os
import matplotlib.pyplot as plt
import numpy as np

def create_grain_images_with_neighbors(input_directory, output_directory):
    print("Processing grain images...")

    for folder, subfolders, files in os.walk(input_directory):
        dat_files = [f for f in files if f.startswith("2d_grain") and f.endswith(".dat")]

        for file in dat_files:
            input_file = os.path.join(folder, file)
            print(f"Processing: {input_file}")

            relative_path = os.path.relpath(folder, input_directory)
            save_path = os.path.join(output_directory, relative_path)
            os.makedirs(save_path, exist_ok=True)
            output_image = os.path.join(save_path, file.replace(".dat", ".png"))

            try:
                data = []
                with open(input_file, "r") as infile:
                    for line in infile:
                        columns = line.strip().split()
                        try:
                            x = int(columns[0])
                            y = int(columns[1])
                            grain_id = int(float(columns[-4]))  
                            value = int(columns[-1])  
                            data.append((x, y, grain_id, value))
                        except ValueError:
                            continue  

                max_x = max(row[0] for row in data) + 1
                max_y = max(row[1] for row in data) + 1

                grain_grid = np.full((max_x, max_y), -1, dtype=int)
                image_grid = np.zeros((max_x, max_y, 3), dtype=np.uint8)
                grain_neighbors = {}

                for x, y, grain_id, value in data:
                    grain_grid[x, y] = grain_id
                    if grain_id not in grain_neighbors:
                        grain_neighbors[grain_id] = set()

                    if value == 1:
                        for i in range(-2, 3):
                            for j in range(-2, 3):
                                if 0 <= x + i < max_x and 0 <= y + j < max_y:
                                    image_grid[x + i, y + j] = [255, 0, 0]  
                    else:
                        image_grid[x, y] = [135, 206, 235]  # Sky Blue

                grain_centers = {}

                for x in range(max_x):
                    for y in range(max_y):
                        if grain_grid[x, y] == -1:
                            continue
                        grain_id = grain_grid[x, y]
                        if grain_id not in grain_centers:
                            grain_centers[grain_id] = []
                        grain_centers[grain_id].append((x, y))

                        neighbors = [
                            (x + dx, y + dy)
                            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]
                            if 0 <= x + dx < max_x and 0 <= y + dy < max_y
                        ]
                        for nx, ny in neighbors:
                            if grain_grid[nx, ny] != grain_id and grain_grid[nx, ny] != -1:
                                image_grid[x, y] = [0, 0, 0]  # Black boundary
                                grain_neighbors[grain_id].add(grain_grid[nx, ny])

                # Find center of each grain and display neighbor count
                plt.figure(figsize=(10, 10))
                plt.imshow(image_grid, origin="upper")

                for grain_id, pixels in grain_centers.items():
                    if pixels:
                        cx, cy = np.mean(pixels, axis=0).astype(int)
                        neighbor_count = len(grain_neighbors.get(grain_id, []))
                        plt.text(cy, cx, str(neighbor_count), color="black", fontsize=6, ha="center", va="center", fontweight="bold")

                plt.axis("off")
                plt.savefig(output_image, bbox_inches="tight", pad_inches=0, dpi=300)
                plt.close()

                print(f"Saved: {output_image}")

            except Exception as e:
                print(f"Error processing {input_file}: {e}")
                continue

    print("Processing completed.")

# Update paths
input_directory = "./processed"  
output_directory = "./grain_images_with_neighbors"  

create_grain_images_with_neighbors(input_directory, output_directory)
