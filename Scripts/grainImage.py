import os
import matplotlib.pyplot as plt
import numpy as np

def create_grain_images_with_neighbors_and_vertices(input_directory, output_directory):
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
            vertices_image = os.path.join(save_path, file.replace(".dat", "_vertices.png"))

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
                vertices_grid = np.zeros((max_x, max_y, 3), dtype=np.uint8)
                vertices_grid.fill(255)  # White background
                
                # Track vulnerable zones
                vulnerable_zones = np.zeros((max_x, max_y), dtype=bool)
                grain_neighbors = {}

                # First pass: Fill grids and identify vulnerable zones
                for x, y, grain_id, value in data:
                    grain_grid[x, y] = grain_id
                    if grain_id not in grain_neighbors:
                        grain_neighbors[grain_id] = set()

                    if value == 1:
                        for i in range(-2, 3):
                            for j in range(-2, 3):
                                if 0 <= x + i < max_x and 0 <= y + j < max_y:
                                    image_grid[x + i, y + j] = [255, 0, 0]  # Red vulnerable zone
                                    vulnerable_zones[x + i, y + j] = True
                    else:
                        image_grid[x, y] = [135, 206, 235]  # Sky Blue

                grain_centers = {}
                vertices = set()  # To store vertices

                # Second pass: Identify grain boundaries and vertices
                for x in range(max_x):
                    for y in range(max_y):
                        if grain_grid[x, y] == -1:
                            continue
                        
                        grain_id = grain_grid[x, y]
                        if grain_id not in grain_centers:
                            grain_centers[grain_id] = []
                        grain_centers[grain_id].append((x, y))

                        # Check for grain boundaries
                        neighbors = [
                            (x + dx, y + dy)
                            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]
                            if 0 <= x + dx < max_x and 0 <= y + dy < max_y
                        ]
                        
                        is_boundary = False
                        for nx, ny in neighbors:
                            if grain_grid[nx, ny] != grain_id and grain_grid[nx, ny] != -1:
                                image_grid[x, y] = [0, 0, 0]  # Black boundary
                                is_boundary = True
                                grain_neighbors[grain_id].add(grain_grid[nx, ny])
                        
                        # If this is a boundary pixel, check if it's a vertex by counting different neighboring grains
                        if is_boundary:
                            adjacent_grains = set()
                            for dx in [-1, 0, 1]:
                                for dy in [-1, 0, 1]:
                                    nx, ny = x + dx, y + dy
                                    if 0 <= nx < max_x and 0 <= ny < max_y and grain_grid[nx, ny] != -1:
                                        adjacent_grains.add(grain_grid[nx, ny])
                            
                            # If the pixel touches 3 or more grains, it's a vertex
                            if len(adjacent_grains) >= 3:
                                vertices.add((x, y))
                                if vulnerable_zones[x, y]:
                                    vertices_grid[x, y] = [255, 0, 0]  # Red vertices in vulnerable zones
                                else:
                                    vertices_grid[x, y] = [0, 0, 255]  # Blue vertices outside vulnerable zones

                # Generate main image with neighbor counts
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

                # Generate vertices-only image
                plt.figure(figsize=(10, 10))
                plt.imshow(vertices_grid, origin="upper")
                
                # Enlarge the vertices for better visibility
                for x, y in vertices:
                    circle = plt.Circle((y, x), radius=2, color="blue" if not vulnerable_zones[x, y] else "red")
                    plt.gca().add_patch(circle)
                
                plt.axis("off")
                plt.savefig(vertices_image, bbox_inches="tight", pad_inches=0, dpi=300)
                plt.close()

                print(f"Saved: {output_image} and {vertices_image}")

            except Exception as e:
                print(f"Error processing {input_file}: {e}")
                continue

    print("Processing completed.")

# Update paths
input_directory = "./processed"  
output_directory = "./grain_images_with_neighbors"  

create_grain_images_with_neighbors_and_vertices(input_directory, output_directory)