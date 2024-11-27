import os
from PIL import Image

def create_gif_from_images(image_folder, output_gif_path, duration=200):
    """
    Create a GIF from a sequence of images in a folder.

    Parameters:
        image_folder (str): Path to the folder containing the images.
        output_gif_path (str): Path to save the output GIF.
        duration (int): Duration (in milliseconds) each frame will be displayed in the GIF.

    Returns:
        None
    """
    print(f"Loading images from folder: {image_folder}")
    print(f"Saving GIF to: {output_gif_path}")

    # List all image files in the folder
    image_files = sorted(
        [os.path.join(image_folder, f) for f in os.listdir(image_folder) if f.endswith(('.png', '.jpg', '.jpeg'))]
    )
    
    if not image_files:
        print("No image files found in the folder. Exiting.")
        return

    print(f"Found {len(image_files)} images to process.")

    # Open images and store them in a list
    images = []
    for image_file in image_files:
        try:
            img = Image.open(image_file)
            images.append(img)
        except Exception as e:
            print(f"Error loading image {image_file}: {e}")
            continue

    # Save the images as a GIF
    try:
        images[0].save(
            output_gif_path,
            save_all=True,
            append_images=images[1:],
            duration=duration,
            loop=0  # Loop infinitely
        )
        print(f"GIF created successfully: {output_gif_path}")
    except Exception as e:
        print(f"Error creating GIF: {e}")

# Folder containing images
image_folder = r"C:\sem7\btp\BTP\grain_images_with_holes\20_100\output\grain"

# Path to save the GIF
output_gif_path = r"C:\sem7\btp\BTP\grain_animation.gif"

# Create GIF
create_gif_from_images(image_folder, output_gif_path, duration=200)  # Duration of 200ms per frame
