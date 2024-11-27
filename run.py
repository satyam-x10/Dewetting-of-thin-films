import subprocess
import os

def run_script(script_name):
    # Define the path to the script
    script_path = os.path.join("Scripts", script_name)

    # Check if the script exists
    if not os.path.exists(script_path):
        print(f"Error: {script_path} does not exist.")
        return False

    # Run the script
    try:
        print(f"Running {script_path}...")
        subprocess.run(["python", script_path], check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error while executing {script_path}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
    return False

def main():
    # Run the first script
    if not run_script("isVulnerable.py"):
        print("Failed to execute isVulnerable.py. Exiting...")
        return

    # Run the second script
    if not run_script("grainImage.py"):
        print("Failed to execute grainImage.py. Exiting...")
        return

    print("All scripts executed successfully.")

if __name__ == "__main__":
    main()
