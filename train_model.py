import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import joblib

# Function to load a single file
def load_data(file_path, is_test=False):
    """
    Load a .dat file and assign proper column names.
    For test data, the `vulnerable` column is absent.
    """
    data = pd.read_csv(file_path, delim_whitespace=True, header=None)
    
    # Define column names
    if is_test:
        col_names = ['x', 'y', 'comp', 'grain_id', 'eta^2', 'eta']
    else:
        col_names = ['x', 'y', 'comp', 'grain_id', 'eta^2', 'eta', 'vulnerable']
    
    data.columns = col_names
    return data

# Function to list all files starting with "2d" in a directory
def list_files(directory, prefix="2d"):
    """
    List all files in a directory starting with a specific prefix.
    """
    files = [os.path.join(directory, file) for file in os.listdir(directory) if file.startswith(prefix)]
    print(f"Found {len(files)} files starting with '{prefix}' in {directory}.")
    return files

# Function to load and combine multiple files
def load_multiple_files(file_paths):
    """
    Load and combine data from multiple .dat files.
    """
    combined_data = pd.DataFrame()
    
    for file_path in file_paths:
        data = load_data(file_path, is_test=False)  # Load training data
        combined_data = pd.concat([combined_data, data], ignore_index=True)
    
    print(f"Combined data shape: {combined_data.shape}")
    return combined_data

# Function to prepare data for training
def prepare_data(data):
    """
    Prepare the dataset by excluding `x`, `y`, and `grain_id` as features and separating the target (`vulnerable`).
    """
    # Exclude x, y, and grain_id from the features
    feature_cols = ['comp', 'eta^2', 'eta']
    X = data[feature_cols]  # Features
    y = data['vulnerable']  # Target column
    return train_test_split(X, y, test_size=0.2, random_state=42)  # 80-20 split

# Function to train and evaluate the model
def train_and_evaluate(X_train, X_test, y_train, y_test):
    """
    Train a RandomForest model and evaluate its accuracy.
    """
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Model Accuracy: {accuracy * 100:.2f}%")
    return model

# Function to save the model to a file
def save_model(model, file_path):
    """
    Save the trained model using joblib.
    """
    joblib.dump(model, file_path)
    print(f"Model saved to {file_path}")

# Main function to load, combine, and train the model
def train_on_multiple_files(directory):
    """
    Train the model using data from multiple .dat files in a directory.
    """
    # Get list of files starting with "2d"
    file_paths = list_files(directory, prefix="2d")
    
    # Load and combine data
    combined_data = load_multiple_files(file_paths)
    
    # Prepare the data for training
    X_train, X_test, y_train, y_test = prepare_data(combined_data)
    
    # Train the model
    trained_model = train_and_evaluate(X_train, X_test, y_train, y_test)
    
    # Save the trained model
    model_file = "trained_model_multiple.pkl"
    save_model(trained_model, model_file)
    return trained_model

# Example usage
if __name__ == "__main__":
    # Directory containing the .dat files
    directory = r"C:\sem7\btp\BTP\processed\20_100\output\grain"
    
    # Train the model on files starting with "2d"
    trained_model = train_on_multiple_files(directory)
