import pandas as pd
import joblib

def load_test_data(file_path):
    """
    Load a .dat test file with 6 columns and assign appropriate column names.
    """
    data = pd.read_csv(file_path, delim_whitespace=True, header=None)
    
    # Assign column names for test data
    col_names = ['x', 'y', 'comp', 'grain_id', 'eta^2', 'eta']
    if data.shape[1] != 6:
        raise ValueError(f"Expected 6 columns, but found {data.shape[1]} in {file_path}")
    
    data.columns = col_names
    return data

def test_model(model, test_data):
    """
    Test the loaded model on the test data and make predictions.
    """
    # Extract relevant features
    feature_cols = ['comp', 'eta^2', 'eta']
    X_test = test_data[feature_cols]
    
    # Make predictions
    predictions = model.predict(X_test)
    
    # Attach predictions to the test data
    test_data['vulnerable_prediction'] = predictions
    
    # Display sample predictions
    print("Sample predictions:")
    print(test_data[['x', 'y', 'comp', 'grain_id', 'eta^2', 'eta', 'vulnerable_prediction']].head())
    
    return test_data

# Main code to load the model and test on new data
if __name__ == "__main__":
    # Path to the saved model
    model_file_path = r"C:\sem7\btp\BTP\trained_model_multiple.pkl"
    
    # Load the pre-trained model
    try:
        loaded_model = joblib.load(model_file_path)
        print(f"Model loaded successfully from {model_file_path}")
    except FileNotFoundError:
        print(f"Model file not found at {model_file_path}. Ensure the model is trained and saved.")
        exit()
    
    # Path to the test .dat file
    test_file_path = r"C:\sem7\btp\data\run\ms_100_24\output\grain\2d_grain142000.dat"
    
    # Load and test the data
    try:
        test_data = load_test_data(test_file_path)
        predictions = test_model(loaded_model, test_data)
        
        # Save predictions to a CSV file for further analysis
        predictions.to_csv("test_predictions.csv", index=False, sep='\t')
        print("Predictions saved to 'test_predictions.csv'")
    except Exception as e:
        print(f"Error: {e}")
