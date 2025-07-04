import json
import numpy as np
import matplotlib.pyplot as plt

# Path to your JSON file
file_path = 'outputs/val/BM_He/raw/He_norm.json'

try:
    with open(file_path, 'r') as f:
        # Load the JSON data into a Python list of lists
        data_list = json.load(f)

    # Convert the list of lists into a 3D NumPy array
    # Ensure the structure of data_list matches the expected 3D shape

    print("3D array loaded successfully:")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Example data
    _x = np.arange(len(data_list))
    _y = np.arange(len(data_list[0]))
    print(len(data_list[0]))
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()

    height = np.zeros_like(x)

    bottom = np.zeros_like(x)  # Base of the bars at z=0
    for i in range(len(data_list)):
        for j in range(len(data_list[0])):
            print(i*len(data_list[0])+j)
            height[i*len(data_list[0])+j]=len(data_list[i][j])

    # height = x + y  # Example height based on x and y values
    width = depth = 0.8 # Width and depth of the bars
    ax.bar3d(x, y, bottom, width, depth, height, shade=True)
    plt.show()

except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found.")
except json.JSONDecodeError:
    print(f"Error: Could not decode JSON from '{file_path}'. Check if the file contains valid JSON.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")