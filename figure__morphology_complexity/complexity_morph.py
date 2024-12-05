import os
import pandas as pd
import matplotlib.pyplot as plt

# Directories and file paths
complexity_folder = "full_complexity"
swc_folder = "swc"
soma_fire_rate_file = "soma_fire_rate.csv"

# List all files in the complexity folder and sort by full identifier (cell type and name)
complexity_files = sorted([f for f in os.listdir(complexity_folder) if f.endswith(".csv")])
cell_names_full = sorted(set(f.split(".csv")[0] for f in complexity_files))  # Includes cell type
cell_names_full = [name for name in cell_names_full if "seed_1" in name]

# Read soma fire rate data
soma_fire_rate_data = pd.read_csv(soma_fire_rate_file)

# Function to set axis limits for equal scaling
def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = x_limits[1] - x_limits[0]
    y_range = y_limits[1] - y_limits[0]
    z_range = z_limits[1] - z_limits[0]
    max_range = max(x_range, y_range, z_range)
    x_middle = (x_limits[0] + x_limits[1]) / 2
    y_middle = (y_limits[0] + y_limits[1]) / 2
    z_middle = (z_limits[0] + z_limits[1]) / 2
    ax.set_xlim(x_middle - max_range / 2, x_middle + max_range / 2)
    ax.set_ylim(y_middle - max_range / 2, y_middle + max_range / 2)
    ax.set_zlim(z_middle - max_range / 2, z_middle + max_range / 2)

# Group cells by cell type
cell_types = sorted(set(f.split("__")[0] for f in complexity_files))

# Loop through each cell name
plots_per_file = 20
fig, axes = plt.subplots(4, 5, figsize=(20, 16), subplot_kw={'projection': '3d'})  # 4x5 grid
axes = axes.flatten()  # Flatten to access each subplot
current_plot = 0

for idx, cell_name_full in enumerate(cell_names_full):
    # Get all files for the cell name and load the seed 1 data
    seed_1_file = next(f for f in complexity_files if cell_name_full in f)
    cell_type = seed_1_file.split("__")[0]  # Extract cell type from file name
    print(cell_name_full)

    table_1 = pd.read_csv(os.path.join(complexity_folder, seed_1_file))

    # Find corresponding .swc file
    swc_file = next(f for f in os.listdir(swc_folder) if cell_name_full.split("__")[1].split("_")[0] in f)
    table_2 = pd.read_csv(os.path.join(swc_folder, swc_file), header=None, delim_whitespace=True)

    # Process soma fire rate
    table_3 = soma_fire_rate_data[soma_fire_rate_data.iloc[:, 0].str.contains(cell_name_full.split("__")[1].split("_")[0])]
    table_3 = table_3[soma_fire_rate_data.iloc[:, 0].str.contains("seed_1")]
    if len(table_3.iloc[:, 2:].max(axis=1)) > 0:
        max_rate = table_3.iloc[:, 2:].max(axis=1).values[0]
        max_rate_string = table_3.iloc[:, 2:].idxmax(axis=1).values[0]
        matching_column_indices = [i for i, col in enumerate(table_1.columns) if max_rate_string in col]

        # Prepare scatter plot data
        scatter_data = table_1.iloc[:, matching_column_indices]
        start_points = table_1.iloc[:, 0].astype(int)
        end_points = table_1.iloc[:, 1].astype(int)

        all_x_coords = []
        all_y_coords = []
        all_z_coords = []
        all_sizes = []
        all_colors = []

        for start, end, data_value in zip(start_points, end_points, scatter_data.values):
            # Get all points between start and end (inclusive)
            range_indices = table_2[(table_2.iloc[:, 0] >= start) & (table_2.iloc[:, 0] <= end)].index
            all_x_coords.extend(table_2.iloc[range_indices, 2].values)
            all_y_coords.extend(table_2.iloc[range_indices, 3].values)
            all_z_coords.extend(table_2.iloc[range_indices, 4].values)
            all_sizes.extend(table_2.iloc[range_indices, 5].values * 10)  # Make scatters 2x larger
            all_colors.extend([data_value * 0.001] * len(range_indices))  # Normalize color as in the working version

        # Create scatter plot in the next subplot
        ax = axes[current_plot]
        sc = ax.scatter(all_x_coords, all_y_coords, all_z_coords, c=all_colors, s=all_sizes, cmap="viridis", alpha=0.8)
        ax.set_facecolor('white')  # Remove grey background
        ax.view_init(elev=90, azim=-90)  # Rotate to XY plane
        ax._axis3don = False  # Remove the 3D box lines
        set_axes_equal(ax)  # Equal axis scaling
        ax.set_title(f"{cell_type} - {cell_name_full.split('__')[1].split('_')[0]}", fontsize=20)  # Larger font for title

        current_plot += 1

        # Save the figure and reset after 20 plots
        if current_plot == plots_per_file or (idx + 1) == len(cell_names_full):
            plt.tight_layout()
            plt.savefig(f"scatter_plot_batch_{(idx // plots_per_file) + 1}.png", dpi=300)
            plt.close()
            fig, axes = plt.subplots(4, 5, figsize=(20, 16), subplot_kw={'projection': '3d'})
            axes = axes.flatten()
            current_plot = 0
# Create summary plot with one cell per type
fig, axes = plt.subplots(1, len(cell_types), figsize=(5 * len(cell_types), 6), subplot_kw={'projection': '3d'})
if len(cell_types) == 1:
    axes = [axes]  # Ensure axes is iterable for a single subplot

for ax, cell_type in zip(axes, cell_types):
    # Find the longest SWC for this cell type
    cells_of_type = [name for name in cell_names_full if name.startswith(cell_type)]
    longest_cell = None
    longest_swc_length = -1

    for cell_name_full in cells_of_type:
        swc_file = next(f for f in os.listdir(swc_folder) if cell_name_full.split("__")[1].split("_")[0] in f)
        swc_path = os.path.join(swc_folder, swc_file)
        swc_length = sum(1 for _ in open(swc_path))  # Count lines in the SWC file
        if swc_length > longest_swc_length:
            longest_swc_length = swc_length
            longest_cell = cell_name_full

    if longest_cell:
        seed_1_file = next(f for f in complexity_files if longest_cell in f)
        table_1 = pd.read_csv(os.path.join(complexity_folder, seed_1_file))

        swc_file = next(f for f in os.listdir(swc_folder) if longest_cell.split("__")[1].split("_")[0] in f)
        table_2 = pd.read_csv(os.path.join(swc_folder, swc_file), header=None, delim_whitespace=True)

        scatter_data = table_1.iloc[:, matching_column_indices]
        start_points = table_1.iloc[:, 0].astype(int)
        end_points = table_1.iloc[:, 1].astype(int)

        all_x_coords = []
        all_y_coords = []
        all_z_coords = []
        all_sizes = []
        all_colors = []

        for start, end, data_value in zip(start_points, end_points, scatter_data.values):
            # Get all points between start and end (inclusive)
            range_indices = table_2[(table_2.iloc[:, 0] >= start) & (table_2.iloc[:, 0] <= end)].index
            all_x_coords.extend(table_2.iloc[range_indices, 2].values)
            all_y_coords.extend(table_2.iloc[range_indices, 3].values)
            all_z_coords.extend(table_2.iloc[range_indices, 4].values)
            all_sizes.extend(table_2.iloc[range_indices, 5].values * 10)  # Make scatters 2x larger
            all_colors.extend([data_value * 0.001] * len(range_indices))  # Normalize color as in the working version

        # Create scatter plot for this cell type
        sc = ax.scatter(all_x_coords, all_y_coords, all_z_coords, c=all_colors, s=all_sizes, cmap="viridis", alpha=0.8)
        ax.set_facecolor('white')  # Remove grey background
        ax.view_init(elev=90, azim=-90)  # Rotate to XY plane
        ax._axis3don = False  # Remove the 3D box lines
        set_axes_equal(ax)  # Equal axis scaling
        ax.set_title(f"{cell_type} - {longest_cell.split('__')[1].split('_')[0]}", fontsize=20)  # Larger font for title
        ax.title.set_position([0.5, 0.9])  # Adjust title position

plt.tight_layout()
plt.savefig("summary_plot.png", dpi=300)
plt.close()
