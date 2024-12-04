import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# 1. Data Name Variable and Column Selection
data_name = "Complexity Peak Fitted Distance from Soma (Apical)"
file_path = 'curfitted_local_lzc_max_location.csv'  # Replace with your file path
column_to_use = 4  # Column index to use (3rd column)

# 2. Read the CSV File
data = pd.read_csv(file_path)

# Map cell types to colors
custom_colors = [
    (0.0, 1.0, 1.0),  # Cyan
    (0.0, 0.8, 0.8),  # Teal
    (0.0, 0.6, 0.6),  # Darker teal
    (0.8, 0.8, 0.0),  # Yellow
    (0.95, 0.0, 0.0),  # Red
    (0.0, 0.95, 0.0),  # Green
    (0.0, 0.0, 0.95),  # Blue
]
cell_types = data['cell_type'].unique()
color_map = {cell_type: custom_colors[i % len(custom_colors)] for i, cell_type in enumerate(cell_types)}

# 3. Prepare Data for Violin Plot, Omit Groups with No Useful Data
valid_cell_type_values = {
    cell_type: data[data['cell_type'] == cell_type].iloc[:, column_to_use].dropna().values
    for cell_type in cell_types
}
valid_cell_type_values = {k: v for k, v in valid_cell_type_values.items() if len(v) > 0}

# Update color map to include only valid cell types
valid_cell_types = list(valid_cell_type_values.keys())
valid_colors = [color_map[ct] for ct in valid_cell_types]
box_data = [values for values in valid_cell_type_values.values()]

# 4. Violin Plot with Statistical Annotations
fig, ax = plt.subplots(figsize=(8, 6))

# Violin plot
sns.violinplot(data=box_data, ax=ax, palette=valid_colors, cut=0, inner=None)

# Overlay boxplot for medians and whiskers
sns.boxplot(
    data=box_data,
    ax=ax,
    palette=["none"] * len(valid_cell_types),
    showcaps=False,
    whiskerprops={"linewidth": 2},
    boxprops={"zorder": 3, "facecolor": "none", "edgecolor": "black"},
    medianprops={"color": "black", "linewidth": 2},
    showfliers=False,
)

# Scatter individual points
for i, values in enumerate(box_data):
    ax.scatter([i] * len(values), values, color="black", alpha=0.5, s=30)

# Statistical comparisons and annotations
p_values = [
    (ct1, ct2, ranksums(valid_cell_type_values[ct1], valid_cell_type_values[ct2]).pvalue)
    for ct1, ct2 in combinations(valid_cell_types, 2)
]
adjusted_p_values = multipletests([p for _, _, p in p_values], method='fdr_bh')[1]
significant_p_values = [
    (ct1, ct2, p_val) for (ct1, ct2, _), p_val in zip(p_values, adjusted_p_values) if p_val < 0.05
]

data_max = max([max(values) for values in box_data])
for n, (ct1, ct2, p_val) in enumerate(significant_p_values):
    x1, x2 = valid_cell_types.index(ct1), valid_cell_types.index(ct2)
    y = data_max * (1.05 + 0.05 * n)
    ax.plot([x1, x2], [y, y], color='black', linewidth=1.5)
    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
    ax.text((x1 + x2) / 2, y + data_max * 0.01, f"{significance} p={p_val:.6f}", ha='center', fontsize=8)

ax.set_xticks(range(len(valid_cell_types)))
ax.set_xticklabels(valid_cell_types, fontsize=10)
ax.set_ylabel(data_name, fontsize=12)
plt.title(f"{data_name} with Rank Sum Comparison", fontsize=14)

# Extend Y-axis for significance annotations
plt.ylim(0, data_max * (1.05 + 0.05 * len(significant_p_values)))

plt.tight_layout()
plt.savefig(f'{data_name}_violin_plot_filtered.png', dpi=300, bbox_inches='tight')
plt.show()
