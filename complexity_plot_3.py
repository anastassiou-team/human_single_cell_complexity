import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from itertools import combinations

# 1. Data Name Variable
# data_name = "Soma Spike Rate"
# data_name = "Soma LZ Complexity"
# data_name = "Max LZ Complexity"
data_name = "Whole Cell Complexity per Second"


# 2. Read the CSV File
# file_path = 'soma_fire_rate.csv'  # Replace with your file path
# file_path = 'lzc_soma_hp.csv'  # Replace with your file path
# file_path = 'lzc_max_hp.csv'  # Replace with your file path
file_path = 'global_lzc_hp.csv'  # Replace with your file path

data = pd.read_csv(file_path)

# Define custom colors for cell types
custom_colors = [
    (0.0, 1.0, 1.0),  # Cyan
    (0.0, 0.8, 0.8),  # Teal
    (0.0, 0.6, 0.6),  # Darker teal
    (0.8, 0.8, 0.0),  # Yellow
    (0.95, 0.0, 0.0),  # Red
    (0.0, 0.95, 0.0),  # Green
    (0.0, 0.0, 0.95),  # Blue
]

# Map cell types to colors
cell_types = data['cell_type'].unique()
color_map = {cell_type: custom_colors[i % len(custom_colors)] for i, cell_type in enumerate(cell_types)}



# Scale input frequency by 1/1000
x_values = data.columns[2:].astype(float) / 1000

# Scalling data only for whole cell complexity
# 
if data_name == "Whole Cell Complexity per Second":
    data.iloc[:, 2:] *= 20


# Define the data for each cell type
y_min, y_max = float("inf"), float("-inf")

# 3. Subplots for Each Cell Type
fig, axes = plt.subplots(2, 4, figsize=(12, 6))  # 4x2 layout
axes = axes.flatten()

# Calculate global y-limits
for cell_type in cell_types:
    cell_data = data[data['cell_type'] == cell_type].iloc[:, 2:]
    y_min = min(y_min, cell_data.min().min())
    y_max = max(y_max, cell_data.max().max())

for i, cell_type in enumerate(cell_types):
    ax = axes[i]
    cell_data = data[data['cell_type'] == cell_type].iloc[:, 2:]
    for _, row in cell_data.iterrows():
        ax.plot(x_values, row, color=color_map[cell_type], alpha=0.7, linewidth=1.5)
    ax.set_title(cell_type, fontsize=10)
    ax.set_xlim(x_values.min(), x_values.max())
    ax.set_ylim(y_min, y_max)
    ax.tick_params(axis="both", labelsize=8)
    ax.set_xlabel("Input Frequency (Hz)", fontsize=8)
    ax.set_ylabel(data_name, fontsize=8)
    ax.grid(False)

# Legend in the 8th subplot
legend_ax = axes[-1]
legend_ax.axis('off')
for cell_type, color in color_map.items():
    legend_ax.plot([], [], color=color, label=cell_type)
legend_ax.legend(loc='center', fontsize=10, title='Cell Type Colors')

plt.tight_layout()
plt.savefig(f'{data_name}_subplots_square.png', dpi=300, bbox_inches='tight')
plt.show()


# Whisker Bar Chart with Adjusted P-values
fig, ax = plt.subplots(figsize=(8, 6))
cell_type_max_values = {
    cell_type: data[data['cell_type'] == cell_type].iloc[:, 2:].max(axis=1).dropna().values
    for cell_type in cell_types
}

box_data = [values for values in cell_type_max_values.values()]
sns.boxplot(data=box_data, ax=ax, palette=[color_map[ct] for ct in cell_types])

# Scatter individual points
for i, values in enumerate(box_data):
    ax.scatter([i] * len(values), values, color="black", alpha=0.5, s=30)

# Calculate p-values and adjust for multiple comparisons
p_values = [
    (ct1, ct2, ranksums(cell_type_max_values[ct1], cell_type_max_values[ct2]).pvalue)
    for ct1, ct2 in combinations(cell_types, 2)
]
adjusted_p_values = multipletests([p for _, _, p in p_values], method='fdr_bh')[1]
adjusted_p_value_results = [
    (ct1, ct2, p_val) for (ct1, ct2, _), p_val in zip(p_values, adjusted_p_values)
]

# Annotate significant adjusted p-values
data_max = max([max(values) for values in box_data])  # Current data max
significant_adjusted_p_values = [
    (ct1, ct2, p_val) for (ct1, ct2, p_val) in adjusted_p_value_results if p_val < 0.05
]
for n, (ct1, ct2, p_val) in enumerate(significant_adjusted_p_values):
    x1, x2 = cell_types.tolist().index(ct1), cell_types.tolist().index(ct2)
    y = data_max * (1.0 + 0.05 * (n + 1))  # Dynamic height for each pair of comparisons
    h = 0.000 * data_max  # Minimal spacing for text and bar
    col = "k"
    ax.plot([x1, x2], [y, y], lw=1.5, color=col)  # Horizontal line for significance
    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
    ax.text(
        x=(x1 + x2) * 0.5,
        y=y + h,
        s=f"{significance} p = {p_val:.6f}",  # Stars and text together
        ha="center",
        va="bottom",
        color=col,
        fontsize=10,
    )

ax.set_xticks(range(len(cell_types)))
ax.set_xticklabels(cell_types, fontsize=10)
ax.set_ylabel(data_name, fontsize=12)
plt.title(f"Max {data_name} with Rank Sum Comparison Significance", fontsize=14)
plt.ylim(0, data_max * (1.05 + 0.05 * len(significant_adjusted_p_values)))  # Extend y-axis for annotations
plt.savefig(f'{data_name}_whisker_chart_pval_stars_text.png', dpi=300, bbox_inches='tight')
plt.show()

# Violin Plot with Scatter and Adjusted P-value Annotations
fig, ax = plt.subplots(figsize=(8, 6))

# Violin plot
sns.violinplot(
    data=box_data,
    ax=ax,
    palette=[color_map[ct] for ct in cell_types],
    cut=0,
    inner=None
)

# Overlay boxplot for median and whiskers
sns.boxplot(
    data=box_data,
    ax=ax,
    palette=["none"] * len(cell_types),  # Transparent boxes
    showcaps=False,
    whiskerprops={"linewidth": 2},
    boxprops={"zorder": 3, "facecolor": "none", "edgecolor": "black"},
    medianprops={"color": "black", "linewidth": 2},
    showfliers=False
)

# Scatter individual points
for i, values in enumerate(box_data):
    ax.scatter(
        [i] * len(values),  # Align scatter with violins
        values,
        color="black",
        alpha=0.5,
        s=30
    )

# Annotate significant adjusted p-values
for n, (ct1, ct2, p_val) in enumerate(significant_adjusted_p_values):
    x1, x2 = cell_types.tolist().index(ct1), cell_types.tolist().index(ct2)
    y = data_max * (1.0 + 0.05 * (n + 1))  # Dynamic height for each pair of comparisons
    h = 0.000 * data_max  # Minimal spacing for text and bar
    col = "k"
    ax.plot([x1, x2], [y, y], lw=1.5, color=col)  # Horizontal line for significance
    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
    ax.text(
        x=(x1 + x2) * 0.5,
        y=y + h,
        s=f"{significance} p = {p_val:.6f}",  # Stars and text together
        ha="center",
        va="bottom",
        color=col,
        fontsize=10,
    )

ax.set_xticks(range(len(cell_types)))
ax.set_xticklabels(cell_types, fontsize=10)
ax.set_ylabel(data_name, fontsize=12)
plt.title(f"{data_name} with Rank Sum Comparison Significance", fontsize=14)

# Adjust Y-axis for significance annotations
plt.ylim(0, data_max * (1.05 + 0.05 * len(significant_adjusted_p_values)))
plt.tight_layout()
plt.savefig(f'{data_name}_violin_plot_pval_stars_text.png', dpi=300, bbox_inches='tight')
plt.show()
