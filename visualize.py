import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import json

# Load hierarchical bin data
with open("filtered_bins_hierarchical.json") as f:
    hierarchical_data = json.load(f)

# Create output folder
output_folder = "hierarchical_zslice_maps"
os.makedirs(output_folder, exist_ok=True)

# Dictionary to collect sub-bins by Z-slice
z_slices = {z: [] for z in range(16)}  # assuming 8 top-level bins in Z

# Organize sub-bins by Z-slice
for bin_info in hierarchical_data:
    top_x, top_y, top_z = bin_info["bin_index"]
    for sub in bin_info["sub_bins"]:
        sub_x, sub_y, sub_z = sub["sub_bin_index"]
        z_global = top_z * 2 + sub_z
        z_slices[z_global].append({
            "global_x": top_x * 2 + sub_x,
            "global_y": top_y * 2 + sub_y,
            "fe": sub["fe"]
        })

# Generate one plot per Z-slice
for z in range(16):  # 8 top bins * 2 sub-bins
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title(f'Fe-rich Sub-bins (Z = {z})')
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)
    ax.set_xticks(range(0, 17, 2))
    ax.set_yticks(range(0, 17, 2))
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlabel('X Sub-bin')
    ax.set_ylabel('Y Sub-bin')

    for item in z_slices.get(z, []):
        x, y, fe = item["global_x"], item["global_y"], item["fe"]
        color = 'red' if fe >= 7 else 'gray'
        rect = patches.Rectangle((x, y), 1, 1, edgecolor='black', facecolor=color)
        ax.add_patch(rect)
        ax.text(x + 0.5, y + 0.5, str(fe), ha='center', va='center', fontsize=6, color='white')

    ax.set_aspect('equal')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"{output_folder}/zslice_{z}.png", dpi=150)
    plt.close()


