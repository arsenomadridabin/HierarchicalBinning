import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import argparse
from scipy.ndimage import binary_erosion

# Parse command-line argument
parser = argparse.ArgumentParser()
parser.add_argument("--threshold", type=int, default=6, help="Fe threshold")
args = parser.parse_args()
fe_threshold = args.threshold

# Load JSON
with open("filtered_bins_hierarchical.json") as f:
    hierarchical_data = json.load(f)

# Build atom lookup
atom_lookup = dict()
for bin_info in hierarchical_data:
    tx, ty, tz = bin_info["bin_index"]
    for sub in bin_info["sub_bins"]:
        sx, sy, sz = sub["sub_bin_index"]
        gx, gy, gz = tx * 2 + sx, ty * 2 + sy, tz * 2 + sz
        atom_lookup[(gx, gy, gz)] = sub

# Create 3D mask
grid_shape = (16, 16, 16)
mask = np.zeros(grid_shape, dtype=bool)
for x, y, z in atom_lookup:
    mask[x, y, z] = True

# Morphological erosion
interior_mask = binary_erosion(mask)
boundary_mask = mask & (~interior_mask)

# Classify
strict_retained_bins = []
strict_discarded_bins = []
for coord in atom_lookup:
    x, y, z = coord
    fe = atom_lookup[coord]["fe"]
    if fe < fe_threshold and boundary_mask[x, y, z]:
        strict_discarded_bins.append((x, y, z, fe))
    else:
        strict_retained_bins.append((x, y, z, fe))

# Count atoms
retained_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}
total_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}
for coord, data in atom_lookup.items():
    for el in total_counts:
        total_counts[el] += data[el]
for x, y, z, _ in strict_retained_bins:
    data = atom_lookup[(x, y, z)]
    for el in retained_counts:
        retained_counts[el] += data[el]

# Weight percent
atomic_weights = {"fe": 55.845, "mg": 24.305, "si": 28.085, "o": 15.999}
total_mass_retained = sum(retained_counts[el] * atomic_weights[el] for el in retained_counts)
weight_percent_retained = {
    el: 100 * retained_counts[el] * atomic_weights[el] / total_mass_retained
    for el in retained_counts
}

print("Total Atom Counts:", total_counts)
print("Retained Atom Counts:", retained_counts)
print("Weight Percent in Retained:", weight_percent_retained)

# Add labels and write updated JSON
for bin_info in hierarchical_data:
    tx, ty, tz = bin_info["bin_index"]
    updated_subs = []
    for sub in bin_info["sub_bins"]:
        sx, sy, sz = sub["sub_bin_index"]
        gx, gy, gz = tx * 2 + sx, ty * 2 + sy, tz * 2 + sz
        label = ""
        fe = sub["fe"]
        if fe >= fe_threshold:
            label = "Fe-rich bin"
        elif boundary_mask[gx, gy, gz] and fe < fe_threshold:
            label = "Boundary Bin"
        sub["label"] = label
        updated_subs.append(sub)
    bin_info["sub_bins"] = updated_subs

with open("filtered_bins_labeled.json", "w") as f:
    json.dump(hierarchical_data, f, indent=2)

# Visualization
output_dir = "morph_slices_strict"
os.makedirs(output_dir, exist_ok=True)
strict_grid = np.zeros((16, 16, 16), dtype=int)
for x, y, z, _ in strict_retained_bins:
    strict_grid[x, y, z] = 1
for x, y, z, _ in strict_discarded_bins:
    strict_grid[x, y, z] = 2
color_map = {0: "white", 1: "red", 2: "gray"}

for z in range(16):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title(f'Strict Morphological Boundary (Z = {z})')
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)
    ax.set_xticks(range(0, 17, 2))
    ax.set_yticks(range(0, 17, 2))
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlabel('X Sub-bin')
    ax.set_ylabel('Y Sub-bin')

    for x in range(16):
        for y in range(16):
            code = strict_grid[x, y, z]
            if code > 0:
                color = color_map[code]
                rect = patches.Rectangle((x, y), 1, 1, facecolor=color, edgecolor='black', linewidth=0.5)
                ax.add_patch(rect)
                fe = atom_lookup.get((x, y, z), {}).get("fe", 0)
                ax.text(x + 0.5, y + 0.5, str(fe), ha='center', va='center',
                        fontsize=7, color='white' if code == 1 else 'black')

    ax.set_aspect('equal')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/zslice_{z}.png", dpi=150)
    plt.close()
