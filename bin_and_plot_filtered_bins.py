
import json
import numpy as np
import argparse
import matplotlib.pyplot as plt
import itertools
import os
import matplotlib.patches as patches

"""Command-line Arguments"""

parser = argparse.ArgumentParser(description="Bin atoms and filter by Fe count.")
parser.add_argument("--cell_size", type=float, default=68.0, help="Size of the cell (default: 68.0)")
parser.add_argument("--num_bins", type=int, default=8, help="Number of bins per dimension (default: 8)")
parser.add_argument("--lower_count", type=int, default=0, help="Lower bound for Fe count")
parser.add_argument("--upper_count", type=int, default=6, help="Upper bound for Fe count")
parser.add_argument("--second_lower_count", type=int, default=44, help="Second lower bound for Fe count")
parser.add_argument("--second_upper_count", type=int, default=60, help="Second upper bound for Fe count")
parser.add_argument("--fe_file", type=str, default="fe_last.json", help="Path to Fe JSON file")
parser.add_argument("--mg_file", type=str, default="mg_last.json", help="Path to Mg JSON file")
parser.add_argument("--si_file", type=str, default="si_last.json", help="Path to Si JSON file")
parser.add_argument("--o_file", type=str, default="o_last.json", help="Path to O JSON file")
parser.add_argument("--output", type=str, default="filtered_bins_with_all_counts.json", help="Output JSON file")

args = parser.parse_args()

"""Load JSON Files"""

def load_coordinates(filepath):
    with open(filepath) as f:
        data = json.load(f)
    return np.array([atom["atom_coordinate"] for atom in data])

fe_coords = load_coordinates(args.fe_file)
mg_coords = load_coordinates(args.mg_file)
si_coords = load_coordinates(args.si_file)
o_coords = load_coordinates(args.o_file)

"""Binning and Counting"""

bin_edges = np.linspace(0, args.cell_size, args.num_bins + 1)

def get_bin_indices(coords):
    return np.array([np.digitize(coords[:, i], bin_edges) - 1 for i in range(3)]).T

def count_atoms(bin_indices, shape):
    count = np.zeros(shape, dtype=int)
    for x, y, z in bin_indices:
        if 0 <= x < shape[0] and 0 <= y < shape[1] and 0 <= z < shape[2]:
            count[x, y, z] += 1
    return count

shape = (args.num_bins,) * 3
fe_bins = get_bin_indices(fe_coords)
mg_bins = get_bin_indices(mg_coords)
si_bins = get_bin_indices(si_coords)
o_bins = get_bin_indices(o_coords)

fe_count = count_atoms(fe_bins, shape)
mg_count = count_atoms(mg_bins, shape)
si_count = count_atoms(si_bins, shape)
o_count = count_atoms(o_bins, shape)

"""Filter Based on Fe Count"""

fe_ranges = [(args.lower_count, args.upper_count), (args.second_lower_count, args.second_upper_count)]
filtered_bins = []

for i, j, k in itertools.product(range(args.num_bins), repeat=3):
    count = fe_count[i, j, k]
    if any(lower <= count <= upper for (lower, upper) in fe_ranges):
        filtered_bins.append({
            "bin_index": [i, j, k],
            "fe_count": int(fe_count[i, j, k]),
            "mg_count": int(mg_count[i, j, k]),
            "si_count": int(si_count[i, j, k]),
            "o_count": int(o_count[i, j, k])
        })

"""Save Output JSON"""
with open(args.output, 'w') as f:
    json.dump(filtered_bins, f, indent=4)

print(f"Filtered bins saved to: {args.output}")

"""Colored 2D Slice Plotting (Fixed Alignment)"""
z_slices = {k: [] for k in range(args.num_bins)}
for bin_info in filtered_bins:
    i, j, k = bin_info["bin_index"]
    z_slices[k].append((i, j, bin_info["fe_count"]))

output_basename = os.path.splitext(args.output)[0]

for k in range(args.num_bins):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title(f'Filtered Bins (Z = {k})')
    ax.set_xlabel('X Bin')
    ax.set_ylabel('Y Bin')
    ax.set_xlim(0, args.num_bins)
    ax.set_ylim(0, args.num_bins)
    ax.set_xticks(range(args.num_bins))
    ax.set_yticks(range(args.num_bins))
    ax.grid(True, which='both', color='gray', linestyle='--', linewidth=0.5)

    for i, j, fe_count in z_slices[k]:
        rect = patches.Rectangle((i, j), 1, 1, linewidth=0.5, edgecolor='black', facecolor='red')
        ax.add_patch(rect)
        ax.text(i + 0.5, j + 0.5, str(fe_count), color='white', ha='center', va='center', fontsize=8)

    plt.gca().invert_yaxis()
    plt.tight_layout()
    plot_path = f"{output_basename}_z{k}.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"Saved: {plot_path}")

