import json
import numpy as np
import argparse
import matplotlib.pyplot as plt
import itertools
import os
import matplotlib.patches as patches

# Command-line Arguments
parser = argparse.ArgumentParser(description="Hierarchical binning of filtered Fe bins.")
parser.add_argument("--cell_size", type=float, default=68.0, help="Size of the simulation cell")
parser.add_argument("--num_bins", type=int, default=8, help="Top-level bins per dimension")
parser.add_argument("--sub_bins", type=int, default=2, help="Sub-bins per dimension inside filtered bin")
parser.add_argument("--lower_count", type=int, default=0, help="Fe lower count threshold")
parser.add_argument("--upper_count", type=int, default=6, help="Fe upper count threshold")
parser.add_argument("--second_lower_count", type=int, default=44, help="Second Fe lower count")
parser.add_argument("--second_upper_count", type=int, default=60, help="Second Fe upper count")
parser.add_argument("--fe_file", type=str, required=True)
parser.add_argument("--mg_file", type=str, required=True)
parser.add_argument("--si_file", type=str, required=True)
parser.add_argument("--o_file", type=str, required=True)
parser.add_argument("--output", type=str, default="filtered_bins_hierarchical.json")

args = parser.parse_args()

# Load coordinates
def load_coordinates(filepath):
    with open(filepath) as f:
        data = json.load(f)
    return np.array([atom["atom_coordinate"] for atom in data])

fe_coords = load_coordinates(args.fe_file)
mg_coords = load_coordinates(args.mg_file)
si_coords = load_coordinates(args.si_file)
o_coords = load_coordinates(args.o_file)

atom_dict = {
    "fe": fe_coords,
    "mg": mg_coords,
    "si": si_coords,
    "o": o_coords
}

# Binning
bin_size = args.cell_size / args.num_bins
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
bin_indices = {k: get_bin_indices(v) for k, v in atom_dict.items()}
counts = {k: count_atoms(bin_indices[k], shape) for k in atom_dict}

#Filter top-level bins
fe_ranges = [(args.lower_count, args.upper_count), (args.second_lower_count, args.second_upper_count)]
filtered_bins = []

for i, j, k in itertools.product(range(args.num_bins), repeat=3):
    fe_cnt = counts["fe"][i, j, k]
    if any(lower <= fe_cnt <= upper for (lower, upper) in fe_ranges):
        filtered_bins.append((i, j, k))

#Hierarchical Binning
sub_bin_results = []
sub_bin_size = bin_size / args.sub_bins

for i, j, k in filtered_bins:
    bin_result = {
        "bin_index": [i, j, k],
        "fe_count": int(counts["fe"][i, j, k]),
        "mg_count": int(counts["mg"][i, j, k]),
        "si_count": int(counts["si"][i, j, k]),
        "o_count": int(counts["o"][i, j, k]),
        "sub_bins": []
    }

    # bounding box of top bin
    xmin, xmax = i * bin_size, (i + 1) * bin_size
    ymin, ymax = j * bin_size, (j + 1) * bin_size
    zmin, zmax = k * bin_size, (k + 1) * bin_size

    def get_sub_bin(coord):
        x_sub = int((coord[0] - xmin) // sub_bin_size)
        y_sub = int((coord[1] - ymin) // sub_bin_size)
        z_sub = int((coord[2] - zmin) // sub_bin_size)
        return (x_sub, y_sub, z_sub)

    sub_counts = {}

    for atom_type, coords in atom_dict.items():
        for coord in coords:
            if xmin <= coord[0] < xmax and ymin <= coord[1] < ymax and zmin <= coord[2] < zmax:
                sub_index = get_sub_bin(coord)
                if sub_index not in sub_counts:
                    sub_counts[sub_index] = {"fe": 0, "mg": 0, "si": 0, "o": 0}
                sub_counts[sub_index][atom_type] += 1

    for sub_index in sorted(sub_counts.keys()):
        sub_data = {"sub_bin_index": list(sub_index)}
        sub_data.update(sub_counts[sub_index])
        bin_result["sub_bins"].append(sub_data)

    sub_bin_results.append(bin_result)

#Save JSON
with open(args.output, "w") as f:
    json.dump(sub_bin_results, f, indent=4)
print(f"Saved hierarchical bin data to: {args.output}")

"""
Plot slices of sub-bins
output_basename = os.path.splitext(args.output)[0]
for bin_data in sub_bin_results:
    i, j, k = bin_data["bin_index"]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_title(f'Sub-bins in Bin {i},{j},{k} (Z-slice = 0)')
    ax.set_xlim(0, args.sub_bins)
    ax.set_ylim(0, args.sub_bins)
    ax.set_xticks(range(args.sub_bins))
    ax.set_yticks(range(args.sub_bins))
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlabel('X Sub-bin')
    ax.set_ylabel('Y Sub-bin')

    for sub in bin_data["sub_bins"]:
        sx, sy, sz = sub["sub_bin_index"]
        if sz == 0:  # just slice Z=0 for now
            rect = patches.Rectangle((sx, sy), 1, 1, linewidth=0.5, edgecolor='black', facecolor='red')
            ax.add_patch(rect)
            ax.text(sx + 0.5, sy + 0.5, str(sub["fe"]), color='white', ha='center', va='center', fontsize=8)

    plt.gca().invert_yaxis()
    plt.tight_layout()
    
    os.makedirs("hierarchical_images", exist_ok=True)
    plt.savefig(f"hierarchical_images/{output_basename}_bin_{i}_{j}_{k}_z0.png", dpi=200)

    plt.close()
    print(f"Plot saved for bin ({i},{j},{k})")
"""
