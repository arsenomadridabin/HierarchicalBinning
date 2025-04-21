import json
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
import os

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

def load_coordinates(filepath):
    with open(filepath) as f:
        data = json.load(f)
    return np.array([atom["atom_coordinate"] for atom in data])

fe_coords = load_coordinates(args.fe_file)
mg_coords = load_coordinates(args.mg_file)
si_coords = load_coordinates(args.si_file)
o_coords = load_coordinates(args.o_file)


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


with open(args.output, 'w') as f:
    json.dump(filtered_bins, f, indent=4)

print(f"Filtered bins saved to: {args.output}")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for bin_info in filtered_bins:
    i, j, k = bin_info["bin_index"]
    ax.scatter(i, j, k, c='red', s=50)

ax.set_xlabel('X Bin')
ax.set_ylabel('Y Bin')
ax.set_zlabel('Z Bin')
ax.set_title('Filtered Bins by Fe Count')
plt.tight_layout()

#Save the figure
plot_path = os.path.splitext(args.output)[0] + "_plot.png"
plt.savefig(plot_path, dpi=300)
print(f"Plot saved to: {plot_path}")

