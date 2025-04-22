import json
import numpy as np
import argparse
import itertools

# ----------------------------
# Command-line Arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Hierarchical binning for multiple snapshots.")
parser.add_argument("--cell_size", type=float, default=68.0)
parser.add_argument("--num_bins", type=int, default=8)
parser.add_argument("--sub_bins", type=int, default=2)
parser.add_argument("--lower_count", type=int, default=0)
parser.add_argument("--upper_count", type=int, default=6)
parser.add_argument("--second_lower_count", type=int, default=44)
parser.add_argument("--second_upper_count", type=int, default=60)
parser.add_argument("--fe_file", type=str, required=True)
parser.add_argument("--mg_file", type=str, required=True)
parser.add_argument("--si_file", type=str, required=True)
parser.add_argument("--o_file", type=str, required=True)
parser.add_argument("--output", type=str, default="filtered_bins_hierarchical.json")

args = parser.parse_args()

# ----------------------------
# Load all snapshots
# ----------------------------
def load_coords_list(filepath):
    with open(filepath) as f:
        data = json.load(f)
    return [[atom["atom_coordinate"] for atom in snapshot] for snapshot in data]

fe_snapshots = load_coords_list(args.fe_file)
mg_snapshots = load_coords_list(args.mg_file)
si_snapshots = load_coords_list(args.si_file)
o_snapshots = load_coords_list(args.o_file)

assert len(fe_snapshots) == len(mg_snapshots) == len(si_snapshots) == len(o_snapshots)
n_snapshots = len(fe_snapshots)

# ----------------------------
# Helper Functions
# ----------------------------
bin_size = args.cell_size / args.num_bins
bin_edges = np.linspace(0, args.cell_size, args.num_bins + 1)
shape = (args.num_bins,) * 3
fe_ranges = [(args.lower_count, args.upper_count), (args.second_lower_count, args.second_upper_count)]

def get_bin_indices(coords):
    coords = np.array(coords)
    return np.array([np.digitize(coords[:, i], bin_edges) - 1 for i in range(3)]).T

def count_atoms(bin_indices, shape):
    count = np.zeros(shape, dtype=int)
    for x, y, z in bin_indices:
        if 0 <= x < shape[0] and 0 <= y < shape[1] and 0 <= z < shape[2]:
            count[x, y, z] += 1
    return count

# ----------------------------
# Main loop over snapshots
# ----------------------------
all_snapshot_results = {}

for snap_idx in range(n_snapshots):
    atom_dict = {
        "fe": np.array(fe_snapshots[snap_idx]),
        "mg": np.array(mg_snapshots[snap_idx]),
        "si": np.array(si_snapshots[snap_idx]),
        "o": np.array(o_snapshots[snap_idx])
    }

    bin_indices = {k: get_bin_indices(v) for k, v in atom_dict.items()}
    counts = {k: count_atoms(bin_indices[k], shape) for k in atom_dict}

    filtered_bins = []
    for i, j, k in itertools.product(range(args.num_bins), repeat=3):
        fe_cnt = counts["fe"][i, j, k]
        if any(lower <= fe_cnt <= upper for (lower, upper) in fe_ranges):
            filtered_bins.append((i, j, k))

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

    all_snapshot_results[f"snapshot_{snap_idx}"] = sub_bin_results

# ----------------------------
# Save Output
# ----------------------------
with open(args.output, "w") as f:
    json.dump(all_snapshot_results, f, indent=2)
print(f"Saved hierarchical bin data for {n_snapshots} snapshots to: {args.output}")

