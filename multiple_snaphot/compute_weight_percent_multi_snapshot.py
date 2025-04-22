import json
import argparse
import numpy as np
from scipy.ndimage import binary_erosion

parser = argparse.ArgumentParser()
parser.add_argument("--threshold", type=int, default=6, help="Fe threshold")
parser.add_argument("--input", type=str, default="filtered_bins_hierarchical.json")
parser.add_argument("--output", type=str, default="filtered_bins_labeled.json")
args = parser.parse_args()

fe_threshold = args.threshold

with open(args.input) as f:
    loaded_data = json.load(f)

# If it's a single snapshot (list), wrap into snapshot_0 format
if isinstance(loaded_data, list):
    all_snapshots = {"snapshot_0": loaded_data}
else:
    all_snapshots = loaded_data

results = {}
labeled_snapshots = {}

for snap_key, hierarchical_data in all_snapshots.items():
    atom_lookup = dict()
    for bin_info in hierarchical_data:
        tx, ty, tz = bin_info["bin_index"]
        for sub in bin_info["sub_bins"]:
            sx, sy, sz = sub["sub_bin_index"]
            gx, gy, gz = tx * 2 + sx, ty * 2 + sy, tz * 2 + sz
            atom_lookup[(gx, gy, gz)] = sub

    grid_shape = (16, 16, 16)
    mask = np.zeros(grid_shape, dtype=bool)
    for x, y, z in atom_lookup:
        mask[x, y, z] = True

    interior_mask = binary_erosion(mask)
    boundary_mask = mask & (~interior_mask)

    strict_retained_bins = []
    strict_discarded_bins = []
    for coord in atom_lookup:
        x, y, z = coord
        fe = atom_lookup[coord]["fe"]
        if fe < fe_threshold and boundary_mask[x, y, z]:
            strict_discarded_bins.append((x, y, z, fe))
        else:
            strict_retained_bins.append((x, y, z, fe))

    retained_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}
    total_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}
    for coord, data in atom_lookup.items():
        for el in total_counts:
            total_counts[el] += data[el]
    for x, y, z, _ in strict_retained_bins:
        data = atom_lookup[(x, y, z)]
        for el in retained_counts:
            retained_counts[el] += data[el]

    atomic_weights = {"fe": 55.845, "mg": 24.305, "si": 28.085, "o": 15.999}
    retained_mass = {el: retained_counts[el] * atomic_weights[el] for el in retained_counts}
    total_mass_retained = sum(retained_mass.values())
    weight_percent_retained = {
        el: 100 * retained_mass[el] / total_mass_retained
        for el in retained_mass
    }

    results[snap_key] = {
        "Total Atom Counts": total_counts,
        "Retained Atom Counts": retained_counts,
        "Weight Percent in Retained": weight_percent_retained
    }

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
    labeled_snapshots[snap_key] = hierarchical_data

# Save output
with open(args.output, "w") as f:
    json.dump(labeled_snapshots["snapshot_0"] if "snapshot_0" in labeled_snapshots and len(labeled_snapshots) == 1 else labeled_snapshots, f, indent=2)

# Print result
for snap, stats in results.items():
    print(f"\n=== {snap} ===")
    print("Total Atom Counts:", stats["Total Atom Counts"])
    print("Retained Atom Counts:", stats["Retained Atom Counts"])
    print("Weight Percent in Retained:")
    for el, wt in stats["Weight Percent in Retained"].items():
        print(f"  {el.upper()}: {wt:.2f}%")

# --- Compute AVERAGE across all snapshots ---
total_mass_all = {"fe": 0.0, "mg": 0.0, "si": 0.0, "o": 0.0}
atomic_weights = {"fe": 55.845, "mg": 24.305, "si": 28.085, "o": 15.999}

for snap, stats in results.items():
    for el in total_mass_all:
        total_mass_all[el] += stats["Retained Atom Counts"][el] * atomic_weights[el]

grand_total_mass = sum(total_mass_all.values())
average_weight_percent = {
    el: 100 * total_mass_all[el] / grand_total_mass for el in total_mass_all
}

print("\n=== AVERAGE Weight Percent across all snapshots ===")
for el, wt in average_weight_percent.items():
    print(f"  {el.upper()}: {wt:.2f}%")
