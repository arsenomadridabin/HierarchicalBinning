
import json
import argparse

parser = argparse.ArgumentParser(description="Count atoms in Fe-rich sub-bins from hierarchical binning.")
parser.add_argument("--file", required=True, help="Path to the filtered_bins_hierarchical.json file")
parser.add_argument("--fe-rich", type=int, default=7, help="Fe threshold to define Fe-rich sub-bins")

args = parser.parse_args()

with open(args.file, "r") as f:
    data = json.load(f)

fe_rich_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}
total_counts = {"fe": 0, "mg": 0, "si": 0, "o": 0}

for bin_info in data:
    for sub_bin in bin_info["sub_bins"]:
        total_counts["fe"] += sub_bin["fe"]
        total_counts["mg"] += sub_bin["mg"]
        total_counts["si"] += sub_bin["si"]
        total_counts["o"] += sub_bin["o"]

        if sub_bin["fe"] >= args.fe_rich:
            fe_rich_counts["fe"] += sub_bin["fe"]
            fe_rich_counts["mg"] += sub_bin["mg"]
            fe_rich_counts["si"] += sub_bin["si"]
            fe_rich_counts["o"] += sub_bin["o"]

print(f"Fe-rich threshold: {args.fe_rich}")
print("\nTotal atoms in Fe-rich sub-bins:")
for element, count in fe_rich_counts.items():
    print(f"  {element.upper()}: {count}")

print("\nTotal atoms in all sub-bins:")
for element, count in total_counts.items():
    print(f"  {element.upper()}: {count}")

