import json

input_files = ["fe.json", "mg.json", "si.json", "o.json"]

for infile in input_files:
    with open(infile, 'r') as f:
        data = json.load(f)

    if isinstance(data, list) and len(data) > 0:
        last_entry = data[-2]
    else:
        print(f"Warning: {infile} is empty or not a list.")
        continue

    outfile = infile.replace(".json", "_last.json")
    with open(outfile, 'w') as f:
        json.dump(last_entry, f, indent=2)

    print(f"Saved last element of {infile} to {outfile}")

