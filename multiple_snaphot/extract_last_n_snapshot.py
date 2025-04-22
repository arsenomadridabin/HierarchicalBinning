import json
import argparse

def extract_last_n(infile, n):
    with open(infile, 'r') as f:
        data = json.load(f)

    if isinstance(data, list) and len(data) >= n:
        extracted = data[-n:]
    else:
        print(f"Warning: {infile} does not contain at least {n} snapshots.")
        extracted = data  # fallback: save whatever is available

    outfile = infile.replace(".json", f"_last{n}.json")
    with open(outfile, 'w') as f:
        json.dump(extracted, f, indent=2)

    print(f"Saved last {n} snapshots of {infile} to {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True, help="Number of last snapshots to extract")
    parser.add_argument("--files", nargs="+", required=True, help="Input JSON files (e.g. fe.json mg.json)")
    args = parser.parse_args()

    for infile in args.files:
        extract_last_n(infile, args.n)

