import numpy as np

def unwrap_and_center_by_fe_com(input_file="last_snapshot.dump",
                                output_file="com_centered_unwrapped.dump",
                                box_size=17.0):
    with open(input_file, "r") as f:
        lines = f.readlines()

    num_atoms = int(lines[3].strip())
    header = lines[:9]
    atom_lines = lines[9:9 + num_atoms]
    atom_data = np.genfromtxt(atom_lines)

    ids = atom_data[:, 0].astype(int)
    types = atom_data[:, 1].astype(int)
    positions = atom_data[:, 2:5]

    # Identify Fe atoms
    fe_positions = positions[types == 1]

    # Unwrap Fe atoms based on first atom
    unwrapped = fe_positions.copy()
    ref = unwrapped[0]
    for i in range(1, len(unwrapped)):
        delta = unwrapped[i] - ref
        delta -= box_size * np.round(delta / box_size)
        unwrapped[i] = ref + delta

    # Compute COM in unwrapped space
    fe_com_unwrapped = unwrapped.mean(axis=0)

    # Desired COM
    box_center = np.array([box_size / 2.0] * 3)
    shift = box_center - fe_com_unwrapped

    # Apply shift to all atoms
    shifted_positions = (positions + shift) % box_size

    # Write shifted snapshot
    with open(output_file, "w") as out:
        out.writelines(header)
        for i in range(num_atoms):
            out.write(f"{ids[i]:.0f} {types[i]:.0f} " +
                      f"{shifted_positions[i,0]:.6f} {shifted_positions[i,1]:.6f} {shifted_positions[i,2]:.6f}\n")

    print(f"Fe-cluster centered snapshot written to: {output_file}")

unwrap_and_center_by_fe_com()
