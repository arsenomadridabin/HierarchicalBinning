# Get last snapshot of out.dump
python get_last_snapshot_from_out_dump.py

# Appy shift to center Fe cluster in middle
python apply_pbc_centered_fe.py

# Run Hierarchical Binning first
python hierarchical_binning.py   --fe_file fe_last.json   --mg_file mg_last.json   --si_file si_last.json   --o_file o_last.json   --cell_size 68   --num_bins 8   --sub_bins 2   --lower_count 0   --upper_count -1   --second_lower_count 40   --second_upper_count 60   --output filtered_bins_hierarchical.json

# Atom Count after Hierarchical binning
This is obsolete now
python count_atoms_hierarchial.py --file filtered_bins_hierarchical.json --fe-rich 5


#Visualize after Hierarchical binning
python visualize.py --threshold 7
we use morphological erosion to shrink the solid object inward — removing the outermost layer of voxels.
It’s the same method used in 3D image segmentation and medical imaging
Imagine our 3D voxel shape as a clay model.
Erosion is like shaving off a thin layer of clay all around — the bits that fall off? That’s our boundary layer.


# Plot filtered_bins
python bin_and_plot_filtered_bins.py   --fe_file fe_last.json   --mg_file mg_last.json   --si_file si_last.json   --o_file o_last.json   --num_bins 8   --lower_count 0   --upper_count -1   --second_lower_count 40   --second_upper_count 60   --output filtered_bins_with_all_counts.json

# Multi-snapshot analysis (Averaging)

python hierarchical_binning_multi.py   --fe_file fe_last20.json   --mg_file mg_last20.json   --si_file si_last20.json   --o_file o_last20.json   --cell_size 68   --num_bins 8   --sub_bins 2   --lower_count 0   --upper_count -1   --second_lower_count 40   --second_upper_count 60   --output filtered_bins_hierarchical.json

 python compute_weight_percent_multi_snapshot.py  --input filtered_bins_hierarchical.json   --threshold 6   --output filtered_bins_labeled.json
