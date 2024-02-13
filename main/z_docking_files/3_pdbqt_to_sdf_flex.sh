#!/bin/bash

folder_name="$1"
input_pdbqt="$folder_name/pipeline_files/3_pdbqt_out_flex_m1"

output_mol2="$folder_name/pipeline_files/4_mol2_out"
output_sdf="$folder_name/pipeline_files/5_sdf_out"

mkdir -p "$output_mol2"
mkdir -p "$output_sdf"

# Convert PDBQT to MOL2 in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$input_pdbqt" -maxdepth 1 -type f -name "*.pdbqt" | wc -l)); do
    batch_files=("$input_pdbqt"/*.pdbqt)
    obabel "${batch_files[@]:i:20000}" -omol2 -O "$output_mol2"/.mol2 -m > /dev/null 2>&1
    ((batch_num++))
done

# Convert MOL2 to SDF in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$output_mol2" -maxdepth 1 -type f -name "*.mol2" | wc -l)); do
    batch_files=("$output_mol2"/*.mol2)
    obabel "${batch_files[@]:i:20000}" -Osdf -O "$output_sdf"/.sdf -m > /dev/null 2>&1
    echo "PDBQT to SDF conversion for batch $batch_num is completed."
    ((batch_num++))
done

echo -e "\033[1m\033[34mPDBQT to SDF conversion completed and files saved in folder: \033[91m$output_sdf\033[0m" >&1
