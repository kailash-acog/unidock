#!/bin/bash

folder_name="$1"
input_sdf="$folder_name/pipeline_files/1_sdf"

output_mol2="$folder_name/pipeline_files/2_mol2"
output_pdbqt="$folder_name/pipeline_files/3_pdbqt"

mkdir -p "$output_mol2"
mkdir -p "$output_pdbqt"


# Minimize structures in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$input_sdf" -maxdepth 1 -type f -name "*.sdf" | wc -l)); do
    batch_files=("$input_sdf"/*.sdf)
    for file in "${batch_files[@]:i:20000}"; do
        obminimize -ff MMFF94 -n 1000 "$file" > /dev/null 2>&1
    done
    ((batch_num++))
done

# Convert to MOL2 in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$input_sdf" -maxdepth 1 -type f -name "*.sdf" | wc -l)); do
    batch_files=("$input_sdf"/*.sdf)
    obabel "${batch_files[@]:i:20000}" -omol2 -O "$output_mol2"/.mol2 -m > /dev/null 2>&1
    ((batch_num++))
done

# Convert to PDBQT in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$output_mol2" -maxdepth 1 -type f -name "*.mol2" | wc -l)); do
    batch_files=("$output_mol2"/*.mol2)
    obabel "${batch_files[@]:i:20000}" -opdbqt -O "$output_pdbqt"/.pdbqt -m > /dev/null 2>&1
    echo "SDF to PDBQT File conversion for batch $batch_num completed."
    ((batch_num++))
done

echo -e "\033[1m\033[34mSDF to PDBQT conversion completed and files saved in folder: \033[91m$output_pdbqt\033[0m" >&1
