#!/bin/bash

folder_name="$1"
input_pdbqt="$folder_name/pipeline_files/3_pdbqt"

output_smi="$folder_name/pipeline_files/4_smiles"
mkdir -p "$output_smi"

# Convert PDBQT to SMILES in batches of 20000
batch_num=1
for i in $(seq 1 20000 $(find "$input_pdbqt" -maxdepth 1 -type f -name "*.pdbqt" | wc -l)); do
    batch_files=("$input_pdbqt"/*.pdbqt)
    obabel -ipdbqt "${batch_files[@]:i:20000}" -osmi -O "$output_smi"/.smi -m > /dev/null 2>&1
    echo "PDBQT to SMILES conversion for batch $batch_num completed."
    ((batch_num++))
done

echo -e "\033[1m\033[34mPDBQT to SMILES conversion completed and files saved in folder: \033[91m$output_smi\033[0m" >&1
