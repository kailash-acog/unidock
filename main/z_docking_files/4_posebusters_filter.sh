#!/bin/bash

folder_name="$1"
pdb_filename="$2"

input_sdf="$folder_name/pipeline_files/9_sdf_out"
posebusters_output="$folder_name/pipeline_files/posebusters_output.txt"

touch "$posebusters_output"

echo -e "\033[1m\033[34mCheck PoseBusters Progress... \033[91m$posebusters_output\033[0m"

find "$input_sdf" -name "*.sdf" -print0 | \
    xargs -0 -P "$(nproc)" -I {} bust "{}" -p "$folder_name/$pdb_filename" >> "$posebusters_output" 2>&1

echo -e "\033[1m\033[34mPoseBusters Filtration Completed\033[0m"
