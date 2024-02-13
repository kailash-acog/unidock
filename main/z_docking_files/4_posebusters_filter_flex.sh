#!/bin/bash

folder_name="$1"
pdb_filename="$2"

input_sdf="$folder_name/pipeline_files/5_sdf_out"
posebusters_output="$folder_name/pipeline_files/2_pb_out.csv"

touch "$posebusters_output"

echo -e "\033[1m\033[34mCheck PoseBusters Progress... \033[91m$posebusters_output\033[0m"

find "$input_sdf" -name "*.sdf" -print0 | \
    xargs -0 -P "$(nproc)" -I {} bust "{}" -p "$pdb_filename" --outfmt csv >> "$posebusters_output" 

echo -e "\033[1m\033[34mPoseBusters Filtration Completed\033[0m"
