#!/bin/bash

folder_with_pdbs="../filter_radius_contacts/ns1_denv1/selected_ns1_denv1_rg0.03_ct0.3"

output_dir="outputs/ns1_denv1"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/selected_structures_ns1_denv1.jsonl"
path_for_fixed_positions=$output_dir"/fixed_positions_ns1_denv1_pdbs.jsonl"

python /home/patricia/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains && \
python make_fixed_positions.py && \
python /home/patricia/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
       	--fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 10 \
        --seed 409118984  \
