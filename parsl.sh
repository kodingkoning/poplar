#!/bin/bash
#SBATCH --job-name=tree_from_annotation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p compute

python ../parsl/main.py /home/erkonin/daily_notes/data/Kickxellomycotina/ncbi_dataset/data/dataset_catalog_abridged.json test.out
