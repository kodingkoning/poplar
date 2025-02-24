#!/bin/bash
#SBATCH --job-name=poplar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p compute
#SBATCH -t 2:05:00

source activate ~/.conda/envs/poplar_env
time python parsl/main.py ncbi_dataset/data/dataset_catalog.json -o example.tree

