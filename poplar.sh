#!/bin/bash
#SBATCH --job-name=poplar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p compute
#SBATCH -t 2:05:00

source activate ~/.conda/envs/poplar_env
# time python parsl/main.py /home/erkonin/feeds/data/Rhodophyta/dataset_catalog.json -o trees/rhodophyla.tree4 -t 200 -f True
# time python parsl/main.py /home/erkonin/feeds/data/common_ncbi/ncbi_dataset/data/dataset_catalog.json -o trees/common.tree2 -t 200 -f True
# time python parsl/main.py /home/erkonin/feeds/data/Kickxellomycotina/dataset_catalog.json -o kick_perf_t100_4.tree -t 100
#time python parsl/main.py /home/erkonin/feeds/data/Mortierellomycotina/dataset_catalog.json -o mort_t1000.tree4 -t 1000
# time python parsl/main.py ncbi_dataset/data/dataset_catalog.json -o out.tree
# time python parsl/main.py /home/erkonin/feeds/data/Kickxellomycotina/gff3.json -o trees/kick_annotation.tree4
# time python parsl/main.py /home/erkonin/feeds/data/Kickxellomycotina/cds.json -o trees/kick_cds.tree4
time python parsl/main.py /home/erkonin/feeds/data/pucciniomycotina/dataset_cds.json -t 500 -o puccini_cds_500.tree0
exit
time python parsl/main.py /home/erkonin/feeds/data/pucciniomycotina/dataset_cds.json -t 500 -o puccini_cds_500.tree1
time python parsl/main.py /home/erkonin/feeds/data/pucciniomycotina/dataset_cds.json -t 500 -o puccini_cds_500.tree2
time python parsl/main.py /home/erkonin/feeds/data/pucciniomycotina/dataset_cds.json -t 500 -o puccini_cds_500.tree3
time python parsl/main.py /home/erkonin/feeds/data/pucciniomycotina/dataset_cds.json -t 500 -o puccini_cds_500.tree4
