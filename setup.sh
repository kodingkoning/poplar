#!/bin/bash

if [ -z $1 ]; then
    echo "First argument must be Slurm queue name"
    exit
fi

QUEUE=$1

sed -i "s/#SBATCH -p[^#]*/#SBATCH -p ${QUEUE}/" poplar.sbatch scripts/*.script
