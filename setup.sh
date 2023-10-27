#!/bin/bash

usage() { echo "USAGE: -q to set queue name and -t to sex max time for a single job"; exit; }

while getopts ":q:t:h:" o; do
    case "${o}" in 
        q)
            QUEUE=${OPTARG}
            sed -i "s/#SBATCH -p[^#]*/#SBATCH -p ${QUEUE}/" poplar.sbatch scripts/*.script
            ;;
        t)
            TIME=${OPTARG}
            sed -i "s/#SBATCH --time[^#]*/#SBATCH --time=${TIME}/" poplar.sbatch scripts/*.script
            sed -i "s/#SBATCH -t[^#]*/#SBATCH -t ${TIME}/" poplar.sbatch scripts/*.script
            ;;
        h) usage;;
    esac
done

