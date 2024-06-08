#!/bin/bash

# Run this in the location you would like to download the dependencies

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

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xf ncbi-blast-2.15.0+-x64-linux.tar.gz
echo "Add to PATH: " $PWD/ncbi-blast-2.15.0+/bin
export PATH=$PATH:$PWD/ncbi-blast-2.15.0+/bin

wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
echo "Add to PATH: " $PWD
export PATH=$PATH:$PWD

wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
unzip Linux.zip
# echo "Add to PATH: " $PWD/Aster-Linux/bin
export PATH=$PATH:$PWD/Aster-Linux/bin

echo ""
echo "Add to PATH: " $PWD:$PWD/Aster-Linux/bin:$PATH:$PWD/ncbi-blast-2.15.0+/bin
