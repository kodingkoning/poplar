#!/bin/bash

# Run this in the location you would like to download the dependencies

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
