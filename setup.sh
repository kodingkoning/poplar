#!/bin/bash

# Run this in the location you would like to download the dependencies

# Download into "dependencies" directory
mkdir -p dependencies
cd dependencies

#conda env create -f ../poplar_env.yml
# OR create with
conda config --add channels conda-forge
conda create -n poplar_env python=3.10 numpy=2.1.1 scikit-learn=1.5.1 biopython=1.84 parsl=2024.9.2 bioconda::orfipy=0.0.4 mafft=7.526 ncbi-datasets-cli=16.27.2
conda activate poplar_env
pip install 'parsl[monitoring]'
conda install flask=3.0.3 conda-forge::flask-sqlalchemy=3.1.1 pandas=2.2.2 plotly=5.24.1 networkx=3.3 pydot=3.0.1 # for monitoring

EXTEND_PATH=""

# Find the most recent version at: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -xf ncbi-blast*
cd ncbi-blast*/bin
echo "Add to PATH: " $PWD
export PATH=$PATH:$PWD
EXTEND_PATH=$EXTEND_PATH:$PWD
cd ../..

wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
unzip raxml-ng_v1.2.2_linux_x86_64.zip
echo "Add to PATH: " $PWD
export PATH=$PATH:$PWD
EXTEND_PATH=$EXTEND_PATH:$PWD

# ASTRAL can be downloaded prebuilt or built from the it repo. Prebuilt requires GLIBC 2.34
git clone https://github.com/chaoszhang/ASTER.git
cd ASTER
make
export PATH=$PATH:$PWD/bin
EXTEND_PATH=$EXTEND_PATH:$PWD/bin
cd ..
#wget https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip
#unzip Linux.zip
#export PATH=$PATH:$PWD/ASTER-Linux/bin

# Automatically updates worker_init in config.py
sed -i -E "s%worker_init\s*=\s*'([^']*)'%worker_init='conda activate poplar_env; export PATH=\$PATH:$EXTEND_PATH'%gm;t" ../parsl/config.py
