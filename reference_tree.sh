#!/bin/bash

# Input a TSV table from https://www.ncbi.nlm.nih.gov/datasets/genome

if [ -z "$1" ]; then
    echo Error: input file name is required
    exit 1
fi

table=$1

if [ ! -f ${table} ]; then 
    echo Error: Input file ${catalog} does not exist
    exit 1
fi

# TODO: write comments explaining the steps
# Skip the header line and save accessions from all other lines
ACCESSIONS=$(tail -n +2 ${table} | cut -f1 | tr '\n' ' ') 
mkdir accessions
cd accessions
echo $ACCESSIONS | tr ' ' '\n' > gb_ids.txt
../scripts/datasets download genome accession $ACCESSIONS --include genome,cds
unzip ncbi_dataset.zip
sed -nr 's/.*"taxId":([0-9]*).*/\1/p' ncbi_dataset/data/assembly_data_report.jsonl > ncbi_ids.txt
paste ncbi_ids.txt gb_ids.txt | while read ncbi_id gb_id; do
    lineage=$(efetch -db taxonomy -id $ncbi_id -format xml | sed -nr 's-.*<Lineage>(.*)</Lineage>-\1-p')
    echo $lineage "; " $gb_id >> species.txt
done
python ../scripts/taxonomy_to_newick.py species.txt > species.tree
cd ..
