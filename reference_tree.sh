#!/bin/bash

# Input a TSV table from https://www.ncbi.nlm.nih.gov/datasets/genome or a NCBI taxon ID

if [ -z "$1" ]; then
    echo Error: input file name is required
    exit 1
fi

mkdir accessions
cd accessions

if [ ! -f $1 ]; then 
    taxon=$1
    ../scripts/datasets download genome taxon $taxon --reference --include genome,cds
    unzip ncbi_dataset.zip
    python ../scripts/nav_jsonl.py ncbi_dataset/data/assembly_data_report.jsonl > gb_ids.txt
else
    table=$1
    ACCESSIONS=$(tail -n +2 ${table} | cut -f1 | tr '\n' ' ')
    ../scripts/datasets download genome accession $ACCESSIONS --include genome,cds
    echo $ACCESSIONS | tr ' ' '\n' > gb_ids.txt
    unzip ncbi_dataset.zip
fi

sed -nr 's/.*"taxId":([0-9]*).*/\1/p' ncbi_dataset/data/assembly_data_report.jsonl > ncbi_ids.txt
paste ncbi_ids.txt gb_ids.txt | while read ncbi_id gb_id; do
    lineage=$(efetch -db taxonomy -id $ncbi_id -format xml | sed -nr 's-.*<Lineage>(.*)</Lineage>-\1-p')
    echo $lineage "; " $gb_id >> species.txt
done
python ../scripts/taxonomy_to_newick.py species.txt > species.tree
cd ..
