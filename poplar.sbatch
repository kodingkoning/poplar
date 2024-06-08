#!/bin/bash
#SBATCH --job-name=tree_from_annotation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p QUEUE

# Usage: sbatch ncbi_main.script <dataset_catalog.json> 

SECONDS=0

if [ -z "$1" ]; then
    echo Error: input file name is required
    exit 1
fi

catalog=$1

if [ ! -f ${catalog} ]; then 
    echo Error: Input file ${catalog} does not exist
    exit 1
fi

# check if script is started via SLURM or bash
# if with SLURM: there variable '$SLURM_JOB_ID' will exist
# `if [ -n $SLURM_JOB_ID ]` checks if $SLURM_JOB_ID is not an empty string
if [ -n $SLURM_JOB_ID ];  then
    # check the original location through scontrol and $SLURM_JOB_ID
    SCRIPT_PATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -d' ' -f1)
else
    # otherwise: started with bash. Get the real location.
    SCRIPT_PATH=$(realpath $0)
fi

# getting location of software_name 
SHARED_PATH=$(dirname ${SCRIPT_PATH})/scripts

if [[ "$catalog" != /* ]]; then
    catalog=$PWD/${catalog}
fi

output=$PWD/job_${SLURM_JOB_ID}.tree

CATALOG_PATH=$(dirname $catalog)

WORKING_DIR=$PWD/job_${SLURM_JOB_ID}
mkdir -p $WORKING_DIR
cd $WORKING_DIR

echo Catalog = $catalog
echo Output = $output

# Parse the catalog for available data, create genes.txt and genomes.txt
echo -n "parsing catalog... "
python ${SHARED_PATH}/parse_catalog.py $catalog
echo "done"

# Relabel genes
# iterate over files in genes.txt, taking first element of each line as a species name and the second the file of the genes
# NOTE: if this gets slow, parallelize
echo -n "Relabeling genes... "
QUERY_GENES=query_genes.query_fasta
while read name file; do
    cat ${CATALOG_PATH}/${file} | ${SHARED_PATH}/seqkit replace --f-by-name -p '.*' -r "${name}_gene{nr}" > ${name}.fasta
done < genes.txt
echo "done"

DB_JSD_ID=0
# Find ORFs in unannotated genomes
if (( ($(wc -l < genomes.txt)-1) > 0 )); then
    ORF_ID=$(sbatch --array=0-$(($(wc -l < genomes.txt)-1)) --parsable ${SHARED_PATH}/orf.script ${QUERY_GENES} ${PWD}/genomes.txt ${CATALOG_PATH})
    echo "Finding orfs in job" $ORF_ID
    DB_JSD_ID=$(sbatch -d afterok:$ORF_ID --parsable ${SHARED_PATH}/makeblastdb.script)
else
    DB_JSD_ID=$(sbatch --parsable ${SHARED_PATH}/makeblastdb.script)
fi
echo "Making BLAST database in job" $DB_JSD_ID
BLAST_JSD_ID=$(sbatch --array=0-$(($(wc -l < genes.txt)+$(wc -l < genomes.txt)-1)) -d afterok:$DB_JSD_ID --parsable ${SHARED_PATH}/blastn_tabular.script)
echo "BLAST search in job" $BLAST_JSD_ID

GROUP_JSD_ID=$(sbatch -d afterok:$BLAST_JSD_ID --parsable ${SHARED_PATH}/group.script)
echo "Grouping sequences in job" $GROUP_JSD_ID
SELECT_JSD_ID=$(sbatch -d afterok:$GROUP_JSD_ID --parsable ${SHARED_PATH}/select_seq.script)
echo "Selecting sequences in job" $SELECT_JSD_ID

# Step 4: Run MAFFT and align genes
MAFFT_JSD_ID=$(sbatch -d afterok:$SELECT_JSD_ID --parsable ${SHARED_PATH}/mafft.script)
echo "Running MAFFT in job" $MAFFT_JSD_ID

# Step 5: Create gene trees for groups
RAXML_JSD_ID=$(sbatch -d afterok:$MAFFT_JSD_ID --parsable ${SHARED_PATH}/raxml.script)
echo "Running RAxML in job" $RAXML_JSD_ID

# Step 6: Run ASTRAL-Pro for species tree
ASTRAL_JSD_ID=$(sbatch -d afterok:$RAXML_JSD_ID --parsable ${SHARED_PATH}/astral.script $output)
echo "Running ASTRAL-Pro in job" $ASTRAL_JSD_ID

echo Elapsed Time: $SECONDS seconds
