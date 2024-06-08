# Poplar

## Purpose

Popular is a software pipeline that uses an input of genes and assembled genomes and generates a phylogenetic tree from the input.

## Dependencies

Required dependencies:

- [Python 3.10](https://www.python.org/downloads/)
- [Parsl](https://parsl.readthedocs.io/en/stable)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [numpy](https://numpy.org/)
- [sklearn 1.3.0](https://scikit-learn.org/stable/index.html)
- [biopython](https://biopython.org/docs/1.75/api/Bio.html)
- [orfipy](https://pypi.org/project/orfipy/)
- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [ASTRAL-Pro](https://github.com/chaoszhang/A-pro)
- [Java >= 1.7](https://www.java.com/en/download/), for ASTRAL-Pro

Recommended installation in a conda environment:

```
conda create -n poplar_env python=3.10 numpy scikit-learn biopython parsl bioconda::orfipy mafft
conda activate poplar_env
```

Tools to download and add to PATH in `config.py`:

- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [ASTRAL-Pro](https://github.com/chaoszhang/A-pro)

All of these can be downloaded to the current director using "setup.sh". The locations will be printed for adding to config.py.

## Parsl Configuration

Parsl requires a configuration file with information required to launch the tasks. Information about the scope of options is available in Parsl's documentation, and the file in `parsl/config.py` offers a template to build on for poplar.

The key things that need to be changed when moving to a new machine are:

- conda environment activation
- PATH for tools that are not installed with conda
- partition name
- nodes per block, maximum blocks, walltime

Example `worker_init` value within config, in order to activate the conda environment and set PATH:

```
worker_init='conda activate poplar_env; export PATH=$PATH:$BLAST_PATH:$RAXML_NG_PATH:$ASTRAL_PRO_PATH'
```

Partition name must be set to the slurm partition that should be used, and the details of the slurm job (nodes per block, maximum blocks, walltime) should be chosen based on the machine and job being run.

## Setup

Before running the pipeline, set the name of the queue and the time limit for the jobs as needed for your system. `setup.sh` accepts `-q` for the name of the slurm queue and `-t` for the time limit. If you want different time limits for different pipeline stages, that must be done manually.

## Input

The pipeline accepts a JSON file matching the format of the `dataset_catalog.json` that comes with a download from NCBI's genome datasets. The JSON file specifies the relative locations of the file contianing FASTA files.

### Downloading Data

Visit [NCBI's Genome Datasets page](https://www.ncbi.nlm.nih.gov/datasets/genome) and search for the desired taxa. From there, check the boxes for the assemblies you are interested in. Click "Download" and "Download Package". Select "Genome sequences (FASTA)", "Genomic coding sequences (FASTA)", and "Annotation features (GFF)" and download your file. Unzip the file, and it will contain the relevant FASTA files and `dataset_catalog.json`.

If you have additional data to include alongside NCBI genomes, add an entry in the file for each additional species. If you are not downloading data from NCBI, use `files/dataset_catalog.json` as a template for your input file. `filePath` should be relative to the location of `dataset_catalog.json`, and file types should be `CDS_NUCLEOTIDE_FASTA` for gene files, `GENOMIC_NUCLEOTIDE_FASTA` for assembled genomes, and `GFF3` for GFF annotation files.

The order of preference for input file types is:

1. `CDS_NUCLEOTIDE_FASTA`, using the sequences as genes
2. `GFF3` and `GENOMIC_NUCLEOTIDE_FASTA`, extracting gene sequences based on GFF3 from the genomic sequences file
3. `GENOMIC_NUCLEOTIDE_FASTA`, identifying ORFs to use as gene sequences

The "accession" field in the JSON file will be used for naming intermediate files as well as the placement in the output tree, so the values must be unique for each entry. If downloaded from NCBI, they will be GenBank IDs.

### Running the Pipeline

Setup by running `./setup.sh <queue name>`. This will allow you to select the name of the Slurm queue that should be used by all jobs.

On a machine using Slurm, run `sbatch poplar.sbatch /path/to/dataset_catalog.json`.

## Output

The pipeline will create a new directory within the current directory to store all the temporary files. One output tree will be named `job_{jobnum}.tree` and be in the current directory. The labels of this tree will match the accession IDs. If `dataset_catalog.json` contains a file of type `DATA_TABLE`, then Aasecond output tree named `job_{jobnum}_scinames.tree` will replace the accession IDs with the scientific names. This may be easier to read, but it will not differentiate between multiple instances of a single species, and it will not rename any species without an entry in the data table.

## Limitations

This tool only works on Slurm machines.

If certain pre-compiled executables (such as for seqkit) cannot run on the provided architecture, there will be errors. These can be resolved by replacing the executables in their current locations.

Jobs spawn other jobs and then wait for their completition. If jobs are forced to wait in the queue beyond the time limit of the originating job, then the pipeline might crash.

## Steps

1. Gene and ORF Location

	- The first step is to extract genes and possible genes from each input organism. If coding sequences are included in the input, they are used. If not, then orfipy is used.

2. Grouping Genes and ORFs

	- All genes and ORFs are combined into a BLAST database, available for search. Following the creation of the database, the databaes is queried for each input sequence. This results in a distance between it and the most similar sequences (the match to itself then excluded).
	- Using the distances from BLAST, a distance matrix is created and used with DBSCAN to create groupings.

3. Gene Trees

	- The sequences in each group are aligned using MAFFT.
	- The sequence alignments are given to RAxML to generate gene trees.

4. Species Tree

	- The RAxML gene trees are given to ASTRAL-Pro to generate a single species tree.

## Feature Yet to be Added

- Checkpointing

### Instructions for Generating NCBI Tree

The tests in our paper compare the trees generated by this pipeline to NCBI's taxonomy. Creating that tree requires downloading data from NCBI's website and through their command line tools. The provided script `reference_tree.sh` downloads data from NCBI and creates a Newick tree with the available information on the spcies. In order to use the script, visit [NCBI's Genome Browser](https://www.ncbi.nlm.nih.gov/datasets/genome) and search for the names or IDs of the species of interest. For our tests, we selected only the reference genomes. Then, download the table. This table includes the GenBank reference numbers. The first column of the table will be GenBank IDs, which will be used for searching NCBI's database for the taxonomy.

## Planned Improvements

- Checkpointing using Parsl
- Analysis of optimal gene sequence group size. Currently is required to be between 4 and 99.
