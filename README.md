# Poplar

## Purpose

Popular is a software pipeline that uses an input of genes and assembled genomes and generates a phylogenetic tree from the input.

## Dependencies

Required dependencies:

- [Python 3.10](https://www.python.org/downloads/)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [numpy](https://numpy.org/)
- [sklearn 1.3.0](https://scikit-learn.org/stable/index.html)
- [biopython](https://biopython.org/docs/1.75/api/Bio.html)

Recommended installation:

```
conda create -n poplar_env python=3.10 numpy scikit-learn biopython
conda activate poplar_env
```

Run the pipeline with the environment activated.

## Input

The pipeline accepts a JSON file matching the format of the `dataset_catalog.json` that comes with a download from NCBI's genome datasets. The JSON file specifies the relative locations of the file contianing FASTA files.

### Downloading Data

Visit [NCBI's Genome Datasets page](https://www.ncbi.nlm.nih.gov/datasets/genome) and search for the desired taxa. From there, check the boxes for the assemblies you are interested in. Click "Download" and "Download Package". Select "Genome sequences (FASTA)" and "Genomic coding sequences (FASTA)" and download your file. Unzip the file, and it will contain the relevant FASTA files and `dataset_catalog.json`.

If you have additional data to include alongside NCBI genomes, add an entry in the file for each additional species. If you are not downloading data from NCBI, use `files/dataset_catalog.json` as a template for your input file. `filePath` should be relative to the location of `dataset_catalog.json`, and file types should be `CDS_NUCLEOTIDE_FASTA` for gene files and `GENOMIC_NUCLEOTIDE_FASTA` for assembled genomes.

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

There are currently not any checkpointing implemented. If the pipeline crashes, then the preceeding steps need to be rerun rather than continue from where it was left off. (TODO: add checkpointing)

TODO: figure out why the full Kickxellomycotina set gets an error with one of the inputs (shows as "" query) -- the file in genomes.txt doesn't exist, may have been a problem with the download or something like that, but may want to include a check that a file exists and pass a warning in that case

TODO: add comments to code
