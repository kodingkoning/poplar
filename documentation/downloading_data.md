# Downloading NCBI Data for Poplar

## Example/Test Run

Download the Pleurotus dataset via: `datasets download genome taxon 5320 --include genome,gff3,cds`, and then unzip the downloaded `ncbi_dataset.zip`.

To run, pass `dataset_catalog.json` into Poplar: `python parsl/main.py ncbi_dataset/data/dataset_catalog.json`

## File Formats

Poplar uses NCBI's dataset catalog JSON file format as input, but your data does not need to come from NCBI. If you don't want to download data from NCBI, use `files/dataset_catalog.json` as a template for your input file. `filePath` should be relative to the location of `dataset_catalog.json`, and file types should be `CDS_NUCLEOTIDE_FASTA` for gene files, `GENOMIC_NUCLEOTIDE_FASTA` for assembled genomes, and `GFF3` for GFF annotation files.

The "accession" field in the JSON file will be used for naming intermediate files as well as the placement in the output tree, so the values must be unique for each entry. If downloaded from NCBI, they will be GenBank IDs.

The order of preference for input file types is:

1. `CDS_NUCLEOTIDE_FASTA`, using the sequences as genes
2. `GFF3` and `GENOMIC_NUCLEOTIDE_FASTA`, extracting gene sequences based on GFF3 from the genomic sequences file
3. `GENOMIC_NUCLEOTIDE_FASTA`, identifying ORFs to use as gene sequences

### Command-line Download

Requires [NCBI Command-line Tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/), which can be installed with conda via `conda install -c conda-forge ncbi-datasets-cli`

Use [NCBI's Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html) to find the ID for your desired organisms.

If, for example, you are downloading Pleurotus data, the ID is `5320`, and the download command to download all three types of data usable by Poplar would be: `datasets download genome taxon 5320 --include genome,gff3,cds`

The tool downloads a zip file, and after unzipping the JSON file should be in `ncbi_dataset/data/dataset_catalog.json`.

### Website Download

Visit [NCBI's Genome Datasets page](https://www.ncbi.nlm.nih.gov/datasets/genome) and search for the desired taxa. From there, check the boxes for the assemblies you are interested in. Click "Download" and "Download Package". Select "Genome sequences (FASTA)", "Genomic coding sequences (FASTA)", and "Annotation features (GFF)" and download your file. Unzip the file, and it will contain the relevant data and `dataset_catalog.json`.

