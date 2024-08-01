import parsl
from parsl import python_app, bash_app, join_app
import argparse
import os
import glob
from config import config

parsl.load(config)

def query_genes_func():
	pass

# Throw error if file_path does not exist
def check_file(file_path):
	if not os.path.isfile(file_path):
		raise argparse.ArgumentTypeError(f"File {file_path} does not exist")

# Get annotation sequences
# iterate over file pairs in annotations.txt
@bash_app
def annotations(line: str, CATALOG_PATH: str, SHARED_PATH: str, WORKING_DIR: str, stdout=parsl.AUTO_LOGNAME):
    import os
    os.chdir(WORKING_DIR)
    name, gff_file, fasta_file = line.split()
    return f'''{SHARED_PATH}/bedtools.static getfasta -fi {fasta_file} -bed {gff_file} -name > {fasta_file}.genes && echo {fasta_file}.genes > genes.txt'''

# Relabel genes
# iterate over files in genes.txt, taking first element of each line as a species name and the second the file of the genes
@bash_app
def relabel_genes(line: str, CATALOG_PATH: str, SHARED_PATH: str, WORKING_DIR: str, stdout=parsl.AUTO_LOGNAME):
    import os
    os.chdir(WORKING_DIR)
    name, file = line.split()
    return f'''cat {CATALOG_PATH}/{file} | {SHARED_PATH}/seqkit replace --f-by-name -p '.*' -r "{name}_gene{{nr}}" > {name}.fasta && echo {name}.fasta'''

@python_app
def parse_catalog_func(filename, WORKING_DIR):
    import os
    import json
    os.chdir(WORKING_DIR)
    with open(filename) as f:
        data = json.load(f)
    assemblies = data["assemblies"]
    with open("genes.txt","w") as genes_fout:
        with open("genomes.txt","w") as genomes_fout:
            with open("annotations.txt","w") as annotations_fout:
                for assembly in assemblies:
                    if "accession" in assembly:
                        files = assembly["files"]
                        cds_found = False
                        for file in files:
                            if "fileType" in file and file["fileType"] == "CDS_NUCLEOTIDE_FASTA":
                                genes_fout.write(assembly["accession"] + " " + file["filePath"] + "\n")
                                cds_found = True
                        if not cds_found:
                            for file in files:
                                if "fileType" in file and file["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA":
                                    gff_found = False
                                    for file2 in files:
                                        if "fileType" in file2 and file["fileType"] == "GFF3":
                                            annotations_fout.write(assembly["accession"] + " " + file2["filePath"] + " " + file["filePath"] + "\n")
                                            gff_found = True
                                    if not gff_found:
                                        genomes_fout.write(assembly["accession"] + " " + file["filePath"] + "\n")
    return True

@bash_app
def find_orfs(line: str, CATALOG_PATH: str, SHARED_PATH: str, WORKING_DIR: str, stdout=parsl.AUTO_LOGNAME):
    import os
    os.chdir(WORKING_DIR)
    name, file = line.split()
    subject = f"{CATALOG_PATH}/{file}"
    output = f"blast_{name}.fasta"
    orf_output = f"orf_{name}.fasta"
    PWD = os.getcwd()
    # install orfipy with: conda install bioconda::orfipy                                                                                                                                    
    return f'''orfipy --dna {orf_output} --ignore-case --outdir {PWD} --min 1000 {subject} && cat {orf_output} | {SHARED_PATH}/seqkit replace --f-by-name -p '.*' -r "{name}_gene{{nr}}" > {output} && rm {orf_output} && echo {output}'''

@bash_app
def build_blast_db(WORKING_DIR: str):
	return f'''cd {WORKING_DIR} && rm -f genes.db* && cat *.fasta > all_species.all_fasta && makeblastdb -input_type fasta -in all_species.all_fasta -blastdb_version 5 -title "grouping genes" -dbtype nucl -out genes.db'''

@bash_app
def search_blast(line: str, WORKING_DIR: str, blast_evalue: str):
    import os
    os.chdir(WORKING_DIR)
    line = line.strip()
    name = line.split('.')[0]
    filename = f"blast_{name}.fasta"
    return f'''blastn -query {line} -db genes.db -out blast_results/{filename}.tab -task dc-megablast -evalue {blast_evalue} -outfmt "6 qaccver saccver evalue" -num_threads 2 && echo 0'''

@bash_app
def copy_blast_to_csv(WORKING_DIR: str):
    return f'''cd {WORKING_DIR} && cat blast_results/*.tab > blast_results.csv'''

@python_app
def group(WORKING_DIR, max_group_size):
    import os
    os.chdir(WORKING_DIR)
    import csv
    import numpy as np
    from sklearn.cluster import DBSCAN
    from scipy.sparse import lil_array
    from scipy.sparse import csr_array
    from sklearn.cluster import SpectralClustering
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import sort_graph_by_row_values
    from collections import defaultdict

    # Define the path to your BLAST results file
    blast_results_file = 'blast_results.csv'  # Replace with your file path

    # Define the minimum number of clusters as 1/N (N specified via command line)
    # Replace N with the appropriate value or use a command-line argument parser
    N = 100

    # Initialize dictionaries to map sequence IDs to matrix indices
    sequence_to_index = {}
    index_to_sequence = {}

    # Create a dictionary of dictionaries to store the similarity scores as a sparse matrix
    similarity_scores = {}

    # Read the BLAST results file and populate the dictionaries and similarity_scores
    with open(blast_results_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')  # Adjust delimiter if needed
        for row in reader:
            query_sequence = row[0]  # Replace with the column index of the query sequence
            subject_sequence = row[1]  # Replace with the column index of the subject sequence
            similarity_score = float(row[2])  # Replace with the column index of the similarity score

            # Split query_sequence and subject_sequence on ORF or gene
            # Only store distance if the sequences are on two different species
            if query_sequence.partition("_ORF")[0].partition("_gene")[0] != subject_sequence.partition("_ORF")[0].partition("_gene")[0]:

                # Check if the query sequence is in the dictionary; if not, assign an index
                if query_sequence not in sequence_to_index:
                    index = len(sequence_to_index)
                    sequence_to_index[query_sequence] = index
                    index_to_sequence[index] = query_sequence

                # Check if the subject sequence is in the dictionary; if not, assign an index
                if subject_sequence not in sequence_to_index:
                    index = len(sequence_to_index)
                    sequence_to_index[subject_sequence] = index
                    index_to_sequence[index] = subject_sequence

                # Store the similarity score in the sparse matrix
                query_index = sequence_to_index[query_sequence]
                subject_index = sequence_to_index[subject_sequence]
                similarity_scores[(query_index, subject_index)] = similarity_score
                similarity_scores[(subject_index, query_index)] = similarity_score  # Symmetric matrix

    # Determine the number of sequences
    num_sequences = len(sequence_to_index)
    print("Sequences: " + str(num_sequences))

    # Create a sparse distance matrix
    distance_matrix = lil_array((num_sequences, num_sequences), dtype=float)

    # Fill the sparse distance matrix with the similarity scores
    for (i, j), score in similarity_scores.items():
        distance_matrix[i, j] = score

    distance_matrix = distance_matrix.tocsr()
    sort_graph_by_row_values(distance_matrix)

    print("Nonzeros in distance matrix: " + str(distance_matrix.count_nonzero()))

    # Compute DBSCAN
    # whether eps is 0.1 or 1, the groups are the same (above the threshold of the blast search)
    db = DBSCAN(metric="precomputed", eps = 1, min_samples=3, n_jobs=-1).fit(distance_matrix)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)

    # Print out the groups with suffix .seq_list.fasta
    # for each label in labels, get all db[labels == label] and print its members
    groups = defaultdict(list)
    for i in range(num_sequences):
        if labels[i] != 1:
            groups[labels[i]].append(index_to_sequence[i])

    for label, group in groups.items():
        if len(group) > 3 and len(group) < max_group_size:
            with open(f"group{label}.seq_list", "w") as fout:
                for g in group:
                    print(g, file=fout)
                    
@bash_app
def select_group_seq(input: str, WORKING_DIR: str, SHARED_PATH: str):
    return f'''search_list=$(grep -f {input} -l *.fasta) && {SHARED_PATH}/seqkit grep -f {input} -o {input}.fasta ${{search_list}} && rm {input}'''

@bash_app
def mafft(input: str):
    return f'''mafft --auto --thread -1 {input} > {input}.aln && rm {input}'''

@bash_app
def seq_list_to_gene_tree(input: str, WORKING_DIR: str, SHARED_PATH: str):
    return f'''search_list=$(grep -f {input} -l {WORKING_DIR}/*.fasta) && {SHARED_PATH}/seqkit grep -f {input} -o {input}.fasta ${{search_list}} && rm {input} && mafft --auto --thread -1 {input}.fasta > {input}.aln && sed -i 's/_gene:.*$//g' {input}.aln && sed -i 's/_CDS:.*$//g' {input}.aln && sed -i 's/_ORF:.*$//g' {input}.aln && {SHARED_PATH}/seqkit rmdup -n {input}.aln > {input}.tmp && mv {input}.tmp {input}.aln && raxml-ng --search1 --msa {input}.aln --model GTR+G --prefix {input} && rm {input}.aln'''

@bash_app
def astralpro(output: str):
    return f'''cat *.bestTree > all_trees.out && sed -i 's/_gene[^:]*//g' all_trees.out && sed -i 's/_ORF[^:]*//g' all_trees.out && sed -i 's/_DN[^:]*//g' all_trees.out && astral-pro -i all_trees.out -o {output}'''

@python_app
def make_temp_dir():
    import tempfile
    import os
    return tempfile.mkdtemp(prefix=os.getcwd()+'/')

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='path to input json file', type=str)
parser.add_argument('output_file', help='path to output file', type=str)
parser.add_argument('-t', '--max_trees', help='maximum gene trees', default=50, type=int, required=False)
parser.add_argument('-g', '--max_group_size', help='maximum number of sequences permitted in a gene group', default=100, type=int, required=False)
parser.add_argument('-e', '--blast_evalue', help='evalue used for blastn search in finding related gene sequences', default='1e-20', type=str, required=False)
args = parser.parse_args()

catalog_file_name = os.path.abspath(args.input_file)
output_file_name = os.getcwd() + "/" + args.output_file
if not float(args.blast_evalue) < 1:
    print(f"Error: blast_evalue must be a number less than 1")
    exit
check_file(args.input_file)

print(f"Input file: {catalog_file_name}")
print(f"Output file: {output_file_name}")

CATALOG_PATH = os.path.dirname(catalog_file_name)
SHARED_PATH = os.path.dirname(os.path.abspath(__file__))
WORKING_DIR = make_temp_dir().result()
os.chdir(WORKING_DIR)
print(WORKING_DIR)

print("Parsing catalog...")
parse_catalog_future = parse_catalog_func(catalog_file_name, WORKING_DIR)
if not parse_catalog_future.result():
    print("Error in parsing catalog")
print("Done")

print("Finding sequences from annotations...")
with open("annotations.txt") as fin:
    lines = fin.readlines()
annotations_future = []
for line in lines:
    annotations_future.append(annotations(line, CATALOG_PATH, SHARED_PATH, WORKING_DIR))
for future in annotations_future:
    if future.result() != 0:
        print(f"Error in extracting sequences from annotation: {future.result()}")
        exit()

print("Relabeling genes...")
QUERY_GENES="query_genes.query_fasta"
with open("genes.txt") as fin:
	lines = fin.readlines()
relabel_genes_future = []
for line in lines:
	relabel_genes_future.append(relabel_genes(line, CATALOG_PATH, SHARED_PATH, WORKING_DIR))

print("Finding ORFs in unannotated genomes...")
blast_search_files = []
with open("genomes.txt") as fin:
	lines = fin.readlines()
find_orfs_future = []
for line in lines:
    find_orfs_future.append(find_orfs(line, CATALOG_PATH, SHARED_PATH, WORKING_DIR))
for future in relabel_genes_future:
    if future.result() != 0:
        print(f"Error in relabeling genes: {future.result()}")
        exit()
    else:
        with open(future.stdout) as f:
            for line in f:
                pass
            blast_search_files.append(line.strip())
for future in find_orfs_future:
    if future.result() != 0:
        print(f"Error in finding ORFs: {future.result()}")
        exit()
    else:
        with open(future.stdout) as f:
            for line in f:
                pass
            blast_search_files.append(line.strip())
print("Done")

print("Building BLAST DB of potential gene sequences")
build_blast_db_future = build_blast_db(WORKING_DIR)
if(build_blast_db_future.result() != 0):
    print(f"Error in building BLAST DB: {build_blast_db_future.result()}")
print("Done")

print("Searching BLAST DB for matching sequences")
os.mkdir("blast_results")
search_blast_future = []
#print(blast_search_files)
for line in blast_search_files:
    search_blast_future.append(search_blast(line, WORKING_DIR, args.blast_evalue))
for future in search_blast_future:
    if future.result() != 0:
        print(f"Error in searching BLAST DB: {future.result()}")
        exit()
print("Done")

copy_future = copy_blast_to_csv(WORKING_DIR)
if copy_future.result() != 0:
    print(f"Error in copying BLAST search results to csv file")

print("Grouping the sequences")
grouping_future = group(WORKING_DIR, args.max_group_size)
if grouping_future.result() != 0:
    print(f"Error in grouping sequences: {grouping_future.result()}")
print("Done")

print("Generate gene trees")
gene_tree_futures = []
count = 0
for filename in glob.glob(WORKING_DIR+"/*.seq_list"):
    count += 1
    gene_tree_futures.append(seq_list_to_gene_tree(filename, WORKING_DIR, SHARED_PATH))
    if count >= args.max_trees:
        break
for future in gene_tree_futures:
    if future.result() != 0:
        print(f"Error in generating gene trees: {future.result()}")
print("Done")

print("Generate species tree")
species_tree_future = astralpro(output_file_name)
if species_tree_future.result() != 0:
    print("Error in ASTRAL-Pro")
print("Done")
