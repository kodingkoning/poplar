import parsl
from parsl import python_app, bash_app, join_app
import argparse
import os
import tempfile
import glob
from config import config

# TODO: draw the dependency graph very clearly - use the tool I used for the graph for the thesis schedule, and make sure to indicate how the number of tasks is related to previous numbers of tasks

parsl.load(config)

def query_genes_func():
	pass

# TODO: add the ability to fail gracefully when any of the function calls exit

# TODO: add the ability to continue after failure

# Throw error if file_path does not exist
def check_file(file_path):
	if not os.path.isfile(file_path):
		raise argparse.ArgumentTypeError(f"File {file_path} does not exist")

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
def search_blast(line: str, WORKING_DIR: str):
    import os
    os.chdir(WORKING_DIR)
    line = line.strip()
    name = line.split('.')[0]
    filename = f"blast_{name}.fasta"
    # TODO: allow user to choose evalue and number of threds
    return f'''blastn -query {line} -db genes.db -out blast_results/{filename}.tab -task dc-megablast -evalue 1e-20 -outfmt "6 qaccver saccver evalue" -num_threads 32 && echo 0'''

@bash_app
def copy_blast_to_csv(WORKING_DIR: str):
    return f'''cd {WORKING_DIR} && cat blast_results/*.tab > blast_results.csv'''

@python_app
def group(WORKING_DIR):
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
        if len(group) > 3 and len(group) < 100: # TODO: Select maximum group size carefully
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
def raxml(input: str):
    # sed is used to remove illegal parts of gene names
    return f'''sed -i 's/_gene:.*$//g' {input} && sed -i 's/_CDS:.*$//g' {input} && sed -i 's/_ORF:.*$//g' {input} && seqkit rmdup -n {input} > {input}.tmp && mv {input}.tmp {input} && raxml-ng --search1 --msa {input} --model GTR+G --prefix {input} && rm {input}'''

@bash_app
def seq_list_to_gene_tree(input: str, WORKING_DIR: str, SHARED_PATH: str):
    return f'''search_list=$(grep -f {input} -l *.fasta) && {SHARED_PATH}/seqkit grep -f {input} -o {input}.fasta ${{search_list}} && rm {input} && mafft --auto --thread -1 {input}.fasta > {input}.aln && sed -i 's/_gene:.*$//g' {input}.aln && sed -i 's/_CDS:.*$//g' {input}.aln && sed -i 's/_ORF:.*$//g' {input}.aln && seqkit rmdup -n {input}.aln > {input}.tmp && mv {input}.tmp {input}.aln && raxml-ng --search1 --msa {input}.aln --model GTR+G --prefix {input} && rm {input}.aln'''

@bash_app
def astralpro(output: str):
    return f'''cat *.bestTree > all_trees.out && sed -i 's/_gene[^:]*//g' all_trees.out && sed -i 's/_ORF[^:]*//g' all_trees.out && sed -i 's/_DN[^:]*//g' all_trees.out && astral-pro -i all_trees.out -o {output}'''

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='path to input json file', type=str)
parser.add_argument('output_file', help='path to output file', type=str)
args = parser.parse_args()

catalog_file_name = os.path.abspath(args.input_file)
output_file_name = os.getcwd() + "/" + args.output_file
check_file(args.input_file)

print(f"Input file: {catalog_file_name}")
print(f"Output file: {output_file_name}")

CATALOG_PATH = os.path.dirname(catalog_file_name)
SHARED_PATH = os.path.dirname(os.path.abspath(__file__))
WORKING_DIR = tempfile.mkdtemp(prefix=os.getcwd()+'/')
os.chdir(WORKING_DIR)
print(WORKING_DIR)

print("Parsing catalog...")
parse_catalog_future = parse_catalog_func(catalog_file_name, WORKING_DIR)
if not parse_catalog_future.result():
    print("Error in parsing catalog")
print("Done")

print("Relabeling genes...")
QUERY_GENES="query_genes.query_fasta"
with open("genes.txt") as fin:
	lines = fin.readlines()
relabel_genes_future = []
for line in lines:
	relabel_genes_future.append(relabel_genes(line, CATALOG_PATH, SHARED_PATH, WORKING_DIR))
    # TODO: take all the relabel_genes output files and send them right to search_blast

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
print(blast_search_files)
for line in blast_search_files:
    search_blast_future.append(search_blast(line, WORKING_DIR))
for future in search_blast_future:
    if future.result() != 0:
        print(f"Error in searching BLAST DB: {future.result()}")
        exit()
print("Done")

copy_future = copy_blast_to_csv(WORKING_DIR)
if copy_future.result() != 0:
    print(f"Error in copying BLAST search results to csv file")

print("Grouping the sequences")
grouping_future = group(WORKING_DIR)
if grouping_future.result() != 0:
    print(f"Error in grouping sequences: {grouping_future.result()}")
print("Done")

print("Generate gene trees")
gene_tree_futures = []
# TODO: add a input for the number of gene trees to run
count = 0
for filename in glob.glob(WORKING_DIR+"/*.seq_list"):
    count += 1
    gene_tree_futures.append(seq_list_to_gene_tree(filename, WORKING_DIR, SHARED_PATH))
    if count > 49:
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

# TODO: make sure this is all correct (we get a resulting tree) and add comments and make it look clean

# TODO: test that this is working, will need to increase the job time limit
# TODO: get group to return the number of seq
# TODO: or should group() be a join_app and spawn jobs?

# TODO: figure out if there is a way to structure this that we specify the dependencies so the main script can die -- use the dynamic workflows where some of the jobs will spawn more jobs
# This should be done with join_app, which calls python_apps or bash_apps. See Fibonacci example
