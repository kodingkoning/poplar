from parsl import join_app

@join_app
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

    group_fastas = []
    print("Writing groups to files")
    for label, group in groups.items():
        if len(group) > 3 and len(group) < 100: # TODO: Select maximum group size carefully
            filename = f"group{label}.seq_list"
            with open(filename, "w") as fout:
                for g in group:
                    print(g, file=fout)
