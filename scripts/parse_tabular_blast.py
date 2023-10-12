import argparse
import os
from Bio.Blast import NCBIXML

def check_file(file_path):
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"{file_path} does not exist")
    return file_path

def merge_sets(list_of_sets):
    parent = list(range(len(list_of_sets)))
    rank = [0] * len(list_of_sets)

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        x_root = find(x)
        y_root = find(y)
        if x_root == y_root:
            return
        if rank[x_root] < rank[y_root]:
            parent[x_root] = y_root
        elif rank[x_root] > rank[y_root]:
            parent[y_root] = x_root
        else:
            parent[y_root] = x_root
            rank[x_root] += 1

    for i in range(len(list_of_sets)):
        for j in range(i + 1, len(list_of_sets)):
            if list_of_sets[i].intersection(list_of_sets[j]):
                # TODO: only merge if all species are unique
                union(i, j)

    merged_sets = {}
    for i in range(len(parent)):
        root = find(i)
        if root not in merged_sets:
            merged_sets[root] = set()
        merged_sets[root].update(list_of_sets[i])

    return list(merged_sets.values())


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_files', nargs='+', help='input file paths')
    parser.add_argument('output_prefix', help='output prefix', type=str)
    args = parser.parse_args()

    counter = 0
    genes = []

    for fin in args.input_files:
        with open(fin) as f:
            next_gene = set()
            gene_name = ""
            for line in f:
                A, B = line.split()
                if A != B:
                    if A == gene_name:
                        # split on "_gene" or "_ORF" to find the species name
                        species_A = ""
                        species_B = ""
                        if "_gene" in A:
                            species_A = A.split("_gene")[0]
                        else:
                            species_A = A.split("_ORF")[0]
                        if "_gene" in B:
                            species_B = B.split("_gene")[0]
                        else:
                            species_B = B.split("_ORF")[0]
                        if species_A != species_B:
                            # check if species_B is unique in the set to ensure a single gene from each species
                            unique_species = True
                            for val in next_gene:
                                species = ""
                                if "_gene" in val:
                                    species = val.split("_gene")[0]
                                else:
                                    species = val.split("_ORF")[0]
                                if species == species_B:
                                    unique_species = False
                                    break
                            if unique_species:
                                next_gene.add(B)
                    else:
                        if len(next_gene) > 10:
                            genes.append(next_gene)
                        gene_name = A
                        next_gene = {A,B}
                    # genes.append({A,B})
        print(f"read records from {fin}, genes = {len(genes)}")

    #print(f"skipping merge, {len(genes)} sets")
    print(f"before merge, {len(genes)} sets")
    genes = merge_sets(genes)

    print(f"completed merge, {len(genes)} sets")
    return
    for gene in genes:
        if len(gene) > 10:
            output_name = args.output_prefix + str(counter) + ".seq_list"
            with open(output_name, 'w') as of:
                for seq_name in gene:
                    print(seq_name,file=of)
            counter += 1

if __name__ == '__main__':
    main()
