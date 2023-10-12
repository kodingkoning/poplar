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

	for fin in args.input_files:
		with open(fin) as f:
			blast_records = NCBIXML.parse(f)
			print("opened records")

			genes = []
			for blast_record in blast_records:
				if len(blast_record.alignments) > 1:
					genes.append( set([aln.hit_def for aln in blast_record.alignments]) )
		print("read records from", fin)

	genes = merge_sets(genes)

	print(f"completed merge, {len(genes)} sets")
	for gene in genes:
		if len(gene) > 3:
			output_name = args.output_prefix + str(counter) + ".seq_list"
			with open(output_name, 'w') as of:
				for seq_name in gene:
					print(seq_name,file=of)
			counter += 1

if __name__ == '__main__':
	main()

