import re
import sys

class Node:
    def __init__(self, name):
        self.name = name
        self.children = []

    def add_child(self, child):
        self.children.append(child)

    def to_newick(self):
        if self.children:
            child_strings = [child.to_newick() for child in self.children]
            return f'({",".join(child_strings)})'
        else:
            return self.name

def build_tree(taxonomic_data):
    root = Node('Root')

    for line in taxonomic_data:
        classifications = line.split(';')
        parent = root

        for classification in classifications:
            classification = classification.strip()
            child = next((node for node in parent.children if node.name == classification), None)

            if not child:
                child = Node(classification)
                parent.add_child(child)

            parent = child

    return root

def read_file(file_name):
    try:
        with open(file_name, 'r') as file:
            return file.readlines()
    except IOError:
        print(f"Error: File '{file_name}' not found.")
        sys.exit(1)

def main(file_name):
    taxonomic_data = read_file(file_name)
    tree = build_tree(taxonomic_data)
    newick_tree = tree.to_newick()
    print(newick_tree,";")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python taxonomy_to_newick.py <file_name>")
        sys.exit(1)

    file_name = sys.argv[1]
    main(file_name)

