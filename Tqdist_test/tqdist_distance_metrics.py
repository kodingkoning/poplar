import sys
import tqdist
from tabulate import tabulate

def compute_distance(tree1, tree2, measure):
    """
    Computes the triplet or quartet distance between two trees.
    
    :param tree1: Newick format string of first tree.
    :param tree2: Newick format string of second tree.
    :param measure: Either "triplet" or "quartet".
    :return: The computed distance.
    """
    try:
        if measure == "triplet":
            return tqdist.triplet_distance(tree1, tree2)
        elif measure == "quartet":
            return tqdist.quartet_distance(tree1, tree2)
        else:
            raise ValueError("Invalid distance measure. Use 'triplet' or 'quartet'.")
    except Exception as e:
        return f"ERROR: {e}"

def load_trees(tree_files):
    """
    Reads tree files into a dictionary.

    :param tree_files: List of file paths.
    :return: Dictionary {filename: tree_string}
    """
    trees = {}
    for file in tree_files:
        try:
            with open(file, "r") as f:
                trees[file] = f.read().strip()
        except Exception as e:
            print(f"Error reading {file}: {e}")
    return trees

def create_distance_matrix(trees, measure):
    """
    Creates a distance matrix for triplet or quartet distances.

    :param trees: Dictionary of tree filenames and their Newick strings.
    :param measure: Either "triplet" or "quartet".
    :return: Distance matrix as a 2D list.
    """
    files = list(trees.keys())
    n = len(files)
    results = [[None for _ in range(n+1)] for _ in range(n)]

    for i in range(n):
        results[i][0] = files[i]  # First column: tree names
        for j in range(n):
            if i == j:
                results[i][j+1] = 0.0  # Identical trees = distance of 0
            else:
                results[i][j+1] = compute_distance(trees[files[i]], trees[files[j]], measure)

    return results, files

def save_results(output_file, triplet_results, quartet_results, files):
    """
    Saves the triplet and quartet distance matrices to a file.

    :param output_file: Path to output file.
    :param triplet_results: Triplet distance matrix.
    :param quartet_results: Quartet distance matrix.
    :param files: List of tree file names.
    """
    with open(output_file, "w") as f:
        f.write("\nTriplet Distance Matrix:\n")
        f.write(tabulate(triplet_results, headers=["Tree"] + files, tablefmt="github", floatfmt=".3f"))
        f.write("\n\nQuartet Distance Matrix:\n")
        f.write(tabulate(quartet_results, headers=["Tree"] + files, tablefmt="github", floatfmt=".3f"))

def main():
    """
    Main function to load trees, compute triplet and quartet distances, and display results.
    """
    # Get tree files from command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python tqdist_distance_metrics.py <tree1> <tree2> ... <treeN>")
        sys.exit(1)

    tree_files = sys.argv[1:]
    trees = load_trees(tree_files)

    # Compute triplet distance matrix
    triplet_results, files = create_distance_matrix(trees, "triplet")
    print("\nTriplet Distance Matrix:")
    print(tabulate(triplet_results, headers=["Tree"] + files, tablefmt="github", floatfmt=".3f"))

    # Compute quartet distance matrix
    quartet_results, _ = create_distance_matrix(trees, "quartet")
    print("\nQuartet Distance Matrix:")
    print(tabulate(quartet_results, headers=["Tree"] + files, tablefmt="github", floatfmt=".3f"))

    # Save results to a file
    output_file = "tqdist_results.txt"
    save_results(output_file, triplet_results, quartet_results, files)
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    main()
