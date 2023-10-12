import json
import sys

def main():
    filename = sys.argv[1]
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

if __name__ == "__main__":
    main()
