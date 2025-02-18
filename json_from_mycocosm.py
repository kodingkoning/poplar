import json
import sys

# Usage: python json_from_mycocosm.py [json input] [file listing fastas]
# Output to stdout
# file listing fastas must follow Mycocosm standard of having NAME/file_name, where NAME is how the species will appear in the tree
# fasta file paths must be relative and only include directory name and the file name

#FILE_TYPE='GENOMIC_NUCLEOTIDE_FASTA'
FILE_TYPE='CDS_NUCLEOTIDE_FASTA'
FILE_TYPE2='GFF3'

if len(sys.argv) < 3:
    raise ValueError("missing input argument")

with open(sys.argv[1]) as fin:
    data = json.load(fin)

if len(sys.argv) > 3:
    with open(sys.argv[2]) as fin:
        lines = fin.readlines()
    with open(sys.argv[3]) as fin:
        lines2 = fin.readlines()
    for line1, line2 in zip(lines, lines2):
        line1 = line1.strip()
        line2 = line2.strip()
        name = line1.split('/')[0]
        data['assemblies'].append({'accession':name, 'files':[{'filePath':line1, 'fileType':FILE_TYPE}, {'filePath':line2, 'fileType':FILE_TYPE2}]})
else:
    with open(sys.argv[2]) as fin:
        for line in fin.readlines():
            line = line.strip()
            name = line.split('/')[0]
            data['assemblies'].append({'accession':name, 'files':[{'filePath':line, 'fileType':FILE_TYPE}]})

print(json.dumps(data, sort_keys=True, indent=4))
