import json
import sys

# Specify the path to your JSONL file
jsonl_file_path = sys.argv[1]

# Open the JSONL file and read it line by line
with open(jsonl_file_path, 'r') as file:
    for line in file:
        # Parse the JSON from each line
        data = json.loads(line)
        
        # Check if the "accession" key is present in the JSON data
        if 'accession' in data:
            # Print the "accession" value
            print(data['accession'])

