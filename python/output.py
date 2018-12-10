import json


def json_for_dashboard(input_fasta, input_json, tree, output, wildcards):
    with open(input_json, 'r') as input_file:
         indices = json.load(input_file)

    with open(input_fasta) as file:
        fasta = file.read()
    with open(tree) as file:
        newick = file.read()

    output_dict = {
        'fasta': fasta,
        'newick': newick,
        'CDR3': indices['CDR3'],
        'FR3': indices['FR3']
    }

    with open(output, 'w') as file:
        json.dump(output_dict, file, indent=4)

