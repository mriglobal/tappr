import argparse
from Bio import SeqIO
import pandas as pd

def create_accession_label_df(fasta_files, labels):

    df = pd.DataFrame(columns=['accession', 'label'])
    for file, label in zip(fasta_files, labels):
        records = SeqIO.parse(file, 'fasta')
        temp_df = pd.DataFrame({
            'accession': [record.id for record in records],
            'label': [label] * len([record.id for record in SeqIO.parse(file, 'fasta')])  # Repeat label for each accession
        })
        # Append to the main DataFrame
        df = pd.concat([df, temp_df], ignore_index=True)

    return df

parser = argparse.ArgumentParser(description="Map accessions to labels to make a label table input for qpcr inclusive/exclusive like analysis")
parser.add_argument('--groups',nargs='*',help="any number of multifasta files. must be the length of labels args (space separated strings)")
parser.add_argument('--labels',nargs='*',help="any number of labels for the records in the mutlifasta files. must be the length of groups args (space separated strings)")
parser.add_argument('-o', type=str, required=False, default='metadata.tsv', help='name of ouptut file. default is metadata.tsv')

args = parser.parse_args()
if len(args.groups) != len(args.labels):
    print("Error: The number of groups and labels must be the same.")
else:
    result_df = create_accession_label_df(args.groups, args.labels)
    result_df.to_csv(args.o, index=False, sep='\t')
    
