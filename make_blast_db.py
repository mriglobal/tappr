# need to make blast db from some number of input fastas (just cat them together)
# OR need to make it from an input taxid, make the pull, then build
import argparse
import subprocess
import zipfile
from Bio import SeqIO
from io import StringIO
import os

def fasta_check(filename):
    if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
        return True
    else:
        return False

parser = argparse.ArgumentParser(description='automate building of blast database for taxid and or specified organism input')

parser.add_argument('-t', '--taxid', required=False, type=int, default=None, help='taxid to construct db for')
parser.add_argument('-s', '--seqs', required=False, type=str, default=None, help='file containing sequences to build a database for')
parser.add_argument('-o', '--outname', required=False, type=str, default='tappr_db/tappr_db', help='name for output blast database. default: tappr_db/tappr_db')

args = parser.parse_args()

if args.seqs and args.taxid:
    exit("Pass only one of --seqs or --taxid")

if args.seqs:
    print("Making Blastdb")
    subprocess.run(['makeblastdb', '-in', args.seqs, '-dbtype', 'nucl', '-parse_seqids', '-out', args.o])
    
elif args.taxid:
    subprocess.run(['datasets', 'download', 'genome', 'taxon', str(args.taxid), '--assembly-level', 'chromosome,complete', '--filename', 'tappr_datasets_pull.zip'])
    with zipfile.ZipFile('tappr_datasets_pull.zip','r') as myzip:
        filelist = [f for f in myzip.namelist() if fasta_check(f)]
        records = [list(SeqIO.parse(StringIO(myzip.read(file).decode()),format="fasta")) for file in filelist]
        records = [item for sublist in records for item in sublist]
    
    print("Writing exclusivity multifasta")

    SeqIO.write(records, "tappr_multi.fasta", "fasta")
    
    print("Making Blastdb")
    subprocess.run(['makeblastdb', '-in', "tappr_multi.fasta", '-dbtype', 'nucl', '-parse_seqids', '-out', args.o])
    
    print("Cleaning Temp Files")
    try:
        os.remove("tappr_multi.fasta")
        print("Removed tappr_multi.fasta")
    except Exception as e:
        print(f"Failed to remove tappr_multi.fasta: {e}")
    try:
        os.remove("tappr_datasets_pull.zip")
        print("Removed tappr_datasets_pull.zip")
    except Exception as e:
        print(f"Failed to remove tappr_datasets_pull.zip: {e}")