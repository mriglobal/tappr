from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Perform set difference operation between two Kmer sets: A - B.")

#command line arguments
parser.add_argument('-k', nargs='?',type=int, default=18,help="kmer size")
parser.add_argument('-A',required=True, type=str,help='A K-mer set. Assumed to be Kmer fasta file. Use --A_set flag if file is Kmer set pickle file.')
parser.add_argument('-B', required=True, type=str, help="B K-mer set. Assumed to be Kmer set pickle file.")
parser.add_argument('--A_set', action='store_true',default=False,help="Flag for A set file type.")
parser.add_argument('-o',type=str,default=None,help="Output filename prefix. Default is inherit A set file name with '_exclusive' added.")
myargs=parser.parse_args()

#variables assigned by arguments
if myargs.o:
    out_prefix = myargs.o
else:
    out_prefix = myargs.A.split(".")[0]

ksize=myargs.k
if not myargs.A_set:
    A_set = set([str(s.seq) for s in SeqIO.parse(myargs.A,'fasta')])
else:
    A_set = pd.read_pickle(myargs.A)

B_set = pd.read_pickle(myargs.B)

A_sniff = A_set.pop()
B_sniff = B_set.pop()

if len(A_sniff) != len(B_sniff) or len(A_sniff) != myargs.k or len(B_sniff) != myargs.k:
    raise ValueError("K-mer size mismatch.")
A_set.add(A_sniff)
B_set.add(B_sniff)

outset = A_set - B_set

output = []
for i,k in enumerate(outset):
    output.append(SeqRecord(Seq(k),id='marker_'+str(i),description=''))

SeqIO.write(output,out_prefix+"_exclusive.fasta",'fasta')
