from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
import pickle
import argparse
import os
import time
from Bio.SeqRecord import SeqRecord


def get_chunks(x, num_chunks):
     chunks=[]
     position=0
     for r in range(num_chunks):
             if r != num_chunks-1:
                     chunks.append(x[position:position+len(x)//num_chunks])
             else:
                     chunks.append(x[position:])
             position+=len(x)//num_chunks
     return chunks

def kmercount(replicon, k):
    kd = set()
    for j in range(len(replicon)-k):
        kmer = Seq(replicon[j:j+k])
        kd.add(str(kmer).upper())
        kd.add(str(kmer.reverse_complement()).upper())
    return kd

def wrapper(replicon, ksize):
    return kmercount(replicon, ksize)

def process_chunk(chunk, ksize, cores):
    with Pool(cores) as pool:
        results = pool.starmap(wrapper, [(str(x), ksize) for x in chunk])
        return results

parser = argparse.ArgumentParser(description="Union Kmer Set Evaluation")
parser.add_argument('--seqs', required=True,help='Folder for Sequences.')
parser.add_argument('-k', nargs='?', default=18,help='K size for counting. Default: 18')
parser.add_argument('-l',default=0, type=int, help='Minimum length filter for inclusion in kmer counting. (Default: 0)')
parser.add_argument('-o',default="",help='Output pickled kmer set file name prefix. (Default: None')
parser.add_argument('-c', type=int, default=1, help='Number of chunks to split the sequences into before kmercounting.')
parser.add_argument('-t',type=int,default=os.cpu_count()-2,help='Number of threads for multiprocessing. Default: All.')
parser.add_argument('-i', type=str, required=False, default=None, help='Pass the fasta output of kmercountinner.py to do set difference (returns additional file)')
parser.add_argument('--fastq', action='store_true',default=False,help="Flag for fastq input instead of fasta.")

myargs = parser.parse_args()
cores=myargs.t
ksize = int(myargs.k)
min_len = myargs.l
outname = myargs.o

os.chdir(os.getcwd())

refdirectory = myargs.seqs
os.chdir(refdirectory)
filelist = os.listdir()
recdata = []
for file in filelist:
    if not myargs.fastq:
        file_format = 'fasta'
    else:
        file_format = 'fastq'
    if ".fna" in file or ".fasta" in file or ".fastq" in file:
        records = [str(s.seq) for s in SeqIO.parse(file,format=file_format) if len(s.seq) > min_len]
        #if len(records) <= len(refstats):
        recdata.extend(records)

kmerset = set()
seq_chunks = get_chunks(recdata,myargs.c)


print("Processing chunks. {}".format(time.asctime()))

num_processed = 0
for seq in seq_chunks:
    kmers = process_chunk(seq, ksize, cores)
    kmerset.update(set.union(*kmers))
    num_processed += len(seq)
    print("{} records processed.".format(num_processed))
    print("{} unique kmers.".format(len(kmerset)))

with open(outname+str(ksize)+"Kmer_set.pickle",mode='wb') as outfile:
    pickle.dump(kmerset,outfile)

if myargs.i:
    outname = myargs.i
    outname = outname.split('.fasta')[0]
    print('Reading in in group kmers to check vs exclusive')
    inner_set = set()
    for record in SeqIO.parse(myargs.i, "fasta"):
        inner_set.update(set(record.seq.upper()))
    
    result_set = inner_set - kmerset
    print('Found '+str(len(result_set))+' kmers in inner set not found in outer set')
    kmerlist = []
    for n,k in enumerate(result_set):
        kmerlist.append(SeqRecord(Seq(k),id=str(n),description=str(ksize)+"mers"))
    os.chdir(myargs.directory)
    SeqIO.write(kmerlist,outname+"_exclusive.fasta",'fasta')




