import os
import pandas as pd
from Bio import SeqIO
#from Bio.SeqUtils import nt_search
from multiprocessing import Pool
import time
from skbio import DNA
import argparse
import re


def nt_search(s, p_seq):
    '''Takes a target sequence and re.compile object and returns a list of start positions for matches'''
    pat = make_regex_pattern(p_seq)
    results = [p_seq]
    for m in re.finditer(pat,s):
        results.append(m.start(1))
    return results

def make_regex_pattern(subseq):
    ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
    pattern = ""
    for nt in subseq:
        value = ambiguous_dna_values[nt]
        if len(value) == 1:
            pattern += value
        else:
            pattern += f"[{value}]"
    pattern_read_ahead = '(?=('+pattern+'))'
    return pattern_read_ahead

def pscore(primer,primerset):
    '''function that populates bed column number 5 currently this is set to number of records the motif appears in divided by the number of total records'''
    setcompletion = len(primerset)/len(totalset)
    return {'setcompletion':setcompletion,'total_score':int(round((setcompletion* 1000),0))}


def make_splits(x, num_splits):
    splits=[]
    position=0
    for r in range(num_splits):
        if r != num_splits-1:
            splits.append(x[position:position+len(x)//num_splits])
        else:
            splits.append(x[position:])
        position+=len(x)//num_splits
    return splits


def find_alignments(kmer_list):
    my_align = {}
    for k in kmer_list:
        my_align[str(k.seq)] = {r.id:nt_search(str(r.seq),str(k.seq))[1:] for r in myseqs}
    return my_align

def multi_map(func, data):
    with Pool(cores) as pool:
        kmer_splits = make_splits(kmers,cores)
        results = pool.map(func, kmer_splits)
        return results
    

parser = argparse.ArgumentParser(description="Mapping Degenerate Kmers to Reference Sequences")
parser.add_argument('-r',required=True, help="Designated Reference Fasta")
parser.add_argument('--seqs',default=None,help="Complete set of sequences to provide scoring.")
parser.add_argument('-k',required=True, help='Degenerate Kmers Fasta')
parser.add_argument('-s',default=0,type=int,help="Motif conservation cutoff for inclusion in bed file a number between 0 and 1000 (corresponding to 0 and 100 percent coverage). Default: 0")
parser.add_argument('-t',default=os.cpu_count()-2,type=int,help="Number of threads. Default all.")

myargs=parser.parse_args()

binaries = [key for key,values in DNA.degenerate_map.items() if len(values) ==2]
trinaries = [key for key,values in DNA.degenerate_map.items() if len(values) > 2]


rfile = myargs.r
refseq= SeqIO.read(rfile,'fasta')
kfile = myargs.k
kmers=list(SeqIO.parse(kfile,'fasta'))

if myargs.seqs:
    sfile = myargs.seqs
    kfile = myargs.k
    os.chdir(os.getcwd())
    if myargs.t > len(kmers):
        cores = len(kmers)
        print("scaling down threads to match # sequences")
    else:
        cores = myargs.t
    myseqs = list(SeqIO.parse(sfile,'fasta'))

    totalset = set([r.id for r in myseqs])

    alignments = {}

    print("Mapping motifs with {} cores. {}".format(cores,time.asctime()))
    mapped_kmers = multi_map(find_alignments,kmers)

    for m in mapped_kmers:
        alignments.update(m)

    #there may be a more effecient way to do this ¯\_(ツ)_/¯
    del(mapped_kmers)

    print("Building score table and bed files. {}".format(time.asctime()))
    score_table = {}
    for primer in alignments.keys():
        score = pscore(primer,set([k for k in alignments[primer].keys() if alignments[primer][k]]))
        score_table[primer] = score


bed = []
my_align = {}
for k in kmers:
    my_align[str(k.seq)] = nt_search(str(refseq.seq),str(k.seq))[1:]
for primer in my_align:
    if my_align[primer]:
        for pos in my_align[primer]:
                bed.append({'chr':refseq.id,'start':pos,'end':pos+len(primer),'name':primer,'score':score_table[primer]['total_score']})

bed_df = pd.DataFrame(bed)
bed_df = bed_df[['chr','start','end','name','score']].sort_values('start')
bed_df = bed_df[bed_df['score'] >= myargs.s]
bed_df.to_csv(refseq.id.replace('|','_')+'_'+str(myargs.s)+'_primers.bed',sep='\t',header=False,index=False)

