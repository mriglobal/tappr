from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.SeqUtils import GC
from skbio import DNA, TabularMSA
import annoy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import average, fcluster
from random import randint
import multiprocessing as mp
import pandas as pd
import numpy as np
import zipfile
from io import StringIO
import time
import os
import argparse
import re


parser = argparse.ArgumentParser(description="Count kmers across a multifasta and cluster these kmers to form degenerate motifs..")
    
    #command line arguments
parser.add_argument('-k', nargs='?',type=int, default=18,help="kmer size")
parser.add_argument('--directory', default=os.getcwd(),help="Output file directory")
parser.add_argument('--seqs', required=True,help="Multifasta file containing sequence records.")
parser.add_argument('--cluster_distance',default=1.0,type=float,help="Hierarchically cluster kmers and produce degenerate motifs at a flat distance.")
parser.add_argument('--singletons',action='store_true',default=False,help='Remove K-mers that only appear once in the data before clustering.')
#parser.add_argument('-f',default=None,help='Optional fasta file containing a restricted set of k-mers to be used for set complete. Default: None')
parser.add_argument('-o',default=None,help="Output filename prefix. Default is to use reference genome sequence ids,flanking sequence identity, K size. Must be specified if no reference is designated.")

myargs=parser.parse_args()

k_size=myargs.k

if myargs.o:
    outfile=myargs.o
else:
    outfile="{}.{}mer".format(myargs.seqs.split(".",maxsplit=1)[0],myargs.k)

seqs = list(SeqIO.parse(myargs.seqs,'fasta'))

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
encode_dict = {'A':0,
             'T':1,
             'C':2,
             'G':3}

decode_dict = {0:'A',
              1:'T',
              2:'C',
              3:'G'}


dim = k_size*4
n_trees = 60 #Hard coded right now
n_neighbors = 20 #Hard coded right now
max_dist = round((myargs.cluster_distance/k_size), 2)
basekey = {tuple(sorted(values)):keys for keys,values in DNA.degenerate_map.items()}
basekey.update({('A',):'A',('G',):'G',('C',):'C',('T',):'T'})
motifs = []

def _kmer_decode(v,decode_dict,ksize):
    return ''.join([decode_dict[np.where(c)[0][0]] for c in np.lib.stride_tricks.as_strided(v,shape=[ksize,4])])

def _kmer_encode(encode_dict, kmer):
    vec = np.zeros(len(kmer)*4,dtype=np.int8)
    for n,nuc in enumerate(kmer):
        vec[4*n+encode_dict[nuc]] = 1
    return vec

def _count_kmers_dict(sequences, k):
    kmers_counter = Counter()
    for sequence in sequences:
        kmers_set = set()
        for i in range(len(sequence) - k + 1):
            kmers_set.update([sequence[i:i + k]])
        kmers_counter.update(kmers_set)
    return kmers_counter


kmers = _count_kmers_dict([str(s.seq) for s in seqs],k_size)
unambig = set('ATCG')
kmers_filtered = {k:c for k,c in kmers.items() if set(k).issubset(unambig)} #in case singletons want to be filtered: and c >= self.min_count}
kmers_ranked = pd.Series(kmers_filtered).sort_values(ascending=False)
if myargs.singletons:
    kmers_ranked = kmers_ranked[kmers_ranked > 1].copy()
ann = annoy.AnnoyIndex(dim,'hamming')
for i,t in enumerate(kmers_ranked.index):
    vec = _kmer_encode(encode_dict,t)
    ann.add_item(i,vec)
ann.build(n_trees)
total_items = ann.get_n_items()
print("{} total kmers in index.".format(total_items))
subtract_set = set()

i = 0
while len(subtract_set) < total_items:
    if i not in subtract_set:
        result = ann.get_nns_by_item(i,n_neighbors,include_distances=True)
        indices = np.array([index for index, distance in zip(result[0],result[1]) if index not in subtract_set])
        subtract_set.update([i])
        i+=1
            #print(len(indices))
        if len(indices) > 1:
            motif_array = np.zeros((len(indices),k_size*4))
            for n,index in enumerate(indices):
                motif_array[n] = ann.get_item_vector(index)
            kmerdist = pdist(motif_array,'hamming')*2
            Z = average(kmerdist)
                    #this could probably be much faster using boolean vector sliding window
            clusters = fcluster(Z,max_dist,criterion='distance')
            #cluster number:list of cluster indices
            myclusters = {key:[] for key in set(clusters)}
            for index, clust in enumerate(clusters):
                myclusters[clust].append(index)
            for clust in myclusters.keys():
                if len(myclusters[clust]) > 1:
                    position_vars = [tuple(set(str(x))) for x in TabularMSA([DNA(_kmer_decode(s,decode_dict,k_size)) for s in motif_array[myclusters[clust]]]).iter_positions()]
                    motif = ''.join([basekey[tuple(sorted(p))] for p in position_vars])
                    motifs.append(motif)
                    subtract_set.update(indices[myclusters[clust]])
                    #Combine accession sets for children kmers in kmer_accessions
                else:
                    subtract_set.update(indices[myclusters[clust]])
                    motif = _kmer_decode(motif_array[myclusters[clust][0]],decode_dict,k_size)
                    motifs.append(motif)
        elif indices.size > 0:
            motif = _kmer_decode(ann.get_item_vector(indices[0]),decode_dict,k_size)
            motifs.append(motif)
    else:
        i+=1
    if i%100000 == 0:
        print("Processed {} kmers.".format(i))


kmerlist = []
for n,k in enumerate(motifs):
    kmerlist.append(SeqRecord(Seq(k),id=str(n),description=str(k_size)+"mers"))
os.chdir(myargs.directory)
SeqIO.write(kmerlist,outfile+".fasta",'fasta')
