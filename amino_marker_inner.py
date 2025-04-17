import zipfile
from io import StringIO
from Bio import SeqIO
from skbio import Protein
import multiprocessing as mp
import pickle
import argparse
import time
import os

def make_splits(x, num_splits):
    '''generic chunking recipe'''
    splits = []
    position = 0
    for r in range(num_splits):
        if r!= num_splits-1:
            splits.append(x[position:position+len(x)//num_splits])
        else:
            splits.append(x[position:])
        position+=len(x)//num_splits
    return splits

def get_kmers(seqs):
    kmer_list = []
    for genome in seqs:
        new = set()
        for seq in genome:
            for feature in seq.features:
                if feature.type == 'CDS':
                    try:
                        new.update(Protein(str(feature.extract(seq).translate().seq)).kmer_frequencies(myargs.k).keys())
                    except ValueError:
                        pass
        kmer_list.append(new)
    inner_kmers = set.intersection(*kmer_list)
    print("Found {} inner kmers from current chunk.".format(len(inner_kmers)))
    return inner_kmers

parser = argparse.ArgumentParser(description='Amino K-mer Inner Join for Marker Discovery')

#ncbi datasets data structure: ncbi_dataset/data/<Assembly Accession>/genomic.gbff
parser.add_argument("-k",type=int,default=8,help="K size for K-mer counting.")
parser.add_argument("--seqs",type=str,help='Either a multifasta containing CDS sequences or a Zip file containing gbff files from NCBI Dataset')
parser.add_argument("-t",type=int,default=os.cpu_count()-2,help="Count K-mers in memory and use T threads.")
parser.add_argument("-c",type=int,default=1,help="Number of chunks to break zip file into when loading into memory. Default: 1")
parser.add_argument('-o',type=str,default=None,help="Output prefix. Default is to inherit from zip file.")
#myargs = parser.parse_args(['-k','5','-z','/home/pdavis/Documents/ASSAY_DESIGN/test_data/Ecoli_O157.zip'])
myargs = parser.parse_args()

if not myargs.o:
    out_prefix = myargs.seqs.split(".")[0]
else:
    out_prefix = myargs.o

if myargs.seqs.split('.')[-1] in set(['fna','fasta','a']):
    seqs = [seq for seq in SeqIO.parse(myargs.seqs,'fasta')]
    seq =  seqs.pop()
    output = set(Protein(str(seq.seq.translate())).kmer_frequencies(myargs.k).keys())
    for s in seqs:
        output = output.intersection(Protein(str(s.seq.translate())).kmer_frequencies(myargs.k).keys())
        print("Remaining Kmers: {}".format(len(output)))

elif myargs.seqs.endswith('.zip'):
    if not myargs.t:
        with zipfile.ZipFile(myargs.seqs,'r') as myzip:
            files = [f for f in myzip.namelist() if f.endswith(".gbff")]
            output = set()
            file = files.pop()
            for seq in SeqIO.parse(StringIO(myzip.read(file).decode()),'genbank'):
                for feature in seq.features:
                    if feature.type == 'CDS':
                        #protein sequences below K size will throw ValueError numpy error
                        try:
                            output.update(Protein(str(feature.extract(seq).translate().seq)).kmer_frequencies(myargs.k).keys())
                        except ValueError:
                            pass
            print("Remaining Kmers: {}".format(len(output)))
            for file in files:
                new = set()
                for seq in SeqIO.parse(StringIO(myzip.read(file).decode()),'genbank'):
                    for feature in seq.features:
                        if feature.type == 'CDS':
                            try:
                                new.update(Protein(str(feature.extract(seq).translate().seq)).kmer_frequencies(myargs.k).keys())
                            except ValueError:
                                pass
                output = output.intersection(new)
                print("Remaining Kmers: {}".format(len(output)))

    else:
        with zipfile.ZipFile(myargs.seqs,'r') as myzip:
            files = [f for f in myzip.namelist() if f.endswith(".gbff")]
            output = set()
            file = files.pop()
            for seq in SeqIO.parse(StringIO(myzip.read(file).decode()),'genbank'):
                for feature in seq.features:
                    if feature.type == 'CDS':
                        #protein sequences below K size will throw ValueError numpy error
                        try:
                            output.update(Protein(str(feature.extract(seq).translate().seq)).kmer_frequencies(myargs.k).keys())
                        except ValueError:
                            pass
            print("Remaining Kmers: {}".format(len(output)))
            print("Reading in Sequence Data in {} chunks. {}".format(myargs.c,time.asctime()))
            file_splits = make_splits(files,myargs.c)
            for file_chunk in file_splits:
                seqs = [list(SeqIO.parse(StringIO(myzip.read(file).decode()),'genbank')) for file in file_chunk]
                if myargs.t > len(seqs):
                    threads = len(seqs)
                    print("scaling down threads to match # sequences")
                else:
                    threads = myargs.t
                print("Counting Kmers with {} processes.".format(threads))
                with mp.Pool(threads) as pool:
                    kmer_splits = make_splits(seqs,threads)
                    kmer_sets = pool.map(get_kmers,kmer_splits)
                output = output.intersection(set.intersection(*kmer_sets))
                print("Remaining Kmers: {} {}".format(len(output),time.asctime()))
            

#TODO: Maybe change this to json output
with open(out_prefix+'_inner_{}mers.pickle'.format(myargs.k),'wb') as outfile:
    pickle.dump(output,outfile)


