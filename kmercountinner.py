from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import pandas as pd
import os
import argparse
import multiprocessing

#function to bin replicons in the references based on length
#In the future this could be expanded to bin using GC content if examples arise where replicon lengths are within 10% of one another
def replicontest(ref, example):
    refindex=0
    for i in range(len(ref)):
        if len(example) < ref[refindex][1]+(ref[refindex][1]*percent) and len(example) > ref[refindex][1]-(ref[refindex][1]*percent):
            return refindex
        refindex+=1
    return -1

#Kmer counting function. Canonical kmers (rev comp) are also counted
def kmercount(replicon, k):
    kd = {}
    for j in range(len(replicon)-k):
        kmer = Seq(replicon[j:j+k])
        if str(kmer) in kd.keys():
            kd[str(kmer)]+=1
            kd[str(kmer.reverse_complement())]+=1
        else:
            kd[str(kmer)]=1
            kd[str(kmer.reverse_complement())]=1
    return kd


parser = argparse.ArgumentParser(description="Count Kmers on records by performing inner join.")

#command line arguments
parser.add_argument('-k', nargs='?',type=int, default=18,help="kmer size")
parser.add_argument('--directory', default=os.getcwd(),help="output file directory")
parser.add_argument('-r',default=None,help="Reference genome fasta for sequence binning.") #reference genome
parser.add_argument('--seqs', required=True,help="Directory containing sequence files")
parser.add_argument('--assembly_level', default=False, action='store_true', )
parser.add_argument('-o',default=None,help="Output filename prefix.")
parser.add_argument('-p',default=.1, type=float, help='Percent Variance in reference length for replicon binning')
myargs=parser.parse_args()

#variables assigned by arguments
ksize=myargs.k
percent=myargs.p
if not os.path.exists(myargs.directory):
    os.makedirs(myargs.directory)
os.chdir(myargs.directory)
reference = myargs.r
refdirectory = myargs.seqs
if myargs.o:
    outfile=myargs.o
elif myargs.r:
    #reference genome given
    outfile="{}.{}mer_".format(myargs.r.rsplit(".",maxsplit=1)[0],myargs.k)
else:
    #no reference genome and no outfile specified
    outfile = "{}.{}mer".format(os.path.split(os.path.abspath(myargs.seqs))[-1],myargs.k)
    

if reference:
    refrec = list(SeqIO.parse(reference,"fasta"))
    refstats = [(GC(f.seq),len(f.seq)) for f in refrec]
    refreplicon=[]
    for s in refrec:
        refreplicon.append(str(s.seq))

    kmerdicts = {key:{} for key in range(len(refreplicon))}
    #kmerpos = {key:{} for key in range(len(refreplicon))}
    for x,replicon in enumerate(refreplicon):
        kmerdicts[x]=kmercount(replicon,ksize)
    
    os.chdir(refdirectory)
    filelist =os.listdir()
    recdata = {key:[] for key in range(len(refstats))}
    filtered = False
    for file in filelist:
        if ".fna" in file or ".fasta" in file:
            recordname = file.rsplit(sep="_",maxsplit=1)[0]
            records = list(SeqIO.parse(file,format="fasta"))
        #if len(records) <= len(refstats):
            for r in records:
                i = replicontest(refstats, r)
                print(i)
                if i >= 0:
                    recdata[i].append(r)
                else:
                    filtered = True

    kmerseries = {key:pd.Series(kmerdicts[key]) for key in kmerdicts}

    for repliconexample in recdata:
        for sequence in recdata[repliconexample]:
            newkmers = pd.Series(kmercount(str(sequence.seq),ksize))
            kmerseries[repliconexample] = pd.concat([kmerseries[repliconexample],newkmers],axis=1,join="inner",ignore_index=True)
            print("Number of conserved kmers: {}".format(kmerseries[repliconexample].shape))
            print("Length of current sequence: {}".format(len(sequence)))
            print("Sequence ID of current sequence: {}".format(sequence.id))

    for kmers in kmerseries.keys():
        kmerlist=[]
        for n,k in enumerate(list(kmerseries[kmers].index)):
            kmerlist.append(SeqRecord(Seq(k),id=str(n),description=str(ksize)+"mers"))
        os.chdir(myargs.directory)
        SeqIO.write(kmerlist,outfile+str(kmers)+'.fasta','fasta')

    if filtered:
        for r in recdata.keys():
            SeqIO.write(recdata[r],refrec[r].name+"_filtered_seqs.fasta",'fasta')

elif myargs.assembly_level:
    # assembly level kmer counting
    os.chdir(refdirectory)
    filelist = os.listdir()
    recdata = {}
    for file in filelist:
        if file.endswith('.fna') or file.endswith('.fasta'):
            recdata[file] = list(SeqIO.parse(file,'fasta'))
    # rec data is now a dictionary of file names by the records in the file
    kmerdict = dict()
    list_of_lists = [[key, value] for key, value in recdata.items()]
    first_sequences = list_of_lists[0][1]
    first_sequences_key = list_of_lists[0][0]
    kmerseries = [pd.Series(kmercount(str(V.seq),ksize)) for V in first_sequences]
    kmerseries = pd.concat(kmerseries, axis=1, join='outer', ignore_index=True)
    kmerseries['collapsed'] = kmerseries.max(axis=1)
    kmerseries = kmerseries[['collapsed']]
    print("Number of conserved kmers: {}".format(kmerseries.shape))
    print("File name of current sequences: {}".format(first_sequences_key))
    for i in range(1, len(list_of_lists)):
        new_seqs = list_of_lists[i][1]
        new_key = list_of_lists[i][0]
        newkmers = [pd.Series(kmercount(str(V.seq),ksize)) for V in new_seqs]
        newkmers = pd.concat(newkmers, axis=1, join='outer', ignore_index=True)
        newkmers['collapsed'] = newkmers.max(axis=1)
        newkmers = newkmers[['collapsed']]
        kmerseries = pd.concat([kmerseries,newkmers],axis=1,join="inner",ignore_index=True)
        print("Number of conserved kmers: {}".format(kmerseries.shape))
        print("File name of current sequences: {}".format(new_key))
        
    kmerlist = []
    for n,k in enumerate(list(kmerseries.index)):
        kmerlist.append(SeqRecord(Seq(k),id=str(n),description=str(ksize)+"mers"))
    os.chdir(myargs.directory)
    SeqIO.write(kmerlist,outfile+".fasta",'fasta')

else:
    os.chdir(refdirectory)
    filelist = os.listdir()
    recdata = []
    for file in filelist:
        if file.endswith('.fna') or file.endswith('.fasta'):
            recdata.extend(list(SeqIO.parse(file,'fasta')))

    kmerdict = dict()
    first_sequence = recdata.pop()
    kmerseries = pd.Series(kmercount(str(first_sequence.seq),ksize))
    print("Number of conserved kmers: {}".format(kmerseries.shape))
    print("Length of current sequence: {}".format(len(first_sequence)))
    print("Sequence ID of current sequence: {}".format(first_sequence.id))
    for sequence in recdata:
        newkmers = pd.Series(kmercount(str(sequence.seq),ksize))
        kmerseries = pd.concat([kmerseries,newkmers],axis=1,join="inner",ignore_index=True)
        print("Number of conserved kmers: {}".format(kmerseries.shape))
        print("Length of current sequence: {}".format(len(sequence)))
        print("Sequence ID of current sequence: {}".format(sequence.id))

    kmerlist = []
    for n,k in enumerate(list(kmerseries.index)):
        kmerlist.append(SeqRecord(Seq(k),id=str(n),description=str(ksize)+"mers"))
    os.chdir(myargs.directory)
    SeqIO.write(kmerlist,outfile+".fasta",'fasta')
    
        
