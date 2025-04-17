import os
import random
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
import argparse as ap
from Bio import SeqIO
import sourmash


def fasta_check(filename):
    if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
        return True
    else:
        return False

parser = ap.ArgumentParser(description="Clustering method to bin sequences into subgroups based on their jaccard similarity.")

#command line arguments
#TODO: add ability to file walk across directories to hash assemblies instead of just multifasta
parser.add_argument('--seqs',required=True,help='Either a single multifasta containing sequences or a directory containing fasta files to cluster.')
parser.add_argument('-o',default=False,help='Directory for all output files.')
parser.add_argument('-l',default=None,type=str,help='Label prefix for clusters. This prefix will also be applied to cluster names if label tsv output is specified.')
parser.add_argument('-k',default=11,type=int,help='K length for sourmash kmer sketching. Default: 11')
parser.add_argument('--abundance',action='store_true',default=False,help='Whether to track kmer abundances. Useful when binning sequences of widely varying sizes. Default: False')
parser.add_argument('--criterion',choices=['maxclust','maxdist'],default='maxclust',help="Criterion for clustering: 'maxclust' specifies desired number of clusters while 'maxdist' specifies maximum distance to cut tree. Default: maxclust")
parser.add_argument('--reference',action='store_true',default=False,help='Whether to randomly sample a sequence representative for each cluster.')
parser.add_argument('-c',default=None,help="Clustering criterion value. Default is set to 2 is clustering criterion is set to 'maxclust' and .2 if criterion is set to 'maxdist'")
parser.add_argument('--tsv',action='store_true',default=False,help='Flag for output of a TSV file with sequence ids and corresponding cluster membership.')
parser.add_argument('-s',default=1,type=int,help='Scaled argument for sourmash sketch. 1 over S sampled kmers. Default: 1')
myargs=parser.parse_args()

if myargs.c:
   cluster_t = myargs.c
else:
   if myargs.criterion == 'maxclust':
       cluster_t = 2
   else:
       cluster_t = .2
    

#inputs assigned by arguments
if myargs.o:
    output_dir = myargs.o
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        raise OSError('Output directory already exists.')
else:
    output_dir = '{}_clusters_{}_{}'.format(myargs.seqs.rsplit('.',maxsplit=1)[0],myargs.k,myargs.c)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        raise OSError('Output directory already exists.')

klen = myargs.k

if myargs.l:
    label_prefix = myargs.l + '_'
else:
    label_prefix = ''

if fasta_check(myargs.seqs):
    seqs = list(SeqIO.parse(myargs.seqs,'fasta'))

    seq_ids = [s.id for s in seqs]

    sketches = []
    for s in seqs:
        mh = sourmash.MinHash(0,ksize=klen,scaled=myargs.s,track_abundance=myargs.abundance)
        mh.add_sequence(str(s.seq),force=True)
        sketches.append(mh)

elif os.path.isdir(myargs.seqs):
    filelist = []
    filewalk = list(os.walk(myargs.seqs))
    for f in filewalk:
        if f[2]:
            for file in f[2]:
                if fasta_check(file):
                    filelist.append(os.path.join(f[0],file))
    sketches = []
    seq_ids = []
    seqs = []
    for fasta in filelist:
        for s in SeqIO.parse(fasta,'fasta'):
            seqs.append(s)
            seq_ids.append(s.id)
            mh = sourmash.MinHash(0,ksize=klen,scaled=myargs.s,track_abundance=myargs.abundance)
            mh.add_sequence(str(s.seq),force=True)
            sketches.append(mh)

else:
    raise SystemExit("Invalid Sequence Input.")

ignore_abund = not myargs.abundance
sim_matrix = sourmash.compare.compare_all_pairs(sketches,ignore_abund)

#dis_matrix = 1 - sim_matrix
#condensed = scipy.spatial.distance.squareform(dis_matrix)

#Z = linkage(condensed, 'average')
Z = linkage(sim_matrix,'ward')

if myargs.criterion == 'maxclust':
    clusters = fcluster(Z,cluster_t,criterion='maxclust')
else:
    clusters = fcluster(Z,cluster_t,criterion='distance')

myclusters = {key:[] for key in set(clusters)}
for index, clust in enumerate(clusters):
    myclusters[clust].append(index)
    
mydatabase = {}
for c in myclusters.keys():
    seq_clust= [seqs[i] for i in myclusters[c]]
    mydatabase[c] = seq_clust

if myargs.tsv:
    table_data = []
    for c in mydatabase.keys():
        table_data.extend([[s.id,'{}cluster_{}'.format(label_prefix,c)] for s in mydatabase[c]])
    out_table = pd.DataFrame(table_data,columns=['accession','label'])
    print("Writing Labels Tabular File.")
    out_table.to_csv('{}cluster_labels.tsv'.format(label_prefix),sep='\t',index=False)

os.chdir(output_dir)
for c in mydatabase.keys():
    print("Cluster {} contains {} records.".format(c,len(mydatabase[c])))
    os.mkdir('cluster_{}'.format(c))
    SeqIO.write(mydatabase[c],'cluster_{0}/{1}cluster_{0}.fasta'.format(c,label_prefix),'fasta')
    if myargs.reference:
        SeqIO.write(random.choice(mydatabase[c]),'{}cluster_{}_reference.fasta'.format(label_prefix,c),'fasta')

