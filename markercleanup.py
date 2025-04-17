import os
import argparse
import subprocess
import pandas as pd
import pybedtools
from Bio import SeqIO
from Bio.Blast.Applications import NcbitblastnCommandline

def formatmarkername(markername):
    if "prot" in markername:
        if len(markername.split(sep="prot")) == 2:
            _, formattedname = markername.split(sep="prot")
            return formattedname
        else:
            formattedname = markername.split(sep="prot")[1]
            return formattedname
    else:
        return markername

parser = argparse.ArgumentParser()
parser.add_argument('-c', required=True,help='Marker fasta file')
parser.add_argument('-r', required=True,help='Reference genome to map to')
parser.add_argument('-e', default=1.0,type=float,help='Evalue cutoff for blast reporting. Default: 1.0')
parser.add_argument('-t',default=os.cpu_count()-2,help='Number of threads for blast. Default: All')
parser.add_argument('-o', required=True,help='Output directory')
myargs = parser.parse_args()


if not os.path.exists(myargs.o):
    os.makedirs(myargs.o)

subprocess.run(["makeblastdb","-in", myargs.r, "-dbtype", "nucl"])
if "/" in myargs.r:
    _, rfile = myargs.r.rsplit(sep="/",maxsplit=1)
else:
    rfile = myargs.r
tblastn_cline = NcbitblastnCommandline(query=myargs.c, db=myargs.r, evalue=myargs.e, max_hsps=1, outfmt=6, out=rfile+".tsv",num_threads=myargs.t)

tblastn_cline()

blast_results_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
df=pd.read_table(rfile+".tsv",names=blast_results_columns)
df = df[(df["pident"] == 100.0) & (df['gapopen'] == 0)]

for i in df.index:
    if df.loc[i,'sstart']>df.loc[i,'send']:
        start = df.loc[i,'send']
        df.loc[i,'send']=df.loc[i,'sstart']
        df.loc[i,'sstart']=start
df["qseqid"]=df["qseqid"].apply(formatmarkername)
grouped = df.groupby('sseqid')
refgroup = {g[0]:g[1] for g in grouped}

myrefs = {r.name:r for r in SeqIO.parse(myargs.r, "fasta")}

blastseqs = {key:[] for key in myrefs.keys()}
for seqid in refgroup.keys():
    for n, s in enumerate(refgroup[seqid][['sstart','send']].itertuples()):
        blastseqs[seqid].append(myrefs[seqid][s.sstart-1:s.send])
        blastseqs[seqid][-1].description=blastseqs[seqid][-1].description+' '+str(n)

for o in blastseqs.keys():
    SeqIO.write(blastseqs[o], '{0}/{1}.faa'.format(myargs.o,o), "fasta")

#for b in blastseqs.keys():
#    blastn_cline = NcbiblastnCommandline(query=b+".faa", db="/home/share/NCBI/blastdb/nt/nt", evalue=0.001, max_hsps=1, outfmt=6, out=b+"_nt.tsv")
#    blastn_cline()


for k in refgroup.keys():
    # Increment all numbers in the 'send' column by 1 because bed files end coordinate is not inclusive
    refgroup[k]['send'] += 1
    bed_tmp = refgroup[k][["sseqid","sstart","send"]].copy()
    bed_tmp.columns = ['chrom','start','end']
    bedout = pybedtools.BedTool().from_dataframe(bed_tmp).sort().merge().to_dataframe()
    bedout.to_csv("{0}/{1}_markers.bed".format(myargs.o,k),sep='\t',header=False, index=False)
