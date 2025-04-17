import os
import argparse
import pandas as pd
from Bio import SeqIO
from collections import Counter
from skbio.alignment import local_pairwise_align_ssw
from Bio import pairwise2
from skbio import DNA, TabularMSA
import numpy as np
from sklearn.preprocessing import Binarizer
from multiprocessing import Pool


parser = argparse.ArgumentParser(description="Take amplicons file from simulate_PCR results and realigns primers to give acccurate mismatch data in the presence of degenerate bases.")

#command line arguments
parser.add_argument('-a', required=True,help="Amplicons file from simulate_PCR")
parser.add_argument('-m',default=3,type=int,help="Maximum number of mismatches allowed for alignment.")
parser.add_argument('-o',default='',help="Prefix for output files.")
parser.add_argument('-t',default=os.cpu_count(),type=int,help="Number of threads for multiprocessing. (Default: All)")
parser.add_argument('-n',default=0,type=int,help="Maximum number of allowed mismatches in the last two 3' basepairs in the case of the primers.")

#myargs=parser.parse_args(['-a','/home/pdavis/Documents/DBPAO/filos/batch2/marburg/MARV1_primers_only.fasta.mux.probe_aligned.amplicons'])
myargs=parser.parse_args()

# match, mismatch, open, extend
scoring_scheme = [4,0,-12,-12]
if not myargs.o:
    out_prefix = os.path.basename(myargs.a.rsplit('.',maxsplit=1)[0])
else:
    out_prefix = myargs.o

assert(myargs.n<=2)

data = pd.read_table(myargs.a,low_memory=False)
data = data[~data['AmpliconSeq'].isna()]
data = data[~data['AmpliconSeq'].str.contains("N")]


primer_dict = dict(data[['FP_ID','FP_seq']].drop_duplicates().values)
primer_dict.update(dict(data[['RP_ID','RP_seq']].drop_duplicates().values))

temp_dict = primer_dict
for k in temp_dict.keys():
    primer_dict[k] = [primer_seq for primer_seq in DNA(temp_dict[k]).expand_degenerates()]

#primer_scores = {key:((len(temp_dict[key])-myargs.m)*scoring_scheme[0]) for key in primer_dict.keys()}

#TODO: Produce multiprocessed version of this routine
#forward primer alignment
def new_alignment(x):
    best_plus = 0
    best_minus = 0
    for pr in primer_dict[x[1]['FP_ID']]:
        alignment_plus = pairwise2.align.localms(str(pr),x[1]['AmpliconSeq'],scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
        #mismatch_count = sum(1 for a, b in zip(alignment_plus[0][0][alignment_plus[0][4]-2:alignment_plus[0][4]], alignment_plus[0][1][alignment_plus[0][4]-2:alignment_plus[0][4]]) if a != b)
        #print(mismatch_count,mismatch_count>0)
        #if mismatch_count>0:assert(False)
        #summary_dict['filterdueton']=mismatch_count

        if alignment_plus[0][2] > best_plus:
            best_plus_alignment = alignment_plus
            best_plus = best_plus_alignment[0][2]
            primer_score = int(((len(pr)*scoring_scheme[0]) - alignment_plus[0][2])/scoring_scheme[0])
    mismatch_count = sum(1 for a, b in zip(alignment_plus[0][0][alignment_plus[0][3]:alignment_plus[0][3]+2], alignment_plus[0][1][alignment_plus[0][3]:alignment_plus[0][3]+2]) if a != b)
    x[1]['filterdueton']= mismatch_count>myargs.n
    x[1]['FP_mismatches'] = primer_score
    #reverse primer alignment
    for pr in primer_dict[x[1]['RP_ID']]:
        alignment_minus = pairwise2.align.localms(str(pr),str(DNA(x[1]['AmpliconSeq']).reverse_complement()),scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
        if alignment_minus[0][2] > best_minus:
            best_minus_alignment = alignment_minus
            best_minus = best_minus_alignment[0][2]
            primer_score = int(((len(pr)*scoring_scheme[0]) - alignment_minus[0][2])/scoring_scheme[0])
    x[1]['RP_mismatches'] = primer_score
    #mismatch_count = sum(1 for a, b in zip(alignment_plus[0][0][alignment_plus[0][3]:alignment_plus[0][3]+2], alignment_plus[0][1][alignment_plus[0][3]:alignment_plus[0][3]+2]) if a != b)
        #print(mismatch_count,mismatch_count>0)
        #if mismatch_count>0:assert(False)
        #summary_dict['filterdueton']=mismatch_count
    #x[1]['filterdueton']= mismatch_count>myargs.n 
    return x[1]

print("Refining alignment statistics with {} processors.".format(myargs.t))

with Pool(myargs.t) as pool:
    new_data = pool.map(new_alignment,data.iterrows())

output = pd.concat([pd.DataFrame(n).T for n in new_data])
output = output[~output['filterdueton']]

if 'ProbeID' in output.columns:
    output_format = ['amplicon_len','HitName','FP_ID','FP_seq','FP_degeneracies','FP_mismatches','RP_ID','RP_seq','RP_degeneracies','RP_mismatches','RevcompRP','ProbeID','Probe_seq','RevcompProbe','Probe_degeneracies','Probe_mismatches','Probe_startOnAmplicon','ProbeStrand','Start','End','Full_Hit_ID','AmpliconSeq']
else:
    output_format = ['amplicon_len','HitName','FP_ID','FP_seq','FP_degeneracies','FP_mismatches','RP_ID','RP_seq','RP_degeneracies','RP_mismatches','RevcompRP','Start','End','Full_Hit_ID','AmpliconSeq']

output = output[(output['FP_mismatches'] <= myargs.m) & (output['RP_mismatches'] <= myargs.m)]
output[output_format].to_csv(out_prefix+'.primer_realigned.amplicons',sep='\t',index=False)
