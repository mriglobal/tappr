import os
import argparse
import pandas as pd
from Bio import SeqIO
from collections import Counter
from Bio import pairwise2
from skbio import DNA
from multiprocessing import Pool

def get_align_stats(x):
    summary_dict = {'index':x[0],'Probe_mismatches':'N/A','ProbeID':'N/A','RevcompProbe':'N/A','Probe_degeneracies':'N/A',
    'Probe_mismatches':'N/A','Probe_startOnAmplicon':'N/A','ProbeStrand':'N/A'}
    probe_id = x[1]['FP_ID'].split('|')[0]
    for pr in probe_dict[probe_id]:
        alignment_minus = pairwise2.align.localms(str(pr),str(DNA(x[1]['AmpliconSeq'],lowercase=True).reverse_complement()),scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
        alignment_plus = pairwise2.align.localms(str(pr),x[1]['AmpliconSeq'],scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
        if max([alignment_plus[0][2],alignment_minus[0][2]]) >= probe_scores[probe_id]:
            summary_dict['ProbeID'] = probe_id
            if alignment_plus[0][2] > alignment_minus[0][2]:
                summary_dict['ProbeStrand'] = 'Plus'
                alignment=alignment_plus
            else:
                summary_dict['ProbeStrand'] = 'Minus'
                alignment=alignment_minus
            #TODO for now degenerate probe sequences are handled, but degenerate target sequences are not
            # msa = TabularMSA([DNA(alignment[0][0][alignment[0][3]:alignment[0][3]+len(pr)]),DNA(alignment[0][1][alignment[0][3]:alignment[0][3]+len(pr)])])
            # for a in msa[1].expand_degenerates():
            #     seq_alignment = local_pairwise_align_ssw(pr,a)
            #     if seq_alignment[1] > probe_max:
            #         best_pr = pr
            #         best_candidate = a
            #         probe_max = seq_alignment[1]
            # best_alignment = TabularMSA([best_pr,best_candidate])
            probe_score = int(((len(pr)*scoring_scheme[0]) - alignment[0][2])/scoring_scheme[0])
            #keeps best scoring alignment so far
            if summary_dict['Probe_mismatches'] == 'N/A' or summary_dict['Probe_mismatches'] > probe_score:
                summary_dict['Probe_mismatches'] = probe_score
                summary_dict['Probe_seq'] = temp_dict[probe_id]
                summary_dict['RevcompProbe'] = str(DNA(temp_dict[probe_id],lowercase=True).reverse_complement())
                summary_dict['Probe_degeneracies'] = Counter(DNA(temp_dict[probe_id],lowercase=True).degenerates())[True]
                summary_dict['Probe_startOnAmplicon'] = alignment[0][3]
        #find best alignment if target sequence has degenerates
    return summary_dict

#TODO right now this will only accept the highest scoring probe alignment per amplicon record no multiple probe hits per amplicon
def get_mux_stats(x):
    summary_dict = {'index':x[0],'Probe_mismatches':'N/A','ProbeID':'N/A','RevcompProbe':'N/A','Probe_degeneracies':'N/A',
    'Probe_mismatches':'N/A','Probe_startOnAmplicon':'N/A','ProbeStrand':'N/A'}
    for probe_id in probe_dict.keys():
        for pr in probe_dict[probe_id]:
            alignment_minus = pairwise2.align.localms(str(pr),str(DNA(x[1]['AmpliconSeq'],lowercase=True).reverse_complement()),scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
            alignment_plus = pairwise2.align.localms(str(pr),x[1]['AmpliconSeq'],scoring_scheme[0],scoring_scheme[1],scoring_scheme[2],scoring_scheme[3],one_alignment_only=True)
            if max([alignment_plus[0][2],alignment_minus[0][2]]) >= probe_scores[probe_id]:
                #summary_dict['ProbeID'] = probe_id
                if alignment_plus[0][2] > alignment_minus[0][2]:
                    summary_dict['ProbeStrand'] = 'Plus'
                    alignment=alignment_plus
                else:
                    summary_dict['ProbeStrand'] = 'Minus'
                    alignment=alignment_minus
                #TODO for now degenerate probe sequences are handled, but degenerate target sequences are not
                # msa = TabularMSA([DNA(alignment[0][0][alignment[0][3]:alignment[0][3]+len(pr)]),DNA(alignment[0][1][alignment[0][3]:alignment[0][3]+len(pr)])])
                # for a in msa[1].expand_degenerates():
                #     seq_alignment = local_pairwise_align_ssw(pr,a)
                #     if seq_alignment[1] > probe_max:
                #         best_pr = pr
                #         best_candidate = a
                #         probe_max = seq_alignment[1]
                # best_alignment = TabularMSA([best_pr,best_candidate])
                probe_score = int(((len(pr)*scoring_scheme[0]) - alignment[0][2])/scoring_scheme[0])
                if summary_dict['Probe_mismatches'] == 'N/A' or summary_dict['Probe_mismatches'] > probe_score:
                    summary_dict['Probe_mismatches'] = probe_score
                    summary_dict['ProbeID'] = probe_id
                    summary_dict['Probe_seq'] = temp_dict[probe_id]
                    summary_dict['RevcompProbe'] = str(DNA(temp_dict[probe_id],lowercase=True).reverse_complement())
                    summary_dict['Probe_degeneracies'] = Counter(DNA(temp_dict[probe_id],lowercase=True).degenerates())[True]
                    summary_dict['Probe_startOnAmplicon'] = alignment[0][3]
            #find best alignment if target sequence has degenerates
    return summary_dict

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Take amplicons file from simulate_PCR results and add alignment data for probe sequences.")

    #command line arguments
    parser.add_argument('--mux',action='store_true',default=False,help="Flag for Multiplex Analysis. This will consider probe mappings with all possible amplicon products.")
    parser.add_argument('-a', required=True,help="Amplicons file from simulate_PCR")
    parser.add_argument('-m',default=3,type=int,help="Maximum number of mismatches allowed for alignment.")
    parser.add_argument('-p', required=True,help="Oligos file containing primer and probe sequences. (Probes indicated with '|IO' tag)")
    parser.add_argument('-o',default='',help="Prefix for output files.")
    parser.add_argument('-t',default=os.cpu_count()-2,type=int,help="Number of threads for multiprocessing. (Default: All)")
    myargs=parser.parse_args()

    # match, mismatch, open, extend
    scoring_scheme = [4,0,-12,-12]
    if not myargs.o:
        out_prefix = os.path.basename(myargs.a.rsplit('.',maxsplit=1)[0])
    else:
        out_prefix = myargs.o

    if 'mux' in myargs.a:
        mux = True
    else:
        mux = False

    oligos = list(SeqIO.parse(myargs.p,'fasta'))

    primers = []
    probes = []
    for o in oligos:
        if o.id.endswith("IO"):
            probes.append(o)
        elif not mux:
            if o.id.endswith("|R") or o.id.endswith("|F"):
                primers.append(o)
            else:
                raise Exception('Oligos File Not Formatted Correctly.')
        else:
            primers.append(o)

    if myargs.t == 1:
        multip = False
    else:
        multip = True
    #data import
    data = pd.read_table(myargs.a,low_memory=False)
    data = data[~data['AmpliconSeq'].isna()]
    data = data[~data['AmpliconSeq'].str.contains("N")]
    #if probe flag specified then it is assumed alignments with probe are the only ones worth considering



    global probe_dict
    probe_dict = {}
    temp_dict ={s.id.split('|')[0]:str(s.seq) for s in probes}
    for k in temp_dict.keys():
        probe_dict[k] = [primer_seq for primer_seq in DNA(temp_dict[k],lowercase=True).expand_degenerates()]

    probe_scores = {key:((len(temp_dict[key])-myargs.m)*scoring_scheme[0]) for key in probe_dict.keys()}

    #multiprocess pool formula
    if multip:
        if myargs.t > data.shape[0]:
            threads = data.shape[0]
            print("scaling down threads to match input data size")
        else:
            threads = myargs.t
        if not myargs.mux:
            with Pool(threads) as p:
                results = p.map(get_align_stats,data.iterrows())
        else:
            with Pool(threads) as p:
                results = p.map(get_mux_stats,data.iterrows())
    else:
        if not myargs.mux:
            results = map(get_align_stats,data.iterrows())
        else:
            results = map(get_mux_stats,data.iterrows())

    results_df = pd.DataFrame(results)
    results_df.set_index('index',inplace=True)
    output = pd.merge(data,results_df,left_index=True,right_index=True)
    #amplicon_len	HitName	FP_ID	FP_seq	FP_degeneracies	FP_mismatches	RP_ID	RP_seq	RP_degeneracies	RP_mismatches	RevcompRP	ProbeID	Probe_seq	RevcompProbe	Probe_degeneracies	Probe_mismatches	Probe_startOnAmplicon	ProbeStrand	Start	End	Full_Hit_ID	SubjectFullLength	AmpliconSeq	Primary_tag	Gene	Product	Protein_id	Note
    output[['amplicon_len','HitName','FP_ID','FP_seq','FP_degeneracies','FP_mismatches','RP_ID','RP_seq','RP_degeneracies','RP_mismatches','RevcompRP','ProbeID','Probe_seq','RevcompProbe','Probe_degeneracies','Probe_mismatches','Probe_startOnAmplicon','ProbeStrand','Start','End','Full_Hit_ID','AmpliconSeq']].to_csv(out_prefix+'.probe_aligned.amplicons',sep='\t',index=False)
    #TODO add ability to output alignments
