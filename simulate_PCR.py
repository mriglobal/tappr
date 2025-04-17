import os
import pandas as pd
import subprocess
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import chain
from io import StringIO
from collections import Counter
from skbio import DNA
import argparse

def check_mismatch(x):
	if oligo_mismatch_dict[x['qseqid']] >= x['n_mismatches']:
		return True
	else:
		return False

def strand_eval(x):
	'''Labeling function to add strand column to results DataFrame'''
	if x['sstart'] > x['send']:
		return 'minus'
	else:
		return 'plus'

def header_correction(x):
	if x == 'sseqid':
		return 'HitName'
	elif x == 'qseqid':
		return 'FP_ID'
	elif x == 'slen':
		return 'SubjectFullLength'
	elif x == 'stitle':
		return 'Full_Hit_ID'
	else:
		return x

#scan_hits
#this covers the mux situation
def scan_hits(x):
	hits = []
	forward = x.groupby('strand').get_group('plus')
	reverse = x.groupby('strand').get_group('minus')
	#considers all combinations of forward and reverse
	for f in forward.iterrows():
		for r in reverse.iterrows():
			if (r[1]['sstart'] - f[1]['sstart'] <= myargs.max) and (r[1]['sstart'] - f[1]['sstart'] >= myargs.min):
				hits.append({'amplicon_len':(r[1]['sstart'] - f[1]['sstart'])+1,'HitName':f[1]['sseqid'],'FP_ID':f[1]['qseqid'],'FP_seq':str(oligo_dict[f[1]['qseqid']]),'FP_degeneracies':Counter(DNA(str(oligo_dict[f[1]['qseqid']]),lowercase=True).degenerates())[True],
			'FP_mismatches': f[1]['n_mismatches'],'RP_ID':r[1]['qseqid'],'RP_seq':str(oligo_dict[r[1]['qseqid']]),'RP_degeneracies':Counter(DNA(str(oligo_dict[r[1]['qseqid']]),lowercase=True).degenerates())[True],
			'RP_mismatches': r[1]['n_mismatches'],'RevcompRP':str(oligo_dict[r[1]['qseqid']].reverse_complement()),'Start':f[1]['sstart'],'End':r[1]['sstart'],'Full_Hit_ID':f[1]['stitle']})
	return hits

parser = argparse.ArgumentParser(description="Simulate_PCR redux using python and a few performance improvements.")
parser.add_argument('-p',required=True, help="Primer Fasta File")
parser.add_argument('-r',default=None, help="Blastout results table.")
parser.add_argument('--mismatches',default=3.0,type=float,help="If integer, number of mismathces allowed in primers for match. If float < 1, fraction of mismatches allowed based on primer length. (Default:3)")
parser.add_argument('--min',default=40,type=int,help="Minimum product length for consideration. (Defaulte:40)")
parser.add_argument('--max',default=1000,type=int,help="Maximum produce length for consideration. (Default:1000)")
parser.add_argument('--mux',action='store_true',default=False,help="Multiplex Flag. (Default: False)")
parser.add_argument('-t',default=os.cpu_count()-2,help="Number of threads for BLAST. (Default: All)")
parser.add_argument('--max_target_seqs',default=1000,type=int,help="Maximum number of target seqs to return from BLAST. (Default:1000)")
parser.add_argument('--evalue',default=10000,type=float,help="Evalue cutoff for BLAST results reporting.")
parser.add_argument('--word_size',default=4,type=int,help="Word size for BLAST. (Default: 4)")
#parser.add_argument('--3prime',default=3,help="Number of perfect matches required in the 3' end of the primer. (Default:3)")
parser.add_argument('--db',required=True,help="Path to blast database.")
parser.add_argument('--outname', default=None, help="Specify output name")

#myargs = parser.parse_args(['-p','/home/pdavis/Documents/CDC/stx_gene/stx_panel_primers.fasta','-r','/home/pdavis/Documents/CDC/stx_gene/stx_panel_primers.fasta.mux.O157_genoms.blastout','--mismatches','3','--extract_amp','--db', '/home/pdavis/Documents/CDC/stx_gene/stx_genes'])

myargs = parser.parse_args()

if myargs.mux:
	file_name_insert = 'mux'
else:
	file_name_insert = 'pair'

db_base_name = os.path.basename(myargs.db).split(".")[0]
out_prefix = '.'.join([os.path.basename(myargs.p),file_name_insert,db_base_name])

class FileFormatError(Exception):
	pass

oligos = list(SeqIO.parse(myargs.p,'fasta'))
oligo_dict = {o.id:o.seq for o in oligos}
oligo_lens = {s.id:len(s.seq) for s in oligos}

if myargs.mismatches.is_integer():
	oligo_mismatch_dict = {o.id:myargs.mismatches for o in oligos}
elif myargs.mismatches < 1.0:
	oligo_mismatch_dict = {o.id:myargs.mismatches*len(o) for o in oligos}


primers = []
probes = []
for o in oligos:
	if o.id.endswith("IO"):
		probes.append(o)
	elif not myargs.mux:
		if o.id.endswith("|R") or o.id.endswith("|F"):
			primers.append(o)
		else:
			raise FileFormatError('Oligos File Not Formatted Correctly.')
	else:
		primers.append(o)

#blast_primers
#shell=True necessary?
if not myargs.r:
	print("Blasting Oligos.")
	subprocess.run(['blastn', '-num_threads',str(myargs.t) , '-db', myargs.db, '-query', myargs.p, '-task', 'blastn-short', '-dust', 'no', '-evalue', str(myargs.evalue), '-penalty', '-2', '-reward', '3', '-gapopen','8', '-gapextend', '5', '-max_target_seqs',
	str(myargs.max_target_seqs), '-word_size', str(myargs.word_size), '-out', out_prefix+'.blastout', '-outfmt', '6 qseqid sseqid nident qstart qend sstart send sseq slen stitle'])


#read_blastout
print("Reading in Blast Results.")
if myargs.r:
	results = pd.read_table(myargs.r,names=['qseqid','sseqid','nident','qstart','qend','sstart','send','sseq','slen','stitle'])
else:
	results = pd.read_table(out_prefix+'.blastout',names=['qseqid','sseqid','nident','qstart','qend','sstart','send','sseq','slen','stitle'])

#database record sniffer
if "|" in results['sseqid'].loc[0]:
	nt_flag = True
else:
	nt_flag = False

#TODO: Handle degeneracies
results['qseq_len'] = [oligo_lens[q] for q in results['qseqid']]
results['n_mismatches'] = (results['qseq_len'] - (results['qend'] - results['qstart'] +1)) + ((results['qend'] - results['qstart'] +1) - results['nident'])
#filter by mismatch criteria
filtered = results[results.apply(check_mismatch,axis=1)].copy()
#data strucure for .amplicons file:
'''['amplicon_len', 'HitName', 'FP_ID', 'FP_seq', 'FP_degeneracies',
       'FP_mismatches', 'RP_ID', 'RP_seq', 'RP_degeneracies', 'RP_mismatches',
       'RevcompRP', 'ProbeID', 'Probe_seq', 'RevcompProbe',
       'Probe_degeneracies', 'Probe_mismatches', 'Probe_startOnAmplicon',
       'ProbeStrand', 'Start', 'End', 'Full_Hit_ID', 'SubjectFullLength',
       'AmpliconSeq', 'Primary_tag', 'Gene', 'Product', 'Protein_id', 'Note']'''

#results = pd.merge(results,primer_stats)
filtered['strand'] = filtered.apply(strand_eval,axis=1)
#only consider hits that have oligos mapping on both strands
if set(filtered['strand'].unique()) != set(['minus','plus']):
	raise Exception("No amplicon predicted.")

if myargs.mux:
	try:
		candidate_set = filtered.groupby('strand')['sseqid'].apply(set)['minus'].intersection(filtered.groupby('strand')['sseqid'].apply(set)['plus'])
	except KeyError:
		raise Exception("No candidates found.")
	filtered = filtered[filtered['sseqid'].isin(candidate_set)]

	print("Scanning for candidate Amplicons.")
#https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby
	output = [scan_hits(g) for _,g in filtered.groupby('sseqid')]

#paired results
else:
	filtered['primer_groups'] = filtered['qseqid'].str.split("|",expand=True)[0]
	filtered_groups = filtered.groupby('primer_groups')
	output = []
	for primer_group in filtered_groups.groups.keys():
		print("Scanning for candidate Amplicons for {} primer pair.".format(primer_group))
		try:
			candidate_set = filtered_groups.get_group(primer_group).groupby('strand')['sseqid'].apply(set)['minus'].intersection(filtered_groups.get_group(primer_group).groupby('strand')['sseqid'].apply(set)['plus'])
		except KeyError:
			continue
		filtered_candidates = filtered_groups.get_group(primer_group)[filtered_groups.get_group(primer_group)['sseqid'].isin(candidate_set)]
		output.extend([scan_hits(g) for _,g in filtered_candidates.groupby('sseqid')])

amplicons = pd.DataFrame(list(chain(*output)))

if not amplicons.empty:
	hitlist = amplicons[['HitName','Start','End']].copy()
else:
	raise Exception("No Predicted Amplicons.")
if nt_flag:
	accession = re.compile(r"[A-Z]{1,2}\_{0,1}[A-Z]*[0-9]+\.[0-9]+")
	hitlist['HitName'] = [re.search(accession,a).group() for a in hitlist['HitName']]

hitlist['range'] = hitlist['Start'].astype(str) + '-' + hitlist['End'].astype(str)
hitlist['strand'] = 'plus'
entry_plus = StringIO(hitlist[['HitName','range','strand']].to_string(index=False,header=False))
#extract_amplicons
#blastdbcmd -db <db name> -entry_batch <input_file> ([seq id|range|strand|mask_algo_id] separated by spaces) -outfmt %s
#shell=True may or may not be necessary still testing
proc_plus = subprocess.Popen(["blastdbcmd","-db",myargs.db,"-entry_batch","-"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,encoding='ascii')
output, output_err = proc_plus.communicate(entry_plus.read())
output_handle = StringIO(output)
plus_amplicons = list(SeqIO.parse(output_handle,'fasta'))
print("Writing Amplicons Table.")
amplicons['AmpliconSeq'] = [str(a.seq) for a in plus_amplicons]
amplicons.to_csv(out_prefix+".amplicons",sep='\t',index=False)

amplicon_seqs = []
for r in amplicons.iterrows():
    amplicon_seqs.append(SeqRecord(Seq(r[1]['AmpliconSeq']),id=r[1]['HitName'],description=' '.join([r[1]['FP_ID'],r[1]['RP_ID']])))
print("Writing Amplicons Fasta File.")
if not myargs.outname:
    SeqIO.write(amplicon_seqs,out_prefix+'.amplicons.fasta','fasta')
else:
    SeqIO.write(amplicon_seqs,myargs.outname+'.amplicons.fasta','fasta')