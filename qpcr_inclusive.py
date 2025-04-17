import pandas as pd
from Bio import SeqIO
import argparse

def get_taxid(x):
	if rank == 'species':
		try:
			return taxid.getnodeatrank(accession2taxid.get_taxid(x),rank)
		except KeyError:
			return 0
	else:
		try:
			return accession2taxid.get_taxid(x)
		except KeyError:
			return 0

def get_accession(x):
     if "|" in x:
         return x.split("|")[split_col]
     else:
         return x
         
def get_org(x):
	try:
		return taxid2name.get_name(x)
	except KeyError:
		return "Unknown"
    
def get_coverage(x):
    return len(pd.Series.unique(x))

parser = argparse.ArgumentParser(description="Taxonomic inclusivity analysis for each primer set.")
parser.add_argument('-a',required=True, help="Amplicons file from simulate PCR.")
parser.add_argument('--subspecies', action='store_true',default=False,help="Flag for subspecies level exclusivity analysis.")
parser.add_argument('--probe', action='store_true',default=False,help='This assay depends on an interior oligo.')
#Could add argument to only evaluate specific taxids but currently will evaluate all taxids in input seqs
#parser.add_argument('--orgs',nargs='*',required=True)
parser.add_argument('-m',type=float,default=None,help='Fraction of matches as a percentage of oligo length to cutoff alignment at. Will filter results table.')
parser.add_argument('-f',required=True,help="Fasta file of sequences evaluated for inclusivity.")
parser.add_argument('--name', action='store_true',default=False,help='Return results summary by organism name rather than taxid. (Default: False)')
parser.add_argument('-t',default=None,help='Tab delimited metadata table for custom inclusivity labeling instead of taxonomy.')
parser.add_argument('--by_primer_pair',action='store_true',default=False,help='Split out coverage results by primer pairs observed in amplicons table.')
parser.add_argument('--blacklist',nargs='*',default=[],help="Accession numbers to exclude from consideration.")
parser.add_argument('--min',type=int,default=40,help="Minimum amplicon length to filter from analysis.")
#parser.add_argument('--max',type=int,default=250,help="Maximum amplicon length to filter from analysis.")
parser.add_argument('--split_col',type=int,default=1,help="Column to split on to excise accession number from HitName columns (Default: 1 for genbank sequences)")
parser.add_argument('-o',default=None,help="Output prefix (Default is to inherit from input.)")
myargs=parser.parse_args()
if myargs.t and myargs.name:
	parser.error("Organism names can't be returned if custom inclusivity labels have been provided.")
if myargs.subspecies and myargs.t:
	parser.error("Subspecies flag has no meaning if custom inclusivity labels have been provided.")

if not myargs.t:
    from MRItaxonomy import accession2taxid
    from MRItaxonomy import taxid
if myargs.name:
	from MRItaxonomy import taxid2name


file = myargs.a
min_length = myargs.min
#max_length = myargs.max
#TODO: add argument for custom groupings table to be provided for inclusivity analysis outside of taxonomy
split_col = myargs.split_col

if myargs.subspecies:
    rank = 'subspecies'
else:
    rank = 'species'

if myargs.o:
    outprefix = myargs.o
else:
    outprefix = file

df = pd.read_table(file,low_memory=False)




if not myargs.t:
	meta = {seq.id:get_taxid(seq.id) for seq in SeqIO.parse(myargs.f,'fasta')}
	meta_table = pd.DataFrame(list(meta.items()),columns=['accession','label'])

	if myargs.name:
		meta_table['species_name'] = [get_org(t) for t in meta_table['label']]

	if myargs.name:
		summary_grouping = 'species_name'
	else:
		summary_grouping = 'label'

	meta_dict = meta_table[summary_grouping].value_counts().to_dict()

	meta_sets = meta_table.groupby(summary_grouping)['accession'].apply(set).to_dict()

else:
	meta_table = pd.read_table(myargs.t)
	meta_table = meta_table[meta_table['accession'].isin(set([seq.id for seq in SeqIO.parse(myargs.f,'fasta')]))]
	summary_grouping = 'label'
	meta_dict = meta_table[summary_grouping].value_counts().to_dict()
	meta_sets = meta_table.groupby(summary_grouping)['accession'].apply(set).to_dict()

df['accession'] = [get_accession(h) for h in df['HitName']]

#TODO: Implement primer grouping for multiplex output
df['primer_group'] = df['FP_ID'].str.split("|",expand=True)[0]

#filter amplicon table on min and max amplicon length parameters
df = df[df['amplicon_len'] > min_length]

df = pd.merge(df,meta_table,left_on='accession',right_on='accession')

df = df[df['label'] != 0]

#restrict analysis to sequences incorporating probes according to simulate_PCR
if myargs.probe:
	no_probe = df[df['ProbeID'].isna()]
	hits = df[~df['ProbeID'].isna()]
	if myargs.m:
		probe_mismatch = len(hits['Probe_seq'].iloc[0]) - round(len(hits['Probe_seq'].iloc[0]) * myargs.m)
		hits = hits[hits['Probe_mismatches'] <= probe_mismatch]
	no_probe.to_csv(outprefix +"_amplicon_but_no_probe.tsv",sep='\t',index=False)
else:
	hits = df.copy()

if myargs.m:
	f_mismatch = len(hits['FP_seq'].iloc[0]) - round(len(hits['FP_seq'].iloc[0]) * myargs.m)
	r_mismatch = len(hits['RP_seq'].iloc[0]) - round(len(hits['RP_seq'].iloc[0]) * myargs.m)
	#TODO: This filter statement needs to be rewritten to consider complement strand records
	hits = hits[(hits['FP_mismatches'] <= f_mismatch) & (hits['RP_mismatches'] <= r_mismatch)]

hit_groups = set(hits[summary_grouping])

hits_summary = hits[[summary_grouping,'primer_group']].pivot_table(index='primer_group',columns=summary_grouping,aggfunc=len).fillna(0).T



if myargs.by_primer_pair:
	panel_coverage = hits.pivot_table(index=[summary_grouping,'primer_group'],values='accession',aggfunc=get_coverage)
	panel_coverage['coverage'] = [panel_coverage.loc[s]['accession']/meta_dict[s[0]] for s in panel_coverage.index]
else:
	panel_coverage = hits.pivot_table(index=summary_grouping,values='accession',aggfunc=get_coverage)
	panel_coverage['coverage'] = [panel_coverage.loc[s]['accession']/meta_dict[s] for s in panel_coverage.index]

no_coverage = pd.DataFrame({'accession':[0 for s in meta_dict.keys() if s not in hit_groups]},index=[s for s in meta_dict.keys() if s not in hit_groups])

panel_coverage = pd.concat([panel_coverage,no_coverage])
#print(panel_coverage)

if panel_coverage.empty:
	print('empty coverage!')

panel_coverage.columns = ['Number of Hits','Coverage']

hits_summary.to_csv(outprefix + "_hits_summary_table.tsv",sep='\t')
panel_coverage.to_csv(outprefix + "_panel_coverage.tsv",sep='\t')

for p in panel_coverage.index:
	if myargs.by_primer_pair:
		g = p[0]
	else:
		g = p
	if len(meta_sets[g].difference(set(hits[hits[summary_grouping] == g]['accession']))) > 0:
		with open(outprefix + "_" + str(g).replace(' ','_') + "_missing.txt", 'w') as outfile:
			for l in meta_sets[g].difference(set(hits[hits[summary_grouping] == g]['accession'])):
				outfile.write(l+'\n')
	
