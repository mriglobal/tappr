import time
import json
import pandas as pd
import argparse
from MRItaxonomy import accession2taxid as accession
from MRItaxonomy import taxid


def taxidfunc(x):
    if rank=='species':
        return taxid.getnodeatrank(accession.get_taxid(x),rank)
    elif rank=='subspecies':
        if taxid.getnodeatrank(accession.get_taxid(x),rank) == 0:
            return taxid.getnodeatrank(accession.get_taxid(x),'species')
        else:
            return taxid.getnodeatrank(accession.get_taxid(x),rank)

def get_accession(x):
     if "|" in x:
         return x.split("|")[split_col]
     else:
         return x

def makeset(x):
    fp = x['FP_ID'].split('|')[0]
    rp = x['RP_ID'].split('|')[0]
    return tuple(set([fp,rp]))

def makeset_probes(x):
    fp = x['FP_ID'].split('|')[0]
    rp = x['RP_ID'].split('|')[0]
    probe = x['ProbeID'].split("|")[0]
    return tuple(set([fp,rp,probe]))

def package(r,primer):
    '''function for packaging of data produced by simulate_PCR'''
    amplicon = r[1]['AmpliconSeq']
    return_dict = {'primer_pair':r[1]['RP_ID'].split('|')[0],
                'fp_mismatches':r[1]['FP_mismatches'],
                'rp_mismatches':r[1]['RP_mismatches'],
                'name':r[1]['Full_Hit_ID'],
                'primer':r[1]['FP_ID'].split('|')[0],"label":r[1]["label"],
                'AmpliconSeq':amplicon,
                'accession':r[1]['accession'],
                'primer_set':primer}
    if myargs.probe:
        return_dict['ProbeID'] = r[1]['ProbeID']
    return return_dict


parser = argparse.ArgumentParser(description="Taxonomic exclusivity analysis for each primer set.")
parser.add_argument('-a',required=True, help="Amplicons file from simulate PCR.")
parser.add_argument('--subspecies', action='store_true',default=False,help="Flag for subspecies level exclusivity analysis.")
parser.add_argument('-t',default='',help="Optional TSV to assign labels to each record for comparison (Optional. Overrides subspecies flag)")
parser.add_argument('--probe', action='store_true',default=False,help='This assay depends on an interior oligo.')
parser.add_argument('--orgs',nargs='*',required=True)
parser.add_argument('--blacklist',nargs='*',default=[],help="Accession numbers to exclude from consideration.")
parser.add_argument('--min',type=int,default=100,help="Minimum amplicon length to filter from analysis.")
parser.add_argument('--max',type=int,default=250,help="Maximum amplicon length to filter from analysis.")
parser.add_argument('--split_col',type=int,default=1,help="Column to split on to excise accession number from HitName columns (Default: 1)")
parser.add_argument('-o',default=None,help="Output prefix (Default is to inherit from input.)")
myargs=parser.parse_args()



file = myargs.a
min_length = myargs.min
max_length = myargs.max
split_col = myargs.split_col
#TODO Change this behavior to allow assignment of tax rank at run time for exclusivity consideration
if myargs.subspecies:
    rank = 'subspecies'
else:
    rank = 'species'

if myargs.t:
    labels = pd.read_table(myargs.t)
    orgs = list(map(str,myargs.orgs))
    print('Using .tsv for labels')
else:
    orgs = list(map(int,myargs.orgs))
    print('No .tsv specified, assuming --orgs to be taxids')

if myargs.o:
    outprefix = myargs.o
else:
    outprefix = file
df = pd.read_table(file)
#filter amplicon table on min and max amplicon length parameters
df = df[df['amplicon_len'] > min_length]
#restrict analysis to sequences incorporating probes according to simulate_PCR
if myargs.probe:
    df = df[~df['ProbeID'].isna()]

# taxidfunc = lambda x: taxid.getnodeatrank(accession.get_taxid(x),rank)

df['accession'] = [get_accession(h) for h in df['HitName']]
if not myargs.t:
    df['label'] = df.accession.apply(taxidfunc)
else:
    df = pd.merge(df,labels,left_on='accession',right_on='accession')

#Removes accession numbers identified as being problematic
if myargs.blacklist:
    print("Removing accession numbers from consideration: {}:".format(time.asctime()))
    for b in myargs.blacklist:
        print(b)
        df = df[~df['accession'].str.contains(b)]

if myargs.probe:
    df['primer_set'] = df.apply(makeset_probes,axis=1)
else:
    df['primer_set'] = df.apply(makeset,axis=1)

primergroups = df.groupby('primer_set')

primernames = set([i for i in primergroups.groups.keys()])

amplicongroups = {key:[] for key in primernames}
print("Evaluating exclusivity for the following taxids at the {} level {}:".format(rank,time.asctime()))
for t in orgs:
    print(str(t))

#THIS ASSUMES PROBES WILL ONLY INTEGRATE IN THEIR CORRESPONDING PRIMER PRODUCTS
for primer in amplicongroups.keys():
    for row in primergroups.get_group(primer).iterrows():
        amplicongroups[primer].append(package(row,primer))

output = pd.concat([pd.DataFrame(amplicongroups[k]) for k in amplicongroups.keys()],ignore_index=True)
print(output)
#TFilter synthetic construct records e.g. 32630
output = output[output["label"] != 0]
output = output[output["label"] != 32630]
#drop products produced by interactions between multiple primer sets
#output = output[~(output['primer'] != output['primer_pair'])]
#Create column in dataframe that is a SET object of 'primer' and 'primer_pair' to group on

finalgroups = output.groupby('primer_set')
finalprimers = {key:[] for key in orgs}
failures = {key:[] for key in orgs}
print("Comparing amplicons... {}".format(time.asctime()))

#FOR NOW THIS ASSUMES EACH PRIMER SET SHOULD TARGET A SINGLE TAXID, BUT THAT MIGHT NOT ALWAYS BE TRUE
for tid in orgs:
    for hits in finalgroups.groups.keys():
        if not finalgroups.get_group(hits)[finalgroups.get_group(hits)["label"] != tid].empty:
            off_target = finalgroups.get_group(hits)[finalgroups.get_group(hits)["label"] != tid]
            failures[tid].append({'primers':hits,'failure_acc':list(off_target['accession'].unique())})
        else:
            finalprimers[tid].append(hits)
    
print("Writing output files. {}".format(time.asctime()))

for o in orgs:
    pd.DataFrame(failures[o]).to_csv(outprefix + str(o)+"_failure_table.tsv",sep='\t')
with open(outprefix+"_exclusive_primers.json",'w') as outfile:
    json.dump(finalprimers,outfile,indent=0)
output.to_csv(outprefix+"_normalized_table.amplicons",sep='\t',index=False)
