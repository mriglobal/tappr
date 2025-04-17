import time
import json
import pandas as pd
import argparse
from skbio import DNA

# taxidfunc = lambda x: taxid.getnodeatrank(accession.get_taxid(x),rank)
def taxidfunc(x):
    if rank=='species':
        return taxid.getnodeatrank(accession.get_taxid(x),rank)
    elif rank=='subspecies':
        if taxid.getnodeatrank(accession.get_taxid(x),rank) == 0:
            return taxid.getnodeatrank(accession.get_taxid(x),'species')
        else:
            return taxid.getnodeatrank(accession.get_taxid(x),rank)

def labelfunc(x):
    try:
        return label_dict[x]
    except KeyError:
        return "Off Target"

def makeset(x):
    return tuple(set([x['primer'],x['primer_pair']]))

def makeset_spe(x):
    return tuple(set([x['primer'],x['primer']]))

def package(r,primer,trim=False,rev=False):
    '''function for packaging of data produced by simulate_PCR'''
    if not trim:
        if rev:
            amplicon = str(DNA(r[1]['AmpliconSeq']).reverse_complement())
        else:
            amplicon = r[1]['AmpliconSeq']
        return_dict = {'primer_pair':r[1]['RP_ID'].split('|')[0],
                'fp_mismatches':r[1]['FP_mismatches'],
                'rp_mismatches':r[1]['RP_mismatches'],
                'name':r[1]['Full_Hit_ID'],
                'primer':primer,"label":r[1]["label"],
                'AmpliconSeq':amplicon,
                'accession_num':r[1]['accession_num']}
        return return_dict
    else:
        if rev:
            fragment = str(DNA(r[1]['AmpliconSeq']).reverse_complement())[:read_length]
        else:
            fragment = r[1]['AmpliconSeq'][:read_length]
        return_dict = {'primer_pair':r[1]['RP_ID'].split('|')[0],
                'fp_mismatches':r[1]['FP_mismatches'],
                'rp_mismatches':r[1]['RP_mismatches'],
                'name':r[1]['Full_Hit_ID'],
                'primer':primer,"label":r[1]["label"],
                'AmpliconSeq':fragment,
                'accession_num':r[1]['accession_num']}
        return return_dict

def package_spe(r,primer,trim=False):
    '''function for packaging of data produced by simulate_SPE'''
    if not trim:
        return_dict = {'fp_mismatches':r[1]['FP_mismatches'],
                'name':r[1]['Full_Hit_ID'],
                'primer':primer,"label":r[1]["label"],
                'AmpliconSeq':str(DNA(r[1]['AmpliconSeq'])),
                'accession_num':r[1]['accession_num']}
        return return_dict

parser = argparse.ArgumentParser(description="Taxonomic exclusivity analysis for each primer set.")
parser.add_argument('-a',required=True, help="Amplicons file from simulate PCR.")
parser.add_argument('--read_length',type=int,default=150,help="Expected read length for sequencing experiment. (0 will indicate variable read length i.e. minion)")
parser.add_argument('-t',default=None,help="Provide custom labels for ingroup identifiers in exclusivity analysis. Labels specified in --targets must match label in provided table. Not recommended for use with NT database.")
parser.add_argument('--subspecies', action='store_true',default=False,help="Flag for subspecies level exclusivity analysis.")
parser.add_argument('--spe', action='store_true',default=False,help="Flag for single primer extension amplicon analysis.")
parser.add_argument('--targets',nargs='*',required=True)
parser.add_argument('--blacklist',nargs='*',default=[],help="Accession numbers to exclude from consideration.")
parser.add_argument('--min',type=int,default=100,help="Minimum amplicon length to filter from analysis.")
parser.add_argument('--max',type=int,default=250,help="Maximum amplicon length to filter from analysis.")
parser.add_argument('--split_col',type=int,default=0,help="Column to split on to excise accession number from HitName columns (Default: 0 for genbank sequences)")
myargs=parser.parse_args()

if not myargs.t:
    from MRItaxonomy import accession2taxid as accession
    from MRItaxonomy import taxid
    orgs = list(map(int,myargs.targets))
else:
    labels = pd.read_table(myargs.t,index_col=0)
    labels.columns = ['label']
    label_dict = labels.to_dict()['label']
    orgs = myargs.targets

file = myargs.a
min_length = myargs.min
max_length = myargs.max
read_length = myargs.read_length

if myargs.subspecies:
    rank = 'subspecies'
else:
    rank = 'species'

df = pd.read_table(file)
#filter amplicon table on min and max amplicon length parameters
df = df[df['amplicon_len'] > min_length]

#TODO: Figure out column splitting strategy
split_col = myargs.split_col

df['accession_num'] = df.HitName.str.split('|',expand=True)[split_col]
if not myargs.t:
    df['label'] = df['accession_num'].apply(taxidfunc)
else:
    df['label'] = df['accession_num'].apply(labelfunc)

#Removes accession numbers identified as being problematic
if myargs.blacklist:
    print("Removing accession numbers from consideration: {}:".format(time.asctime()))
    for b in myargs.blacklist:
        print(b)
        df = df[~df['accession_num'].str.contains(b)]
primergroups = df.groupby('FP_ID')

primernames = set([i.split(sep='|')[0] for i in primergroups.groups.keys()])

amplicongroups = {key:[] for key in primernames}
print("Evaluating exclusivity for the following Labels at the {} level {}:".format(rank,time.asctime()))
for t in orgs:
    print(str(t))

if not myargs.spe:
    print("Normalizing amplicon strand. {}".format(time.asctime()))
    #groups forward and reverse primer products together and chops up those outside of max length
    #the products where the reverse read is in the forward primer position is reverse complemented
    for primer in amplicongroups.keys():
        if primer+'|F' in primergroups.groups.keys():
            for row in primergroups.get_group(primer+"|F").iterrows():
                if read_length and len(row[1]['AmpliconSeq']) > max_length:
                    #read length is non-zero and the amplicon is greater than maximum merge length
                    #therefore, trim amplicon into its 5 primer and 3 prime ends
                    amplicongroups[primer].append(package(row,primer,True))
                    amplicongroups[primer].append(package(row,primer,True,True))
                else:
                    amplicongroups[primer].append(package(row,primer))
                    
        if primer+'|R' in primergroups.groups.keys():
            for row in primergroups.get_group(primer+"|R").iterrows():
                if read_length and len(row[1]['AmpliconSeq']) > max_length:
                    amplicongroups[primer].append(package(row,primer,True))
                    amplicongroups[primer].append(package(row,primer,True,True))
                else:
                    amplicongroups[primer].append(package(row,primer,False,True))
else:
    print("Packaging SPE amplicon products. {}".format(time.asctime()))
    for primer in amplicongroups.keys():
        for row in primergroups.get_group(primer+"|F").iterrows():
            amplicongroups[primer].append(package_spe(row,primer))
            
output = pd.concat([pd.DataFrame(amplicongroups[k]) for k in amplicongroups.keys()],ignore_index=True)
print(output)
#TFilter synthetic construct records e.g. 32630
if not myargs.t:
    output = output[output["label"] != 0]
    output = output[output["label"] != 32630]
#drop products produced by interactions between multiple primer sets
#output = output[~(output['primer'] != output['primer_pair'])]
#Create column in dataframe that is a SET object of 'primer' and 'primer_pair' to group on
if not myargs.spe:
    output['primer_set'] = output.apply(makeset,axis=1)
else:
    output['primer_set'] = output.apply(makeset_spe,axis=1)
finalgroups = output.groupby('primer_set')
finalprimers = {key:[] for key in orgs}
failures = {key:[] for key in orgs}
print("Comparing amplicons... {}".format(time.asctime()))
for tid in orgs:
    for hits in finalgroups.groups.keys():
        query = finalgroups.get_group(hits)[finalgroups.get_group(hits)["label"] == tid]
        db = finalgroups.get_group(hits)[finalgroups.get_group(hits)["label"] != tid]
        hit_results = []
        if not query.empty:
            if db.empty:
                finalprimers[tid].append(hits)
            else:
                for i in query.iterrows():
                    #TODO: Replace this rediculous code with a straight string comparison
                    results = [bool(DNA(i[1].AmpliconSeq).kmer_frequencies(3) == DNA(r[1].AmpliconSeq).kmer_frequencies(3)) for r in db.iterrows()]
                    hit_results+=results
                    if True in results:
                        failures[tid].append({'primers':hits,'q_accession':i[1].accession_num,'failure_acc':db[results].name.unique()})
                if len(hit_results) > 1 and True not in hit_results:
                    finalprimers[tid].append(hits)
    
print("Writing output files. {}".format(time.asctime()))
for o in orgs:
    pd.DataFrame(failures[o]).to_csv(str(o)+"_failure_table.tsv",sep='\t')
with open(file+"_exclusive_primers.json",'w') as outfile:
    json.dump(finalprimers,outfile,indent=0)
output.to_csv(file+"_normalized_table.amplicons",sep='\t',index=False)

