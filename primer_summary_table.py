import pandas as pd
import argparse
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO

def calc_mt(oligo):
    return mt.Tm_NN(oligo, saltcorr=4, Mg=1.5, dNTPs=.6) 

def filter_on_ids(df, id_series):
# Create a mask for FP_ID where the part before the '|' is in the pandas Series
    mask_fp = df['FP_ID'].apply(lambda x: x.split('|')[0] in id_series.values)
    # Create a mask for RP_ID where the part before the '|' is in the pandas Series
    mask_rp = df['RP_ID'].apply(lambda x: x.split('|')[0] in id_series.values)
    
    # Filter the DataFrame where either condition is True
    return df[mask_fp | mask_rp]

parser = argparse.ArgumentParser(description="Read in a primer pair file and produce a .tsv summarizing the details")
parser.add_argument('-o', required=False, default=None, help="Optionally specify output prefix")
parser.add_argument('-a', required=True, help='Probe Aligned Amplicons File')
parser.add_argument('-i', required=False, default=None, type=str, help='optional inclusivity .tsv of hits summary table, containing counts')
parser.add_argument('-e', required=False, default=None, type=str, help='optional exclusivity json output to filter by')
parser.add_argument('-r', required=False, default=None, type=str, help='optional reference genome accession identifier to report positional information')

myargs = parser.parse_args()

probe_aligned_amplicons = pd.read_table(myargs.a)

if myargs.r:
    probe_aligned_amplicons = probe_aligned_amplicons[probe_aligned_amplicons['HitName'].str.contains(myargs.r, na=False, case=False)]

out_df = probe_aligned_amplicons[['FP_ID','FP_seq','RP_ID','RP_seq','ProbeID','Probe_seq','ProbeStrand']]
out_df.dropna(inplace=True)
out_df.drop_duplicates(inplace=True)
out_df = out_df.join(probe_aligned_amplicons[['Start', 'End']])

out_df['FP_mt'] = out_df['FP_seq'].apply(calc_mt)
out_df['RP_mt'] = out_df['RP_seq'].apply(calc_mt)
out_df['Probe_mt'] = out_df['Probe_seq'].apply(calc_mt)
out_df['FP_len'] = out_df['FP_seq'].str.len()
out_df['RP_len'] = out_df['RP_seq'].str.len()
out_df['Probe_len'] = out_df['Probe_seq'].str.len()

out_df['Primer_mt_diff'] = abs(out_df['FP_mt'] - out_df['RP_mt'])
out_df.index = out_df['FP_ID'].apply(lambda x: x.split('|')[0])


if myargs.e:
    exclusive_df = pd.read_json(myargs.e)
    exclusive_series = exclusive_df[exclusive_df.columns[0]].apply(lambda x: x[0])
    out_df = filter_on_ids(out_df, exclusive_series)


if myargs.i:
    inclusive_df = pd.read_table(myargs.i)
    inclusive_df = inclusive_df.T
    inclusive_df.columns = inclusive_df.iloc[0]
    inclusive_df = inclusive_df.drop(inclusive_df.index[0])
    out_df = out_df.merge(inclusive_df, left_index=True, right_index=True, how='outer')


if myargs.o:
    out_df.to_csv(myargs.o, sep='\t', index=False)
else:
    out_df.to_csv(myargs.a+'.stats_summary.tsv', sep='\t', index=False)




























