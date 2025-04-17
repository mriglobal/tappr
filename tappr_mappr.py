import pysam
from Bio import SeqIO
import pandas as pd
import argparse
import multiprocessing
import pybedtools
import os

def process_reference(r, k_size, kmer_set, kmer_ids):
    alignments = []
    print(f"Finding Alignments for {r.id}")
    # Process the direct sequence
    for i in range(len(r) - k_size + 1):
        candidate_kmer = str(r[i:i + k_size].seq)
        if candidate_kmer in kmer_set:
            alignments.append({
                'query_id': kmer_ids[candidate_kmer],
                'query_seq': candidate_kmer,
                'flag': 0,
                'reference_id': r.id,
                'reference_start': i + 1,  # Adjust for 1-based index
                'qual': 42,
                'cigar': [(0, k_size)],
                'next_reference': 0,
                'next_ref_start': 0,
                'template_len': 0,
                'query_qual': 'I' * k_size
            })

    # Process the reverse complement
    rev_comp = r.reverse_complement()
    for i in range(len(rev_comp) - k_size + 1):
        candidate_kmer = str(rev_comp[i:i + k_size].seq)
        if candidate_kmer in kmer_set:
            alignments.append({
                'query_id': kmer_ids[candidate_kmer],
                'query_seq': candidate_kmer,
                'flag': 16,
                'reference_id': r.id,
                'reference_start': len(r) - (i + k_size - 1),  # Correcting the index for reverse
                'qual': 42,
                'cigar': [(0, k_size)],
                'next_reference': 0,
                'next_ref_start': 0,
                'template_len': 0,
                'query_qual': 'I' * k_size
            })

    return alignments


parser = argparse.ArgumentParser(description="Script for aligning k-mers to a single reference genome..")
    
#command line arguments
parser.add_argument('--kmers', required=True, type=str, help="K-mer fasta file for mapping. All K-mers are assumed to be the same length.")
parser.add_argument('-r', required=True ,help="Reference genome fasta") #reference genome
parser.add_argument('--bed',action="store_true",default=False,help="Flag for ben output instead of SAM.")
parser.add_argument('-t',type=int,default=os.cpu_count()-2,help='Number of threads to use for concurrent kmer counting.')
parser.add_argument('-o',type=str,default=None,help='Output file prefix. Default it to inherit from reference file.')

myargs = parser.parse_args()

reference = list(SeqIO.parse(myargs.r,'fasta'))
kmers = list(SeqIO.parse(myargs.kmers,'fasta'))

kmer_set = set([str(k.seq) for k in kmers])

kmer_ids = {str(k.seq):k.id for k in kmers}

print("Total kmers to map: {}".format(len(kmer_set)))

alignments = []

#sniff ksize
k_size = len(kmers[0])


if not myargs.o:
    outfile = myargs.r.rsplit('.',maxsplit=1)[0] + '.{}mers'.format(k_size)
else:
    outfile = myargs.o

#prepare data structures

#if sam output

if not myargs.bed:
    sam_header = {'HD':{'VN':'1.0'},
                  'SQ':[{'LN':len(r),'SN':r.id} for r in reference]}
    ref_index = {r.id:i for i,r in enumerate(reference)}

if myargs.t > len(reference):
    threads = len(reference)
    print("Using "+str(threads)+" threads due to num records in reference")
else:
    threads=myargs.t

with multiprocessing.Pool(processes=threads) as pool:
        results = pool.starmap(process_reference, [(r, k_size, kmer_set, kmer_ids) for r in reference])

# Flatten the list of lists
alignments = [item for sublist in results for item in sublist]
print("Total alignments found: ", len(alignments))

if not myargs.bed:
    print("Writing SAM alignment file.")
    with pysam.AlignmentFile(outfile+'.sam','wh',header=sam_header) as outf:
        for align in alignments:
            a = pysam.AlignedSegment()
            a.query_name = align['query_id']
            a.query_sequence=align['query_seq']
            a.flag = align['flag']
            a.reference_id = ref_index[align['reference_id']]
            a.reference_start = align['reference_start']
            a.mapping_quality = align['qual']
            a.cigar = align['cigar']
            a.next_reference_id = align['next_reference']
            a.next_reference_start=align['next_ref_start']
            a.template_length=align['template_len']
            a.query_qualities = pysam.qualitystring_to_array(align['query_qual'])
            outf.write(a)

else:
    print("Writing BED file.")
    #BED files are 0 indexed
    outf = pd.DataFrame({'chr':[align['reference_id'] for align in alignments], 'start':[align['reference_start']-1 for align in alignments]})
    outf['end'] = outf['start'] + k_size
    bed = pybedtools.BedTool.from_dataframe(outf)
    merged_bed = bed.sort().merge()
    merged_bed.saveas(outfile + ".bed")
    #outf.to_csv(outfile+".bed",sep='\t',index=False,header=False)

