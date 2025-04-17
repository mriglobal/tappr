import os
import pybedtools
from Bio import SeqIO
import pysam
import argparse
from pybedtools import BedTool

def sam_to_bed(sam_file, bed_file):
    # First convert SAM to BAM using pysam
    bam_file_path = sam_file.replace('.sam', '.bam')

    # Open the SAM file and convert it to BAM
    with pysam.AlignmentFile(sam_file, "r") as samfile:
        with pysam.AlignmentFile(bam_file_path, "wb", template=samfile) as bamfile:
            for read in samfile:
                bamfile.write(read)

    # Sort the BAM file
    sorted_bam_file_path = bam_file_path.replace('.bam', '.sorted.bam')
    pysam.sort("-o", sorted_bam_file_path, bam_file_path)

    # Now use pybedtools to convert sorted BAM to BED
    bed = BedTool(sorted_bam_file_path).bam_to_bed().merge()
    bed.saveas(bed_file)

    # Optionally, clean up intermediate BAM files
    os.remove(bam_file_path)
    os.remove(sorted_bam_file_path)


parser = argparse.ArgumentParser(description="Return complement Marker locations based on a sam or bed file input for primer regions.")

parser.add_argument("-i",required=True,help="Primer regions bed or sam file.")
parser.add_argument("-f",required=True,help="Fasta file for reference sequence.")
parser.add_argument("-o",default=None,help="Output file prefix. Default: Inherit from bed file.")

myargs = parser.parse_args()

if not myargs.o:
    outprefix = myargs.i.rsplit(".",maxsplit=1)[0]
else:
    outprefix = myargs.o

seq_stats = {s.id:len(s) for s in SeqIO.parse(myargs.f,'fasta')}

for s in seq_stats:
    seq_stats[s] = tuple([0,seq_stats[s]])


genome_file = pybedtools.helpers.chromsizes_to_file(seq_stats)

if '.sam' in myargs.i:
    sam_to_bed(myargs.i, myargs.i.rsplit(".",maxsplit=1)[0]+'.bed')
    primer_bed = pybedtools.bedtool.BedTool(myargs.i.rsplit(".",maxsplit=1)[0]+'.bed').sort()

else:
    primer_bed = pybedtools.bedtool.BedTool(myargs.i).sort()

marker_bed = primer_bed.complement(g=genome_file)

marker_bed.saveas(outprefix+'_complement_markers.bed')







