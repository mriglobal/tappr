#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 --exc EXC_DIR --seqs SEQS_DIR \
          --ref REF  \
          -k KMER_LENGTH --ampmax MAX --target TARGET_NAME 

  --exc                exclusive sequences directory
  --seqs               inclusive sequences directory
  --ref                Reference Sequence Fasta File
  -k KMER_LENGTH       k-mer length for all k-mer based scripts
  --ampmax MAX         maximum amplicon size for make_primers.py
  --target TARGET_NAME Target name for db naming and inclusivity tagging
EOF
  exit 1
}

# --- parse CLI ---
EXC_DIR=
SEQS_DIR=
REF=
KMER_LENGTH=
AMPMAX=
TARGET_NAME=

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exc) EXC_DIR="$2"; shift 2;;
    --seqs)             SEQS_DIR="$2";             shift 2;;
    --ref)              REF="$2";              shift 2;;
    -k)                 KMER_LENGTH="$2";      shift 2;;
    --ampmax)           AMPMAX="$2";           shift 2;;
    --target)           TARGET_NAME="$2";           shift 2;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

# --- validate ---
if [[ -z $EXC_DIR || -z $SEQS_DIR || -z $REF || -z $KMER_LENGTH || -z $AMPMAX || -z $TARGET_NAME ]]; then
  usage
fi

REF_PREFIX="${REF%%.*}"
REF_ID=$(grep -m1 '^>' $REF | cut -d' ' -f1 | sed 's/^>//')

# --- path to tappr repository ---
TAPPR_DIR="${HOME}/tappr"

echo "Running pipeline with:
  exc=$EXC_DIR
  seqs=$SEQS_DIR
  ref=$REF
  kmer_length=$KMER_LENGTH
  ampmax=$AMPMAX
  target_name=$TARGET_NAME
"

# 1. get inner kmers
python "${TAPPR_DIR}/kmercountinner.py" \
  -k "$KMER_LENGTH" \
  --seqs "$SEQS_DIR" \
  -r "$REF" 

# 2. tappr mapper
python "${TAPPR_DIR}/tappr_mappr.py" \
  -r "$REF" \
  --kmers "${REF_PREFIX}.${KMER_LENGTH}mer_0.fasta" \


# 3. get exclusive kmers
python "${TAPPR_DIR}/kmercountouter_set.py" \
  --seqs "$EXC_DIR" \
  -k "$KMER_LENGTH" \
  -c "4"

# 4. get kmer set difference
python "${TAPPR_DIR}/kmer_set_difference.py" \
  -A "${REF_PREFIX}.${KMER_LENGTH}mer_0.fasta" \
  -B "${EXC_DIR}${KMER_LENGTH}Kmer_set.pickle"

# 5. map markers
python "${TAPPR_DIR}/tappr_mappr.py" \
  --bed \
  -r "$REF" \
  --kmers "${REF_PREFIX}_exclusive.fasta"

# 6. make probes
"${TAPPR_DIR}/make_probes" \
  -m "${REF_PREFIX}.${KMER_LENGTH}mers.bed" \
  -r "$REF" \
  -k "$KMER_LENGTH" \

# 7. make primers
"${TAPPR_DIR}/make_primers" \
  -m "${REF_PREFIX}.${KMER_LENGTH}mers.probes.bed" \
  -i "${REF_PREFIX}.${KMER_LENGTH}mers.sam" \
  -r "$REF" \
  --ampmax "$AMPMAX" \
  -k "$KMER_LENGTH"

# 8. build inclusivity BLAST DB
makeblastdb \
  -in "${SEQS_DIR}"*fasta \
  -dbtype nucl \
  -parse_seqids \
  -out "$TARGET_NAME"

# 9. simulate PCR
python "${TAPPR_DIR}/simulate_PCR.py" \
  -p "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta" \
  --db "$TARGET_NAME" \
  --evalue 10000 \
  --max_target_seqs 10000

# 10. probe alignment
python "${TAPPR_DIR}/probe_alignment.py" \
  -a "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.${TARGET_NAME}.amplicons" \
  -p "${REF_PREFIX}.${KMER_LENGTH}mers.probes.fasta"

# 10. make data labels
python "${TAPPR_DIR}/make_label_table.py" \
  --groups "${SEQS_DIR}"*fasta "${EXC_DIR}"*fasta \
  --labels "$TARGET_NAME" outgroup \

# 11. inclusivity analysis
python "${TAPPR_DIR}/qpcr_inclusive.py" \
  -a "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.${TARGET_NAME}.probe_aligned.amplicons" \
  -f "${SEQS_DIR}"*fasta \
  --probe \
  -t metadata.tsv \
  --by_primer_pair

# 12. build exclusivity BLAST DB
makeblastdb \
  -in "${EXC_DIR}"*fasta \
  -dbtype nucl \
  -parse_seqids \
  -out "outgroup"

# 12. simulate PCR on exclusive data
python "${TAPPR_DIR}/simulate_PCR.py" \
  -p "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta" \
  --db outgroup \
  --evalue 10000 \
  --max_target_seqs 10000

# 13. probe alignment for outgroup
python "${TAPPR_DIR}/probe_alignment.py" \
  -a "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.outgroup.amplicons" \
  -p "${REF_PREFIX}.${KMER_LENGTH}mers.probes.fasta"

# 14. compile results
python "${TAPPR_DIR}/compile_results.py" \
  -i "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.${TARGET_NAME}.probe_aligned.amplicons" \
  -e "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.outgroup.probe_aligned.amplicons" \

# 15. exclusivity analysis
python "${TAPPR_DIR}/qpcr_exclusive.py" \
  -a "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.outgroup.probe_aligned.compiled.amplicons" \
  -t metadata.tsv \
  --probe \
  --orgs "$TARGET_NAME"

# 16. results summary
python "${TAPPR_DIR}/primer_summary_table.py" \
  -a "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.${TARGET_NAME}.probe_aligned.amplicons" \
  -i "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.${TARGET_NAME}.probe_aligned.amplicons_panel_coverage.tsv" \
  -e "${REF_PREFIX}.${KMER_LENGTH}mers.sam.primers.fasta.pair.outgroup.probe_aligned.compiled.amplicons_exclusive_primers.json" \
  -r "$REF_ID"
