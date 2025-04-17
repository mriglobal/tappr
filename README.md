# TAPPR: Targeted Amplicon Primer/Probe Retrieval

A pipeline for designing targeted PCR primers and probes.

## Installation

### Conda Environment Setup

```
conda env create -f tappr.yml

conda activate tappr

pip install MRITaxonomy
```

## Workflow Overview

TAPPR provides several paths for primer design depending on the input data and goals:

1.  **Conserved Kmer Search:** Start by searching for conserved kmers within your target sequence group (`kmercountinner.py`).
2.  **Clustering (Optional):** If conserved kmers are insufficient, cluster sequences into subgroups (`clustering.py`) and repeat the kmer search on clusters.
3.  **Degenerate Primer Approach (Alternative):** If direct conservation fails, use a degenerate kmer approach (`kmercountouter_sparse.py`, `hammingclusters_fast.py`, `motifmapper.py`, `collapseoverlapingbedentries.py`).
4.  **Marker Identification:**
    * **Large Genomes (CDS-based):** Find unique amino acid markers (`amino_marker_outer.py`, `amino_marker_inner.py`, `kmer_set_difference.py`, `markercleanup.py`).
    * **Small Genomes/Genes (Sybr Green):** Map kmers to a reference (`bowtie2`) for `make_primers_gene.py`.
    * **Small Genomes (Internal Probes):** Map kmers (`bowtie2`) and find complement regions (`complement_markers.py`).
5.  **Primer Design:** Generate primers using appropriate scripts (`make_primers.py` or `make_primers_gene.py`).
6.  **Validation:** Perform in silico PCR (`simulate_PCR.py`).
7.  **Probe Alignment (Optional):** Align internal probes if designed (`probe_alignment.py`).
8.  **Evaluation:** Assess inclusivity and exclusivity (`qpcr_inclusive.py`, `qpcr_exclusive.py`, `amplicon_exclusive.py`).

## Detailed Workflow Steps & Usage

### 1. Search for Conserved Kmers in Target Group

Use `kmercountinner.py` to find kmers conserved across a set of sequences.

```
python kmercountinner.py --seqs /path/to/sequence/directory -k 18
```

**Arguments:**

* `--seqs`: (Required) Directory containing sequence files.
* `-k`: (Required) Kmer size to search for conservation.
* `--directory`: Output directory (optional).
* `-r`: Reference genome FASTA for sequence binning (optional).
* `-o`: Output filename prefix (optional).
* `-p`: Allowed percent variance in reference length for replicon binning (default: 0.1, useful for large genomes) (optional).
* `--assembly_level`: Flag to count kmers by assembly rather than per record (optional).

**Next Steps:**

* **If successful:** Proceed to Marker Identification (Step 4) or Primer Design (Step 5).
* **If fails:** Proceed to Clustering (Step 2) or Degenerate Primer Approach (Step 3).

### 2. Cluster Target Group into Subgroups (Optional)

If `kmercountinner.py` doesn't find enough conserved kmers, cluster the sequences.

```
python clustering.py --seqs clustering_records.fasta -k 11 --criterion maxclust -c 3 -o outprefix
```

**Arguments:**

* `--seqs`: (Required) Single multifasta file or directory of FASTA files to cluster.
* `-k`: Kmer length for Sourmash sketching (default: 11).
* `--criterion`: Clustering criterion: `maxclust` (target number of clusters) or `maxdist` (maximum distance cutoff) (default: `maxclust`).
* `-c`: Clustering criterion value (default: 2 if `maxclust`, 0.2 if `maxdist`).
* `-o`: Output directory prefix (optional).
* `-l`: Label prefix for clusters and TSV output (optional).
* `--abundance`: Track kmer abundance (useful for varying sequence sizes) (optional flag).
* `--reference`: Randomly sample a sequence representative for each cluster (optional flag).
* `--tsv`: Output a TSV file mapping sequence IDs to cluster membership (optional flag).
* `-s`: Scaled argument for Sourmash sketch (1/S kmers sampled) (default: 1).

**Next Steps:** Feed the output clusters back into `kmercountinner.py` (Step 1). Repeat clustering/kmer finding as needed.

### 3. Vorpal Degenerate Kmer Clustering (Alternative)

An alternative approach if direct kmer conservation fails.

**3.1. Generate Sparse Kmer Counts:**

```
python kmercountouter_sparse.py -r ref_genome.fasta --seqs /sequence_folder/ -k 17
```

**Arguments:**

* `-r`: (Required) Reference genome FASTA.
* `--seqs`: (Required) Folder of sequences to count from.
* `-k`: Size of kmers to count (default: 17).
* `-n`: Number of samples for resampling (optional, default: 0).
* `-p`: Allowed percent variance in reference length for replicon binning (default: 0.1) (optional).

**3.2. Cluster Sparse Kmers by Hamming Distance:**

```
python hammingclusters_fast.py -p kmers.pkl -n 1.0 -q 0.2 -c 3
```

**Arguments:**

* `-p`: (Required) `.pkl` output file from `kmercountouter_sparse.py`.
* `-n`: Average allowed degenerate bases (default: 4.0).
* `-q`: Quantile of clustered kmers to use (optional, default: all kmers).
* `-c`: Number of chunks for matrix processing (default: 0) (optional).
* `--temp`: Location of temp directory for memory-mapped matrix (optional).
* `--mem`: Memory allocation for processing distance matrix chunks (optional).

**3.3. Map Degenerate Primers to Sequences:**

```
python motifmapper.py -r reference.fasta --seqs seqs.fasta -k degenerate_primers.fasta -s 800
```

**Arguments:**

* `-r`: (Required) Reference FASTA.
* `--seqs`: (Required) Set of sequences to score degenerate primers against.
* `-k`: (Required) Degenerate primers FASTA (`degenerate_primers.fasta` from previous step).
* `-s`: Motif conservation cutoff (0-1000, corresponds to 0-100%) (default: 800) (optional).
* `-t`: Number of threads to use (optional).

**3.4. Collapse Overlapping BED Entries:**

```
python collapseoverlapingbedentries.py primers.bed
```

**Arguments:**

* Positional: Input `.bed` file (`primers.bed` from previous step).

**Next Steps:** The output `primers.bed` can potentially be used as input for primer design scripts (Step 5), though the primary path involves marker identification.

### 4. Marker Identification

Prepare marker regions for primer design based on genome size and assay type.

**4.1. For Small Genomes / Gene Targets (Sybr Green or Internal Probes):**

* **Map Kmers to Reference (Example using Bowtie2):**
    ```
    # Build index
    bowtie2-build reference.fasta ref_prefix

    # Align kmers (e.g., from kmercountinner.py)
    # Ensure your kmer file is in FASTA format
    bowtie2 --end-to-end -x ref_prefix -f -U conserved_kmers.fasta -S mapped_kmers.sam
    ```
    * The output `mapped_kmers.sam` (or converted `.bed`) is used as input `-i` for `make_primers_gene.py` (Step 5.2).

* **Find Complement Regions for Internal Probes (Optional):**
    ```
    python complement_markers.py -i mapped_kmers.sam -f reference.fasta -o complement_out
    ```
    **Arguments:**
    * `-i`: (Required) Input file (`.sam` or `.bed`) of kmers mapped to reference.
    * `-f`: (Required) Reference FASTA used for mapping.
    * `-o`: Output prefix (optional).
    * The output `.bed` file can potentially inform internal oligo design.

**4.2. For Larger Genomes (Using Amino Acid Markers):**

* **Find Kmers Unique to Outer Group (Exclusivity):**
    ```
    python amino_marker_outer.py -k 8 --seqs outer_group_sequences.zip -t 4 -o outer_group
    ```
    **Arguments:**
    * `-k`: (Required) Kmer size to count.
    * `--seqs`: (Required) Outer/exclusive group CDS (multifasta or NCBI datasets zip).
    * `-t`: Number of threads (optional).
    * `-c`: Use memory chunks (optional flag).
    * `-o`: Output prefix for `outer_kmers.pickle` (optional).

* **Find Kmers Unique to Inner Group (Inclusivity):**
    ```
    python amino_marker_inner.py -k 8 --seqs inner_group_sequences.zip -t 4 -o inner_group
    ```
    **Arguments:** (Identical to `amino_marker_outer.py` but for the inner group)
    * `-k`: (Required) Kmer size.
    * `--seqs`: (Required) Inner/inclusive group CDS.
    * `-t`: Threads (optional).
    * `-c`: Memory chunks (optional flag).
    * `-o`: Output prefix for `inner_kmers.pickle` (optional).

* **Consolidate Inner/Outer Kmers:**
    ```
    python kmer_set_difference.py -k 8 -A inner_group_inner_kmers.pickle -B outer_group_outer_kmers.pickle --A_set -o unique_inner_markers
    ```
    **Arguments:**
    * `-k`: (Required) Kmer size used in previous steps.
    * `-A`: (Required) Output from `amino_marker_inner.py` (`.pickle` or `.fasta`).
    * `-B`: (Required) Output from `amino_marker_outer.py` (`.pickle` or `.fasta`).
    * `--A_set`: Pass this flag if inputs are kmer set `.pickle` files.
    * `-o`: Output prefix (optional). Output is `protein.fasta`.

* **Clean Up Markers Against Reference Genome:**
    ```
    python markercleanup.py -c unique_inner_markers_protein.fasta -r ref_genome.fasta -t 4 -o cleaned_markers_dir
    ```
    **Arguments:**
    * `-c`: (Required) Protein marker file (`protein.fasta` from previous step).
    * `-r`: (Required) Reference genome FASTA to map markers against.
    * `-o`: (Required) Output directory.
    * `-e`: E-value cutoff for BLAST (optional).
    * `-t`: Threads for BLAST (optional).
    * The output `markers.bed` file is used as input `-m` for `make_primers.py` (Step 5.1).
 
* **Generate Probes from Cleaned Markers (make_probes):**
   Use the cleaned marker regions to design potential internal probes (oligos).
   ```
   ./make_probes -m cleaned_markers_dir/markers.bed -r reference.fasta --output_name potential_probes -t 4
   ```
   **Arguments:**
   
   * `-m M`: (Required) Marker file post-marker cleanup in BED format (markers.bed from previous step).
   
   * `-r R`: (Required) Reference genome of interest (FASTA).
   
   * `--iomin IOMIN`: Interior oligo minimum size (default: 18).
   
   * `--iomax IOMAX`: Interior oligo maximum size (default: 28).
   
   * `--io_opt_tm IO_OPT_TM`: Interior oligo optimal melting temp (default: 65 C).
   
   * `--io_5prime IO_5PRIME`: 5' sequence requirement for Interior Oligo (e.g., HNNNN) (default: None).
   
   * `-g G`: GFF annotation file for reference (optional, default: None).
   
   * `--cds_limit`: Flag to limit candidate amplicons to regions entirely within a CDS (optional flag).
   
   * `-k K`: Kmer size (default: 18) (Note: Ensure consistency if using kmer info, though primary input is markers).
   
   * `--output_name OUTPUT_NAME`: Output name prefix for probe files (optional).
   
   * `-t THREADS`, --threads THREADS: Number of threads for multiprocessing (optional).

This step generates candidate probe sequences based on the marker locations and specified parameters.

### 5. Make Primers

Generate candidate primers based on identified kmer or marker regions.

**`5.1. Using Kmers and Amino Markers (make_primers):`**

```
./make_primers -m cleaned_markers_dir/markers.bed -i mapped_kmers.sam -r reference.fasta
```

**Arguments:**

* `-m`: (Required) Marker file (`.bed` from `markercleanup.py`).
* `-i`: (Required) Input kmer alignment file (`.bed` or `.sam`).
* `-r`: (Required) Reference genome FASTA.
* `--ampmax`: Maximum amplicon length (default: 250).
* `--ampmin`: Minimum amplicon length (default: 150).
* `--io`: Flag to pick an interior oligo from marker regions.
* `--iomin`: Interior oligo minimum size (default: 18).
* `--io_optimum`: Interior oligo optimum size (default: 20).
* `--iomin_temp`: Interior oligo min melting temp (default: 57 C).
* `--iomax_temp`: Interior oligo max melting temp (default: 63 C).
* `--io_5prime`: 5' sequence requirement for Interior Oligo (e.g., `HNNNN`) (default: None).
* `-k`: Kmer size used for `-i` file (default: 18).
* `--directory`: Output directory (optional).
* `-g`: GFF annotation file for reference (optional).
* `--primermin`: Minimum primer length (default: 15).
* `--primermax`: Maximum primer length (default: 25).
* `--primermin_temp`: Minimum primer melting temp (default: 56 C).
* `--primermax_temp`: Maximum primer melting temp (default: 64 C).
* `--primer_optimum`: Optimum primer length (default: 18).
* `--cds_limit`: Flag to limit candidate amplicons to regions entirely within a CDS.
* 

**Next Steps:** If candidate primers are produced (`primers.fasta`), proceed to In Silico PCR (Step 6).

### 6. In Silico PCR

Simulate PCR using the designed primers against a BLAST database.

```
python simulate_PCR.py -p primers.fasta --db /path/to/blastdb --extract_amp
```

**Arguments:**

* `-p`: (Required) Primer FASTA file (output from Step 5).
* `--db`: (Required) Path to BLAST database.
* `--extract_amp`: Flag to extract amplicon sequences (outputs `.amplicons` file).
* `-r`: Pre-computed BLAST output table (optional).
* `--mismatches`: Allowed mismatches (integer count or float fraction < 1) (default: 3).
* `--min`: Minimum product length (default: 40).
* `--max`: Maximum product length (default: 1000).
* `--mux`: Multiplex flag (default: False).
* `-t`: Number of threads for BLAST (default: all).
* `--max_target_seqs`: Max target sequences from BLAST (default: 1000).
* `--evalue`: E-value cutoff for BLAST (optional).
* `--word_size`: Word size for BLAST (default: 4).

**Next Steps:** Analyze the `.amplicons` file. If internal probes were designed, proceed to Probe Alignment (Step 7). Then evaluate results (Step 8).

### 7. Probe Alignment (Optional)

Align internal probes to the simulated amplicons.

```
python probe_alignment.py -a simulated_pcr.amplicons -p primers_and_probes.fasta --mux -o probe_aln_out -t 4
```

**Arguments:**

* `-a`: (Required) Amplicons file from `simulate_PCR.py` (with `--extract_amp`).
* `-p`: (Required) Oligos file (FASTA) containing primers and probes (probes tagged with `|IO`).
* `--mux`: (Required) Flag for Multiplex Analysis (considers all probe/amplicon combinations).
* `-o`: (Required) Prefix for output files.
* `-t`: (Required) Number of threads for multiprocessing.
* `-m`: Maximum number of mismatches allowed (default: 3).

### 8. Evaluate Primer Inclusivity and Exclusivity

Assess primer performance against target (inclusive) and non-target (exclusive) sequences.

**8.1. Inclusivity Check:**

```
python qpcr_inclusive.py -a simulated_pcr.amplicons -f inclusivity_sequences.fasta
```

**Arguments:**

* `-a`: (Required) Amplicons file from `simulate_PCR.py`.
* `-f`: (Required) FASTA file of sequences for inclusivity evaluation.
* `--subspecies`: Flag for subspecies level analysis.
* `--probe`: Flag if the assay uses an interior oligo.
* `-m`: Minimum match fraction (percentage of oligo length) for alignment cutoff (optional).
* `--name`: Return summary by organism name instead of taxid (default: False).
* `-t`: Tab-delimited metadata table for custom labeling (optional).
* `--by_primer_pair`: Split results by primer pairs (optional flag).
* `--blacklist`: List of accession numbers to exclude (optional).
* `--min`: Minimum amplicon length to filter (optional).
* `--split_col`: Column index (0-based) to split HitName for accession (default: 0 for GenBank).
* `-o`: Output prefix (optional, defaults to input prefix).

**8.2. Exclusivity Check (Taxonomy-based):**

```
python qpcr_exclusive.py -a simulated_pcr.amplicons --orgs 9606 10090 # Example taxids
```

**Arguments:**

* `-a`: (Required) Amplicons file from `simulate_PCR.py`.
* `--orgs`: (Required) List of target taxids (space-separated). Exclusivity is checked against hits *not* matching these taxids.
* `--subspecies`: Flag for subspecies level analysis.
* `-t`: TSV file for custom labels (optional, overrides `--subspecies`).
* `--probe`: Flag if the assay uses an interior oligo.
* `--blacklist`: List of accession numbers to exclude (optional).
* `--min`: Minimum amplicon length to filter (optional).
* `--max`: Maximum amplicon length to filter (optional).
* `--split_col`: Column index (0-based) to split HitName for accession (default: 1).
* `-o`: Output prefix (optional, defaults to input prefix).

**8.3. Exclusivity Check (Amplicon Sequence Level):**

```
python amplicon_exclusive.py -a simulated_pcr.amplicons --targets 9606 10090 # Example taxids
```

**Arguments:**

* `-a`: (Required) Amplicons file from `simulate_PCR.py`.
* `--targets`: (Required) List of target taxids (space-separated).
* `--read_length`: Expected read length for sequencing (0 for variable length) (optional).
* `-t`: TSV file for custom labels matching `--targets` (optional).
* `--subspecies`: Flag for subspecies level analysis.
* `--spe`: Flag for single primer extension amplicon analysis.
* `--blacklist`: List of accession numbers to exclude (optional).
* `--min`: Minimum amplicon length to filter (optional).
* `--max`: Maximum amplicon length to filter (optional).
* `--split_col`: Column index (0-based) to split HitName for accession (default: 0 for GenBank).
