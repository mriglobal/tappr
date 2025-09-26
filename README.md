# TAPPR: Targeted Amplicon Primer/Probe Retrieval

A pipeline for designing targeted PCR primers and probes.

## Installation

### Conda Environment Setup

```
conda env create -f tappr.yml
```

## Workflow Overview

TAPPR provides several paths for primer design depending on the input data and goals:

### Inclusive Intervals
1.  **Conserved Kmer Search (recommended to start):** Start by searching for conserved kmers within your target sequence group (`kmercountinner.py`).
2.  **Clustering (Optional):** If conserved kmers are insufficient, cluster sequences into subgroups (`clustering.py`) and repeat the kmer search on clusters.
3.  **Degenerate Primer Approach (Alternative):** If direct conservation fails, use a degenerate kmer approach (`kmercountouter_sparse.py`, `hammingclusters_fast.py`, `motifmapper.py`, `collapseoverlapingbedentries.py`).
### Exclusive Intervals
1. **Kmer Outer Set (recommended):** Find the corpus of all kmers that occur in the outer group (`kmercountouter_set.py`), to then set difference with inner kmers (`kmer_set_difference.py`) to find exclusive kmers.
2.  **Alternative Methods:**
    * **Large Genomes (CDS-based):** Find unique amino acid markers (`amino_marker_outer.py`, `amino_marker_inner.py`, `kmer_set_difference.py`, `markercleanup.py`).
    * **Small Genomes (Internal Probes):** Map kmers (`tappr_mappr.py`) and find complement regions (`complement_markers.py`).
### Designing Primers & Probes
1.  **Probe Design:** Generate probes using (`make_probes`).
2.  **Primer Design:** Generate primers using (`make_primers`).
### Validating Primer & Probe Designs
1.  **Validation:** Perform in silico PCR (`simulate_PCR.py`).
2.  **Probe Alignment (if applicable):** Align internal probes if designed (`probe_alignment.py`).
3.  **Evaluation:** Assess inclusivity and exclusivity depending on PCR type (`qpcr_inclusive.py`, `qpcr_exclusive.py`, `amplicon_exclusive.py`).

## Detailed Workflow Steps & Usage

[📄 View the full flowchart figure (PDF)](img/Fig_1_TAPPR_flow.pdf)

## Example workflow can be found in tappr_pipeline_example.sh

### 1. Search for Conserved Kmers in Target Group (Inner-join kmer sets)

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


### 2. Cluster Target Group into Subgroups (Subcluster)

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

**Next Steps:** Feed the output clusters back into `kmercountinner.py` (Step 1). Repeat clustering/kmer finding as needed to find a set of groups that can be used to in aggregate cover the desired in group. 

### 3. Vorpal Degenerate Kmer Clustering (Degenerate bases)

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

### 4. Marker Identification

Prepare marker regions for primer design based on genome size and assay type.

**4.1. Find outer kmer set of exclusive records (if exclusivity is desired)**
* **Find Outer Set of Outer Group Kmers (Exclusivity set):**
   ```
   python kmercountouter_set.py -k 18 --seqs outer_group_sequences.zip
   ```
   **Arguments:**
   * `--seqs`: (Required) Folder containing FASTA/FASTQ sequences (or NCBI datasets zip).  
   * `-k`: K-mer size (default: `18`).  
   * `-l`: Minimum length filter for sequences to be included (default: `0`).  
   * `-c`: Number of chunks to split input into for memory management.  
   * `-t`: Number of threads (default: all available).  
   * `-o`: Output prefix for pickled k-mer set file (default: `None`).  
   * `-i`: Pass in FASTA output from `kmercountinner.py` to do set difference (creates additional output).  
   * `--fastq`: Flag if input is in FASTQ format instead of FASTA.
     
* **Find Set Difference of Outer Group Kmer Set and Inner Kmers:**
   ```
   kmer_set_difference.py -A inner_kmers -B outer_kmers
   ```
   **Arguments:**
   * `-k`: K-mer size.  
   * `-A`: K-mer set A (FASTA by default).  
   * `-B`: K-mer set B (pickle file).  
   * `--A_set`: Flag if `-A` is also a pickled k-mer set (instead of FASTA).  
   * `-o`: Output prefix (default: `<A>_exclusive`).
   
**4.2. (Alternative method) Use Coding Domain Sequences Data to find Amino Markers:**

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
    * The output `markers.bed` file is used as input `-m` for `make_primers` (Step 5.1).

**4.2. (Alternative method) Use the complement of inner kmers as marker regions:**
This method is useful if exclusivity isn't a requirement, or for small genomes.
   ```
   python complement_markers.py -i inner_kmers.bed -f reference.fasta -o complement_markers.bed
   ```
   **Arguments:**
   * `-i`: (Required) Primer regions bed or sam file.
   * `-f`: (Required) Fasta file for reference sequence.
   * `-o`: utput file prefix. Default: Inherit from bed/sam file.

### 5. Probe and Primer Design

Generate candidate primers based on identified kmer and/or marker regions.

### A note on marker regions and primer regions
If amplification is the primary reporter, ensure that the input primer regions are exclusive through one of the above methods.
If a probe is used as the reporter, ensure that the input marker regions are exclusive through one of the above methods.
If the input markers file to the make_primers step is ok to draw primers from, use the --marker-primers flag. Note that this flag should not be used if the complement_markers method was used to generate markers.

* **Generate Probes from Cleaned Markers (make_probes):**
   Use the cleaned marker regions to design potential internal probes (oligos).
   ```
   ./make_probes -m cleaned_markers_dir/markers.bed -r reference.fasta --output_name potential_probes -t 4
   ```
   **Arguments:**
   
   * `-m M`: (Required) Marker file
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
 
**`5.1. Using Kmers and Amino Markers (make_primers):`**

```
./make_primers -m cleaned_markers_dir/markers.bed -i mapped_kmers.sam -r reference.fasta
```

**Arguments:**

* `-m`: (Required) Marker file
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

### 9. Evaluate Primer Inclusivity and Exclusivity

```
python dataset_operators.py {subtract,add} 
```
If add:
**Arguments:**
* `-z`: List of zip archives to aggregate.
* `-o`: Output archive name. Default is to inherit name from first archive in list.
If subtract:
**Arguments:**
* `-a`: Zip archive A to be filtered.
* `-b`: Zip archive B whose members will be filtered from A.
* `-o`: Output archive name. Default is to inherit names from A and B archive file names.

```
python make_blast_db.py
```
**Arguments:**
* `-t`: taxid to construct db for
* `-s`: file containing sequences to build a database for
* `-o`: name for output blast database. default: tappr_db/tappr_db

```
python make_label_table.py
```
**Arguments:**
* `--groups`: any number of multifasta files. must be the length of labels args (space separated strings)
* `--labels`: any number of labels for the records in the mutlifasta files. must be the length of groups args (space separated strings)
* `--o`: name of ouptut file. default is metadata.tsv

```
python primer_summary_table.py
```
**Arguments:**
* `-o`: Optionally specify output prefix
* `-a`: Probe Aligned Amplicons File
* `-i`: optional inclusivity .tsv coverage table
* `-e`: optional exclusivity json output to filter by
* `-r`: optional reference genome accession identifier to report positional information


