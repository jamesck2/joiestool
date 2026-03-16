# Joie's Tool

**Gene-level comparison of highly similar bacterial genomes.**
 
Joie's Tool compares complete or near-complete bacterial genomes within a species,
reporting nucleotide-level variants (SNPs, indels), gene presence/absence, and
structural variants (inversions, duplications, gene-order rearrangements) for each
genome relative to a chosen reference.
 
Two comparison modes are supported:
 
- **Single-reference** (default): all genomes compared against one reference.
- **All-vs-all** (`--all-vs-all`): every unique genome pair compared once.
  Outputs a pan-genome presence/absence matrix and pairwise distance summary.


## Installation

### 1. Clone this repository

```bash
git clone https://github.com/jamesck2/joiestool.git
```

### 2. Set up a conda environment to contain the installed packages:

```bash
conda create -n joiestool 'python>=3.8' pip
```

### 3. Activate the conda environment:

```bash
conda activate joiestool
```

### 4. Install the tool and its dependencies:

```bash
pip install ./joiestool
```


## Quick Start

```bash
# Compare all genomes against an auto-selected reference
joiestool genomes/ -o results/ -t 4
 
# Specify a reference genome by name (file stem)
joiestool genomes/ -r reference_strain -o results/
 
# All-vs-all: every unique pair compared, pan-genome matrix produced
joiestool genomes/ --all-vs-all -o results/ -t 4
 
# Use GenBank files with existing annotation
joiestool genbank_genomes/ -f genbank -o results/
 
# Use pre-predicted genes in Prodigal FASTA format
joiestool gene_fastas/ -f genes -o results/
```

> **Note:** If you are providing raw genome sequences (`.fa`/`.fna`/`.fasta` files) rather than pre-annotated GenBank files, Pyrodigal will predict genes automatically.


## Options

```
positional arguments:
  INPUT_DIR             Directory containing one genome file per strain.
 
options:
  -o, --output-dir DIR  Output directory (default: joiestool_out/)
  -r, --reference NAME  Reference genome name (file stem, no extension).
                        Default: genome with the most predicted genes.
                        Ignored when --all-vs-all is set.
  --all-vs-all          Compare every unique genome pair. The more gene-rich
                        genome is the reference for each pair, the most-complete
                        genome anchors the pan-genome matrix.
  -f, --format {auto,fasta,genbank,genes}
                        Input format. Default: auto-detect per file.
  -t, --threads N       Worker threads (default: 1)
  --min-identity FLOAT  Minimum nucleotide identity for ortholog assignment
                        (default: 0.70).
  --kmer-size K         K-mer size for fast candidate screening (default: 11).
  --meta                Use metagenomic Prodigal models for gene prediction.
                        Auto-applied for genomes < 100 kbp.
  -v, --verbose         Debug-level logging.
  --version             Show version and exit.
```


## Input Formats

Place one file per genome (or strain) in a single directory. Joie's Tool will detect the format of each file automatically unless you specify it with `--format`.

| Extension | What it contains | How genes are obtained |
|---|---|---|
| `.fa` `.fna` `.fasta` | Full nucleotide genome sequence | Predicted automatically using Prodigal (via Pyrodigal) |
| `.gbk` `.gb` `.genbank` `.gbff` | GenBank-format annotated genome | Loaded from existing CDS (coding sequence) features |
| `.ffn` | Gene nucleotide sequences in Prodigal header format | Loaded directly (no prediction step) |

You can mix formats within a single run by leaving `--format` set to `auto` (the default).


## Notes
Joie's Tool can handle fragmented assemblies. Genes predicted on short or incomplete contigs may be flagged as `is_partial=yes` in the diversity table. Structural variant detection (inversions, large deletions) works per-contig, so it is less sensitive for highly fragmented assemblies.

Joie's Tool is optimised for within-species comparisons, where most genes are expected to be shared and highly similar. For more divergent comparisons (e.g., genus-level), many genes will fall below the identity threshold and be reported as novel or absent. Lowering `--min-identity` can help, but results should be interpreted with care.


## How It Works

### Overview

1. **Gene prediction** (if needed): If you supply raw genome FASTA files, genes are predicted using [Prodigal](https://github.com/hyattpd/Prodigal).

2. **Reference selection**: User-specified or auto-selected (most genes). In all-vs-all mode the more gene-rich genome is reference for each pair.

3. **K-mer indexing**: *k*-mers (default k=11) of all reference genes indexed for O(1) candidate lookup.  In all-vs-all mode each unique reference genome is indexed once and reused across all pairs it anchors.

4. **Ortholog detection**: Each query gene is matched to its best reference gene by k-mer Jaccard similarity, then verified by global pairwise alignment (BioPython PairwiseAligner).

5. **Variant calling**: For each matched gene pair (query vs. reference ortholog), the alignment is examined base by base to identify point mutations (SNPs), insertions, and deletions. The position and nature of every change is recorded.

6. **Structural variants**: Gene-order profiles compared to detect inversions, duplications, and large-scale deletions.

7. **Non-coding regions**: Intergenic sequences aligned when full genome sequences are available (genome FASTA or GenBank input).

### All-vs-all pairing rules
 
- N genomes produce N*(N-1)/2 unique comparisons.
- For each pair the genome with more genes is the reference.
- Ties in gene count are broken by lexicographic name order.
- The pan-genome anchor (most genes overall) is guaranteed to be the reference in all pairs it participates in, so a single index build covers those pairs.

### Single-reference mode (`-o OUTPUT_DIR/`)
 
#### `gene_table.tsv`: Wide-format gene presence/variant table
 
Rows = genes (reference-ordered, novel query genes appended at end).
Two columns per genome: `{genome}_status` and `{genome}_variants`.
 
| gene_id | ref_contig | ref_start | ref_end | ref_strand | strain_B_status | strain_B_variants | ... |
|---|---|---|---|---|---|---|---|
| locus_0001 | chr1 | 100 | 1000 | + | identical | | ... |
| locus_0002 | chr1 | 1050 | 2100 | + | variant | SNP:A42T; DEL:200-3bp | ... |
| locus_0003 | chr1 | 2200 | 3100 | - | absent | | ... |
| novel:strain_B_99 | chr1 | 5000 | 5900 | + | novel | | ... |
 
#### `gene_diversity.tsv`: Long-format per-gene statistics
 
One row per genome x gene combination.
 
Columns: `genome`, `gene_id`, `is_present`, `is_novel`, `is_partial`,
`is_truncated`, `identity`, `num_snps`, `num_indels`, `num_variants`,
`ref_gene_length`, `query_gene_length`, `notes`
 
#### `variant_summary.tsv`: Variant event summary
 
One row per variant event (nucleotide or structural).
 
| ref_genome | query_genome | variant_type | gene_or_region | ref_start | ref_end | details |
|---|---|---|---|---|---|---|
| ref | strain_B | SNP | locus_0002 | 142 | 142 | SNP:A42T |
| ref | strain_B | DEL | locus_0002 | 200 | 202 | DEL:200-3bp |
| ref | strain_B | GENE_LOSS | locus_0003 | | | Reference gene absent from query genome |
| ref | strain_B | INVERSION | locus_0010,locus_0011 | | | Block of 2 genes appears inverted |
 
Non-coding variants (when full genome FASTA is provided) are prefixed `NONCODING_`.

 
### All-vs-all mode (`--all-vs-all`)
 
#### `gene_presence_absence.tsv`: Pan-genome presence/absence matrix
 
Anchored to the most-complete genome (most genes), its column appears first.
One status/variants column pair per genome.
 
| gene_id | ref_contig | ref_start | ref_end | ref_strand | anchor_status | anchor_variants | strain_B_status | ... |
|---|---|---|---|---|---|---|---|---|
| locus_0001 | chr1 | 100 | 1000 | + | identical | | identical | ... |
| locus_0002 | chr1 | 1050 | 2100 | + | identical | | variant | SNP:A42T |
| locus_0003 | chr1 | 2200 | 3100 | - | identical | | absent | | |
| novel:strain_B_99 | chr1 | 5000 | 5900 | + | absent | | novel | | |
 
#### `pairwise_summary.tsv`: Per-pair aggregate statistics
 
One row per unique genome pair.  All counts are relative to the reference
genome used for that pair (the more gene-rich genome).
 
| ref_genome | query_genome | shared_genes | identical_genes | variant_genes | truncated_genes | novel_query_genes | deleted_ref_genes | total_snps | total_indels | total_variants | structural_variants | mean_identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| strain_A | strain_B | 3800 | 3640 | 150 | 10 | 5 | 8 | 320 | 45 | 365 | 2 | 0.9993 |
 
#### `variant_summary.tsv`: Variant events across all pairs
 
Same format as single-reference mode. `ref_genome` and `query_genome` columns identify which pair each row belongs to.


## References

- **Pyrodigal**: Larralde M. 2022. *JOSS* doi:10.21105/joss.04296
- **Prodigal**: Hyatt D et al. 2010. *BMC Bioinformatics* doi:10.1186/1471-2105-11-119
- **Biopython**: Cock PJA et al. 2009. *Bioinformatics* doi:10.1093/bioinformatics/btp163