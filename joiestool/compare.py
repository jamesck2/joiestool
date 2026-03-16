"""
Core genome comparison logic.

Pipeline per query genome (relative to a chosen reference):

  1. KmerIndex: hash all reference gene sequences for fast candidate lookup.
  2. find_orthologs(): map each query gene to its best-matching reference gene
     using k-mer Jaccard similarity, verified with pairwise alignment.
  3. call_gene_variants(): align ortholog pairs and call SNPs + indels.
  4. detect_structural_variants(): compare gene-order profiles to find
     inversions, duplications, large-scale deletions.
  5. compare_noncoding(): if full sequences are present, align corresponding
     intergenic regions.

All-vs-all mode: compare_all_vs_all() runs every unique genome pair once,
using the gene-richer genome as reference for each pair.
"""
from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from joiestool.genome import Gene, Genome, NoncodingRegion

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Variant data classes
# ---------------------------------------------------------------------------


@dataclass
class NucleotideVariant:
    kind: str           # 'SNP', 'INS', 'DEL', 'MNP'
    ref_pos: int        # 0-based position in the reference gene/region
    ref_allele: str     # reference allele ('' for pure insertion)
    alt_allele: str     # alternate allele   ('' for pure deletion)

    def __str__(self) -> str:
        if self.kind == "SNP":
            return f"SNP:{self.ref_allele}{self.ref_pos + 1}{self.alt_allele}"
        if self.kind == "INS":
            return f"INS:{self.ref_pos + 1}+{self.alt_allele}"
        if self.kind == "DEL":
            return f"DEL:{self.ref_pos + 1}-{len(self.ref_allele)}bp"
        return f"{self.kind}:{self.ref_pos + 1}"


@dataclass
class GeneComparison:
    """Result of comparing one query gene to its reference ortholog."""
    query_gene: Gene
    ref_gene: Optional[Gene]          # None -> novel gene (no ref ortholog)
    jaccard: float = 0.0              # k-mer Jaccard before alignment
    identity: float = 0.0            # nucleotide identity from alignment
    alignment_length: int = 0
    variants: List[NucleotideVariant] = field(default_factory=list)
    is_novel: bool = False            # in query but not in reference
    is_truncated: bool = False        # large length difference

    @property
    def num_snps(self) -> int:
        return sum(1 for v in self.variants if v.kind == "SNP")

    @property
    def num_indels(self) -> int:
        return sum(1 for v in self.variants if v.kind in ("INS", "DEL"))

    @property
    def variant_summary(self) -> str:
        if not self.variants:
            return ""
        return "; ".join(str(v) for v in self.variants[:20]) + (
            f" (+{len(self.variants) - 20} more)" if len(self.variants) > 20 else ""
        )


@dataclass
class StructuralVariant:
    kind: str                   # 'INVERSION', 'DUPLICATION', 'DELETION', 'INSERTION', 'REARRANGEMENT'
    genes_ref: List[str]        # reference gene names involved
    genes_query: List[str]      # query gene names involved
    description: str = ""

    def __str__(self) -> str:
        genes = ",".join(self.genes_ref) if self.genes_ref else "(none)"
        return f"{self.kind}[{genes}]: {self.description}"


@dataclass
class NoncodingComparison:
    ref_region: NoncodingRegion
    query_region: Optional[NoncodingRegion]
    identity: float = 0.0
    alignment_length: int = 0
    variants: List[NucleotideVariant] = field(default_factory=list)

    @property
    def num_snps(self) -> int:
        return sum(1 for v in self.variants if v.kind == "SNP")

    @property
    def num_indels(self) -> int:
        return sum(1 for v in self.variants if v.kind in ("INS", "DEL"))


@dataclass
class GenomeComparisonResult:
    """All comparison results for one query genome vs the reference."""
    query_name: str
    ref_name: str
    # Gene-level results keyed by reference gene name (or query gene name for novels)
    gene_comparisons: Dict[str, GeneComparison] = field(default_factory=dict)
    # Reference genes with no orthologs in query (deleted)
    deleted_ref_genes: List[str] = field(default_factory=list)
    # Structural variants detected
    structural_variants: List[StructuralVariant] = field(default_factory=list)
    # Non-coding comparisons keyed by region name
    noncoding_comparisons: Dict[str, NoncodingComparison] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# K-mer index
# ---------------------------------------------------------------------------

class KmerIndex:
    """Inverted k-mer index over a list of genes for fast similarity lookups."""

    def __init__(self, genes: List[Gene], k: int = 11) -> None:
        self.k = k
        self.genes = genes
        # gene_kmers[i] = frozenset of k-mers in genes[i]
        self.gene_kmers: List[frozenset] = []
        # inverted_index[kmer] = list of gene indices
        self.inverted_index: Dict[str, List[int]] = defaultdict(list)
        self._build(genes)

    def _build(self, genes: List[Gene]) -> None:
        for i, gene in enumerate(genes):
            # Use a frozenset so each k-mer is counted once per gene (for Jaccard)
            kmers = frozenset(_extract_kmers(gene.sequence, self.k))
            self.gene_kmers.append(kmers)
            for km in kmers:
                self.inverted_index[km].append(i)

    def best_match(
        self,
        query_seq: str,
        min_jaccard: float = 0.3,
    ) -> Tuple[Optional[int], float]:
        """
        Return the index of the most similar gene in the index, and its
        Jaccard similarity, or (None, 0.0) if none exceeds min_jaccard.
        """
        query_kmers = frozenset(_extract_kmers(query_seq, self.k))
        if not query_kmers:
            return None, 0.0

        # Count shared k-mers per candidate gene
        shared_counts: Dict[int, int] = defaultdict(int)
        for km in query_kmers:
            for idx in self.inverted_index.get(km, []):
                shared_counts[idx] += 1

        best_idx: Optional[int] = None
        best_j = 0.0
        for idx, shared in shared_counts.items():
            ref_kmers = self.gene_kmers[idx]
            union = len(query_kmers) + len(ref_kmers) - shared
            j = shared / union if union > 0 else 0.0
            if j > best_j:
                best_j = j
                best_idx = idx

        if best_j < min_jaccard:
            return None, 0.0
        return best_idx, best_j


def _extract_kmers(seq: str, k: int) -> List[str]:
    seq = seq.upper()
    return [seq[i: i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i: i + k]]


# ---------------------------------------------------------------------------
# Pairwise alignment and variant calling
# ---------------------------------------------------------------------------

def _align_biopython(
    ref_seq: str,
    query_seq: str,
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    Use BioPython PairwiseAligner to globally align two sequences.
    Returns (ref_blocks, query_blocks) each a list of (start, end) pairs
    covering matching runs.
    """
    from Bio.Align import PairwiseAligner
    aln = PairwiseAligner()
    aln.mode = "global"
    aln.match_score = 2
    aln.mismatch_score = -1
    aln.open_gap_score = -3
    aln.extend_gap_score = -0.5
    alignment = aln.align(ref_seq, query_seq)[0]
    return list(alignment.aligned[0]), list(alignment.aligned[1])


def _align_nw(ref_seq: str, query_seq: str) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    Pure-Python affine-gap Needleman-Wunsch global alignment.

    Optimised for the common case of highly similar sequences (same-species
    strains) where most genes are >95% identical.  Uses a simple O(n*m)
    DP with linear space traceback for sequences up to ~5 kbp.  For longer
    sequences the comparison is chunked at conserved anchor points.

    Returns (ref_blocks, query_blocks): lists of (start, end) pairs for
    consecutive equal-sequence blocks in the alignment.
    """
    n, m = len(ref_seq), len(query_seq)

    # For very long sequences use anchor-chunked comparison
    if n > 5000 or m > 5000:
        return _align_nw_chunked(ref_seq, query_seq)

    MATCH = 2
    MISMATCH = -1
    GAP_OPEN = -3
    GAP_EXTEND = -0.5

    # Affine gap: three matrices M (match/mismatch), X (gap in ref), Y (gap in query)
    NEG_INF = float("-inf")
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0.0
    for i in range(1, n + 1):
        X[i][0] = GAP_OPEN + (i - 1) * GAP_EXTEND
    for j in range(1, m + 1):
        Y[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND

    for i in range(1, n + 1):
        ri = ref_seq[i - 1]
        for j in range(1, m + 1):
            score = MATCH if ri == query_seq[j - 1] else MISMATCH
            M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + score
            X[i][j] = max(M[i-1][j] + GAP_OPEN, X[i-1][j] + GAP_EXTEND,
                          Y[i-1][j] + GAP_OPEN)
            Y[i][j] = max(M[i][j-1] + GAP_OPEN, X[i][j-1] + GAP_OPEN,
                          Y[i][j-1] + GAP_EXTEND)

    # Traceback -- determine which matrix to start in
    best = max(M[n][m], X[n][m], Y[n][m])
    if best == M[n][m]:
        state = "M"
    elif best == X[n][m]:
        state = "X"
    else:
        state = "Y"

    cigar: List[str] = []   # 'M'=match/mismatch, 'D'=del in query, 'I'=ins in query
    i, j = n, m
    while i > 0 or j > 0:
        if state == "M":
            cigar.append("M")
            score = MATCH if ref_seq[i - 1] == query_seq[j - 1] else MISMATCH
            prev_m = M[i-1][j-1] + score
            prev_x = X[i-1][j-1] + score
            prev_y = Y[i-1][j-1] + score
            best_prev = max(prev_m, prev_x, prev_y)
            if best_prev == prev_m:
                state = "M"
            elif best_prev == prev_x:
                state = "X"
            else:
                state = "Y"
            i -= 1; j -= 1
        elif state == "X":
            cigar.append("D")
            if X[i-1][j] + GAP_EXTEND == X[i][j]:
                state = "X"
            else:
                state = "M"
            i -= 1
        else:  # Y
            cigar.append("I")
            if Y[i][j-1] + GAP_EXTEND == Y[i][j]:
                state = "Y"
            else:
                state = "M"
            j -= 1

    cigar.reverse()

    # Convert cigar to block lists
    ref_blocks: List[Tuple[int, int]] = []
    qry_blocks: List[Tuple[int, int]] = []
    ri = qi = 0
    k = 0
    while k < len(cigar):
        op = cigar[k]
        run = 1
        while k + run < len(cigar) and cigar[k + run] == op:
            run += 1
        if op == "M":
            ref_blocks.append((ri, ri + run))
            qry_blocks.append((qi, qi + run))
            ri += run; qi += run
        elif op == "D":
            ri += run
        else:  # I
            qi += run
        k += run

    return ref_blocks, qry_blocks


def _align_nw_chunked(
    ref_seq: str,
    query_seq: str,
    chunk: int = 500,
    kmer_anchor: int = 15,
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """
    For long sequences: find shared k-mer anchors, then NW-align between anchors.
    Falls back to treating the whole region as a single (unresolved) block if no
    anchors are found.
    """
    # Find unique anchors
    k = kmer_anchor
    ref_kmers: Dict[str, int] = {}
    for i in range(len(ref_seq) - k + 1):
        km = ref_seq[i: i + k]
        if km not in ref_kmers:
            ref_kmers[km] = i
        else:
            ref_kmers[km] = -1  # mark duplicate

    anchors: List[Tuple[int, int]] = []  # (ref_pos, qry_pos)
    for j in range(len(query_seq) - k + 1):
        km = query_seq[j: j + k]
        ri = ref_kmers.get(km, -1)
        if ri >= 0:
            anchors.append((ri, j))
            ref_kmers[km] = -1   # consume anchor (use each once)

    anchors.sort()

    if not anchors:
        # No shared anchors: return a single block for the shorter sequence
        min_len = min(len(ref_seq), len(query_seq))
        return [(0, min_len)], [(0, min_len)]

    all_ref_blocks: List[Tuple[int, int]] = []
    all_qry_blocks: List[Tuple[int, int]] = []

    prev_r, prev_q = 0, 0
    for (r_anchor, q_anchor) in anchors:
        if r_anchor < prev_r or q_anchor < prev_q:
            continue
        # Align the gap before this anchor
        gap_r = ref_seq[prev_r: r_anchor]
        gap_q = query_seq[prev_q: q_anchor]
        if gap_r or gap_q:
            sub_r, sub_q = _align_nw(gap_r, gap_q)
            all_ref_blocks.extend((b + prev_r, e + prev_r) for b, e in sub_r)
            all_qry_blocks.extend((b + prev_q, e + prev_q) for b, e in sub_q)
        # Add the anchor as a matching block
        all_ref_blocks.append((r_anchor, r_anchor + k))
        all_qry_blocks.append((q_anchor, q_anchor + k))
        prev_r = r_anchor + k
        prev_q = q_anchor + k

    # Align the tail
    tail_r = ref_seq[prev_r:]
    tail_q = query_seq[prev_q:]
    if tail_r or tail_q:
        sub_r, sub_q = _align_nw(tail_r, tail_q)
        all_ref_blocks.extend((b + prev_r, e + prev_r) for b, e in sub_r)
        all_qry_blocks.extend((b + prev_q, e + prev_q) for b, e in sub_q)

    # Merge adjacent/overlapping blocks
    return _merge_blocks(all_ref_blocks), _merge_blocks(all_qry_blocks)


def _merge_blocks(blocks: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not blocks:
        return []
    merged = [blocks[0]]
    for b, e in blocks[1:]:
        if b <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((b, e))
    return merged


def align_and_call_variants(
    ref_seq: str,
    query_seq: str,
    ref_offset: int = 0,
) -> Tuple[float, int, List[NucleotideVariant]]:
    """
    Globally align ref_seq to query_seq and extract SNPs and indels.

    Tries BioPython PairwiseAligner first; falls back to pure-Python NW for
    environments where BioPython is not installed.  For same-species strains
    (>90% identity) the two approaches produce equivalent results.

    Returns
    -------
    identity:         fraction of aligned positions that match
    alignment_length: total alignment length (including gaps)
    variants:         list of NucleotideVariant
    """
    ref_seq_u = ref_seq.upper()
    qry_seq_u = query_seq.upper()

    # Fast path: identical sequences (most common for same-species strains)
    if ref_seq_u == qry_seq_u:
        return 1.0, len(ref_seq), []

    try:
        ref_blocks, qry_blocks = _align_biopython(ref_seq_u, qry_seq_u)
    except ImportError:
        ref_blocks, qry_blocks = _align_nw(ref_seq_u, qry_seq_u)
    except Exception as exc:
        logger.debug("BioPython alignment failed (%s), falling back to NW", exc)
        ref_blocks, qry_blocks = _align_nw(ref_seq_u, qry_seq_u)

    if not ref_blocks:
        return 0.0, 0, []

    variants: List[NucleotideVariant] = []
    matches = 0
    aln_len = 0

    # Within each matching block: walk base-by-base to catch SNPs.
    # BioPython blocks may contain mismatches; difflib "equal" blocks do not.
    for (r_s, r_e), (q_s, q_e) in zip(ref_blocks, qry_blocks):
        block_len = r_e - r_s
        aln_len += block_len
        for i in range(block_len):
            r_base = ref_seq_u[r_s + i]
            q_base = qry_seq_u[q_s + i]
            if r_base == q_base:
                matches += 1
            else:
                variants.append(NucleotideVariant(
                    kind="SNP",
                    ref_pos=ref_offset + r_s + i,
                    ref_allele=r_base,
                    alt_allele=q_base,
                ))

    # Between consecutive matching blocks: extract indels / complex variants
    for i in range(1, len(ref_blocks)):
        prev_r_end = ref_blocks[i - 1][1]
        curr_r_start = ref_blocks[i][0]
        prev_q_end = qry_blocks[i - 1][1]
        curr_q_start = qry_blocks[i][0]

        r_gap = curr_r_start - prev_r_end   # ref bases absent from query -> deletion
        q_gap = curr_q_start - prev_q_end   # query bases absent from ref -> insertion

        aln_len += max(r_gap, q_gap)

        if r_gap > 0 and q_gap == 0:
            variants.append(NucleotideVariant(
                kind="DEL",
                ref_pos=ref_offset + prev_r_end,
                ref_allele=ref_seq_u[prev_r_end: curr_r_start],
                alt_allele="",
            ))
        elif q_gap > 0 and r_gap == 0:
            variants.append(NucleotideVariant(
                kind="INS",
                ref_pos=ref_offset + prev_r_end,
                ref_allele="",
                alt_allele=qry_seq_u[prev_q_end: curr_q_start],
            ))
        elif r_gap > 0 and q_gap > 0:
            if r_gap == q_gap:
                # Equal-length replacement: call individual SNPs
                for k in range(r_gap):
                    r_base = ref_seq_u[prev_r_end + k]
                    q_base = qry_seq_u[prev_q_end + k]
                    if r_base != q_base:
                        variants.append(NucleotideVariant(
                            kind="SNP",
                            ref_pos=ref_offset + prev_r_end + k,
                            ref_allele=r_base,
                            alt_allele=q_base,
                        ))
            else:
                # Unequal-length replacement: report as MNP/complex
                variants.append(NucleotideVariant(
                    kind="MNP",
                    ref_pos=ref_offset + prev_r_end,
                    ref_allele=ref_seq_u[prev_r_end: curr_r_start],
                    alt_allele=qry_seq_u[prev_q_end: curr_q_start],
                ))

    identity = matches / aln_len if aln_len > 0 else 0.0
    return identity, aln_len, variants


# ---------------------------------------------------------------------------
# Ortholog detection + gene-level comparison
# ---------------------------------------------------------------------------

def compare_genome_to_reference(
    query: Genome,
    ref: Genome,
    ref_index: KmerIndex,
    min_jaccard: float = 0.3,
    min_identity: float = 0.7,
    kmer_size: int = 11,
) -> GenomeComparisonResult:
    """
    Compare query against ref using the pre-built ref_index.

    Returns a GenomeComparisonResult with per-gene comparisons, deleted genes,
    structural variants, and non-coding comparisons (if sequences are available).
    """
    result = GenomeComparisonResult(query_name=query.name, ref_name=ref.name)

    # Track which ref genes have been claimed
    ref_assigned: Dict[int, str] = {}      # ref_gene_index -> query_gene_name

    query_gene_to_ref: Dict[str, Optional[int]] = {}  # query_gene_name -> ref_gene_idx

    # ------------------------------------------------------------------ #
    # 1. Map query genes to reference genes                               #
    # ------------------------------------------------------------------ #
    for qgene in query.genes:
        best_ref_idx, jaccard = ref_index.best_match(qgene.sequence, min_jaccard=min_jaccard)

        if best_ref_idx is None:
            # Novel gene: no reference ortholog
            gc = GeneComparison(
                query_gene=qgene,
                ref_gene=None,
                jaccard=0.0,
                identity=0.0,
                is_novel=True,
            )
            result.gene_comparisons[qgene.name] = gc
            query_gene_to_ref[qgene.name] = None
            continue

        ref_gene = ref.genes[best_ref_idx]

        # Alignment for verification and variant calling
        identity, aln_len, variants = align_and_call_variants(
            ref_gene.sequence, qgene.sequence
        )

        if identity < min_identity:
            # Similarity too low -- treat as novel
            gc = GeneComparison(
                query_gene=qgene,
                ref_gene=None,
                jaccard=jaccard,
                identity=identity,
                alignment_length=aln_len,
                is_novel=True,
            )
            result.gene_comparisons[qgene.name] = gc
            query_gene_to_ref[qgene.name] = None
            continue

        is_truncated = (
            abs(len(qgene.sequence) - len(ref_gene.sequence)) / len(ref_gene.sequence) > 0.1
        )

        gc = GeneComparison(
            query_gene=qgene,
            ref_gene=ref_gene,
            jaccard=jaccard,
            identity=identity,
            alignment_length=aln_len,
            variants=variants,
            is_novel=False,
            is_truncated=is_truncated,
        )
        # Key comparisons by the reference gene name so the table aligns
        key = ref_gene.name

        # Handle paralogs: if ref gene already assigned, pick the better match
        if key in result.gene_comparisons and not result.gene_comparisons[key].is_novel:
            existing = result.gene_comparisons[key]
            if identity <= existing.identity:
                # Current match is worse; store as novel (potential paralog)
                gc.is_novel = True
                gc.ref_gene = None
                result.gene_comparisons[qgene.name] = gc
                query_gene_to_ref[qgene.name] = None
                continue

        result.gene_comparisons[key] = gc
        ref_assigned[best_ref_idx] = qgene.name
        query_gene_to_ref[qgene.name] = best_ref_idx

    # ------------------------------------------------------------------ #
    # 2. Identify reference genes absent in query (deletions)             #
    # ------------------------------------------------------------------ #
    assigned_ref_indices = set(ref_assigned.keys())
    for i, rgene in enumerate(ref.genes):
        if i not in assigned_ref_indices:
            result.deleted_ref_genes.append(rgene.name)

    # ------------------------------------------------------------------ #
    # 3. Structural variant detection                                      #
    # ------------------------------------------------------------------ #
    result.structural_variants = _detect_structural_variants(
        query=query,
        ref=ref,
        query_gene_to_ref=query_gene_to_ref,
        ref_index=ref_index,
    )

    # ------------------------------------------------------------------ #
    # 4. Non-coding region comparison (if full sequences available)        #
    # ------------------------------------------------------------------ #
    if query.contigs and ref.contigs:
        result.noncoding_comparisons = _compare_noncoding(query, ref)

    return result


# ---------------------------------------------------------------------------
# Structural variant detection
# ---------------------------------------------------------------------------

def _detect_structural_variants(
    query: Genome,
    ref: Genome,
    query_gene_to_ref: Dict[str, Optional[int]],
    ref_index: KmerIndex,
) -> List[StructuralVariant]:
    """
    Detect structural variants by comparing gene-order profiles.

    Approach: build a sequence of reference gene indices corresponding to the
    ordered query genes on each contig, then look for:
      - Inversions:    consecutive decreasing runs of ref indices
      - Large gaps:    many consecutive ref genes skipped
      - Duplications:  same ref gene index appearing multiple times
    """
    svs: List[StructuralVariant] = []

    from collections import defaultdict

    # Group query genes by contig, ordered by position
    contig_genes: Dict[str, List[Gene]] = defaultdict(list)
    for g in query.genes:
        contig_genes[g.contig].append(g)
    for glist in contig_genes.values():
        glist.sort(key=lambda g: g.start)

    ref_name_to_idx = {g.name: i for i, g in enumerate(ref_index.genes)}

    for contig_id, genes_on_contig in contig_genes.items():
        profile: List[Tuple[Optional[int], str, int]] = []
        for qg in genes_on_contig:
            ref_idx = query_gene_to_ref.get(qg.name)
            if ref_idx is None:
                profile.append((None, qg.name, qg.strand))
            else:
                ref_gene = ref_index.genes[ref_idx]
                strand_match = 1 if qg.strand == ref_gene.strand else -1
                profile.append((ref_idx, qg.name, strand_match))

        # Filter to only mapped positions
        mapped = [(idx, name, sm) for idx, name, sm in profile if idx is not None]
        if len(mapped) < 2:
            continue

        ref_indices = [idx for idx, _, _ in mapped]

        # ---- Duplication detection: same ref gene appears multiple times ----
        seen: Dict[int, List[str]] = defaultdict(list)
        for idx, name, _ in mapped:
            seen[idx].append(name)
        for ref_idx, qnames in seen.items():
            if len(qnames) > 1:
                ref_gene = ref_index.genes[ref_idx]
                svs.append(StructuralVariant(
                    kind="DUPLICATION",
                    genes_ref=[ref_gene.name],
                    genes_query=qnames,
                    description=(
                        f"Reference gene '{ref_gene.name}' matched by "
                        f"{len(qnames)} query genes: {', '.join(qnames)}"
                    ),
                ))

        # ---- Inversion detection: runs of reversed synteny -------
        _detect_inversions(mapped, svs, ref_index.genes)

        # ---- Gene-block deletion detection: large consecutive gaps in ref -----
        _detect_block_deletions(ref_indices, svs, ref_index.genes)

    return svs


def _detect_inversions(
    mapped: List[Tuple[int, str, int]],
    svs: List[StructuralVariant],
    ref_genes: List[Gene],
    min_block: int = 2,
) -> None:
    """Identify runs where ref indices are decreasing (candidate inversions)."""
    i = 0
    while i < len(mapped) - 1:
        inv_block = [mapped[i]]
        j = i + 1
        while j < len(mapped):
            prev_idx = inv_block[-1][0]
            curr_idx, curr_name, curr_sm = mapped[j]
            if curr_idx < prev_idx and curr_sm == -1:
                inv_block.append(mapped[j])
                j += 1
            else:
                break
        if len(inv_block) >= min_block:
            ref_names = [ref_genes[b[0]].name for b in inv_block]
            query_names = [b[1] for b in inv_block]
            svs.append(StructuralVariant(
                kind="INVERSION",
                genes_ref=sorted(ref_names),
                genes_query=query_names,
                description=(
                    f"Block of {len(inv_block)} genes appears inverted "
                    f"(ref: {ref_names[0]}...{ref_names[-1]})"
                ),
            ))
            i = j
        else:
            i += 1


def _detect_block_deletions(
    ref_indices: List[int],
    svs: List[StructuralVariant],
    ref_genes: List[Gene],
    min_gap: int = 3,
) -> None:
    """
    Find large consecutive gaps in the reference index sequence.  Gaps of
    min_gap or more indicate a block of consecutive reference genes absent
    from the query contig.
    """
    for i in range(1, len(ref_indices)):
        gap = ref_indices[i] - ref_indices[i - 1] - 1
        if gap >= min_gap:
            missing = [ref_genes[k].name for k in range(ref_indices[i - 1] + 1, ref_indices[i])]
            svs.append(StructuralVariant(
                kind="DELETION",
                genes_ref=missing,
                genes_query=[],
                description=f"Block of {gap} consecutive reference genes absent from query contig",
            ))


# ---------------------------------------------------------------------------
# Non-coding comparison
# ---------------------------------------------------------------------------

def _compare_noncoding(
    query: Genome,
    ref: Genome,
    min_length: int = 20,
) -> Dict[str, NoncodingComparison]:
    """
    Compare intergenic regions between query and reference.
    Regions are paired by their flanking reference gene names.
    """
    ref_regions = {r.name: r for r in ref.get_noncoding_regions() if len(r.sequence) >= min_length}
    query_regions_by_flank: Dict[Tuple, NoncodingRegion] = {}
    for r in query.get_noncoding_regions():
        if len(r.sequence) >= min_length:
            query_regions_by_flank[(r.upstream_gene, r.downstream_gene)] = r

    comparisons: Dict[str, NoncodingComparison] = {}

    for rname, rregion in ref_regions.items():
        key = (rregion.upstream_gene, rregion.downstream_gene)
        qregion = query_regions_by_flank.get(key)

        if qregion is None:
            comparisons[rname] = NoncodingComparison(
                ref_region=rregion, query_region=None, identity=0.0
            )
            continue

        # Limit alignment to 5 kbp to avoid excessive memory usage
        max_len = 5000
        r_seq = rregion.sequence[:max_len]
        q_seq = qregion.sequence[:max_len]

        identity, aln_len, variants = align_and_call_variants(r_seq, q_seq)
        comparisons[rname] = NoncodingComparison(
            ref_region=rregion,
            query_region=qregion,
            identity=identity,
            alignment_length=aln_len,
            variants=variants,
        )

    return comparisons


# ---------------------------------------------------------------------------
# Single-reference: compare all query genomes against one reference
# ---------------------------------------------------------------------------

def compare_all(
    genomes: List[Genome],
    ref: Genome,
    min_jaccard: float = 0.3,
    min_identity: float = 0.7,
    kmer_size: int = 11,
    threads: int = 1,
) -> List[GenomeComparisonResult]:
    """
    Compare every genome in genomes (excluding the reference itself) against
    ref, optionally in parallel.

    Returns results in the same order as the input list (minus the ref).
    """
    import concurrent.futures

    query_genomes = [g for g in genomes if g.name != ref.name]
    if not query_genomes:
        logger.warning("No query genomes to compare (all genomes are the reference).")
        return []

    logger.info(
        "Building k-mer index (k=%d) for reference genome '%s' (%d genes)",
        kmer_size, ref.name, ref.num_genes,
    )
    ref_index = KmerIndex(ref.genes, k=kmer_size)

    logger.info(
        "Comparing %d query genome(s) against '%s' using %d thread(s)",
        len(query_genomes), ref.name, threads,
    )

    def _compare_one(q: Genome) -> GenomeComparisonResult:
        logger.debug("  Comparing '%s' vs '%s'", q.name, ref.name)
        return compare_genome_to_reference(
            query=q, ref=ref, ref_index=ref_index,
            min_jaccard=min_jaccard, min_identity=min_identity,
            kmer_size=kmer_size,
        )

    if threads == 1:
        return [_compare_one(q) for q in query_genomes]
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        return list(pool.map(_compare_one, query_genomes))


# ---------------------------------------------------------------------------
# All-vs-all: compare every unique genome pair
# ---------------------------------------------------------------------------

def compare_all_vs_all(
    genomes: List[Genome],
    min_jaccard: float = 0.3,
    min_identity: float = 0.7,
    kmer_size: int = 11,
    threads: int = 1,
) -> List[GenomeComparisonResult]:
    """
    Compare every unique unordered genome pair.

    For each pair the genome with more genes acts as reference; ties are
    broken by lexicographic name order (earlier name is reference).  Returns
    N*(N-1)/2 results for N input genomes.

    K-mer indices are built once per unique reference genome and reused
    across all pairs sharing the same reference, so the most-complete genome
    (which is reference in the most pairs) pays the index-build cost once.
    """
    import concurrent.futures

    if len(genomes) < 2:
        logger.warning("all-vs-all mode requires at least 2 genomes.")
        return []

    # Build unique pairs (i < j ensures each unordered pair is visited once).
    # The genome with more genes is the reference for each pair.
    pairs: List[Tuple[Genome, Genome]] = []
    for i in range(len(genomes)):
        for j in range(i + 1, len(genomes)):
            a, b = genomes[i], genomes[j]
            if a.num_genes > b.num_genes:
                pairs.append((a, b))    # ref=a, query=b
            elif b.num_genes > a.num_genes:
                pairs.append((b, a))    # ref=b, query=a
            else:
                # Tie: lexicographically earlier name is reference
                if a.name <= b.name:
                    pairs.append((a, b))
                else:
                    pairs.append((b, a))

    # Build one KmerIndex per unique reference genome, reused across pairs.
    ref_indices: Dict[str, KmerIndex] = {}
    for ref, _ in pairs:
        if ref.name not in ref_indices:
            logger.info(
                "Building k-mer index (k=%d) for '%s' (%d genes)",
                kmer_size, ref.name, ref.num_genes,
            )
            ref_indices[ref.name] = KmerIndex(ref.genes, k=kmer_size)

    logger.info(
        "all-vs-all: %d genome(s), %d pair(s), %d thread(s)",
        len(genomes), len(pairs), threads,
    )

    def _compare_pair(pair: Tuple[Genome, Genome]) -> GenomeComparisonResult:
        ref, query = pair
        logger.debug("Comparing '%s' (query) vs '%s' (ref)", query.name, ref.name)
        return compare_genome_to_reference(
            query=query,
            ref=ref,
            ref_index=ref_indices[ref.name],
            min_jaccard=min_jaccard,
            min_identity=min_identity,
            kmer_size=kmer_size,
        )

    if threads == 1:
        return [_compare_pair(p) for p in pairs]
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        return list(pool.map(_compare_pair, pairs))