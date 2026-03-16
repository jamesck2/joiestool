"""
Output report generation.

Single-reference mode produces three output files:

  1. gene_table.tsv       Wide-format gene presence/variant table.
                          Rows = genes (reference-ordered, with gap rows for
                          novel query genes).  Two columns per genome:
                          {genome}_status and {genome}_variants.

  2. gene_diversity.tsv   Long-format per-genome per-gene statistics.

  3. variant_summary.tsv  One row per structural or notable nucleotide event.

All-vs-all mode produces three output files:

  1. gene_presence_absence.tsv  Pan-genome presence/absence matrix anchored
                                to the most-complete genome.  One column pair
                                per genome; anchor genome column is first.

  2. pairwise_summary.tsv       One row per genome pair with aggregate counts
                                and mean pairwise identity.

  3. variant_summary.tsv        Same format as single-reference mode; includes
                                ref_genome and query_genome columns.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional

from joiestool.compare import GenomeComparisonResult, GeneComparison
from joiestool.genome import Genome

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_status(gc: Optional[GeneComparison]) -> str:
    if gc is None:
        return "absent"
    if gc.is_novel:
        return "novel"
    if gc.is_truncated:
        return "truncated"
    if not gc.variants:
        return "identical"
    return "variant"


def _fmt_variants(gc: Optional[GeneComparison]) -> str:
    if gc is None or gc.ref_gene is None:
        return ""
    return gc.variant_summary


# ---------------------------------------------------------------------------
# 1a. Wide gene table (single-reference mode)
# ---------------------------------------------------------------------------

def write_gene_table(
    ref: Genome,
    results: List[GenomeComparisonResult],
    outpath: Path,
) -> None:
    """
    Write the wide-format gene table (single-reference mode).

    Columns:
      gene_id | ref_contig | ref_start | ref_end | ref_strand
      | {g1}_status | {g1}_variants | {g2}_status | {g2}_variants | ...

    Rows:
      - All reference genes, sorted by contig then position.
      - After all reference genes: novel query genes (present in >= 1 query
        but absent from reference), inserted in query-genome order.
    """
    genome_names = [r.query_name for r in results]

    # Build a lookup: ref_gene_name -> GeneComparison for each query
    comparisons: Dict[tuple, Optional[GeneComparison]] = {}
    novel_genes: Dict[str, GeneComparison] = {}  # gene_id -> GeneComparison

    for res in results:
        for key, gc in res.gene_comparisons.items():
            if gc.is_novel or gc.ref_gene is None:
                # key is the query gene name
                novel_genes.setdefault(key, gc)
            else:
                comparisons[(res.query_name, key)] = gc

    # Reference genes sorted by position
    ref_gene_order = sorted(
        ref.genes, key=lambda g: (g.contig, g.start)
    )

    header = [
        "gene_id", "ref_contig", "ref_start", "ref_end", "ref_strand",
    ]
    for gname in genome_names:
        header += [f"{gname}_status", f"{gname}_variants"]

    rows = []

    # Reference-anchored rows
    for rgene in ref_gene_order:
        row = [
            rgene.name,
            rgene.contig,
            rgene.start + 1,    # convert to 1-based for users
            rgene.end,
            "+" if rgene.strand == 1 else "-",
        ]
        for gname in genome_names:
            gc = comparisons.get((gname, rgene.name))
            if gc is None:
                for res in results:
                    if res.query_name == gname and rgene.name in res.deleted_ref_genes:
                        break
                else:
                    gc = None  # no info
            row += [_fmt_status(gc), _fmt_variants(gc)]
        rows.append(row)

    # Novel query gene rows
    for novel_name, gc in sorted(novel_genes.items()):
        qg = gc.query_gene
        row = [
            f"novel:{novel_name}",
            qg.contig,
            qg.start + 1,
            qg.end,
            "+" if qg.strand == 1 else "-",
        ]
        for gname in genome_names:
            is_this_genome = any(
                res.query_name == gname and novel_name in res.gene_comparisons
                for res in results
            )
            if is_this_genome:
                row += ["novel", ""]
            else:
                row += ["absent", ""]
        rows.append(row)

    _write_tsv(outpath, header, rows)
    logger.info("Gene table -> %s  (%d rows)", outpath.name, len(rows))


# ---------------------------------------------------------------------------
# 1b. Pan-genome gene presence/absence matrix (all-vs-all mode)
# ---------------------------------------------------------------------------

def write_gene_presence_absence(
    anchor: Genome,
    anchor_results: List[GenomeComparisonResult],
    all_genome_names: List[str],
    outpath: Path,
) -> None:
    """
    Write a pan-genome gene presence/absence matrix (all-vs-all mode).

    The anchor genome (typically the most gene-rich) defines the row set.
    Columns cover every genome: anchor first, then all others sorted
    alphabetically.  The anchor genome's own column is always 'identical'
    for its own genes.

    Rows:
      - Anchor genes sorted by contig then position.
      - Novel genes found in non-anchor genomes (absent from anchor),
        appended at the end.

    Parameters
    ----------
    anchor:
        Genome whose gene set defines the rows.
    anchor_results:
        Comparison results where anchor.name == result.ref_name.
    all_genome_names:
        Names of all genomes in the dataset (determines column order).
    outpath:
        Destination file path.
    """
    # Anchor column first; remaining genomes in sorted order
    other_names = sorted(n for n in all_genome_names if n != anchor.name)
    col_order = [anchor.name] + other_names

    # Build lookup: (query_name, ref_gene_name) -> GeneComparison
    comparisons: Dict[tuple, Optional[GeneComparison]] = {}
    novel_genes: Dict[str, GeneComparison] = {}

    for res in anchor_results:
        for key, gc in res.gene_comparisons.items():
            if gc.is_novel or gc.ref_gene is None:
                novel_genes.setdefault(key, gc)
            else:
                comparisons[(res.query_name, key)] = gc

    ref_gene_order = sorted(anchor.genes, key=lambda g: (g.contig, g.start))

    header = ["gene_id", "ref_contig", "ref_start", "ref_end", "ref_strand"]
    for gname in col_order:
        header += [f"{gname}_status", f"{gname}_variants"]

    rows = []

    # Anchor-gene rows
    for rgene in ref_gene_order:
        row = [
            rgene.name,
            rgene.contig,
            rgene.start + 1,
            rgene.end,
            "+" if rgene.strand == 1 else "-",
        ]
        for gname in col_order:
            if gname == anchor.name:
                row += ["identical", ""]
                continue
            gc = comparisons.get((gname, rgene.name))
            if gc is None:
                for res in anchor_results:
                    if res.query_name == gname and rgene.name in res.deleted_ref_genes:
                        break
                else:
                    gc = None
            row += [_fmt_status(gc), _fmt_variants(gc)]
        rows.append(row)

    # Novel gene rows (present in non-anchor genomes, absent from anchor)
    for novel_name, gc in sorted(novel_genes.items()):
        qg = gc.query_gene
        row = [
            f"novel:{novel_name}",
            qg.contig,
            qg.start + 1,
            qg.end,
            "+" if qg.strand == 1 else "-",
        ]
        for gname in col_order:
            if gname == anchor.name:
                row += ["absent", ""]
                continue
            is_this_genome = any(
                res.query_name == gname and novel_name in res.gene_comparisons
                for res in anchor_results
            )
            row += (["novel", ""] if is_this_genome else ["absent", ""])
        rows.append(row)

    _write_tsv(outpath, header, rows)
    logger.info(
        "Gene presence/absence matrix -> %s  (%d genes, %d genomes)",
        outpath.name, len(rows), len(col_order),
    )


# ---------------------------------------------------------------------------
# 2a. Long-format diversity table (single-reference mode)
# ---------------------------------------------------------------------------

def write_gene_diversity(
    ref: Genome,
    results: List[GenomeComparisonResult],
    outpath: Path,
) -> None:
    """
    Write the long-format gene diversity table (single-reference mode).

    Columns:
      genome | gene_id | is_present | is_novel | is_partial | is_truncated
      | identity | num_snps | num_indels | num_variants
      | ref_gene_length | query_gene_length | notes
    """
    header = [
        "genome", "gene_id",
        "is_present", "is_novel", "is_partial", "is_truncated",
        "identity", "num_snps", "num_indels", "num_variants",
        "ref_gene_length", "query_gene_length",
        "notes",
    ]
    rows = []

    for res in results:
        ref_gene_by_name = {g.name: g for g in ref.genes}
        handled = set()

        for key, gc in res.gene_comparisons.items():
            if gc.is_novel or gc.ref_gene is None:
                gene_id = key
                rows.append([
                    res.query_name, gene_id,
                    "yes", "yes",
                    "yes" if gc.query_gene.partial else "no",
                    "no",
                    "", "", "", "",
                    "",
                    gc.query_gene.length,
                    "novel gene, no reference ortholog",
                ])
                continue

            ref_gene = gc.ref_gene
            handled.add(ref_gene.name)
            rows.append([
                res.query_name,
                ref_gene.name,
                "yes", "no",
                "yes" if gc.query_gene.partial else "no",
                "yes" if gc.is_truncated else "no",
                f"{gc.identity:.4f}",
                gc.num_snps,
                gc.num_indels,
                len(gc.variants),
                ref_gene.length,
                gc.query_gene.length,
                "",
            ])

        # Reference genes absent in query
        for del_name in res.deleted_ref_genes:
            if del_name in handled:
                continue
            rg = ref_gene_by_name.get(del_name)
            rows.append([
                res.query_name,
                del_name,
                "no", "no", "no", "no",
                "", "", "", "",
                rg.length if rg else "",
                "",
                "absent from query genome",
            ])

    _write_tsv(outpath, header, rows)
    logger.info("Gene diversity table -> %s  (%d rows)", outpath.name, len(rows))


# ---------------------------------------------------------------------------
# 2b. Pairwise summary table (all-vs-all mode)
# ---------------------------------------------------------------------------

def write_pairwise_summary(
    results: List[GenomeComparisonResult],
    outpath: Path,
) -> None:
    """
    Write a per-pair aggregate summary table (all-vs-all mode).

    One row per GenomeComparisonResult.  All counts and identity values are
    relative to the reference genome used for that comparison.

    Columns:
      ref_genome | query_genome
      | shared_genes | identical_genes | variant_genes | truncated_genes
      | novel_query_genes | deleted_ref_genes
      | total_snps | total_indels | total_variants
      | structural_variants | mean_identity
    """
    header = [
        "ref_genome", "query_genome",
        "shared_genes", "identical_genes", "variant_genes", "truncated_genes",
        "novel_query_genes", "deleted_ref_genes",
        "total_snps", "total_indels", "total_variants",
        "structural_variants", "mean_identity",
    ]
    rows = []

    for res in results:
        assigned = [
            gc for gc in res.gene_comparisons.values()
            if not gc.is_novel and gc.ref_gene is not None
        ]
        novel = sum(1 for gc in res.gene_comparisons.values() if gc.is_novel)
        identical = sum(1 for gc in assigned if not gc.variants and not gc.is_truncated)
        variant_genes = sum(1 for gc in assigned if gc.variants)
        truncated = sum(1 for gc in assigned if gc.is_truncated)
        total_snps = sum(gc.num_snps for gc in assigned)
        total_indels = sum(gc.num_indels for gc in assigned)
        total_variants = sum(len(gc.variants) for gc in assigned)
        mean_id = (
            sum(gc.identity for gc in assigned) / len(assigned)
            if assigned else 0.0
        )
        rows.append([
            res.ref_name,
            res.query_name,
            len(assigned),
            identical,
            variant_genes,
            truncated,
            novel,
            len(res.deleted_ref_genes),
            total_snps,
            total_indels,
            total_variants,
            len(res.structural_variants),
            f"{mean_id:.4f}",
        ])

    _write_tsv(outpath, header, rows)
    logger.info("Pairwise summary -> %s  (%d rows)", outpath.name, len(rows))


# ---------------------------------------------------------------------------
# 3. Variant summary (both modes)
# ---------------------------------------------------------------------------

def write_variant_summary(
    results: List[GenomeComparisonResult],
    outpath: Path,
) -> None:
    """
    Write the per-comparison variant summary.

    Columns:
      ref_genome | query_genome | variant_type | gene_or_region
      | ref_start | ref_end | details
    """
    header = [
        "ref_genome", "query_genome", "variant_type", "gene_or_region",
        "ref_start", "ref_end", "details",
    ]
    rows = []

    for res in results:
        # Gene-level nucleotide variants
        for key, gc in res.gene_comparisons.items():
            if gc.is_novel or gc.ref_gene is None:
                rows.append([
                    res.ref_name, res.query_name, "GENE_GAIN",
                    gc.query_gene.name,
                    gc.query_gene.start + 1,
                    gc.query_gene.end,
                    f"Novel gene in query (len={gc.query_gene.length}bp)",
                ])
                continue

            for var in gc.variants:
                rows.append([
                    res.ref_name, res.query_name, var.kind,
                    gc.ref_gene.name,
                    var.ref_pos + 1,
                    var.ref_pos + max(1, len(var.ref_allele)),
                    str(var),
                ])

        # Deleted reference genes
        for del_name in res.deleted_ref_genes:
            rows.append([
                res.ref_name, res.query_name, "GENE_LOSS",
                del_name,
                "", "",
                "Reference gene absent from query genome",
            ])

        # Structural variants
        for sv in res.structural_variants:
            genes_str = ",".join(sv.genes_ref) if sv.genes_ref else "(none)"
            rows.append([
                res.ref_name, res.query_name, sv.kind,
                genes_str,
                "", "",
                sv.description,
            ])

        # Non-coding variants (SNPs and indels only; skip identical regions)
        for rname, nc in res.noncoding_comparisons.items():
            if nc.query_region is None:
                rows.append([
                    res.ref_name, res.query_name, "NONCODING_LOSS",
                    rname,
                    nc.ref_region.start + 1,
                    nc.ref_region.end,
                    "Intergenic region absent from query",
                ])
                continue
            for var in nc.variants:
                rows.append([
                    res.ref_name, res.query_name,
                    f"NONCODING_{var.kind}",
                    rname,
                    nc.ref_region.start + var.ref_pos + 1,
                    nc.ref_region.start + var.ref_pos + max(1, len(var.ref_allele)),
                    str(var),
                ])

    _write_tsv(outpath, header, rows)
    logger.info("Variant summary -> %s  (%d rows)", outpath.name, len(rows))


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _write_tsv(path: Path, header: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)