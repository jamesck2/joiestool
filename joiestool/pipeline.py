"""
Main pipeline: load -> predict -> compare -> report.

Supports two comparison modes:

  single-reference (default):
    One genome acts as reference; all others are compared against it.
    Outputs: gene_table.tsv, gene_diversity.tsv, variant_summary.tsv

  all-vs-all (--all-vs-all):
    Every unique genome pair is compared once.  The most-complete genome
    (most genes) acts as the anchor for the pan-genome presence/absence
    matrix and as reference for each pair it participates in.
    Outputs: gene_presence_absence.tsv, pairwise_summary.tsv, variant_summary.tsv
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def run_pipeline(args) -> None:
    """
    Execute the full pipeline from parsed CLI arguments.

    Steps
    -----
    1. Load all genomes from the input directory.
    2. Predict genes for any genome_fasta inputs (pyrodigal).
    3. Run comparisons (single-reference or all-vs-all).
    4. Write output tables.
    """
    from joiestool.genome import load_genomes_from_directory
    from joiestool.predict import predict_genes_parallel

    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir
    reference_name: Optional[str] = args.reference
    fmt: str = args.format
    threads: int = args.threads
    min_identity: float = args.min_identity
    kmer_size: int = args.kmer_size
    meta: bool = args.meta
    all_vs_all: bool = args.all_vs_all

    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path '{input_dir}' is not a directory.")

    if all_vs_all and reference_name:
        logger.warning(
            "--reference is ignored in all-vs-all mode; "
            "the most-complete genome is used as the pan-genome anchor."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Load genomes                                                      #
    # ------------------------------------------------------------------ #
    format_hint = None if fmt == "auto" else {
        "fasta": "genome_fasta",
        "genbank": "genbank",
        "genes": "genes_fasta",
    }[fmt]

    genomes = load_genomes_from_directory(input_dir, format_hint=format_hint)

    if len(genomes) < 2:
        raise ValueError(
            f"Need at least 2 genome files; found {len(genomes)} in '{input_dir}'."
        )

    # ------------------------------------------------------------------ #
    # 2. Predict genes where needed                                        #
    # ------------------------------------------------------------------ #
    genomes = predict_genes_parallel(genomes, meta=meta, threads=threads)

    empty = [g.name for g in genomes if g.num_genes == 0]
    if empty:
        logger.warning(
            "Genomes with 0 genes (skipped): %s", ", ".join(empty)
        )
        genomes = [g for g in genomes if g.num_genes > 0]

    logger.info(
        "Loaded %d genomes; gene counts: %s",
        len(genomes),
        ", ".join(f"{g.name}={g.num_genes}" for g in genomes),
    )

    # ------------------------------------------------------------------ #
    # 3. Compare + 4. Write reports                                        #
    # ------------------------------------------------------------------ #
    logger.info("Writing output reports to '%s/'", output_dir)

    if all_vs_all:
        _run_all_vs_all(genomes, output_dir, min_identity, kmer_size, threads)
    else:
        _run_single_reference(
            genomes, output_dir, reference_name, min_identity, kmer_size, threads
        )

    logger.info("Done. Output in '%s/'", output_dir)


# ---------------------------------------------------------------------------
# Single-reference mode
# ---------------------------------------------------------------------------

def _run_single_reference(genomes, output_dir, reference_name, min_identity, kmer_size, threads):
    from joiestool.compare import compare_all
    from joiestool.report import write_gene_table, write_gene_diversity, write_variant_summary

    # Select reference
    if reference_name:
        ref_matches = [g for g in genomes if g.name == reference_name]
        if not ref_matches:
            raise ValueError(
                f"Reference genome '{reference_name}' not found. "
                f"Available: {[g.name for g in genomes]}"
            )
        ref = ref_matches[0]
    else:
        ref = max(genomes, key=lambda g: g.num_genes)
        logger.info("Auto-selected reference: '%s' (%d genes)", ref.name, ref.num_genes)

    results = compare_all(
        genomes=genomes,
        ref=ref,
        min_jaccard=min_identity * 0.5,
        min_identity=min_identity,
        kmer_size=kmer_size,
        threads=threads,
    )

    if not results:
        logger.warning("No comparison results generated.")
        return

    write_gene_table(ref, results, output_dir / "gene_table.tsv")
    write_gene_diversity(ref, results, output_dir / "gene_diversity.tsv")
    write_variant_summary(results, output_dir / "variant_summary.tsv")

    total_novel = sum(
        sum(1 for gc in r.gene_comparisons.values() if gc.is_novel)
        for r in results
    )
    total_deleted = sum(len(r.deleted_ref_genes) for r in results)
    total_svs = sum(len(r.structural_variants) for r in results)
    total_variants = sum(
        sum(len(gc.variants) for gc in r.gene_comparisons.values())
        for r in results
    )
    logger.info(
        "Summary: %d query genome(s), %d nucleotide variants, "
        "%d gene gains, %d gene losses, %d structural variants",
        len(results), total_variants, total_novel, total_deleted, total_svs,
    )


# ---------------------------------------------------------------------------
# All-vs-all mode
# ---------------------------------------------------------------------------

def _run_all_vs_all(genomes, output_dir, min_identity, kmer_size, threads):
    from joiestool.compare import compare_all_vs_all
    from joiestool.report import (
        write_gene_presence_absence,
        write_pairwise_summary,
        write_variant_summary,
    )

    results = compare_all_vs_all(
        genomes=genomes,
        min_jaccard=min_identity * 0.5,
        min_identity=min_identity,
        kmer_size=kmer_size,
        threads=threads,
    )

    if not results:
        logger.warning("No comparison results generated.")
        return

    # Anchor = genome with most genes; the same criterion used to assign
    # reference within pairs, so anchor is ref_name in all its pairs.
    anchor = max(genomes, key=lambda g: g.num_genes)
    logger.info("Pan-genome anchor: '%s' (%d genes)", anchor.name, anchor.num_genes)

    anchor_results = [r for r in results if r.ref_name == anchor.name]
    all_genome_names = [g.name for g in genomes]

    write_gene_presence_absence(
        anchor, anchor_results, all_genome_names,
        output_dir / "gene_presence_absence.tsv",
    )
    write_pairwise_summary(results, output_dir / "pairwise_summary.tsv")
    write_variant_summary(results, output_dir / "variant_summary.tsv")

    # Summary stats across all pairs
    mean_ids = []
    for res in results:
        assigned = [
            gc for gc in res.gene_comparisons.values()
            if not gc.is_novel and gc.ref_gene is not None
        ]
        if assigned:
            mean_ids.append(sum(gc.identity for gc in assigned) / len(assigned))

    logger.info(
        "Summary: %d pair(s) compared; "
        "mean pairwise identity = %.4f (range %.4f-%.4f)",
        len(results),
        sum(mean_ids) / len(mean_ids) if mean_ids else 0.0,
        min(mean_ids) if mean_ids else 0.0,
        max(mean_ids) if mean_ids else 0.0,
    )