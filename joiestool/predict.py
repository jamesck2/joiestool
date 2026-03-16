"""
Gene prediction using Pyrodigal (Cython bindings to Prodigal).

Handles single-contig genomes, draft genomes (multiple contigs), and thread-safe
parallel prediction across multiple genomes.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

from joiestool.genome import Gene, Genome

logger = logging.getLogger(__name__)


def predict_genes(
    genome: Genome,
    meta: bool = False,
    min_gene_len: int = 90,
    translation_table: int = 11,
) -> Genome:
    """
    Run Prodigal gene prediction on a genome whose contigs are loaded but genes
    list is empty.  Populates and returns the genome (mutates in place).

    Parameters
    ----------
    genome:
        A Genome with .contigs populated and .genes == [].
    meta:
        Use metagenomic mode (pre-trained models). Useful for draft genomes
        or when training on the genome itself is unreliable (< 100 kbp).
    min_gene_len:
        Minimum gene length in nucleotides (default 90 = Prodigal default).
    translation_table:
        Genetic code for translation (default 11 = bacteria/archaea).
    """
    try:
        import pyrodigal
    except ImportError as exc:
        raise ImportError(
            "pyrodigal is required for gene prediction from FASTA genomes. "
            "Install it with:  pip install pyrodigal"
        ) from exc

    if genome.contigs is None:
        raise ValueError(f"Genome '{genome.name}' has no contig sequences to predict from.")

    contigs: Dict[str, str] = genome.contigs
    total_len = sum(len(s) for s in contigs.values())

    # Choose mode: single (train on genome) vs meta (use pre-trained models)
    # Auto-switch to meta for small inputs where training is unreliable
    use_meta = meta or total_len < 100_000

    logger.debug(
        "Predicting genes in '%s' (mode=%s, contigs=%d, total_bp=%d)",
        genome.name,
        "meta" if use_meta else "single",
        len(contigs),
        total_len,
    )

    gene_finder = pyrodigal.GeneFinder(
        meta=use_meta,
        min_gene=min_gene_len,
    )

    if not use_meta:
        # Train on all contig sequences concatenated (single mode)
        training_seqs = [bytes(seq, "ascii") for seq in contigs.values()]
        gene_finder.train(*training_seqs, translation_table=translation_table)

    predicted: List[Gene] = []
    for contig_id, contig_seq in contigs.items():
        try:
            preds = gene_finder.find_genes(bytes(contig_seq, "ascii"))
        except Exception as exc:
            logger.warning(
                "Prediction failed on contig '%s' of genome '%s': %s",
                contig_id, genome.name, exc,
            )
            continue

        for i, pred in enumerate(preds, start=1):
            try:
                nuc_seq = pred.sequence()
            except Exception:
                # Fallback: extract manually from contig using coordinates
                b, e = pred.begin - 1, pred.end   # convert to 0-based
                raw = contig_seq[b:e]
                nuc_seq = raw if pred.strand == 1 else _revcomp(raw)

            predicted.append(Gene(
                name=f"{contig_id}_{i}",
                sequence=nuc_seq,
                contig=contig_id,
                start=pred.begin - 1,   # Prodigal is 1-based
                end=pred.end,
                strand=pred.strand,
                partial=pred.partial_begin or pred.partial_end,
            ))

    genome.genes = predicted
    logger.debug("  → %d genes predicted in '%s'", len(predicted), genome.name)
    return genome


def predict_genes_parallel(
    genomes: List[Genome],
    meta: bool = False,
    min_gene_len: int = 90,
    translation_table: int = 11,
    threads: int = 1,
) -> List[Genome]:
    """
    Predict genes for multiple genomes, optionally in parallel.

    Only genomes with contigs but no genes (source_format == 'genome_fasta')
    are processed; others are returned unchanged.
    """
    needs_prediction = [
        g for g in genomes if g.source_format == "genome_fasta" and not g.genes
    ]
    already_annotated = [
        g for g in genomes if g not in needs_prediction
    ]

    if not needs_prediction:
        return genomes

    logger.info(
        "Predicting genes for %d genome(s) using %d thread(s)",
        len(needs_prediction),
        threads,
    )

    import concurrent.futures

    def _predict(g: Genome) -> Genome:
        return predict_genes(g, meta=meta, min_gene_len=min_gene_len,
                             translation_table=translation_table)

    if threads == 1:
        results = [_predict(g) for g in needs_prediction]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
            results = list(pool.map(_predict, needs_prediction))

    # Merge back in original order
    result_map = {g.name: g for g in results}
    return [result_map.get(g.name, g) for g in genomes]


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]