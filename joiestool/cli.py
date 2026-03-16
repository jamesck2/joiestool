"""
Command-line interface for Joie's Tool.

Usage
-----
    joiestool <input_dir> [options]
    python -m joiestool <input_dir> [options]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from joiestool import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="joiestool",
        description=(
            "Gene-level comparison of highly similar bacterial genomes.\n\n"
            "Compares complete or near-complete genomes within a species, reporting\n"
            "nucleotide variants (SNPs/indels), gene presence/absence, and structural\n"
            "variants (inversions, duplications, rearrangements) for each genome\n"
            "relative to a chosen reference. Genomes can be closed (single contigs/\n"
            "scaffolds) or fragmented (multi-contig/scaffold).\n\n"
            "Input files (one per genome) may be:\n"
            " - Nucleotide FASTA (.fa/.fna/.fasta; genes will be predicted using prodigal)\n"
            " - GenBank          (.gbk/.gb/.gbff; existing CDS features are used)\n"
            " - Genes FASTA      (.ffn/.fna/.fasta; genes directly used)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_dir",
        type=Path,
        metavar="INPUT_DIR",
        help="Directory containing one genome file per strain.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("joiestool_out"),
        metavar="DIR",
        help="Output directory (default: joiestool_out/).",
    )
    parser.add_argument(
        "-r", "--reference",
        type=str,
        default=None,
        metavar="NAME",
        help=(
            "Reference genome name (stem of the file, without extension). "
            "If omitted, the genome with the most predicted genes is chosen. "
            "Ignored when --all-vs-all is set."
        ),
    )
    parser.add_argument(
        "--all-vs-all",
        action="store_true",
        default=False,
        help=(
            "Compare every unique genome pair instead of all genomes against a "
            "single reference.  For each pair the more gene-rich genome is used "
            "as reference.  The most-complete genome anchors the pan-genome "
            "presence/absence matrix.  Outputs: gene_presence_absence.tsv, "
            "pairwise_summary.tsv, variant_summary.tsv."
        ),
    )
    parser.add_argument(
        "-f", "--format",
        choices=["auto", "fasta", "genbank", "genes"],
        default="auto",
        help=(
            "Input file format. 'auto' detects per-file (default). "
            "'fasta'=genome FASTA, 'genbank'=GenBank, 'genes'=genes FASTA."
        ),
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        metavar="N",
        help="Number of worker threads for parallel genome comparison (default: 1).",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.70,
        metavar="FLOAT",
        help=(
            "Minimum nucleotide identity (0-1) for ortholog assignment "
            "(default: 0.70).  Lower values increase sensitivity but may "
            "include non-orthologous gene pairs."
        ),
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=11,
        metavar="K",
        help="K-mer size for fast ortholog candidate screening (default: 11).",
    )
    parser.add_argument(
        "--meta",
        action="store_true",
        help=(
            "Use metagenomic (pre-trained) Prodigal models for gene prediction. "
            "Automatically applied to genomes < 100 kbp."
        ),
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug-level logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Configure logging (always to stderr; clean format for pipelines)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)-8s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        level=log_level,
        stream=sys.stderr,
    )
    # Quiet noisy third-party loggers unless verbose
    if not args.verbose:
        for noisy in ("Bio", "urllib3", "requests"):
            logging.getLogger(noisy).setLevel(logging.WARNING)

    logger = logging.getLogger("joiestool")
    logger.info("Joie's Tool %s", __version__)

    try:
        from joiestool.pipeline import run_pipeline
        run_pipeline(args)
    except (FileNotFoundError, NotADirectoryError, ValueError) as exc:
        logger.error("%s", exc)
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted.")
        sys.exit(130)
    except Exception as exc:
        logger.exception("Unexpected error: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()