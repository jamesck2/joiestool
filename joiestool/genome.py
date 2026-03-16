"""
Genome data structures and file I/O.

Supports three input formats:
  - genome_fasta: nucleotide FASTA of the full genome (genes will be predicted)
  - genbank:      GenBank file with CDS features annotated
  - genes_fasta:  nucleotide FASTA of pre-predicted genes in Prodigal header format
"""
from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Regex to detect Prodigal-format FASTA headers:
#   >contig_N # start # end # strand # ID=N;partial=...
_PRODIGAL_HEADER_RE = re.compile(
    r"^.+\s+#\s+\d+\s+#\s+\d+\s+#\s+[+-]?1\s+#\s+ID=\d+"
)


@dataclass
class Gene:
    """A single gene (CDS) with its nucleotide sequence and genomic coordinates."""

    name: str           # unique ID (e.g. contig_1, locus_tag, etc.)
    sequence: str       # nucleotide sequence, coding strand 5'→3'
    contig: str         # parent contig / chromosome ID
    start: int          # 0-based start position on contig
    end: int            # 0-based end position on contig (exclusive)
    strand: int         # +1 or -1
    partial: bool = False  # gene at contig edge (incomplete)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def seq_hash(self) -> str:
        """MD5 digest of the nucleotide sequence for fast exact-match comparison."""
        return hashlib.md5(self.sequence.upper().encode()).hexdigest()


@dataclass
class NoncodingRegion:
    """An intergenic region between two consecutive genes on a contig."""

    name: str           # e.g. "contig_gene3-gene4"
    sequence: str       # nucleotide sequence
    contig: str
    start: int          # 0-based
    end: int            # 0-based, exclusive
    upstream_gene: Optional[str] = None   # gene name to the 5' side
    downstream_gene: Optional[str] = None # gene name to the 3' side


@dataclass
class Genome:
    """A genome with its gene predictions and optionally the full contig sequences."""

    name: str
    genes: List[Gene] = field(default_factory=list)
    contigs: Optional[Dict[str, str]] = None  # contig_id -> sequence
    source_file: str = ""
    source_format: str = ""

    @property
    def num_genes(self) -> int:
        return len(self.genes)

    @property
    def gene_by_name(self) -> Dict[str, Gene]:
        return {g.name: g for g in self.genes}

    def get_noncoding_regions(self) -> List[NoncodingRegion]:
        """
        Extract intergenic regions between consecutive gene predictions.
        Only available when full contig sequences were loaded.
        """
        if not self.contigs:
            return []

        regions: List[NoncodingRegion] = []

        # Group genes by contig, sorted by start position
        from collections import defaultdict
        contig_genes: Dict[str, List[Gene]] = defaultdict(list)
        for gene in self.genes:
            contig_genes[gene.contig].append(gene)

        for contig_id, seq in self.contigs.items():
            local_genes = sorted(contig_genes.get(contig_id, []), key=lambda g: g.start)

            # Region before first gene
            if local_genes and local_genes[0].start > 0:
                regions.append(NoncodingRegion(
                    name=f"{contig_id}_5prime",
                    sequence=seq[: local_genes[0].start],
                    contig=contig_id,
                    start=0,
                    end=local_genes[0].start,
                    downstream_gene=local_genes[0].name,
                ))

            # Regions between genes
            for i in range(len(local_genes) - 1):
                g1, g2 = local_genes[i], local_genes[i + 1]
                if g2.start > g1.end:
                    regions.append(NoncodingRegion(
                        name=f"{contig_id}_{g1.name}-{g2.name}",
                        sequence=seq[g1.end: g2.start],
                        contig=contig_id,
                        start=g1.end,
                        end=g2.start,
                        upstream_gene=g1.name,
                        downstream_gene=g2.name,
                    ))

            # Region after last gene
            if local_genes and local_genes[-1].end < len(seq):
                regions.append(NoncodingRegion(
                    name=f"{contig_id}_3prime",
                    sequence=seq[local_genes[-1].end:],
                    contig=contig_id,
                    start=local_genes[-1].end,
                    end=len(seq),
                    upstream_gene=local_genes[-1].name,
                ))

        return regions


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

def detect_format(filepath: Path) -> str:
    """
    Determine the file format.  Returns one of:
      'genbank', 'genes_fasta', 'genome_fasta'
    """
    ext = filepath.suffix.lower()

    if ext in {".gbk", ".gb", ".genbank", ".gbff"}:
        return "genbank"

    if ext in {".fa", ".fna", ".fasta", ".ffn"}:
        # Peek at the first header to distinguish genes from genome
        try:
            with open(filepath) as fh:
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if _PRODIGAL_HEADER_RE.match(line):
                            return "genes_fasta"
                        # .ffn extension is commonly used for gene nucleotides
                        if ext == ".ffn":
                            return "genes_fasta"
                        return "genome_fasta"
        except Exception:
            pass
        return "genome_fasta"

    # Fallback: try to open as GenBank
    try:
        with open(filepath) as fh:
            if fh.read(5) == "LOCUS":
                return "genbank"
    except Exception:
        pass

    raise ValueError(
        f"Cannot determine format of '{filepath}'. "
        "Expected .fa/.fna/.fasta/.ffn (FASTA) or .gbk/.gb (GenBank)."
    )


# ---------------------------------------------------------------------------
# Format-specific loaders
# ---------------------------------------------------------------------------


def _parse_fasta(filepath) -> list:
    """
    Minimal pure-Python FASTA parser (no BioPython required).
    Returns list of (seq_id, description, sequence) triples.
    """
    records = []
    seq_id = desc = ""
    buf = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id:
                    records.append((seq_id, desc, "".join(buf)))
                header = line[1:]
                parts = header.split(None, 1)
                seq_id = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                buf = []
            else:
                buf.append(line)
    if seq_id:
        records.append((seq_id, desc, "".join(buf)))
    return records

def _load_contigs_from_fasta(filepath: Path) -> Dict[str, str]:
    """Read all sequences from a FASTA file into a contig_id → sequence dict."""
    try:
        from Bio import SeqIO
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(filepath, "fasta")}
    except ImportError:
        return {seq_id: seq for seq_id, _, seq in _parse_fasta(filepath)}


def _load_from_genbank(filepath: Path) -> Tuple[List[Gene], Dict[str, str]]:
    """Load genes and full sequences from a GenBank file."""
    from Bio import SeqIO

    genes: List[Gene] = []
    contigs: Dict[str, str] = {}
    gene_count: Dict[str, int] = {}

    for record in SeqIO.parse(filepath, "genbank"):
        contigs[record.id] = str(record.seq)
        cds_idx = 0
        for feature in record.features:
            if feature.type != "CDS":
                continue
            cds_idx += 1
            # Prefer locus_tag, fall back to protein_id, then positional name
            locus_tags = feature.qualifiers.get("locus_tag", [])
            protein_ids = feature.qualifiers.get("protein_id", [])
            name = (
                locus_tags[0] if locus_tags
                else protein_ids[0] if protein_ids
                else f"{record.id}_{cds_idx}"
            )
            try:
                seq = str(feature.extract(record.seq))
            except Exception as exc:
                logger.warning("Skipping CDS %s in %s: %s", name, record.id, exc)
                continue

            strand = feature.location.strand if feature.location.strand is not None else 1
            genes.append(Gene(
                name=name,
                sequence=seq,
                contig=record.id,
                start=int(feature.location.start),
                end=int(feature.location.end),
                strand=strand,
                partial=bool(
                    feature.qualifiers.get("pseudo")
                    or ">" in str(feature.location)
                    or "<" in str(feature.location)
                ),
            ))

    return genes, contigs


def _load_from_genes_fasta(filepath: Path) -> List[Gene]:
    """
    Load pre-predicted genes from a FASTA file with Prodigal-format headers.

    Prodigal header format:
        >{contig}_{N} # {start} # {end} # {strand} # ID={N};partial={pp};...

    Works with or without BioPython installed.
    """
    genes: List[Gene] = []

    # Build an iterator of (seq_id, description, sequence) — prefer Bio for robustness
    try:
        from Bio import SeqIO
        records = [(r.id, r.description, str(r.seq)) for r in SeqIO.parse(filepath, "fasta")]
    except ImportError:
        records = _parse_fasta(filepath)

    for gene_id, desc, seq_str in records:
        # The full header is seq_id + " " + rest; join them back for parsing
        full_header = f"{gene_id} {desc}" if desc else gene_id
        parts = full_header.split("#")
        if len(parts) >= 4:
            try:
                start = int(parts[1].strip()) - 1  # convert to 0-based
                end = int(parts[2].strip())
                strand = int(parts[3].strip())
                meta = parts[4].strip() if len(parts) > 4 else ""
                partial_match = re.search(r"partial=([01]{2})", meta)
                partial = bool(partial_match and partial_match.group(1) != "00")
            except (ValueError, IndexError):
                start, end, strand, partial = 0, len(seq_str), 1, False
        else:
            start, end, strand, partial = 0, len(seq_str), 1, False

        contig_match = re.match(r"^(.+)_(\d+)$", gene_id)
        contig = contig_match.group(1) if contig_match else gene_id

        genes.append(Gene(
            name=gene_id,
            sequence=seq_str,
            contig=contig,
            start=start,
            end=end,
            strand=strand,
            partial=partial,
        ))

    return genes


# ---------------------------------------------------------------------------
# Public loader
# ---------------------------------------------------------------------------

def load_genome(filepath: Path, format_hint: Optional[str] = None) -> Genome:
    """
    Load a genome from *filepath*, auto-detecting the format unless *format_hint*
    is provided ('genbank', 'genes_fasta', or 'genome_fasta').

    For genome_fasta inputs, genes will be empty — call predict.predict_genes()
    to populate them afterwards.
    """
    fmt = format_hint or detect_format(filepath)
    name = filepath.stem

    logger.debug("Loading '%s' as format '%s'", filepath.name, fmt)

    if fmt == "genbank":
        genes, contigs = _load_from_genbank(filepath)
        logger.debug("  → %d CDS features, %d contigs", len(genes), len(contigs))
        return Genome(
            name=name, genes=genes, contigs=contigs,
            source_file=str(filepath), source_format="genbank",
        )

    if fmt == "genes_fasta":
        genes = _load_from_genes_fasta(filepath)
        logger.debug("  → %d genes", len(genes))
        return Genome(
            name=name, genes=genes, contigs=None,
            source_file=str(filepath), source_format="genes_fasta",
        )

    if fmt == "genome_fasta":
        contigs = _load_contigs_from_fasta(filepath)
        logger.debug("  → %d contigs (genes not yet predicted)", len(contigs))
        return Genome(
            name=name, genes=[], contigs=contigs,
            source_file=str(filepath), source_format="genome_fasta",
        )

    raise ValueError(f"Unknown format '{fmt}'")


def load_genomes_from_directory(
    directory: Path,
    format_hint: Optional[str] = None,
    extensions: Optional[Tuple[str, ...]] = None,
) -> List[Genome]:
    """
    Load all genome files found in *directory*.

    Returns genomes sorted by name for reproducibility.
    """
    if extensions is None:
        extensions = (".fa", ".fna", ".fasta", ".ffn", ".gbk", ".gb", ".genbank", ".gbff")

    paths = sorted(
        p for p in directory.iterdir()
        if p.is_file() and p.suffix.lower() in extensions
    )
    if not paths:
        raise FileNotFoundError(f"No genome files found in '{directory}'.")

    logger.info("Found %d genome files in '%s'", len(paths), directory)
    genomes: List[Genome] = []
    for p in paths:
        try:
            genomes.append(load_genome(p, format_hint=format_hint))
        except Exception as exc:
            logger.warning("Skipping '%s': %s", p.name, exc)

    return genomes