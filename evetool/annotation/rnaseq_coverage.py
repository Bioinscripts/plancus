#!/usr/bin/env python3
"""
rnaseq_coverage.py — RNA-seq coverage over each EVE locus and its flanking region.

Runs when RNA-seq is given, either on the command line (``--rnaseq-r1/--rnaseq-r2`` or
``--rnaseq-bam``) or in the species TSV (FASTQ columns, or any cell ending in ``.bam``).

Two ways in:
  * FASTQ reads: they are aligned with HISAT2 (splice-aware) against the non-redundant
    EVE windows (locus +/- flank, overlapping windows merged so a read is never split
    between two copies of the same sequence), or against the whole genome with
    ``reference="genome"``.
  * BAM: used as is, no alignment. The BAM may be aligned to the whole genome or to
    context windows; the coordinate system is worked out from the header (see
    ``BamReference``), so a BAM made outside the pipeline works too.

Depth, depth-weighted base quality and strand-split depth come from a single pass over
the BAM with pysam, which allows what ``samtools mpileup`` cannot:
  * strand is taken per fragment (read 2 flipped), then turned into transcript strand
    with the library type (fr / rf), inferred from the host genes of the GFF when
    ``library="auto"``;
  * intron gaps (N) do not count as coverage;
  * reads touching an EVE are split into *exonic* (aligned bases inside the EVE, i.e. the
    EVE is in the transcript) and *intron-spanning* (only the N gap crosses it: the EVE sits
    in an intron of an expressed host gene and is spliced out). The two mean opposite
    things, and depth alone cannot tell them apart.

Outputs (in ``10_rnaseq_expressao/``):
  eve_coverage_summary.tsv            one row per sample x EVE locus
  coverage_bins/<sample>/<window>.tsv binned depth, quality and strand depth per window
  figures/<window>_rnaseq_coverage.png/.pdf
  rnaseq_depth_per_base.tsv           per-base depth per Cluster_N, read by the context figures
  bam/<sample>.bam                    only when reads were aligned here
  samples.json                        the BAMs used, their index and library type (read by the curation page)
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

N_BINS = 600            # max plotted points per window
NM_NATIVE = 3           # exonic reads with NM <= this look native, not recruited paralogues
MIN_STRAND_READS = 200  # informative reads needed to infer the library type

# flag bits skipped: unmapped, secondary, QC fail, duplicate, supplementary
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


def is_bam(path: str) -> bool:
    """
    True for a BAM file, whatever its name: BGZF-compressed, or uncompressed (a BAM saved
    through a browser often loses its gzip layer and its extension) - both start with 'BAM\\1'.
    """
    import gzip
    try:
        with open(path, "rb") as fh:
            inicio = fh.read(4)
        if inicio == b"BAM\1":
            return True
        if inicio[:2] == b"\x1f\x8b":
            with gzip.open(path, "rb") as fh:
                return fh.read(4) == b"BAM\1"
    except (OSError, EOFError):
        pass
    return False


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

@dataclass
class EveLocus:
    locus: str          # "contig:start-end", 0-based half-open (as in eves_classified.tsv)
    contig: str
    start: int
    end: int
    strand: str         # "+", "-" or "."
    label: str          # representative ORF (or cluster id)
    group: str = ""
    family: str = ""
    product: str = ""


@dataclass
class Window:
    name: str
    contig: str
    start: int          # 0-based half-open
    end: int
    loci: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Sample:
    name: str
    bam: str | None = None
    r1: str | None = None
    r2: str | None = None


@dataclass
class WindowProfile:
    depth: np.ndarray
    plus: np.ndarray        # read-1 orientation until the library type is applied
    minus: np.ndarray
    qsum: np.ndarray        # summed base quality, divided by depth for the mean
    in_ref: np.ndarray      # position present in the BAM reference
    reads: dict             # locus -> {"exonic", "exonic_nm_ok", "spanning"}


# ---------------------------------------------------------------------------
# EVE loci and windows
# ---------------------------------------------------------------------------

_LOCUS_RE = re.compile(r"^(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+)$")


def _strands_from_bed(results_dir: str) -> dict:
    """ORF id -> strand, from the direct-hit BED written by the anchoring step."""
    strands = {}
    for nome in ("eve_hits_diretos.sorted.bed", "eve_hits_diretos.bed"):
        caminho = os.path.join(results_dir, "06_anotacao_genomica", nome)
        if not os.path.exists(caminho):
            continue
        with open(caminho) as fh:
            for linha in fh:
                c = linha.rstrip("\n").split("\t")
                if len(c) >= 6:
                    strands[c[3]] = c[5]
        break
    return strands


def load_loci(results_dir: str, species: str | None = None,
              classified_tsv: str | None = None) -> list[EveLocus]:
    """
    One EveLocus per genomic locus. Read from eves_classified.tsv (the `locus` column groups
    overlapping ORFs of one insertion); falls back to the merged cluster BED of older runs.
    In multi-species mode the table lives at the run root and is filtered by species.
    """
    candidatos = [classified_tsv] if classified_tsv else [
        os.path.join(results_dir, "eves_classified.tsv"),
        os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(results_dir))), "eves_classified.tsv"),
    ]
    tabela = next((c for c in candidatos if c and os.path.exists(c)), None)
    strands = _strands_from_bed(results_dir)
    loci: dict[str, EveLocus] = {}

    if tabela:
        with open(tabela, newline="") as fh:
            linhas = list(csv.DictReader(fh, delimiter="\t"))
        if linhas and "locus" in linhas[0]:
            for r in linhas:
                if species and r.get("species") and r["species"] != species:
                    continue
                m = _LOCUS_RE.match(r.get("locus", ""))
                if not m:
                    continue
                orf = (r.get("orf") or r.get("query", "")).split("|")[-1]
                representante = str(r.get("locus_representative", "")).lower() == "true"
                if r["locus"] in loci and not representante:
                    continue
                na = lambda v: "" if v in (None, "n/a") else v
                loci[r["locus"]] = EveLocus(
                    locus=r["locus"], contig=m["contig"], start=int(m["start"]), end=int(m["end"]),
                    strand=strands.get(orf, "."), label=orf, group=na(r.get("group")),
                    family=na(r.get("viral_family")), product=na(r.get("viral_product")))

    if not loci:
        bed = os.path.join(results_dir, "06_anotacao_genomica", "eve_hits_diretos.merged.bed")
        if not os.path.exists(bed):
            return []
        print(f"[RNAseq] No locus table; using the cluster BED: {bed}")
        with open(bed) as fh:
            for linha in fh:
                c = linha.rstrip("\n").split("\t")
                if len(c) < 4:
                    continue
                locus = f"{c[0]}:{c[1]}-{c[2]}"
                loci[locus] = EveLocus(locus, c[0], int(c[1]), int(c[2]),
                                       c[5] if len(c) > 5 and c[5] in "+-" else ".", c[3])

    return sorted(loci.values(), key=lambda l: (l.contig, l.start))


def build_windows(loci: list[EveLocus], flank: int, contig_sizes: dict) -> list[Window]:
    """Locus +/- flank, clipped to the contig; overlapping windows are merged into one."""
    janelas: list[Window] = []
    for l in loci:
        ini = max(0, l.start - flank)
        fim = l.end + flank
        if l.contig in contig_sizes:
            fim = min(fim, contig_sizes[l.contig])
        ultima = janelas[-1] if janelas else None
        if ultima and ultima.contig == l.contig and ini <= ultima.end:
            ultima.end = max(ultima.end, fim)
            ultima.loci.append(l)
        else:
            janelas.append(Window("", l.contig, ini, fim, [l]))
    for i, w in enumerate(janelas, 1):
        w.name = f"win{i}_{w.contig}_{w.start + 1}-{w.end}"
    return janelas


def contig_sizes_for(genome: str | None, results_dir: str) -> dict:
    """Contig lengths from the genome .fai, else from the chromosome sizes of the anchoring step."""
    tamanhos = {}
    for caminho in ([f"{genome}.fai"] if genome else []) + [
            os.path.join(results_dir, "06_anotacao_genomica", "cromossomos.sizes")]:
        if caminho and os.path.exists(caminho):
            with open(caminho) as fh:
                for linha in fh:
                    c = linha.split("\t")
                    if len(c) >= 2:
                        tamanhos[c[0]] = int(c[1])
            if tamanhos:
                break
    if not tamanhos and genome and os.path.exists(genome):
        import pysam
        with pysam.FastaFile(genome) as fa:      # builds the .fai when missing
            tamanhos = dict(zip(fa.references, fa.lengths))
    return tamanhos


# ---------------------------------------------------------------------------
# BAM coordinate system
# ---------------------------------------------------------------------------

_RANGE_RE = re.compile(r"(\d+)-(\d+)(?:\([+.-]\))?$")


class BamReference:
    """
    Places each BAM reference sequence on the genome.

      * a reference named like a genome contig -> whole-genome BAM, offset 0;
      * a reference whose name ends in "<contig><sep><start>-<end>" (the pipeline's
        "Cluster_1::NC_074710.1:10045006-10088674", or "win1_NC_074710.1_10045007-10094415")
        -> context window. Whether start is 0-based (bedtools) or 1-based (samtools faidx)
        is decided by the reference length: end-start for the first, end-start+1 for the second.

    References that match neither are ignored and reported.
    """

    def __init__(self, path: str, contigs: set, index: str | None = None):
        import pysam
        self.path = path
        self.index = index
        self.segments: dict[str, list] = defaultdict(list)   # contig -> [(g0, g1, ref)]
        self.ignored: list[str] = []
        with pysam.AlignmentFile(path) as fh:
            refs = list(zip(fh.references, fh.lengths))
        ordenados = sorted(contigs, key=len, reverse=True)
        for nome, comp in refs:
            if nome in contigs:
                self.segments[nome].append((0, comp, nome))
                continue
            m = _RANGE_RE.search(nome)
            prefixo = nome[:m.start()].rstrip(":_|") if m else ""
            contig = next((c for c in ordenados if prefixo == c or
                           (prefixo.endswith(c) and prefixo[-len(c) - 1] in ":_|")), None) if m else None
            if not contig:
                self.ignored.append(nome)
                continue
            s, e = int(m[1]), int(m[2])
            if e - s == comp:
                g0 = s
            elif e - s + 1 == comp:
                g0 = s - 1
            else:
                self.ignored.append(nome)
                continue
            self.segments[contig].append((g0, g0 + comp, nome))

    def overlapping(self, contig: str, start: int, end: int) -> list:
        """Segments overlapping [start, end), largest overlap first."""
        segs = [s for s in self.segments.get(contig, []) if s[0] < end and s[1] > start]
        return sorted(segs, key=lambda s: min(s[1], end) - max(s[0], start), reverse=True)


def prepare_bam(bam: str, pasta: str) -> tuple[str, str]:
    """
    (BAM to read, its index). Uses the index beside the BAM when present (.bam.bai, .bai, .csi).
    Otherwise the index is built in our own output folder, so the user's BAM folder is never
    written to; an uncompressed BAM cannot be indexed, so it is first recompressed (BGZF) there.
    """
    raiz = bam[:-4] if bam.lower().endswith(".bam") else bam
    for c in (bam + ".bai", bam + ".csi", raiz + ".bai", raiz + ".csi"):
        if os.path.exists(c) and os.path.getmtime(c) >= os.path.getmtime(bam):
            return bam, c
    os.makedirs(pasta, exist_ok=True)
    nome = os.path.basename(raiz)
    with open(bam, "rb") as fh:
        comprimido = fh.read(2) == b"\x1f\x8b"
    if not comprimido:
        copia = os.path.join(pasta, nome + ".bam")
        if not (os.path.exists(copia) and os.path.getmtime(copia) >= os.path.getmtime(bam)):
            print(f"[RNAseq] {nome}: uncompressed BAM, recompressing to {copia}")
            subprocess.run(["samtools", "view", "-b", "-o", copia, bam], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        bam = copia
    destino = os.path.join(pasta, nome + ".bam.bai")
    if not (os.path.exists(destino) and os.path.getmtime(destino) >= os.path.getmtime(bam)):
        print(f"[RNAseq] {nome}: indexing -> {destino}")
        r = subprocess.run(["samtools", "index", "-o", destino, bam], capture_output=True, text=True)
        if r.returncode != 0:
            linhas = [l for l in r.stderr.splitlines() if "no version information" not in l]
            raise RuntimeError(linhas[-1] if linhas else "samtools index failed")
    return bam, destino


# ---------------------------------------------------------------------------
# Alignment (FASTQ -> BAM)
# ---------------------------------------------------------------------------

def sample_name(r1: str) -> str:
    nome = os.path.basename(r1)
    nome = re.sub(r"\.(fastq|fq)(\.gz)?$", "", nome, flags=re.I)
    return re.sub(r"([._](R?1))$", "", nome) or nome


def write_window_fasta(genome: str, windows: list[Window], destino: str) -> str:
    """Non-redundant window FASTA; names are "contig:start-end" (0-based), so the BAM places itself."""
    import pysam
    with pysam.FastaFile(genome) as fa, open(destino, "w") as out:
        for w in windows:
            seq = fa.fetch(w.contig, w.start, w.end)
            out.write(f">{w.contig}:{w.start}-{w.end}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    return destino


def _run_pipe(cmd_align: list, bam: str, threads: int, log: str) -> None:
    with open(log, "w") as err:
        alin = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=err)
        sort = subprocess.run(["samtools", "sort", "-@", str(max(1, threads // 2)), "-o", bam, "-"],
                              stdin=alin.stdout, stderr=err)
        alin.stdout.close()
        alin.wait()
    if alin.returncode != 0 or sort.returncode != 0:
        raise RuntimeError(f"alignment failed; see {log}")
    subprocess.run(["samtools", "index", bam], check=True)


def align_reads(sample: Sample, ref_fasta: str, outdir: str, threads: int = 4,
                read_type: str = "sr", force: bool = False) -> str:
    """
    Align one sample and return a sorted, indexed BAM. Short reads go to HISAT2 (RNA-seq
    reads cross introns), then bowtie2 --local, then minimap2 -ax sr; long reads to
    minimap2 -ax splice / splice:hq. An existing BAM is reused unless `force`.
    """
    pasta_bam = os.path.join(outdir, "bam")
    pasta_idx = os.path.join(outdir, "index")
    os.makedirs(pasta_bam, exist_ok=True)
    os.makedirs(pasta_idx, exist_ok=True)
    bam = os.path.join(pasta_bam, f"{sample.name}.bam")
    log = os.path.join(pasta_bam, f"{sample.name}.align.log")
    if os.path.exists(bam) and os.path.getsize(bam) > 0 and os.path.exists(bam + ".bai") and not force:
        print(f"[RNAseq] {sample.name}: BAM already present, reusing {bam}")
        return bam

    prefixo = os.path.join(pasta_idx, os.path.basename(ref_fasta).rsplit(".", 1)[0])
    pareado = bool(sample.r2)

    if read_type in ("ont", "pb", "hifi"):
        if not shutil.which("minimap2"):
            raise RuntimeError("minimap2 is needed for long reads (conda install -c bioconda minimap2)")
        preset = "splice" if read_type == "ont" else "splice:hq"
        cmd = ["minimap2", "-ax", preset, "-t", str(threads), ref_fasta, sample.r1]
        aligner = f"minimap2 -ax {preset}"
    elif shutil.which("hisat2") and shutil.which("hisat2-build"):
        if not os.path.exists(prefixo + ".1.ht2"):
            print(f"[RNAseq] Building the HISAT2 index for {os.path.basename(ref_fasta)}...")
            subprocess.run(["hisat2-build", "-p", str(threads), ref_fasta, prefixo], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cmd = ["hisat2", "-p", str(threads), "-x", prefixo, "--no-unal"]
        cmd += ["-1", sample.r1, "-2", sample.r2] if pareado else ["-U", sample.r1]
        aligner = "HISAT2"
    elif shutil.which("bowtie2") and shutil.which("bowtie2-build"):
        print("[RNAseq] [warning] HISAT2 not found; bowtie2 is not splice-aware, so reads "
              "crossing introns will be soft-clipped or lost.")
        if not os.path.exists(prefixo + ".1.bt2"):
            subprocess.run(["bowtie2-build", "--threads", str(threads), ref_fasta, prefixo], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cmd = ["bowtie2", "-p", str(threads), "-x", prefixo, "--very-sensitive-local", "--no-unal"]
        cmd += ["-1", sample.r1, "-2", sample.r2] if pareado else ["-U", sample.r1]
        aligner = "bowtie2 --local"
    elif shutil.which("minimap2"):
        cmd = ["minimap2", "-ax", "sr", "-t", str(threads), ref_fasta, sample.r1] + ([sample.r2] if pareado else [])
        aligner = "minimap2 -ax sr"
    else:
        raise RuntimeError("no aligner found (install hisat2: conda install -c bioconda hisat2 samtools)")

    print(f"[RNAseq] {sample.name}: aligning with {aligner} ({'paired' if pareado else 'single'}-end)...")
    _run_pipe(cmd, bam, threads, log)
    return bam


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

def _aligned_blocks(read):
    """(ref_start, ref_end, query_start) of each aligned block; deletions and N gaps excluded."""
    r, q = read.reference_start, 0
    for op, n in read.cigartuples:
        if op in (0, 7, 8):          # M = X
            yield r, r + n, q
            r += n
            q += n
        elif op in (2, 3):           # D N
            r += n
        elif op in (1, 4):           # I S
            q += n


def profile_window(bam_ref: BamReference, window: Window, min_mapq: int = 0,
                   exon_strand: np.ndarray | None = None, strand_votes: list | None = None) -> WindowProfile:
    """
    Per-base depth, quality and strand depth across one window, from whichever BAM
    references cover it. Where two references overlap, each position is taken from one only.
    `strand_votes` ([same, opposite]) collects read-1 vs host-exon strand for the library test.
    """
    import pysam
    n = window.length
    depth = np.zeros(n)
    plus = np.zeros(n)
    minus = np.zeros(n)
    qsum = np.zeros(n)
    dono = np.full(n, -1, dtype=np.int16)        # which segment owns each position
    leituras = {l.locus: {"exonic": 0, "exonic_nm_ok": 0, "spanning": 0} for l in window.loci}
    loci_idx = [(l.locus, l.start - window.start, l.end - window.start) for l in window.loci]

    segs = bam_ref.overlapping(window.contig, window.start, window.end)
    with pysam.AlignmentFile(bam_ref.path, index_filename=bam_ref.index) as fh:
        for k, (g0, g1, ref) in enumerate(segs):
            a, b = max(g0, window.start) - window.start, min(g1, window.end) - window.start
            livre = dono[a:b] == -1
            dono[a:b][livre] = k
            meu = dono == k
            if not meu.any():
                continue
            desloc = g0 - window.start               # BAM position + desloc = window index
            d_tot, d_mais, d_menos = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
            q_k = np.zeros(n)
            r_ini = max(0, window.start - g0)
            r_fim = min(g1, window.end) - g0
            for read in fh.fetch(ref, r_ini, r_fim):
                if read.flag & _SKIP_FLAGS or read.mapping_quality < min_mapq or not read.cigartuples:
                    continue
                reverso = read.is_reverse != (read.is_paired and read.is_read2)   # fragment = read-1 strand
                quals = read.query_qualities
                blocos = []
                for rs, re_, qs in _aligned_blocks(read):
                    x0, x1 = rs + desloc, re_ + desloc
                    c0, c1 = max(x0, 0), min(x1, n)
                    if c0 >= c1:
                        continue
                    blocos.append((c0, c1))
                    d_tot[c0] += 1; d_tot[c1] -= 1
                    alvo = d_menos if reverso else d_mais
                    alvo[c0] += 1; alvo[c1] -= 1
                    if quals is not None:
                        q0 = qs + (c0 - x0)
                        q_k[c0:c1] += np.asarray(quals[q0:q0 + (c1 - c0)], dtype=float)
                if not blocos:
                    continue
                if strand_votes is not None and exon_strand is not None:
                    meio = (blocos[0][0] + blocos[0][1]) // 2
                    fita = exon_strand[meio]
                    if fita in (1, -1):
                        strand_votes[0 if (fita == -1) == reverso else 1] += 1
                ini_read, fim_read = blocos[0][0], blocos[-1][1]
                for locus, ls, le in loci_idx:
                    if not meu[min(max((ls + le) // 2, 0), n - 1)]:
                        continue                       # this locus is counted from another segment
                    if any(c0 < le and c1 > ls for c0, c1 in blocos):
                        leituras[locus]["exonic"] += 1
                        if read.has_tag("NM") and read.get_tag("NM") <= NM_NATIVE:
                            leituras[locus]["exonic_nm_ok"] += 1
                    elif ini_read < le and fim_read > ls:
                        leituras[locus]["spanning"] += 1
            depth[meu] = np.cumsum(d_tot)[:n][meu]
            plus[meu] = np.cumsum(d_mais)[:n][meu]
            minus[meu] = np.cumsum(d_menos)[:n][meu]
            qsum[meu] = q_k[meu]
    return WindowProfile(depth, plus, minus, qsum, dono >= 0, leituras)


def load_host_features(gff: str | None, windows: list[Window]) -> dict:
    """{window name: {"genes": [(start, end, strand, name)], "exons": [(start, end, strand)]}} (0-based)."""
    saida = {w.name: {"genes": [], "exons": []} for w in windows}
    if not gff or not os.path.exists(gff):
        return saida
    por_contig = defaultdict(list)
    for w in windows:
        por_contig[w.contig].append(w)
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#"):
                continue
            c = linha.rstrip("\n").split("\t")
            if len(c) < 9 or c[0] not in por_contig or c[2] not in ("gene", "exon"):
                continue
            s, e = int(c[3]) - 1, int(c[4])
            for w in por_contig[c[0]]:
                if s >= w.end or e <= w.start:
                    continue
                if c[2] == "exon":
                    saida[w.name]["exons"].append((s, e, c[6]))
                else:
                    nome = ""
                    for campo in c[8].split(";"):
                        if campo.startswith(("Name=", "gene=", "gene_name=")):
                            nome = campo.split("=", 1)[1]
                            break
                    saida[w.name]["genes"].append((s, e, c[6], nome))
    return saida


def _exon_strand_array(window: Window, exons: list) -> np.ndarray:
    """+1 / -1 where only exons of that strand lie, 2 where both, 0 elsewhere."""
    mais = np.zeros(window.length, dtype=bool)
    menos = np.zeros(window.length, dtype=bool)
    for s, e, fita in exons:
        a, b = max(s, window.start) - window.start, min(e, window.end) - window.start
        (mais if fita == "+" else menos)[a:b] = True
    arr = np.zeros(window.length, dtype=np.int8)
    arr[mais] = 1
    arr[menos] = -1
    arr[mais & menos] = 2
    return arr


def infer_library(votes: list) -> str:
    """fr if read 1 follows the host exon strand, rf if it opposes it, unstranded otherwise."""
    mesma, oposta = votes
    total = mesma + oposta
    if total < MIN_STRAND_READS:
        return "unstranded"
    frac = mesma / total
    return "fr" if frac >= 0.75 else "rf" if frac <= 0.25 else "unstranded"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def _status(eve_mean, breadth1, flank_mean, exonic, spanning, ref_frac) -> str:
    if ref_frac < 0.5:
        return "Not in BAM reference"
    if eve_mean >= 10 and breadth1 >= 50:
        return "Transcribed (high)"
    if eve_mean >= 1 and breadth1 >= 30:
        return "Transcribed (moderate)"
    if spanning >= 3 and exonic <= max(2, 0.1 * spanning):
        return "Intronic, spliced out (host gene expressed)"
    if eve_mean > 0 and flank_mean >= 5 * max(eve_mean, 0.2):
        return "Flank only (host region active, EVE silent)"
    if eve_mean > 0:
        return "Low / residual"
    if flank_mean >= 2:
        return "Flank only (host region active, EVE silent)"
    return "No expression detected"


def summarise(sample: str, window: Window, prof: WindowProfile, flank: int, library: str) -> list[dict]:
    linhas = []
    ocupado = np.zeros(window.length, dtype=bool)
    for l in window.loci:
        ocupado[l.start - window.start:l.end - window.start] = True
    stranded = library in ("fr", "rf")

    def media(arr, mask):
        return float(arr[mask].mean()) if mask.any() else 0.0

    for l in window.loci:
        a, b = l.start - window.start, l.end - window.start
        m = np.zeros(window.length, dtype=bool)
        m[a:b] = True
        eve = m & prof.in_ref
        up = np.zeros(window.length, dtype=bool)
        down = np.zeros(window.length, dtype=bool)
        up[max(0, a - flank):a] = True
        down[b:min(window.length, b + flank)] = True
        # flanks never include another EVE, so a neighbour cannot inflate the host signal
        up &= prof.in_ref & ~ocupado
        down &= prof.in_ref & ~ocupado
        if l.strand == "-":
            up, down = down, up                  # upstream follows the EVE strand

        d = prof.depth
        eve_mean = media(d, eve)
        breadth1 = 100.0 * float((d[eve] >= 1).mean()) if eve.any() else 0.0
        breadth5 = 100.0 * float((d[eve] >= 5).mean()) if eve.any() else 0.0
        cobertos = eve & (d > 0)
        qual = float(prof.qsum[cobertos].sum() / d[cobertos].sum()) if cobertos.any() else 0.0
        flank_mask = up | down
        flank_mean = media(d, flank_mask)
        r = prof.reads[l.locus]

        sense = antisense = frac_sense = ""
        if stranded and l.strand in "+-" and eve.any():
            s_arr, a_arr = (prof.plus, prof.minus) if l.strand == "+" else (prof.minus, prof.plus)
            sense, antisense = round(media(s_arr, eve), 2), round(media(a_arr, eve), 2)
            frac_sense = round(sense / (sense + antisense), 3) if (sense + antisense) > 0 else ""

        linhas.append({
            "sample": sample, "window": window.name, "locus": l.locus,
            "representative_orf": l.label, "group": l.group, "family": l.family,
            "product": l.product, "strand": l.strand, "eve_length": l.end - l.start,
            "eve_mean_depth": round(eve_mean, 2),
            "eve_median_depth": round(float(np.median(d[eve])), 2) if eve.any() else 0.0,
            "eve_max_depth": int(d[eve].max()) if eve.any() else 0,
            "eve_breadth_1x_pct": round(breadth1, 1), "eve_breadth_5x_pct": round(breadth5, 1),
            "eve_mean_base_quality": round(qual, 1),
            "eve_sense_mean_depth": sense, "eve_antisense_mean_depth": antisense,
            "eve_sense_fraction": frac_sense,
            "upstream_mean_depth": round(media(d, up), 2),
            "downstream_mean_depth": round(media(d, down), 2),
            "flank_mean_depth": round(flank_mean, 2),
            "log2_eve_vs_flank": round(float(np.log2((eve_mean + 1) / (flank_mean + 1))), 3),
            "exonic_reads": r["exonic"], f"exonic_reads_nm_le{NM_NATIVE}": r["exonic_nm_ok"],
            "intron_spanning_reads": r["spanning"],
            "eve_in_reference_pct": round(100.0 * float(prof.in_ref[a:b].mean()), 1),
            "library_type": library,
            "status": _status(eve_mean, breadth1, flank_mean, r["exonic"], r["spanning"],
                              float(prof.in_ref[a:b].mean())),
        })
    return linhas


def bin_profile(window: Window, prof: WindowProfile) -> dict:
    """Binned arrays (<= N_BINS): mean depth, strand depth, depth-weighted quality, in-reference fraction."""
    n = window.length
    nb = min(N_BINS, n)
    cortes = np.array_split(np.arange(n), nb)
    ini = np.array([c[0] for c in cortes])
    fim = np.array([c[-1] + 1 for c in cortes])
    soma = lambda arr: np.add.reduceat(arr, ini)
    tam = fim - ini
    dsum = soma(prof.depth)
    return {
        "start": ini + window.start, "end": fim + window.start,
        "depth": dsum / tam, "plus": soma(prof.plus) / tam, "minus": soma(prof.minus) / tam,
        "quality": np.divide(soma(prof.qsum), dsum, out=np.zeros(nb), where=dsum > 0),
        "in_ref": soma(prof.in_ref.astype(float)) / tam,
    }


def export_bins(bins: dict, destino: str, stranded: bool) -> None:
    with open(destino, "w") as fh:
        fh.write("bin\tstart\tend\tmean_depth\tmean_base_quality\tplus_strand_depth\t"
                 "minus_strand_depth\tin_reference_fraction\n")
        for i in range(len(bins["depth"])):
            mais = f"{bins['plus'][i]:.2f}" if stranded else ""
            menos = f"{bins['minus'][i]:.2f}" if stranded else ""
            fh.write(f"{i + 1}\t{bins['start'][i] + 1}\t{bins['end'][i]}\t{bins['depth'][i]:.2f}\t"
                     f"{bins['quality'][i]:.1f}\t{mais}\t{menos}\t{bins['in_ref'][i]:.2f}\n")


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

# Strand hues: categorical slots 1-2 of the reference palette. Quality is magnitude, so it
# gets one sequential hue (light = low Phred, dark = high); EVE red and host-gene gray are
# kept out of both.
COR_MAIS, COR_MENOS = "#2a78d6", "#eb6834"
COR_EVE, COR_GENE, COR_FORA = "#e34948", "#8b8a86", "#d9d8d3"
TINTA, TINTA_2, GRADE = "#0b0b0b", "#52514e", "#e3e5e0"


def _quality_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("phred", ["#cfe8dc", "#5fb894", "#13795b", "#0a3d2e"])


def _pack_rows(itens, largura_min):
    """Greedy interval packing: list of rows, each a list of items (start, end, ...)."""
    linhas = []
    for it in sorted(itens, key=lambda x: x[0]):
        alvo = next((l for l in linhas if l[-1][1] + largura_min < it[0]), None)
        (alvo.append(it) if alvo is not None else linhas.append([it]))
    return linhas


def plot_window(window: Window, perfis: dict, host: dict, bibliotecas: dict, destino: str) -> str:
    """
    One figure per window. Per sample: total depth (bars coloured by mean base quality) and,
    for stranded libraries, the transcript-strand split (+ above, - mirrored below). The
    bottom track holds the EVE loci and the host genes. EVEs are shaded in every panel.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.patches import Patch, Rectangle

    amostras = list(perfis)
    linhas_g = _pack_rows([(s, e, f, nm) for s, e, f, nm in host.get("genes", [])], window.length * 0.01)
    linhas_e = _pack_rows([(l.start, l.end, l) for l in window.loci], window.length * 0.07)
    altura_trilha = 0.55 + 0.32 * (len(linhas_e) + len(linhas_g))

    razoes, papeis = [], []
    for a in amostras:
        razoes.append(2.6); papeis.append(("depth", a))
        if bibliotecas[a] in ("fr", "rf"):
            razoes.append(1.5); papeis.append(("strand", a))
    razoes.append(altura_trilha); papeis.append(("track", None))

    fig = plt.figure(figsize=(13, sum(razoes) * 0.95 + 1.1))
    gs = fig.add_gridspec(len(razoes), 2, height_ratios=razoes, width_ratios=[60, 1],
                          hspace=0.12, wspace=0.02)
    eixos = []
    for i in range(len(razoes)):
        eixos.append(fig.add_subplot(gs[i, 0], sharex=eixos[0] if eixos else None))

    cmap, norm = _quality_cmap(), Normalize(vmin=10, vmax=40, clip=True)
    teto = max((b["depth"].max() for b in perfis.values()), default=1)
    teto_fita = max((max(b["plus"].max(), b["minus"].max()) for a, b in perfis.items()
                     if bibliotecas[a] in ("fr", "rf")), default=1)

    def sombrear(ax):
        for l in window.loci:
            ax.axvspan(l.start, l.end, color=COR_EVE, alpha=0.13, lw=0, zorder=0)

    def fora_da_referencia(ax, b):
        for s, e, f in zip(b["start"], b["end"], b["in_ref"]):
            if f < 0.5:
                ax.axvspan(s, e, color=COR_FORA, alpha=0.6, lw=0, zorder=0.5)

    def estilo(ax):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(TINTA_2)
        ax.tick_params(colors=TINTA_2, labelsize=8)
        ax.grid(axis="y", color=GRADE, lw=0.7, zorder=0)
        ax.set_axisbelow(True)

    primeiro_depth = None
    for ax, (papel, a) in zip(eixos, papeis):
        if papel == "track":
            continue
        b = perfis[a]
        largura = b["end"] - b["start"]
        sombrear(ax)
        fora_da_referencia(ax, b)
        estilo(ax)
        if papel == "depth":
            ax.bar(b["start"], b["depth"], width=largura, align="edge", lw=0, zorder=2,
                   color=cmap(norm(b["quality"])))
            ax.set_yscale("symlog", linthresh=1, linscale=0.4)
            ax.set_ylim(0, teto * 1.6 + 1)
            ax.set_ylabel("depth", fontsize=9, color=TINTA_2)
            ax.text(0.005, 0.96, a, transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
                    fontweight="bold", color=TINTA)
            if primeiro_depth is None:
                primeiro_depth = ax
                cax = fig.add_subplot(gs[papeis.index(("depth", a)), 1])
                cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
                cb.set_label("mean base\nquality (Phred)", fontsize=8, color=TINTA_2)
                cb.ax.tick_params(labelsize=7, colors=TINTA_2)
                cb.outline.set_visible(False)
        else:
            ax.bar(b["start"], b["plus"], width=largura, align="edge", lw=0, color=COR_MAIS, zorder=2)
            ax.bar(b["start"], -b["minus"], width=largura, align="edge", lw=0, color=COR_MENOS, zorder=2)
            ax.axhline(0, color=TINTA_2, lw=0.6, zorder=3)
            ax.set_yscale("symlog", linthresh=1, linscale=0.4)
            lim = teto_fita * 1.6 + 1
            ax.set_ylim(-lim, lim)
            # mirrored axis: a few decades, labelled as magnitudes on both sides of zero
            dec = max(1, int(np.floor(np.log10(lim))))
            marcas = sorted({10 ** dec, 10 ** max(1, (dec + 1) // 2)})
            ax.set_yticks([-m for m in reversed(marcas)] + [0] + marcas)
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(int(v)):,}"))
            ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.set_ylabel(f"strand\n({bibliotecas[a]})", fontsize=8.5, color=TINTA_2)
            ax.legend(handles=[Patch(color=COR_MAIS, label="+ strand transcripts"),
                               Patch(color=COR_MENOS, label="− strand transcripts")],
                      loc="upper right", fontsize=7.5, frameon=False, ncol=2)
        ax.tick_params(labelbottom=False)

    # Track: EVE lane(s) on top, host genes underneath
    ax_t = eixos[-1]
    sombrear(ax_t)
    total = len(linhas_e) + len(linhas_g)
    passo = 1.0 / max(total, 1)
    y = 1.0
    for linha in linhas_e:
        y -= passo
        for s, e, l in linha:
            ax_t.add_patch(Rectangle((s, y + passo * 0.12), e - s, passo * 0.4, color=COR_EVE, lw=0, zorder=3))
            seta = {"+": " ▸", "-": " ◂"}.get(l.strand, "")
            rotulo = re.sub(rf"^{re.escape(l.contig)}_", "", l.label) + seta
            ax_t.text((s + e) / 2, y + passo * 0.58, rotulo, ha="center", va="bottom", fontsize=7,
                      color="#a32c2b", fontweight="bold", clip_on=True)
    for linha in linhas_g:
        y -= passo
        for s, e, fita, nome in linha:
            s_v, e_v = max(s, window.start), min(e, window.end)
            ax_t.add_patch(Rectangle((s_v, y + passo * 0.12), e_v - s_v, passo * 0.38, color=COR_GENE,
                                     alpha=0.85, lw=0, zorder=2))
            if (e_v - s_v) > window.length * 0.06 and nome:
                ax_t.text((s_v + e_v) / 2, y + passo * 0.56, f"{nome} {fita}", ha="center", va="bottom",
                          fontsize=6.5, color=TINTA_2, clip_on=True)
    if not linhas_g:
        ax_t.text(0.995, 0.05, "no host annotation (GFF) for this window", transform=ax_t.transAxes,
                  ha="right", va="bottom", fontsize=7, color=TINTA_2)
    ax_t.set_ylim(0, 1)
    ax_t.set_yticks([])
    ax_t.set_ylabel("EVEs /\nhost genes", fontsize=8.5, color=TINTA_2)
    ax_t.spines[["top", "right", "left"]].set_visible(False)
    ax_t.tick_params(colors=TINTA_2, labelsize=8)
    ax_t.set_xlim(window.start, window.end)
    ax_t.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax_t.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_t.set_xlabel(f"position on {window.contig} (bp)", fontsize=9, color=TINTA_2)

    extra = ""
    if any(b["in_ref"].min() < 0.5 for b in perfis.values()):
        extra = "  ·  gray = outside the BAM reference"
    eixos[0].set_title(f"RNA-seq coverage · {window.contig}:{window.start + 1:,}-{window.end:,} "
                       f"({window.length / 1000:.1f} kb, {len(window.loci)} EVE locus/loci shaded){extra}",
                       fontsize=10.5, fontweight="bold", color=TINTA, loc="left", pad=8)

    fig.savefig(destino, dpi=170, bbox_inches="tight")
    fig.savefig(destino.rsplit(".", 1)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return destino


# ---------------------------------------------------------------------------
# Legacy per-base table (read by plot_genomic_context through curadoria)
# ---------------------------------------------------------------------------

def write_legacy_depth(results_dir: str, windows: list[Window], perfis_por_amostra: dict, destino: str) -> str | None:
    """Per-base depth summed over samples, per Cluster_N of eve_clusters_flanked.bed (pos 1..N)."""
    bed = os.path.join(results_dir, "06_anotacao_genomica", "eve_clusters_flanked.bed")
    if not os.path.exists(bed) or not perfis_por_amostra:
        return None
    total = {w.name: sum(p[w.name].depth for p in perfis_por_amostra.values()) for w in windows}
    with open(bed) as fh, open(destino, "w") as out:
        for linha in fh:
            c = linha.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            contig, s, e = c[0], int(c[1]), int(c[2])
            for w in windows:
                if w.contig != contig or e <= w.start or s >= w.end:
                    continue
                a, b = max(s, w.start), min(e, w.end)
                vals = total[w.name][a - w.start:b - w.start]
                for i, v in enumerate(vals):
                    out.write(f"{c[3]}\t{a - s + i + 1}\t{int(v)}\n")
                break
    return destino


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_eve_coverage(
    results_dir: str,
    genome: str | None = None,
    gff: str | None = None,
    flank: int = 20000,
    r1: str | None = None,
    r2: str | None = None,
    bams: list | None = None,
    species: str | None = None,
    classified_tsv: str | None = None,
    threads: int = 4,
    library: str = "auto",
    min_mapq: int = 0,
    read_type: str = "sr",
    reference: str = "context",
    outdir: str | None = None,
    plots: bool = True,
    force: bool = False,
) -> dict:
    """
    Coverage of every EVE locus of one results folder (one species) for every sample given.
    Returns {"summary": path, "figures": [...], "legacy_depth_tsv": path} (empty dict when
    there is nothing to do).
    """
    bams = [b for b in (bams or []) if b]
    if not r1 and not bams:
        return {}
    outdir = outdir or os.path.join(results_dir, "10_rnaseq_expressao")
    os.makedirs(outdir, exist_ok=True)
    tag = f"[RNAseq]{f' [{species}]' if species else ''}"

    loci = load_loci(results_dir, species, classified_tsv)
    if not loci:
        print(f"{tag} No EVE locus found in {results_dir}; skipping coverage.")
        return {}
    tamanhos = contig_sizes_for(genome, results_dir)
    janelas = build_windows(loci, flank, tamanhos)
    print(f"{tag} {len(loci)} EVE loci in {len(janelas)} non-redundant window(s) (flank {flank:,} bp).")

    amostras = [Sample(re.sub(r"\.bam$", "", os.path.basename(b), flags=re.I), bam=b) for b in bams]
    if r1:
        if not genome or not os.path.exists(genome):
            print(f"{tag} [warning] reads given but the genome FASTA is missing; cannot align.")
        else:
            amostra = Sample(sample_name(r1), r1=r1, r2=r2)
            if reference == "genome":
                ref = genome
            else:
                ref = write_window_fasta(genome, janelas, os.path.join(outdir, "eve_windows.fasta"))
            try:
                amostra.bam = align_reads(amostra, ref, outdir, threads, read_type, force)
                amostras.append(amostra)
            except Exception as erro:
                print(f"{tag} [warning] {amostra.name}: {erro}")

    contigs = {w.contig for w in janelas}
    host = load_host_features(gff, janelas)
    fitas_exon = {w.name: _exon_strand_array(w, host[w.name]["exons"]) for w in janelas} if gff else {}

    resumo, perfis, bibliotecas, registro = [], {}, {}, []
    for am in amostras:
        if not am.bam or not os.path.exists(am.bam):
            print(f"{tag} [warning] BAM not found: {am.bam}")
            continue
        try:
            caminho, indice = prepare_bam(am.bam, os.path.join(outdir, "bam_index"))
        except Exception as erro:
            print(f"{tag} [warning] {am.name}: could not index {am.bam} ({erro}); "
                  f"it must be coordinate-sorted (samtools sort).")
            continue
        ref = BamReference(caminho, contigs, indice)
        if not ref.segments:
            print(f"{tag} [warning] {am.name}: no BAM reference matches the EVE contigs "
                  f"(BAM: {', '.join(ref.ignored[:3])}...; EVEs: {', '.join(sorted(contigs)[:3])}). Skipped.")
            continue
        votos = [0, 0]
        p_am = {w.name: profile_window(ref, w, min_mapq, fitas_exon.get(w.name), votos) for w in janelas}

        lib = library
        if lib == "auto":
            lib = infer_library(votos)
            detalhe = f"{votos[0]} read-1 on the exon strand, {votos[1]} opposite"
            if sum(votos) < MIN_STRAND_READS:
                detalhe = ("too few reads on annotated exons" if gff else "no GFF to test against") + \
                          "; strand split not reported (set --rnaseq-library fr/rf if known)"
            print(f"{tag} {am.name}: library {lib} ({detalhe}).")
        if lib == "rf":                       # read 1 is antisense to the transcript
            for p in p_am.values():
                p.plus, p.minus = p.minus, p.plus
        bibliotecas[am.name] = lib
        registro.append({"name": am.name, "bam": os.path.abspath(caminho), "index": indice,
                         "library": lib, "source": am.bam})

        cobertas = sum(p.in_ref.any() for p in p_am.values())
        print(f"{tag} {am.name}: {cobertas}/{len(janelas)} window(s) present in the BAM reference.")
        perfis[am.name] = p_am
        for w in janelas:
            resumo += summarise(am.name, w, p_am[w.name], flank, lib)

    if registro:
        # The curation page reads the BAMs back from here to draw coverage over any region
        with open(os.path.join(outdir, "samples.json"), "w", encoding="utf-8") as fh:
            json.dump(registro, fh, indent=1)
    if not resumo:
        print(f"{tag} No coverage computed.")
        return {}

    caminho_resumo = os.path.join(outdir, "eve_coverage_summary.tsv")
    with open(caminho_resumo, "w", newline="") as fh:
        escritor = csv.DictWriter(fh, fieldnames=list(resumo[0].keys()), delimiter="\t")
        escritor.writeheader()
        escritor.writerows(resumo)

    figuras = []
    pasta_fig = os.path.join(outdir, "figures")
    os.makedirs(pasta_fig, exist_ok=True)
    for w in janelas:
        # A window no BAM covers (e.g. a BAM made for one chromosome) would only give an
        # empty gray panel; its loci are already reported as "Not in BAM reference".
        if not any(p_am[w.name].in_ref.any() for p_am in perfis.values()):
            continue
        bins = {}
        for a, p_am in perfis.items():
            bins[a] = bin_profile(w, p_am[w.name])
            pasta_bins = os.path.join(outdir, "coverage_bins", a)
            os.makedirs(pasta_bins, exist_ok=True)
            export_bins(bins[a], os.path.join(pasta_bins, f"{w.name}.tsv"), bibliotecas[a] in ("fr", "rf"))
        if plots:
            try:
                figuras.append(plot_window(w, bins, host[w.name], bibliotecas,
                                           os.path.join(pasta_fig, f"{w.name}_rnaseq_coverage.png")))
            except Exception as erro:
                print(f"{tag} [warning] figure for {w.name} failed: {erro}")

    legado = write_legacy_depth(results_dir, janelas, perfis, os.path.join(outdir, "rnaseq_depth_per_base.tsv"))

    contagem = defaultdict(int)
    for r in resumo:
        contagem[r["status"]] += 1
    print(f"{tag} Coverage summary: {caminho_resumo}")
    for st, n in sorted(contagem.items(), key=lambda x: -x[1]):
        print(f"{tag}   {n:4d}  {st}")
    if figuras:
        print(f"{tag} {len(figuras)} figure(s) in {pasta_fig}")
    return {"summary": caminho_resumo, "figures": figuras, "legacy_depth_tsv": legado}


def main():
    ap = argparse.ArgumentParser(
        description="RNA-seq coverage over each EVE and its flanks, for an existing results folder.")
    ap.add_argument("-r", "--results", required=True, help="Results folder of one species (main.py -o, or <outdir>/especies/<sp>)")
    ap.add_argument("--bam", nargs="+", default=[], help="BAM file(s), aligned to the genome or to context windows")
    ap.add_argument("--r1", help="RNA-seq FASTQ (R1 or single-end); aligned here")
    ap.add_argument("--r2", help="RNA-seq FASTQ R2 (paired-end)")
    ap.add_argument("--genome", help="Host genome FASTA (default: pipeline_info.json)")
    ap.add_argument("--gff", help="Host GFF3 (default: pipeline_info.json); used for gene tracks and library type")
    ap.add_argument("--flank", type=int, help="Flank size in bp (default: pipeline_info.json, else 20000)")
    ap.add_argument("--species", help="Species name, to filter eves_classified.tsv of a multi-species run")
    ap.add_argument("--library", default="auto", choices=["auto", "fr", "rf", "unstranded"],
                    help="Library strandedness; auto infers it from the GFF host exons (default)")
    ap.add_argument("--min-mapq", type=int, default=0, help="Ignore alignments below this MAPQ (default 0)")
    ap.add_argument("--read-type", default="sr", choices=["sr", "ont", "pb", "hifi"])
    ap.add_argument("--reference", default="context", choices=["context", "genome"],
                    help="Align reads to the EVE windows (fast) or the whole genome (no paralogue recruitment)")
    ap.add_argument("-o", "--outdir", help="Output folder (default: <results>/10_rnaseq_expressao)")
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--force", action="store_true", help="Re-align even when the BAM already exists")
    a = ap.parse_args()

    info = {}
    caminho_info = os.path.join(a.results, "pipeline_info.json")
    if os.path.exists(caminho_info):
        with open(caminho_info, encoding="utf-8") as fh:
            info = json.load(fh)
    res = run_eve_coverage(
        results_dir=a.results, genome=a.genome or info.get("genome"), gff=a.gff or info.get("gff"),
        flank=a.flank or info.get("flank") or 20000, r1=a.r1, r2=a.r2, bams=a.bam,
        species=a.species, threads=a.threads, library=a.library, min_mapq=a.min_mapq,
        read_type=a.read_type, reference=a.reference, outdir=a.outdir, plots=not a.no_plots, force=a.force)
    if not res:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
