#!/usr/bin/env python3
"""
Small RNA reads over the EVEs, with the methodology of sRNAProfiler
(https://github.com/LymF/srnaprofiler): its mapping, its 1/NH-weighted features and its rule
set for siRNA and piRNA.

Mapping (sRNAProfiler `mapping`). The reads must come with the adapters already removed.
  * U is turned into T when the file is RNA;
  * optional structural-RNA filter: reads matching a Bowtie index of rRNA/tRNA/snRNA/snoRNA
    with at most one mismatch are dropped before the genome;
  * Bowtie1 against the whole genome, `-v 1 -k 20 --best --strata`: at most ONE mismatch,
    always, and only the best stratum, so a perfect match hides the 1-mismatch placements;
  * unmapped reads dropped, coordinate-sorted and indexed BAM.

Counting (sRNAProfiler `classify`). A read placed at N positions weighs 1/N at each one
(Bowtie1 has no NH tag: its XM holds the number of reported alignments), so every count is a
weighted sum, and the library size used for RPM is the sum of the weights, not the number of
alignment lines.

For a span (an EVE, a cluster, a tblastn segment) the reads with at least half their length
inside it are kept - a read that only touches the edge belongs to the flank - and described
as sRNAProfiler describes a locus: size profile (18-35 nt), strand balance, 5' U and A,
A at position 10, the 10-nt 5'-5' overlap of ping-pong (fraction and z-score), phasing, the
share of uniquely mapped reads, and the class given by its rule set (siRNA 20-22 nt, both
strands, no ping-pong; piRNA 24-31 nt with 5' U or ping-pong), with the softmax confidence
and the `review` flag.

Files, in 13_curadoria/srna/: libraries.json (BAMs, library sizes), bam/<species>/, the
Bowtie index of each genome in index/<species>/ and the mapping logs.

  python -m evetool.annotation.srna -o <results> --species Sp --reads trimmed.fastq.gz
  python -m evetool.annotation.srna -o <results> --species Sp --region contig:1000-5000
"""
import argparse
import gzip
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import threading
from collections import Counter, OrderedDict, defaultdict

import numpy as np

SUBPASTA = os.path.join("13_curadoria", "srna")
MISMATCHES = 1                 # fixed: never more than one mismatch
MAX_ALINHAMENTOS = 20          # bowtie -k
MAX_REGIAO = 2_000_000
MAX_BINS = 1200

# sRNAProfiler rules.yaml (animal). The class windows also colour the track.
REGRAS = {
    "extractor": {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8},
    "decision": {"review_confidence": 0.6, "min_margin": 0.15, "min_reads": 10},
    "classes": OrderedDict([
        ("siRNA", {"size": {"lo": 20, "hi": 22}, "w_size": 1.0, "w_fractop_mid": 1.0,
                   "w_pingpong_penalty": 1.0, "w_phase": 0.5, "defining_evidence": ["fractop_bid", "size"]}),
        ("piRNA", {"size": {"lo": 24, "hi": 31}, "w_size": 1.0, "w_pingpong": 1.5, "w_1u": 1.0,
                   "w_phase": 0.5, "defining_evidence": [["1u", "pingpong"], "size"]}),
    ]),
}
CLASSES_TRILHA = ("siRNA", "piRNA", "other")
PHASE_WINDOW_CYCLES = 11
NORM_LIB_SIZE = 20_000_000

_trava = threading.Lock()
_cache = OrderedDict()
_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------

def _pasta(outdir: str) -> str:
    return os.path.join(outdir, SUBPASTA)


def _arq_registro(outdir: str) -> str:
    return os.path.join(_pasta(outdir), "libraries.json")


def _ler_registro(outdir: str) -> dict:
    try:
        with open(_arq_registro(outdir), encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def _gravar_registro(outdir: str, dados: dict) -> None:
    os.makedirs(_pasta(outdir), exist_ok=True)
    with open(_arq_registro(outdir), "w", encoding="utf-8") as fh:
        json.dump(dados, fh, indent=1)
    _cache.clear()


def bibliotecas(outdir: str) -> dict:
    """{species: [{name, bam, lib_size, ...}]} of the small RNA libraries whose BAM still exists."""
    return {sp: [b for b in l if os.path.exists(b.get("bam", ""))]
            for sp, l in _ler_registro(outdir).items() if l}


def resumo(outdir: str) -> dict:
    """What the page needs up front: the libraries of every species."""
    campos = ("name", "lib_size", "input_reads", "aligned_reads", "aligner", "source", "ncrna_removed")
    return {sp: [{k: b.get(k) for k in campos} for b in l] for sp, l in bibliotecas(outdir).items() if l}


def _usa_xm(bam: str) -> bool:
    """True when the @PG header shows Bowtie1, whose XM is the number of reported alignments."""
    import pysam
    with pysam.AlignmentFile(bam) as fh:
        for pg in fh.header.get("PG", []):
            prog = " ".join(str(pg.get(k, "")) for k in ("PN", "ID", "CL")).lower()
            if "bowtie2" in prog:
                return False
            if "bowtie" in prog:
                return True
    return False


def _nh(rec, xm: bool) -> int:
    if rec.has_tag("NH"):
        return max(int(rec.get_tag("NH")), 1)
    if xm and rec.has_tag("XM"):
        return max(int(rec.get_tag("XM")), 1)
    return 1


def tamanho_ponderado(bam: str, xm: bool) -> float:
    """Library size on the 1/NH scale (one full pass): the number of mapped reads, not of lines."""
    import pysam
    total = 0.0
    with pysam.AlignmentFile(bam) as fh:
        for rec in fh.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            total += 1.0 / _nh(rec, xm)
    return round(total, 2)


def _registrar(outdir: str, species: str, item: dict) -> dict:
    dados = _ler_registro(outdir)
    dados[species] = [b for b in dados.get(species, []) if b["name"] != item["name"]] + [item]
    _gravar_registro(outdir, dados)
    return item


def adicionar_bam(outdir: str, species: str, bam: str, nome: str = "", passo=None) -> dict:
    """Register a small RNA BAM made elsewhere (e.g. by sRNAProfiler `mapping`); it is not copied."""
    from evetool.annotation.rnaseq_coverage import is_bam, prepare_bam
    from evetool.curation.curation import pastas_especies
    passo = passo or (lambda t: None)
    bam = os.path.abspath(os.path.expanduser(str(bam or "").strip()))
    if not os.path.exists(bam):
        raise RuntimeError(f"file not found: {bam}")
    if not is_bam(bam):
        raise RuntimeError(f"not a BAM file: {bam}")
    if species not in pastas_especies(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    passo("indexing the BAM")
    caminho, indice = prepare_bam(bam, os.path.join(_pasta(outdir), "bam_index"))
    xm = _usa_xm(caminho)
    passo("measuring the library size (sum of 1/NH)")
    nome = nome or re.sub(r"\.bam$", "", os.path.basename(bam), flags=re.I)
    return _registrar(outdir, species, {
        "name": nome, "bam": os.path.abspath(caminho), "index": indice, "xm_as_nh": xm,
        "lib_size": tamanho_ponderado(caminho, xm), "source": bam,
        "aligner": "bowtie1 (XM = alignments)" if xm else "external BAM (NH tag, else unique)"})


def remover(outdir: str, species: str, nome: str) -> dict:
    """Forget a library (its BAM stays on disk)."""
    dados = _ler_registro(outdir)
    dados[species] = [b for b in dados.get(species, []) if b["name"] != nome]
    _gravar_registro(outdir, dados)
    return {"libraries": resumo(outdir)}


# ---------------------------------------------------------------------------
# Mapping (sRNAProfiler `mapping`)
# ---------------------------------------------------------------------------

def _abrir(caminho: str):
    return gzip.open(caminho, "rt") if str(caminho).endswith(".gz") else open(caminho, "rt")


def _formato(caminho: str) -> str:
    with _abrir(caminho) as fh:
        primeiro = fh.readline()
    if primeiro.startswith("@"):
        return "fastq"
    if primeiro.startswith(">"):
        return "fasta"
    raise RuntimeError(f"{os.path.basename(caminho)} is neither FASTQ nor FASTA")


def _tem_uracila(caminho: str, fmt: str, n: int = 2000) -> bool:
    passo = 4 if fmt == "fastq" else 2
    with _abrir(caminho) as fh:
        for k, linha in enumerate(fh):
            if k >= n * passo:
                break
            if k % passo == 1 and ("U" in linha or "u" in linha):
                return True
    return False


def _converter_u(origem: str, destino: str, fmt: str) -> str:
    tr = str.maketrans("Uu", "Tt")
    with _abrir(origem) as fi, gzip.open(destino, "wt", compresslevel=3) as fo:
        if fmt == "fastq":
            for k, linha in enumerate(fi):
                fo.write(linha.translate(tr) if k % 4 == 1 else linha)
        else:
            for linha in fi:
                fo.write(linha if linha.startswith(">") else linha.translate(tr))
    return destino


def _indice_ok(prefixo: str) -> bool:
    return all(os.path.exists(prefixo + e) or os.path.exists(prefixo + e + "l")
               for e in (".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"))


def _genoma(outdir: str, species: str) -> str:
    """The genome the reads are mapped on: the whole genome when the run was made on segments."""
    from evetool.curation.curation import genomas_para_sintenia
    g = genomas_para_sintenia(outdir).get(species) or {}
    if not g.get("genome"):
        raise RuntimeError(f"genome of '{species or 'this run'}' not found")
    return g["genome"]


def _estatisticas_bowtie(log: str) -> dict:
    est = {}
    try:
        with open(log, errors="replace") as fh:
            for linha in fh:
                if "reads processed:" in linha:
                    est["input_reads"] = int(linha.split(":")[1].strip())
                elif "reads with at least one" in linha and "alignment:" in linha:
                    est["aligned_reads"] = int(linha.split(":")[1].split()[0])
    except (OSError, ValueError):
        pass
    return est


def _memoria_sort(threads: int) -> int:
    """MB per thread for samtools sort: half of the free memory, at most 4 GB, split by thread."""
    try:
        livre = (os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")) // (1024 ** 2)
    except (ValueError, OSError, AttributeError):
        livre = 0
    total = min(4096, max(768, livre // 2)) if livre else 768
    return max(64, total // max(1, threads))


def mapear(outdir: str, species: str, reads: str, nome: str = "", threads: int = 4,
           ncrna_index: str = "", passo=None, forcar: bool = False) -> dict:
    """
    Map adapter-trimmed small RNA reads (FASTQ/FASTA, .gz or not) on the genome of `species`
    with Bowtie1 (-v 1 -k 20 --best --strata) and register the sorted BAM as a library.
    """
    from evetool.curation.curation import pastas_especies
    passo = passo or (lambda t: print(f"[sRNA] {t}", flush=True))
    reads = os.path.abspath(os.path.expanduser(str(reads or "").strip()))
    if not os.path.exists(reads):
        raise RuntimeError(f"file not found: {reads}")
    if species not in pastas_especies(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    for ferramenta in ("bowtie", "bowtie-build", "samtools"):
        if not shutil.which(ferramenta):
            raise RuntimeError(f"{ferramenta} not found on the PATH (conda install -c bioconda bowtie samtools)")
    genoma = _genoma(outdir, species)
    fmt = _formato(reads)
    nome = nome or re.sub(r"\.(fastq|fq|fasta|fa|fna)(\.gz)?$", "", os.path.basename(reads), flags=re.I)
    nome = re.sub(r"[^\w.\-]+", "_", nome)
    sp_dir = species or "genome"
    pasta_bam = os.path.join(_pasta(outdir), "bam", sp_dir)
    pasta_idx = os.path.join(_pasta(outdir), "index", sp_dir)
    os.makedirs(pasta_bam, exist_ok=True)
    os.makedirs(pasta_idx, exist_ok=True)
    bam = os.path.join(pasta_bam, nome + ".bam")
    log = os.path.join(pasta_bam, nome + ".bowtie.log")

    prefixo = os.path.join(pasta_idx, "genome")
    if not _indice_ok(prefixo):
        passo(f"building the Bowtie index of {os.path.basename(genoma)} (once per genome)")
        r = subprocess.run(["bowtie-build", "--threads", str(threads), genoma, prefixo],
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0 or not _indice_ok(prefixo):
            raise RuntimeError(f"bowtie-build failed: {r.stderr.strip()[-300:]}")

    temporarios, est = [], {}
    try:
        entrada = reads
        if _tem_uracila(reads, fmt):
            passo("turning U into T")
            entrada = _converter_u(reads, os.path.join(pasta_bam, f"{nome}.DNA.{fmt}.gz"), fmt)
            temporarios.append(entrada)
        flag = "-f" if fmt == "fasta" else "-q"
        if ncrna_index:
            ncrna_index = os.path.expanduser(ncrna_index.strip())
            if not _indice_ok(ncrna_index):
                raise RuntimeError(f"no Bowtie index with the prefix {ncrna_index}")
            passo("dropping structural RNA (rRNA, tRNA, snRNA, snoRNA)")
            filtrado = os.path.join(pasta_bam, f"{nome}.ncrna_filtered.{fmt}")
            log_nc = os.path.join(pasta_bam, nome + ".ncrna.log")
            with open(log_nc, "w") as err:
                r = subprocess.run(["bowtie", "-v", str(MISMATCHES), flag, "-p", str(threads), "-x", ncrna_index,
                                    entrada, "--un", filtrado, "-S", "/dev/null"], stderr=err)
            if r.returncode != 0:
                raise RuntimeError(f"structural RNA filter failed; see {log_nc}")
            nc = _estatisticas_bowtie(log_nc)
            est["ncrna_removed"] = nc.get("aligned_reads")
            est["ncrna_index"] = ncrna_index
            temporarios.append(filtrado)
            entrada = filtrado

        passo(f"Bowtie1 -v {MISMATCHES} -k {MAX_ALINHAMENTOS} --best --strata on the genome")
        cmd = ["bowtie", flag, "-S", "-p", str(threads), "-k", str(MAX_ALINHAMENTOS), "-v", str(MISMATCHES),
               "--best", "--strata", "-x", prefixo, entrada]
        tmp = bam + ".tmp"
        with open(log, "w") as err:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
            p2 = subprocess.Popen(["samtools", "view", "-h", "-F", "4", "-"], stdin=p1.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            p1.stdout.close()
            p3 = subprocess.run(["samtools", "sort", "-@", str(max(1, threads // 2)),
                                 "-m", f"{_memoria_sort(max(1, threads // 2))}M", "-T", bam + ".sort", "-o", tmp, "-"],
                                stdin=p2.stdout, stderr=subprocess.PIPE, text=True)
            p2.stdout.close()
            rc1, rc2 = p1.wait(), p2.wait()
        if rc1 != 0 or rc2 != 0 or p3.returncode != 0:
            raise RuntimeError(f"Bowtie mapping failed; see {log}{(' · ' + p3.stderr.strip()[-200:]) if p3.stderr else ''}")
        os.replace(tmp, bam)
        subprocess.run(["samtools", "index", bam], check=True)
        est.update(_estatisticas_bowtie(log))
        passo("measuring the library size (sum of 1/NH)")
        tamanho = tamanho_ponderado(bam, True)
    finally:
        for t in temporarios + [bam + ".tmp"]:
            if os.path.exists(t):
                os.remove(t)
    return _registrar(outdir, species, {
        "name": nome, "bam": bam, "index": bam + ".bai", "xm_as_nh": True, "lib_size": tamanho,
        "source": reads, "genome": genoma, "mismatches": MISMATCHES, "k": MAX_ALINHAMENTOS,
        "aligner": f"bowtie1 -v {MISMATCHES} -k {MAX_ALINHAMENTOS} --best --strata", **est})


# ---------------------------------------------------------------------------
# Signals and rule set (sRNAProfiler signals.py, evidence.py, decider.py)
# ---------------------------------------------------------------------------

def _espectro(a5: dict, b5: dict, lo: int = 1, hi: int = 20) -> Counter:
    """5'-5' offsets between strands: counts[d] = sum over p of a5[p] * b5[p + d - 1]."""
    c = Counter()
    for p, ac in a5.items():
        for d in range(lo, hi + 1):
            bc = b5.get(p + d - 1, 0)
            if bc:
                c[d] += ac * bc
    return c


def _pingpong_z(counts, alvo: int = 10, lo: int = 1, hi: int = 20, min_bg: int = 3, teto: float = 100.0):
    fundo = [counts.get(o, 0) for o in range(lo, hi + 1) if o != alvo]
    pico = counts.get(alvo, 0)
    if sum(1 for v in fundo if v > 0) < min_bg:
        return 0.0, pico
    sd = statistics.pstdev(fundo)
    if sd == 0:
        return 0.0, pico
    return round(max(-teto, min(teto, (pico - statistics.fmean(fundo)) / sd)), 3), pico


def _fase(posicoes, inicio: int, pl: int, lib_size: float):
    """PhasiHunter-style score at the dominant register: ratio * cycles * ln(mass per 20 M reads)."""
    abun, ciclos, total = defaultdict(float), defaultdict(set), 0.0
    for p5, w in posicoes:
        r = (p5 - inicio) % pl
        abun[r] += w
        ciclos[r].add((p5 - inicio) // pl)
        total += w
    if not total:
        return 0.0, 0, 0.0
    dom = max(abun, key=abun.get)
    massa = abun[dom] * (NORM_LIB_SIZE / lib_size) if lib_size else abun[dom]
    numero = len(ciclos[dom])
    razao = abun[dom] / total
    return (round(razao * numero * math.log(massa), 3) if massa > 0 else 0.0), numero, round(razao, 3)


def _fase_melhor_janela(posicoes, inicio: int, pl: int, lib_size: float):
    """Best window of 11 cycles, chosen by cycles in register and then score (never by p-value)."""
    if not posicoes:
        return 0.0, 0, 0.0
    coords = sorted(p for p, _ in posicoes)
    janela = PHASE_WINDOW_CYCLES * pl
    if coords[-1] - coords[0] + 1 <= janela:
        return _fase(posicoes, inicio, pl, lib_size)
    ordenadas = sorted(posicoes)
    melhor = (0.0, 0, 0.0)
    import bisect
    for ini in range(coords[0], coords[-1] + 1, pl):
        i, j = bisect.bisect_left(coords, ini), bisect.bisect_left(coords, ini + janela)
        if i == j:
            continue
        r = _fase(ordenadas[i:j], inicio, pl, lib_size)
        if (r[1], r[0]) > (melhor[1], melhor[0]):
            melhor = r
    return melhor


def _forca(k: float, n: float, p0: float, conf: float = 0.95) -> float:
    """Wilson lower bound of k/n rescaled by the room the null leaves, in [0, 1]."""
    n, k = int(round(n)), int(round(k))
    if n <= 0 or k <= 0 or k > n or not (0.0 <= p0 < 1.0):
        return 0.0
    z = 1.959963984540054 if conf == 0.95 else statistics.NormalDist().inv_cdf(1 - (1 - conf) / 2)
    ph = k / n
    den = 1.0 + z * z / n
    centro = (ph + z * z / (2 * n)) / den
    meia = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return round(max(0.0, (centro - meia - p0) / (1.0 - p0)), 4)


def _nulo_base(seq: str, base: str, frac_top: float) -> float:
    """Chance that a cut anywhere in the locus DNA gives a read starting with `base`."""
    par = {"A": "T", "T": "A", "C": "G", "G": "C"}
    s = (seq or "").upper()
    n = sum(s.count(b) for b in "ACGT")
    if not n:
        return 0.25
    return frac_top * s.count(base) / n + (1.0 - frac_top) * s.count(par[base]) / n


def _janela(tam, lo, hi, borda=3.0):
    if tam is None:
        return 0.0
    if lo <= tam <= hi:
        return 1.0
    return max(0.0, 1.0 - (lo - tam if tam < lo else tam - hi) / borda)


def _decidir(f: dict) -> None:
    """sRNAProfiler RuleScorer: score every class, softmax confidence, review flag."""
    tam = int(f["dicer_call"]) if f["dicer_call"] != "N" else None
    pp = min(f["pingpong_frac"] / 0.1, 1.0)
    notas, termos = {}, {}
    for nome, spec in REGRAS["classes"].items():
        porta = _janela(tam, spec["size"]["lo"], spec["size"]["hi"])
        t = {}
        if "w_size" in spec:
            t["size"] = spec["w_size"] * porta
        if "w_pingpong" in spec:
            t["pingpong"] = spec["w_pingpong"] * pp * porta
        if "w_1u" in spec:
            t["1u"] = spec["w_1u"] * f["pct_1u"] * porta
        if "w_fractop_mid" in spec:
            t["fractop_mid"] = spec["w_fractop_mid"] * math.exp(-((f["frac_top"] - 0.5) ** 2) / (2 * 0.15 ** 2)) * porta
        if "w_pingpong_penalty" in spec:
            t["pingpong_penalty"] = -spec["w_pingpong_penalty"] * pp * porta
        if "w_phase" in spec:
            t["phase"] = spec["w_phase"] / (1.0 + math.exp(-(f["phase_score"] - 5.0) / 3.0)) * porta
        notas[nome], termos[nome] = sum(t.values()), {k: round(v, 3) for k, v in t.items()}
    ordem = sorted(notas.items(), key=lambda kv: kv[1], reverse=True)
    topo, nota = ordem[0]
    exps = {k: math.exp(v) for k, v in notas.items()}
    conf = exps[topo] / sum(exps.values())
    margem = nota - (ordem[1][1] if len(ordem) > 1 else 0.0)
    spec = REGRAS["classes"][topo]

    def forca(nome_ev):
        if nome_ev == "size":
            total = sum(f["sizes"].values())
            if total <= 0 or f["distinct_seqs"] <= 0:
                return 0.0
            ext = REGRAS["extractor"]
            dentro = sum(c for L, c in f["sizes"].items() if spec["size"]["lo"] <= int(L) <= spec["size"]["hi"])
            p0 = (spec["size"]["hi"] - spec["size"]["lo"] + 1) / (ext["dicermax"] - ext["dicermin"] + 1)
            return _forca(round(dentro / total * f["distinct_seqs"]), f["distinct_seqs"], p0)
        return f.get("e_" + nome_ev, 0.0)

    definidoras = [max(forca(n) for n in e) if isinstance(e, list) else forca(e) for e in spec["defining_evidence"]]
    dec = REGRAS["decision"]
    f["class_scores"] = {k: round(v, 3) for k, v in notas.items()}
    f["cls"] = topo
    f["score"] = round(min(definidoras), 3) if definidoras and tam is not None else 0.0
    f["confidence"] = round(conf, 3)
    f["evidence"] = "; ".join(f"{k}={v}" for k, v in termos[topo].items())
    f["review"] = bool(conf < dec["review_confidence"] or margem < dec["min_margin"]
                       or f["reads"] < dec["min_reads"] or tam is None)
    if tam is None or not (spec["size"]["lo"] <= tam <= spec["size"]["hi"]):
        f["cls"], f["score"], f["review"] = "ambiguous", None, True
        f["evidence"] = f"nearest={topo}; " + f["evidence"]


def caracterizar(alns: dict, inicio: int, fim: int, lib_size: float, seq_locus: str = "",
                 fita_eve: str = ".") -> dict:
    """
    Features of one span from the alignments assigned to it (sRNAProfiler build_record).
    `alns` = {(pos, strand, sequence as aligned, NH): number of identical alignment lines}.
    """
    ext = REGRAS["extractor"]
    dmin, dmax = ext["dicermin"], ext["dicermax"]
    f = {"start": inicio, "end": fim, "reads": 0.0, "alignments": sum(alns.values())}
    if not alns:
        return f
    fitas, seqs, tam = Counter(), Counter(), Counter()
    por_fita = {"+": Counter(), "-": Counter()}
    n1u = n5a = n10a = unicos = gc = bases = 0.0
    multi = 0
    mais5, menos5, fase_pos = Counter(), Counter(), []
    # piRNA-sized reads only, per strand: primary piRNAs carry a 5' U, secondary (ping-pong)
    # piRNAs an A at position 10 on the opposite strand, 10 nt overlapping their primary partner
    plo, phi = REGRAS["classes"]["piRNA"]["size"]["lo"], REGRAS["classes"]["piRNA"]["size"]["hi"]
    pi = {"+": Counter(), "-": Counter()}
    pi5 = {"+": Counter(), "-": Counter()}
    for (pos, fita, seq, nh), n in alns.items():
        w = n / nh
        orientada = seq if fita == "+" else _revcomp(seq)
        L = len(orientada)
        fitas[fita] += w
        seqs[orientada] += w
        chave = "short" if L < dmin else "long" if L > dmax else L
        tam[chave] += w
        if isinstance(chave, int):
            por_fita[fita][L] += w
        gc += w * sum(1 for b in orientada if b in "GCgc")
        bases += w * L
        if orientada[:1] in ("T", "U"):
            n1u += w
        elif orientada[:1] == "A":
            n5a += w
        if L >= 10 and orientada[9] == "A":
            n10a += w
        if nh > 1:
            multi += n
        else:
            unicos += w
        cinco = pos if fita == "+" else pos + L - 1
        (mais5 if fita == "+" else menos5)[cinco] += w
        fase_pos.append((cinco if fita == "+" else cinco - 2, w, L))
        if plo <= L <= phi:
            p = pi[fita]
            p["n"] += w
            if orientada[:1] in ("T", "U"):
                p["u1"] += w
            if orientada[9] == "A":
                p["a10"] += w
            pi5[fita][cinco] += w

    leituras = fitas["+"] + fitas["-"]
    f["reads"] = round(leituras, 3)
    f["reads_rpm"] = round(leituras / lib_size * 1e6, 3) if lib_size else 0.0
    f["distinct_seqs"] = n = len(seqs)
    f["major_rna"], maior = seqs.most_common(1)[0]
    f["major_rna_reads"] = round(maior, 3)
    f["frac_top"] = round(fitas["+"] / leituras, 4)
    f["strand_call"] = "+" if f["frac_top"] >= ext["strand_cutoff"] else "-" if f["frac_top"] <= 1 - ext["strand_cutoff"] else "."
    if fita_eve in "+-":
        f["sense_frac"] = round(fitas[fita_eve] / leituras, 4)       # same strand as the EVE's ORF
    f["short"], f["long"] = round(tam["short"], 3), round(tam["long"], 3)
    f["sizes"] = {L: round(tam.get(L, 0.0), 3) for L in range(dmin, dmax + 1)}
    f["sizes_plus"] = {L: round(por_fita["+"].get(L, 0.0), 3) for L in range(dmin, dmax + 1)}
    f["sizes_minus"] = {L: round(por_fita["-"].get(L, 0.0), 3) for L in range(dmin, dmax + 1)}
    f["frac_in_window"] = round((leituras - tam["short"] - tam["long"]) / leituras, 4)
    predominante = max(range(dmin, dmax + 1), key=lambda L: tam.get(L, 0.0))
    f["dicer_call"] = str(predominante) if f["frac_in_window"] >= 0.8 else "N"
    total_janela = sum(f["sizes"].values())
    if total_janela:
        banda = lambda c: sum(tam.get(L, 0.0) for L in (c - 1, c, c + 1))
        f["size_purity"] = round(banda(predominante) / total_janela, 4)
        fora = {L: c for L, c in f["sizes"].items() if abs(L - predominante) >= 2 and c}
        if fora:
            segundo = max(fora, key=fora.get)
            if banda(segundo) / total_janela >= 0.2:
                f["dicer_call_2"] = str(segundo)
    f["class_mass"] = {nome: round(sum(c for L, c in f["sizes"].items() if s["size"]["lo"] <= L <= s["size"]["hi"]) / leituras, 4)
                       for nome, s in REGRAS["classes"].items()}
    f["pct_1u"], f["pct_5a"], f["pct_10a"] = round(n1u / leituras, 4), round(n5a / leituras, 4), round(n10a / leituras, 4)
    nulo_u, nulo_a = _nulo_base(seq_locus, "T", f["frac_top"]), _nulo_base(seq_locus, "A", f["frac_top"])
    f["null_1u"] = round(nulo_u, 4)
    f["e_1u"] = _forca(f["pct_1u"] * n, n, nulo_u)
    f["e_5a"] = _forca(f["pct_5a"] * n, n, nulo_a)
    f["e_fractop_bid"] = _forca(min(f["frac_top"], 1 - f["frac_top"]) * n, n, 0.10)
    f["gc"] = round(gc / bases, 4) if bases else 0.0
    f["multimap_frac"] = round(multi / f["alignments"], 4)
    f["unique_frac"] = round(unicos / leituras, 4)

    largo = _espectro(mais5, menos5, 1, dmax)
    sobre = Counter({d: c for d, c in largo.items() if d <= 20})
    f["pingpong_z"], pico = _pingpong_z(sobre)
    f["pingpong_pairs10"] = round(pico, 3)
    pares = sum(sobre.values())
    f["pingpong_pairs"] = round(pares, 3)
    f["pingpong_frac"] = round(sobre.get(10, 0) / pares, 4) if pares >= 10 else 0.0
    f["overlaps"] = [round(sobre.get(d, 0.0), 3) for d in range(1, 21)]
    if pares > 0:
        fundo = [sobre.get(o, 0.0) for o in range(1, 21) if o != 10]
        f["e_pingpong"] = _forca(sobre.get(10, 0), pares, (sum(fundo) / len(fundo)) / pares)
    if f["dicer_call"] != "N" and sum(largo.values()):
        f["overlap_dup_frac"] = round(largo.get(int(f["dicer_call"]) - 2, 0) / sum(largo.values()), 4)
    # the same signals restricted to piRNA-sized reads, kept per strand (the page calls
    # primary / secondary / ping-pong from them, with thresholds the curator can change)
    sobre_pi = Counter({d: c for d, c in _espectro(pi5["+"], pi5["-"], 1, dmax).items() if d <= 20})
    z_pi, pico_pi = _pingpong_z(sobre_pi)
    fr = lambda p, k: round(p[k] / p["n"], 4) if p["n"] else None
    f["pirna"] = {"size": [plo, phi],
                  "+": {"reads": round(pi["+"]["n"], 3), "pct_1u": fr(pi["+"], "u1"), "pct_10a": fr(pi["+"], "a10")},
                  "-": {"reads": round(pi["-"]["n"], 3), "pct_1u": fr(pi["-"], "u1"), "pct_10a": fr(pi["-"], "a10")},
                  "pingpong_z": z_pi, "pingpong_pairs10": round(pico_pi, 3),
                  "pingpong_pairs": round(sum(sobre_pi.values()), 3),
                  "overlaps": [round(sobre_pi.get(d, 0.0), 3) for d in range(1, 21)]}
    f["phase_score"], f["phase_number"], f["phase_ratio"] = 0.0, 0, 0.0
    if f["dicer_call"] != "N":
        pl = int(f["dicer_call"])
        f["phase_score"], f["phase_number"], f["phase_ratio"] = _fase_melhor_janela(
            [(p, w) for p, w, L in fase_pos if L == pl], inicio, pl, lib_size)
    _decidir(f)
    return f


# ---------------------------------------------------------------------------
# Profile of a region, for the page
# ---------------------------------------------------------------------------

def _classe_do_tamanho(L: int) -> int:
    for k, nome in enumerate(CLASSES_TRILHA[:-1]):
        s = REGRAS["classes"][nome]["size"]
        if s["lo"] <= L <= s["hi"]:
            return k
    return len(CLASSES_TRILHA) - 1


def perfil(outdir: str, species: str, contig: str, start: int, end: int, bins: int = 600, foco=None) -> dict:
    """
    Small RNA reads of every library of `species` over contig:[start, end) (0-based, half-open):
    1/NH-weighted depth per strand and size class, binned for the track, and the sRNAProfiler
    features of each span of `foco` = [[start, end, strand], ...] (the EVEs) and of the region.
    """
    import pysam
    from evetool.annotation.rnaseq_tracks import _origem
    from evetool.comparative.flank_synteny import mapear_para_genoma
    start, end = int(start), int(end)
    if end <= start:
        raise RuntimeError("empty region")
    if end - start > MAX_REGIAO:
        raise RuntimeError(f"region too large ({end - start:,} nt; at most {MAX_REGIAO:,})")
    bins = max(10, min(int(bins or 600), MAX_BINS, end - start))
    foco = [(int(f[0]), int(f[1]), (f[2] if len(f) > 2 and f[2] in ("+", "-") else ".")) for f in (foco or [])
            if int(f[1]) > int(f[0])]
    libs = bibliotecas(outdir).get(species, [])
    chave = (outdir, species, contig, start, end, bins, tuple(foco), tuple((b["name"], b["lib_size"]) for b in libs))
    if chave in _cache:
        return _cache[chave]
    resposta = {"species": species, "contig": contig, "start": start, "end": end, "bins": bins,
                "libraries": [], "placed_on": None, "classes": list(CLASSES_TRILHA),
                "windows": {n: [s["size"]["lo"], s["size"]["hi"]] for n, s in REGRAS["classes"].items()},
                "size_range": [REGRAS["extractor"]["dicermin"], REGRAS["extractor"]["dicermax"]]}
    if not libs:
        resposta["note"] = "no small RNA library for this species"
        return resposta

    origens = _origem(outdir, species)
    g_contig, g0, g1, invertido = mapear_para_genoma(contig, start, end, origens)
    if g_contig != contig:
        resposta["placed_on"] = {"contig": g_contig, "start": g0, "end": g1, "flipped": invertido}
    troca = {"+": "-", "-": "+", ".": "."}
    # spans in genome coordinates; "region" is the whole window
    alvos = [("region", g0, g1, ".")]
    for k, (a, b, fita) in enumerate(foco):
        _c, x, y, inv = mapear_para_genoma(contig, a, b, origens)
        alvos.append((f"f{k}", x, y, troca[fita] if inv else fita))
    seq_locus = {}
    try:
        with pysam.FastaFile(_genoma(outdir, species)) as fa:
            for nome, x, y, _f in alvos:
                seq_locus[nome] = fa.fetch(g_contig, x, y)
    except Exception:
        pass

    n = g1 - g0
    cortes = np.array_split(np.arange(n), bins)
    ini = np.array([c[0] for c in cortes])
    tam_bin = np.diff(np.append(ini, n))
    nc = len(CLASSES_TRILHA)
    for lib in libs:
        item = {"name": lib["name"], "lib_size": lib.get("lib_size"), "aligner": lib.get("aligner")}
        try:
            xm = bool(lib.get("xm_as_nh"))
            dif = np.zeros((2, nc, n + 1))
            baldes = {nome: Counter() for nome, *_ in alvos}
            with pysam.AlignmentFile(lib["bam"], index_filename=lib.get("index")) as fh:
                if g_contig not in fh.references:
                    raise RuntimeError(f"{g_contig} is not a reference of this BAM")
                for rec in fh.fetch(g_contig, g0, g1):
                    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.query_sequence:
                        continue
                    pos, seq = rec.reference_start, rec.query_sequence
                    L = len(seq)
                    nh = _nh(rec, xm)
                    fita = "-" if rec.is_reverse else "+"
                    a, b = max(pos, g0) - g0, min(pos + L, g1) - g0
                    if a < b:
                        d = dif[1 if rec.is_reverse else 0, _classe_do_tamanho(L)]
                        d[a] += 1.0 / nh
                        d[b] -= 1.0 / nh
                    for nome, x, y, _f in alvos:
                        # at least half of the read inside the span: an edge read belongs to the flank
                        if 2 * (min(pos + L, y) - max(pos, x)) >= L:
                            baldes[nome][(pos, fita, seq, nh)] += 1
            prof = np.cumsum(dif, axis=2)[:, :, :n]
            if invertido:
                prof = prof[::-1, :, ::-1]
            tamanho = lib.get("lib_size") or 0
            item["plus"] = [[round(float(v), 3) for v in np.add.reduceat(prof[0, c], ini) / tam_bin] for c in range(nc)]
            item["minus"] = [[round(float(v), 3) for v in np.add.reduceat(prof[1, c], ini) / tam_bin] for c in range(nc)]
            item["max_depth"] = round(float(max(prof[0].sum(axis=0).max(), prof[1].sum(axis=0).max())), 3)
            feitos = {}
            for nome, x, y, fita in alvos:
                f = caracterizar(baldes[nome], x, y, tamanho, seq_locus.get(nome, ""), fita)
                if invertido and f.get("reads"):          # back to the strands of the segment asked for
                    f["frac_top"] = round(1 - f["frac_top"], 4)
                    f["strand_call"] = troca[f["strand_call"]]
                    f["sizes_plus"], f["sizes_minus"] = f["sizes_minus"], f["sizes_plus"]
                    if f.get("pirna"):
                        f["pirna"]["+"], f["pirna"]["-"] = f["pirna"]["-"], f["pirna"]["+"]
                feitos[nome] = f
            item["region"] = {**feitos["region"], "start": start, "end": end}
            item["focus"] = [{**feitos[f"f{k}"], "start": a, "end": b, "strand": fita} for k, (a, b, fita) in enumerate(foco)]
        except Exception as erro:
            item["error"] = str(erro)
        resposta["libraries"].append(item)
    resposta["max_depth"] = max([l.get("max_depth", 0) for l in resposta["libraries"]] or [0])
    with _trava:
        _cache[chave] = resposta
        while len(_cache) > 64:
            _cache.popitem(last=False)
    return resposta


_cache_matriz = OrderedDict()   # (outdir, species, data stamp, libraries) -> answer (last 8)
CAMPOS_MATRIZ = ("reads", "reads_rpm", "distinct_seqs", "cls", "confidence", "review", "dicer_call", "frac_top",
                 "sense_frac", "pct_1u", "pct_10a", "pingpong_z", "pingpong_frac", "pingpong_pairs10", "unique_frac",
                 "multimap_frac", "class_mass", "sizes_plus", "sizes_minus", "overlaps", "pirna", "major_rna",
                 "major_rna_reads", "evidence")


def matriz(outdir: str, species: str) -> dict:
    """
    Small RNA features of every approved EVE locus of `species` in every library, for the
    heatmap and the piRNA table of the page (the same features as the inspection, plus the
    piRNA-sized reads per strand). Cached until the curation, the data or the libraries change.
    """
    from evetool.curation.curation import estado_dados, loci_aprovados
    libs = bibliotecas(outdir).get(species, [])
    marca = tuple((b["name"], b["bam"], os.path.getmtime(b["bam"]) if os.path.exists(b["bam"]) else 0) for b in libs)
    chave = (outdir, species, estado_dados(outdir), marca)
    if chave in _cache_matriz:
        return _cache_matriz[chave]
    resposta = {"species": species, "loci": [],
                "libraries": [{"name": b["name"], "lib_size": b.get("lib_size")} for b in libs],
                "windows": {n: [s["size"]["lo"], s["size"]["hi"]] for n, s in REGRAS["classes"].items()},
                "size_range": [REGRAS["extractor"]["dicermin"], REGRAS["extractor"]["dicermax"]]}
    if not libs:
        resposta["note"] = "no small RNA library for this species"
        return resposta
    for l in loci_aprovados(outdir, species):
        celulas = {}
        try:
            r = perfil(outdir, species, l["contig"], l["start"], l["end"], 10, [[l["start"], l["end"], l["strand"]]])
            for lib in r["libraries"]:
                if lib.get("error"):
                    celulas[lib["name"]] = {"error": lib["error"]}
                    continue
                f = (lib.get("focus") or [{}])[0]
                celulas[lib["name"]] = {k: f[k] for k in CAMPOS_MATRIZ if k in f}
        except Exception as erro:
            l = {**l, "error": str(erro)}
        resposta["loci"].append({**l, "cells": celulas})
    with _trava:
        _cache_matriz[chave] = resposta
        while len(_cache_matriz) > 8:
            _cache_matriz.popitem(last=False)
    return resposta


COLUNAS_TSV = ["species", "eve", "contig", "start", "end", "strand", "family", "protein_class", "library",
               "reads", "reads_rpm", "distinct_seqs", "cls", "confidence", "review", "score", "dicer_call",
               "dicer_call_2", "size_purity", "frac_in_window", "frac_siRNA", "frac_piRNA", "frac_top",
               "sense_frac", "pct_1u", "pct_5a", "pct_10a", "pingpong_frac", "pingpong_z", "pingpong_pairs10",
               "phase_score", "unique_frac", "multimap_frac", "major_rna", "evidence"]


def tabela_aprovadas(outdir: str, queries=None, passo=None) -> dict:
    """
    One row per approved locus x library with its small RNA features, written to
    13_curadoria/srna/eve_srna_features.tsv (loci without reads are kept, with reads = 0).
    """
    import csv
    from evetool.curation.curation import listar_candidatas
    passo = passo or (lambda t: None)
    itens = [i for i in listar_candidatas(outdir)["items"] if not i["removed_flag"] and i["contig"] and i["end"]]
    if queries is not None:
        quais = set(queries)
        itens = [i for i in itens if i["query"] in quais]
    loci = {}
    for i in itens:
        l = loci.setdefault((i["species"], i["locus"] or i["query"]), {"rep": i, "start": i["start"], "end": i["end"]})
        l["start"], l["end"] = min(l["start"], i["start"]), max(l["end"], i["end"])
        if i["locus_representative"] and not l["rep"]["locus_representative"]:
            l["rep"] = i
    libs = bibliotecas(outdir)
    destino = os.path.join(_pasta(outdir), "eve_srna_features.tsv")
    os.makedirs(_pasta(outdir), exist_ok=True)
    linhas = 0
    with open(destino, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUNAS_TSV)
        lista = [l for l in loci.values() if libs.get(l["rep"]["species"])]
        for k, l in enumerate(sorted(lista, key=lambda x: (x["rep"]["species"], x["rep"]["contig"], x["start"]))):
            i = l["rep"]
            if k % 25 == 0:
                passo(f"small RNA features: locus {k + 1} of {len(lista)}")
            try:
                r = perfil(outdir, i["species"], i["contig"], l["start"], l["end"], 10, [[l["start"], l["end"], i["strand"]]])
            except Exception:
                continue
            for lib in r["libraries"]:
                f = (lib.get("focus") or [{}])[0]
                cm = f.get("class_mass") or {}
                w.writerow([i["species"], i["orf"], i["contig"], l["start"] + 1, l["end"], i["strand"], i["family"],
                            i["protein_class"], lib["name"]] + [
                    f.get("reads", 0), f.get("reads_rpm", 0), f.get("distinct_seqs", 0), f.get("cls", ""),
                    f.get("confidence", ""), "yes" if f.get("review") else ("" if not f.get("reads") else "no"),
                    "" if f.get("score") is None else f.get("score", ""), f.get("dicer_call", ""),
                    f.get("dicer_call_2", ""), f.get("size_purity", ""), f.get("frac_in_window", ""),
                    cm.get("siRNA", ""), cm.get("piRNA", ""), f.get("frac_top", ""), f.get("sense_frac", ""),
                    f.get("pct_1u", ""), f.get("pct_5a", ""), f.get("pct_10a", ""), f.get("pingpong_frac", ""),
                    f.get("pingpong_z", ""), f.get("pingpong_pairs10", ""), f.get("phase_score", ""),
                    f.get("unique_frac", ""), f.get("multimap_frac", ""), f.get("major_rna", ""), f.get("evidence", "")])
                linhas += 1
    return {"file": destino, "loci": len(loci), "rows": linhas,
            "libraries": sum(len(v) for v in libs.values())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--reads", help="adapter-trimmed small RNA reads (FASTQ/FASTA, .gz or not) to map")
    ap.add_argument("--bam", help="register a small RNA BAM made elsewhere instead of mapping")
    ap.add_argument("--name", default="")
    ap.add_argument("--ncrna-index", default="", help="Bowtie index prefix of structural RNA to drop first")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--region", help="contig:start-end (1-based, inclusive): print its features")
    ap.add_argument("--table", action="store_true", help="features of every approved locus, to a TSV")
    a = ap.parse_args()
    if a.reads:
        print(json.dumps(mapear(a.outdir, a.species, a.reads, a.name, a.threads, a.ncrna_index), indent=1))
    if a.bam:
        print(json.dumps(adicionar_bam(a.outdir, a.species, a.bam, a.name), indent=1))
    if a.region:
        contig, faixa = a.region.rsplit(":", 1)
        s, e = (int(x.replace(",", "")) for x in faixa.split("-"))
        r = perfil(a.outdir, a.species, contig, s - 1, e, 50, [[s - 1, e, "."]])
        for lib in r["libraries"]:
            f = (lib.get("focus") or [{}])[0]
            print(lib["name"], lib.get("error") or {k: f.get(k) for k in (
                "reads", "reads_rpm", "cls", "confidence", "review", "dicer_call", "frac_top", "pct_1u",
                "pingpong_frac", "pingpong_z", "unique_frac")})
    if a.table:
        print(json.dumps(tabela_aprovadas(a.outdir), indent=1))
    if not (a.reads or a.bam or a.region or a.table):
        print(json.dumps(resumo(a.outdir), indent=1))


if __name__ == "__main__":
    main()
