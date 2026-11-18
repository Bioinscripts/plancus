#!/usr/bin/env python3
"""
Genomic DNA reads over the EVEs: is the insertion really in the genome, and in which samples?

Two questions, one mapping:

  assembly reads     the DNA reads the genome was assembled from. An integrated EVE is covered
                     like its flanks (EVE/flank ~ 1), like the rest of the genome (EVE/genome ~ 1),
                     and reads cross both EVE/host junctions. A depth far above the genome points
                     to an exogenous virus or a contaminant collapsed into the assembly; no read
                     over a junction points to a chimeric contig.
  population reads   DNA of other individuals or pools (--sample-type individual|pool), to see in
                     which of them the EVE is present. Only the measurements are stored; the
                     presence/absence call is made on the curation page, with thresholds the
                     researcher sets while looking at the matrix.

Reads are mapped on the WHOLE genome with `minimap2 -ax sr|map-ont|map-pb|map-hifi`, so no read
is forced onto an EVE, but only the alignments inside the EVE loci +/- a flank are kept
(`samtools view -L`): a whole-genome BAM is never written. The genome-wide mean depth is taken
from the stream (`samtools stats`) before the reads are dropped.

Per locus and library (`medir`):

  eve_depth, breadth        mean depth over the EVE and fraction of it covered by >= 1 read
  left/right_depth          mean depth of the host flank on each side (`BORDA` bp, skipping
                            the first `FOLGA` bp next to the junction)
  span_left/right           reads aligned across the junction (>= `ANCORA` bp on each side)
  frag_left/right           the same counting the sequenced fragment of a proper pair
  empty_left/right          reads aligned on the flank and soft-clipped at the junction: the
                            signature of a chromosome WITHOUT the insertion (or with another
                            sequence there)
  mapq0_frac                fraction of the EVE reads that map equally well elsewhere

Files (13_curadoria/dna/):  libraries.json, <species>/<name>.eve_windows.bam, eve_dna_reads.tsv
"""

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess

SUBPASTA = os.path.join("13_curadoria", "dna")
PRESET = {"sr": "sr", "ont": "map-ont", "pb": "map-pb", "hifi": "map-hifi"}
TIPOS_AMOSTRA = ("individual", "pool")
PAPEIS = ("assembly", "population")
FLANCO_JANELA = 5000      # kept around each locus in the BAM
BORDA = 1000              # flank used for the depth of the host side
FOLGA = 50                # skipped next to the junction (alignment edge effects)
ANCORA = 10               # bases a read must have on each side of a junction to span it
CLIP_MIN = 10             # soft-clip length that counts as "clipped at the junction"
CLIP_TOL = 5              # distance between the clip and the junction

COLUNAS_TSV = ["species", "orf", "locus", "contig", "start", "end", "length", "removed", "library", "role",
               "sample_type", "group", "genome_depth", "eve_depth", "breadth", "left_depth", "right_depth",
               "flank_depth", "eve_over_flank", "eve_over_genome", "span_left", "span_right", "frag_left",
               "frag_right", "empty_left", "empty_right", "mapq0_frac"]

_cache = {}


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


def _tem_dados(b: dict) -> bool:
    """A library is usable while its BAM exists, or its precomputed measurements (`importar_medidas`)."""
    return os.path.exists(b.get("bam") or "") or os.path.exists(b.get("measurements") or "")


def bibliotecas(outdir: str) -> dict:
    """{species: [library]} of the DNA libraries whose BAM (or measurement table) still exists."""
    return {sp: [b for b in l if _tem_dados(b)]
            for sp, l in _ler_registro(outdir).items() if l}


def resumo(outdir: str) -> dict:
    campos = ("name", "role", "sample_type", "group", "read_type", "genome_depth", "mapped_reads",
              "window_reads", "aligner", "source", "windows_only", "precomputed")
    return {sp: [{k: b.get(k) for k in campos} for b in l] for sp, l in bibliotecas(outdir).items() if l}


def _registrar(outdir: str, species: str, item: dict) -> dict:
    dados = _ler_registro(outdir)
    dados[species] = [b for b in dados.get(species, []) if b["name"] != item["name"]] + [item]
    _gravar_registro(outdir, dados)
    return item


def remover(outdir: str, species: str, nome: str, apagar: bool = True) -> dict:
    """Forget a library; the window BAM made here is deleted, a BAM added from elsewhere is not."""
    dados = _ler_registro(outdir)
    for b in dados.get(species, []):
        if b["name"] == nome and apagar and b.get("windows_only"):
            base = b.get("bam", "")[:-len(".eve_windows.bam")] if b.get("bam", "").endswith(".eve_windows.bam") else ""
            for arq in (b.get("bam", ""), b.get("bam", "") + ".bai") + ((base + ".stats.txt", base + ".windows.bed") if base else ()):
                if arq and os.path.abspath(arq).startswith(os.path.abspath(_pasta(outdir))) and os.path.exists(arq):
                    os.remove(arq)
    dados[species] = [b for b in dados.get(species, []) if b["name"] != nome]
    # a measurement table shared by imported libraries goes when the last of them is forgotten
    em_uso = {b.get("measurements") for l in dados.values() for b in l}
    for b in _ler_registro(outdir).get(species, []):
        arq = b.get("measurements") or ""
        if (b["name"] == nome and apagar and arq and arq not in em_uso and os.path.exists(arq)
                and os.path.abspath(arq).startswith(os.path.abspath(_pasta(outdir)))):
            os.remove(arq)
    _gravar_registro(outdir, dados)
    return {"libraries": resumo(outdir)}


def _genoma(outdir: str, species: str) -> str:
    from evetool.annotation.srna import _genoma as genoma_da_especie
    return genoma_da_especie(outdir, species)


def _tamanhos(genoma: str) -> dict:
    fai = genoma + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", genoma], check=True)
    with open(fai, encoding="utf-8", errors="replace") as fh:
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in fh if l.strip()}


def loci(outdir: str, species: str = None, itens: list = None) -> list:
    """
    Every locus of the run (removed ones included, so a restore needs no new mapping), in
    genome coordinates: runs made on extracted segments are placed back on the whole genome.
    """
    from evetool.annotation.rnaseq_tracks import _origem
    from evetool.comparative.flank_synteny import mapear_para_genoma
    from evetool.curation.curation import listar_candidatas
    grupos = {}
    for i in (itens if itens is not None else listar_candidatas(outdir)["items"]):
        if (species is not None and i["species"] != species) or not i["contig"] or not i["end"]:
            continue
        g = grupos.setdefault((i["species"], i["locus"] or i["query"]),
                              {"rep": i, "start": i["start"], "end": i["end"], "removed": True})
        g["start"], g["end"] = min(g["start"], i["start"]), max(g["end"], i["end"])
        g["removed"] = g["removed"] and i["removed_flag"]
        if i["locus_representative"] and not g["rep"]["locus_representative"]:
            g["rep"] = i
    origens, saida = {}, []
    for (sp, locus), g in grupos.items():
        if sp not in origens:
            try:
                origens[sp] = _origem(outdir, sp)
            except Exception:
                origens[sp] = {}
        i = g["rep"]
        contig, a, b, _inv = mapear_para_genoma(i["contig"], g["start"], g["end"], origens[sp])
        saida.append({"species": sp, "locus": locus, "query": i["query"], "orf": i["orf"], "family": i["family"],
                      "protein_class": i["protein_class"], "verdict": i.get("verdict", ""),
                      "contig": contig, "start": int(a), "end": int(b), "removed": g["removed"]})
    return sorted(saida, key=lambda l: (l["species"], l["contig"], l["start"]))


def _bed_janelas(outdir: str, species: str, destino: str, flanco: int = FLANCO_JANELA) -> int:
    tam = _tamanhos(_genoma(outdir, species))
    faixas = {}
    for l in loci(outdir, species):
        faixas.setdefault(l["contig"], []).append((max(0, l["start"] - flanco),
                                                   min(tam.get(l["contig"], l["end"] + flanco), l["end"] + flanco)))
    n = 0
    with open(destino, "w") as fh:
        for contig, lista in sorted(faixas.items()):
            lista.sort()
            a, b = lista[0]
            for x, y in lista[1:]:
                if x <= b:
                    b = max(b, y)
                else:
                    fh.write(f"{contig}\t{a}\t{b}\n"); n += 1
                    a, b = x, y
            fh.write(f"{contig}\t{a}\t{b}\n"); n += 1
    return n


def _ler_stats(caminho: str) -> dict:
    d = {}
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("SN\t"):
                p = linha.rstrip("\n").split("\t")
                try:
                    d[p[1].rstrip(":")] = float(p[2])
                except (ValueError, IndexError):
                    pass
    return d


def _nome_limpo(texto: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", texto).strip("_") or "dna"


def _validar(sample_type: str, role: str):
    if sample_type not in TIPOS_AMOSTRA:
        raise RuntimeError(f"sample type must be one of {', '.join(TIPOS_AMOSTRA)}")
    if role not in PAPEIS:
        raise RuntimeError(f"role must be one of {', '.join(PAPEIS)}")


def mapear(outdir: str, species: str, r1: str, r2: str = "", read_type: str = "sr", nome: str = "",
           sample_type: str = "individual", role: str = "population", group: str = "",
           threads: int = 4, passo=None) -> dict:
    """
    Map DNA reads (FASTQ, gzip or not) on the whole genome and keep the alignments over the EVE
    loci +/- FLANCO_JANELA. Returns the registered library.
    """
    from evetool.curation.curation import pastas_especies
    passo = passo or (lambda t: None)
    _validar(sample_type, role)
    if read_type not in PRESET:
        raise RuntimeError(f"read type must be one of {', '.join(PRESET)}")
    for prog in ("minimap2", "samtools"):
        if not shutil.which(prog):
            raise RuntimeError(f"{prog} not found in the PATH")
    if species not in pastas_especies(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    arquivos = [os.path.abspath(os.path.expanduser(x.strip())) for x in (r1, r2) if x and x.strip()]
    if not arquivos:
        raise RuntimeError("give the path of the reads")
    for arq in arquivos:
        if not os.path.exists(arq):
            raise RuntimeError(f"file not found: {arq}")
    genoma = _genoma(outdir, species)
    nome = _nome_limpo(nome or re.sub(r"(_R?1)?(\.f(ast)?q)?(\.gz)?$", "", os.path.basename(arquivos[0]), flags=re.I))
    destino = os.path.join(_pasta(outdir), _nome_limpo(species or "genome"))
    os.makedirs(destino, exist_ok=True)
    bed = os.path.join(destino, f"{nome}.windows.bed")
    n_jan = _bed_janelas(outdir, species, bed)
    if not n_jan:
        raise RuntimeError("no locus with coordinates for this species")
    bam = os.path.join(destino, f"{nome}.eve_windows.bam")
    stats = os.path.join(destino, f"{nome}.stats.txt")
    t = max(1, int(threads))
    passo(f"minimap2 -ax {PRESET[read_type]} on the whole genome, keeping {n_jan} EVE window(s) (± {FLANCO_JANELA:,} bp)")
    q = shlex.quote
    comando = (
        f"set -o pipefail; minimap2 -ax {PRESET[read_type]} -t {t} {q(genoma)} {' '.join(q(a) for a in arquivos)} "
        f"2> {q(bam + '.minimap2.log')} | samtools view -u -F 4 - "
        f"| tee >(samtools stats - 2>/dev/null | grep '^SN' > {q(stats)}) "
        f"| samtools view -u -L {q(bed)} - | samtools sort -@ {min(t, 4)} -m 768M -o {q(bam)} -")
    r = subprocess.run(["bash", "-c", comando], capture_output=True, text=True)
    if r.returncode != 0:
        log = ""
        try:
            with open(bam + ".minimap2.log", encoding="utf-8", errors="replace") as fh:
                log = fh.read()[-600:]
        except OSError:
            pass
        raise RuntimeError(f"mapping failed: {(r.stderr or '')[-400:]} {log}".strip())
    passo("indexing")
    subprocess.run(["samtools", "index", bam], check=True)
    # the stats writer may lag a moment behind the pipe
    import time
    for _ in range(50):
        if os.path.exists(stats) and os.path.getsize(stats) > 0:
            break
        time.sleep(0.2)
    st = _ler_stats(stats) if os.path.exists(stats) else {}
    tam_genoma = sum(_tamanhos(genoma).values())
    n_janela = int(subprocess.run(["samtools", "view", "-c", bam], capture_output=True, text=True).stdout.strip() or 0)
    try:
        os.remove(bam + ".minimap2.log")
    except OSError:
        pass
    return _registrar(outdir, species, {
        "name": nome, "bam": os.path.abspath(bam), "index": os.path.abspath(bam) + ".bai",
        "role": role, "sample_type": sample_type, "group": group.strip(), "read_type": read_type,
        "genome_depth": round(st.get("bases mapped (cigar)", 0) / tam_genoma, 3) if tam_genoma else 0,
        "mapped_reads": int(st.get("reads mapped", 0)), "window_reads": n_janela, "windows_only": True,
        "flank": FLANCO_JANELA, "aligner": f"minimap2 -ax {PRESET[read_type]}", "source": "; ".join(arquivos)})


def adicionar_bam(outdir: str, species: str, bam: str, nome: str = "", sample_type: str = "individual",
                  role: str = "population", group: str = "", read_type: str = "sr", passo=None) -> dict:
    """Register a whole-genome DNA BAM made elsewhere (same genome); it is used in place, not copied."""
    import pysam
    from evetool.annotation.rnaseq_coverage import is_bam, prepare_bam
    from evetool.curation.curation import pastas_especies
    passo = passo or (lambda t: None)
    _validar(sample_type, role)
    bam = os.path.abspath(os.path.expanduser(str(bam or "").strip()))
    if not os.path.exists(bam):
        raise RuntimeError(f"file not found: {bam}")
    if not is_bam(bam):
        raise RuntimeError(f"not a BAM file: {bam}")
    if species not in pastas_especies(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    passo("indexing the BAM")
    caminho, indice = prepare_bam(bam, os.path.join(_pasta(outdir), "bam_index"))
    passo("measuring the genome-wide depth")
    with pysam.AlignmentFile(caminho, index_filename=indice) as fh:
        contigs_bam = set(fh.references)
        tam = _tamanhos(_genoma(outdir, species))
        if not contigs_bam & set(tam):
            raise RuntimeError("the BAM has no contig of this genome: it was mapped on another reference")
        mapeados = fh.mapped
        compr, n = 0, 0
        for rec in fh.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            compr += rec.query_alignment_length
            n += 1
            if n >= 20000:
                break
    media = compr / n if n else 0
    nome = _nome_limpo(nome or re.sub(r"\.bam$", "", os.path.basename(bam), flags=re.I))
    return _registrar(outdir, species, {
        "name": nome, "bam": os.path.abspath(caminho), "index": indice, "role": role,
        "sample_type": sample_type, "group": group.strip(), "read_type": read_type,
        "genome_depth": round(mapeados * media / sum(tam.values()), 3) if tam else 0,
        "mapped_reads": int(mapeados), "window_reads": None, "windows_only": False,
        "aligner": "external BAM", "source": bam})


CAMPOS_MEDIDA = {"eve_depth": float, "breadth": float, "left_depth": float, "right_depth": float,
                 "left_bp": int, "right_bp": int, "span_left": int, "span_right": int, "frag_left": int,
                 "frag_right": int, "empty_left": int, "empty_right": int, "mapq0_frac": float}


def _ler_medidas(caminho: str) -> dict:
    """{(species, library, contig, start, end): measurements} of a table written by `importar_medidas`."""
    saida = {}
    with open(caminho, encoding="utf-8") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            m = {k: (None if r.get(k) in ("", None) else f(float(r[k])) if f is int else f(r[k]))
                 for k, f in CAMPOS_MEDIDA.items()}
            for k in ("eve_depth", "breadth"):
                m[k] = m[k] or 0.0
            saida[(r["species"], r["library"], r["contig"], int(r["start"]), int(r["end"]))] = m
    return saida


def importar_medidas(outdir: str, species: str, medidas_tsv: str, bibliotecas_json: str, bams: str = "") -> dict:
    """
    Register DNA libraries mapped elsewhere (e.g. a cluster) from the measurements made there with
    `_medir_locus`, without their reads or BAM: `medidas_tsv` has one row per library x locus
    (species, library, locus, contig, start, end + CAMPOS_MEDIDA) and `bibliotecas_json` the
    libraries (name, role, sample_type, group, read_type, genome_depth, mapped_reads, window_reads,
    aligner, source). A library whose window BAM is in the folder `bams`
    (<name>.eve_windows.bam) is registered with that BAM instead, copied here.
    A precomputed row is used only while the locus keeps the coordinates it was measured on; a
    locus edited later shows no data for that library until it is measured again.
    """
    from evetool.curation.curation import pastas_especies
    if species not in pastas_especies(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    destino = os.path.join(_pasta(outdir), _nome_limpo(species or "genome"))
    os.makedirs(destino, exist_ok=True)
    tabela_local = os.path.join(destino, _nome_limpo(os.path.basename(medidas_tsv)))
    if os.path.abspath(medidas_tsv) != os.path.abspath(tabela_local):
        shutil.copy(medidas_tsv, tabela_local)
    medidas = _ler_medidas(tabela_local)
    atuais = {(species, l["contig"], l["start"], l["end"]) for l in loci(outdir, species)}
    with open(bibliotecas_json, encoding="utf-8") as fh:
        libs = [b for b in json.load(fh) if b.get("kind", "dna") == "dna"]
    feitos = []
    for b in libs:
        nome = _nome_limpo(b["name"])
        _validar(b.get("sample_type", "individual"), b.get("role", "population"))
        item = {"name": nome, "role": b.get("role", "population"), "sample_type": b.get("sample_type", "individual"),
                "group": str(b.get("group", "")).strip(), "read_type": b.get("read_type", "sr"),
                "genome_depth": b.get("genome_depth", 0), "mapped_reads": int(b.get("mapped_reads", 0)),
                "window_reads": b.get("window_reads"), "flank": b.get("flank", FLANCO_JANELA),
                "aligner": b.get("aligner", ""), "source": b.get("source", "")}
        bam_fora = os.path.join(bams, f"{b['name']}.eve_windows.bam") if bams else ""
        if bam_fora and os.path.exists(bam_fora):
            bam = os.path.join(destino, f"{nome}.eve_windows.bam")
            for x, y in ((bam_fora, bam), (bam_fora + ".bai", bam + ".bai")):
                if os.path.exists(x) and os.path.abspath(x) != os.path.abspath(y):
                    shutil.copy(x, y)
            if not os.path.exists(bam + ".bai"):
                subprocess.run(["samtools", "index", bam], check=True)
            item.update(bam=os.path.abspath(bam), index=os.path.abspath(bam) + ".bai", windows_only=True)
        else:
            item.update(bam="", index="", measurements=os.path.abspath(tabela_local), windows_only=False,
                        precomputed=True)
        _registrar(outdir, species, item)
        linhas = [k for k in medidas if k[1] == b["name"]]
        feitos.append({"name": nome, "with_bam": bool(item.get("bam")), "rows": len(linhas),
                       "rows_on_current_loci": sum((k[0], k[2], k[3], k[4]) in atuais for k in linhas)})
    return {"species": species, "measurements": tabela_local, "libraries": feitos, "current_loci": len(atuais)}


def editar(outdir: str, species: str, nome: str, **campos) -> dict:
    """Change the role, sample type or group of a library already mapped."""
    dados = _ler_registro(outdir)
    for b in dados.get(species, []):
        if b["name"] == nome:
            for k in ("role", "sample_type", "group"):
                if campos.get(k) is not None:
                    b[k] = str(campos[k]).strip()
            _validar(b["sample_type"], b["role"])
    _gravar_registro(outdir, dados)
    return {"libraries": resumo(outdir)}


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------

def _medir_locus(fh, contig: str, ini: int, fim: int, tam_contig: int) -> dict:
    """Depth, breadth, junction and clip counts of one locus (0-based, half-open) in one BAM."""
    import numpy as np
    a = max(0, ini - BORDA - FOLGA)
    b = min(tam_contig or fim + BORDA + FOLGA, fim + BORDA + FOLGA)
    prof = np.zeros(b - a, dtype=np.int32)
    span = {"l": 0, "r": 0}
    frag = {"l": 0, "r": 0}
    vazio = {"l": 0, "r": 0}
    n_eve = n_mapq0 = 0
    try:
        registros = fh.fetch(contig, a, b)
    except (ValueError, KeyError):
        return None
    for rec in registros:
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.is_duplicate or rec.is_qcfail:
            continue
        rs, re_ = rec.reference_start, rec.reference_end
        if re_ is None:
            continue
        for x, y in rec.get_blocks():
            x, y = max(x, a), min(y, b)
            if y > x:
                prof[x - a:y - a] += 1
        if rs < fim and re_ > ini:
            n_eve += 1
            if rec.mapping_quality == 0:
                n_mapq0 += 1
        # sequenced fragment: both mates of a proper pair, counted once (on the leftmost mate)
        fs, fe, conta_frag = rs, re_, True
        if rec.is_paired and rec.is_proper_pair and rec.template_length:
            if rec.template_length > 0:
                fe = rs + rec.template_length
            else:
                conta_frag = False
        for lado, j in (("l", ini), ("r", fim)):
            if rs <= j - ANCORA and re_ >= j + ANCORA:
                span[lado] += 1
            if conta_frag and fs <= j - ANCORA and fe >= j + ANCORA:
                frag[lado] += 1
        cig = rec.cigartuples or []
        # aligned on the left flank and clipped where the EVE starts
        if cig and cig[-1][0] == 4 and cig[-1][1] >= CLIP_MIN and abs(re_ - ini) <= CLIP_TOL and rs < ini - ANCORA:
            vazio["l"] += 1
        # aligned on the right flank and clipped where the EVE ends
        if cig and cig[0][0] == 4 and cig[0][1] >= CLIP_MIN and abs(rs - fim) <= CLIP_TOL and re_ > fim + ANCORA:
            vazio["r"] += 1
    eve = prof[ini - a:fim - a]
    esq = prof[max(0, ini - FOLGA - BORDA - a):max(0, ini - FOLGA - a)]
    dir_ = prof[fim + FOLGA - a:fim + FOLGA + BORDA - a]
    media = lambda v: round(float(v.mean()), 2) if len(v) else None
    return {"eve_depth": media(eve) or 0.0,
            "breadth": round(float((eve > 0).mean()), 3) if len(eve) else 0.0,
            "left_depth": media(esq), "right_depth": media(dir_),
            "left_bp": int(len(esq)), "right_bp": int(len(dir_)),
            "span_left": span["l"], "span_right": span["r"], "frag_left": frag["l"], "frag_right": frag["r"],
            "empty_left": vazio["l"], "empty_right": vazio["r"],
            "mapq0_frac": round(n_mapq0 / n_eve, 3) if n_eve else None}


def medir(outdir: str, species: str = None, itens: list = None) -> dict:
    """
    {"libraries": {species: [...]}, "loci": [{...locus, "dna": {library: measurements}}]}.
    The flank depth is the mean of the two sides weighted by the bases available on each.
    """
    import pysam
    libs = bibliotecas(outdir)
    if species is not None:
        libs = {species: libs.get(species, [])}
    mtime = lambda p: os.path.getmtime(p) if p and os.path.exists(p) else 0
    chave = (outdir, species, tuple((sp, b["name"], b.get("bam"), b.get("measurements"), mtime(b.get("measurements")),
                                     b.get("role"), b.get("sample_type"), b.get("group"))
                                    for sp, l in sorted(libs.items()) for b in l))
    todos = loci(outdir, species, itens)
    carimbo = tuple((l["query"], l["start"], l["end"], l["removed"]) for l in todos)
    if _cache.get("chave") == (chave, carimbo):
        return _cache["valor"]
    tamanhos, abertos, prontas = {}, {}, {}
    try:
        for sp, lista in libs.items():
            if lista:
                tamanhos[sp] = _tamanhos(_genoma(outdir, sp))
            for b in lista:
                if os.path.exists(b.get("bam") or ""):
                    abertos[(sp, b["name"])] = pysam.AlignmentFile(b["bam"], index_filename=b.get("index"))
                elif b.get("measurements") and b["measurements"] not in prontas:
                    prontas[b["measurements"]] = _ler_medidas(b["measurements"])
        saida = []
        for l in todos:
            lista = libs.get(l["species"], [])
            if not lista:
                continue
            l = dict(l, dna={})
            for b in lista:
                if (l["species"], b["name"]) in abertos:
                    m = _medir_locus(abertos[(l["species"], b["name"])], l["contig"], l["start"], l["end"],
                                     tamanhos[l["species"]].get(l["contig"], 0))
                else:                 # measured elsewhere, on these same coordinates (or nothing)
                    m = prontas.get(b.get("measurements"), {}).get(
                        (l["species"], b["name"], l["contig"], int(l["start"]), int(l["end"])))
                    m = dict(m) if m else None
                if m is None:
                    continue
                pesos = [(m["left_depth"], m["left_bp"]), (m["right_depth"], m["right_bp"])]
                bp = sum(n for d, n in pesos if d is not None)
                m["flank_depth"] = round(sum(d * n for d, n in pesos if d is not None) / bp, 2) if bp else None
                m["flank_bp"] = bp
                l["dna"][b["name"]] = m
            saida.append(l)
    finally:
        for fh in abertos.values():
            fh.close()
    valor = {"libraries": resumo(outdir) if species is None else {species: resumo(outdir).get(species, [])},
             "loci": saida, "params": {"flank": BORDA, "skip": FOLGA, "anchor": ANCORA}}
    _cache["chave"], _cache["valor"] = (chave, carimbo), valor
    return valor


def tabela(outdir: str) -> dict:
    """Every locus x library, written to 13_curadoria/dna/eve_dna_reads.tsv (raw measurements, no call)."""
    d = medir(outdir)
    destino = os.path.join(_pasta(outdir), "eve_dna_reads.tsv")
    os.makedirs(_pasta(outdir), exist_ok=True)
    info = {(sp, b["name"]): b for sp, l in d["libraries"].items() for b in l}
    linhas = 0
    razao = lambda x, y: round(x / y, 3) if (x is not None and y) else ""
    with open(destino, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUNAS_TSV)
        for l in d["loci"]:
            for nome, m in l["dna"].items():
                b = info.get((l["species"], nome), {})
                vazio = lambda v: "" if v is None else v
                w.writerow([l["species"], l["orf"], l["locus"], l["contig"], l["start"] + 1, l["end"],
                            l["end"] - l["start"], "yes" if l["removed"] else "no", nome, b.get("role", ""),
                            b.get("sample_type", ""), b.get("group", ""), b.get("genome_depth", ""),
                            m["eve_depth"], m["breadth"], vazio(m["left_depth"]), vazio(m["right_depth"]),
                            vazio(m["flank_depth"]), razao(m["eve_depth"], m["flank_depth"]),
                            razao(m["eve_depth"], b.get("genome_depth")), m["span_left"], m["span_right"],
                            m["frag_left"], m["frag_right"], m["empty_left"], m["empty_right"],
                            vazio(m["mapq0_frac"])])
                linhas += 1
    return {"file": destino, "rows": linhas, "loci": len(d["loci"])}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="results folder of the run")
    ap.add_argument("--species", default="", help="species of the run the reads belong to")
    ap.add_argument("--r1", help="DNA reads (FASTQ, gzip or not)")
    ap.add_argument("--r2", default="", help="mate file of a pair")
    ap.add_argument("--bam", help="whole-genome DNA BAM already mapped on the same genome")
    ap.add_argument("--read-type", default="sr", choices=list(PRESET))
    ap.add_argument("--name", default="")
    ap.add_argument("--sample-type", default="individual", choices=TIPOS_AMOSTRA,
                    help="one individual, or a pool of individuals (the page then reports a frequency)")
    ap.add_argument("--role", default="population", choices=PAPEIS,
                    help="assembly = the reads the genome was assembled from (integration evidence)")
    ap.add_argument("--group", default="", help="population or locality label")
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--table", action="store_true", help="only write eve_dna_reads.tsv")
    a = ap.parse_args()
    if a.r1:
        print(json.dumps(mapear(a.outdir, a.species, a.r1, a.r2, a.read_type, a.name, a.sample_type, a.role,
                                a.group, a.threads, print), indent=1))
    elif a.bam:
        print(json.dumps(adicionar_bam(a.outdir, a.species, a.bam, a.name, a.sample_type, a.role, a.group,
                                       a.read_type, print), indent=1))
    print(json.dumps(tabela(a.outdir), indent=1))


if __name__ == "__main__":
    main()
