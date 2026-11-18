"""
RNA-seq read coverage of any region, on demand, for the curation page (inspection, context
clusters and tblastn segments).

The BAMs of a species come from:
  * <species folder>/10_rnaseq_expressao/samples.json, written by the coverage step;
  * otherwise the BAMs left in 10_rnaseq_expressao/bam/ and bam_index/ by older runs;
  * BAMs added on the page, or mapped from reads there (`mapear`, minimap2 as ViralQuest),
    kept in 13_curadoria/rnaseq_bams.json.

A BAM may be aligned to the whole genome or to the pipeline's context windows (see
rnaseq_coverage.BamReference). In runs on extracted segments the region is first placed on the
original contig (calibrated segment origins), so a genome BAM covers it too.

  python -m evetool.annotation.rnaseq_tracks -o <results> --species Sp --region contig:1000-5000
"""
import argparse
import glob
import json
import os
import threading
from collections import OrderedDict

import numpy as np

from evetool.annotation.rnaseq_coverage import (BamReference, EveLocus, Window, is_bam, prepare_bam,
                                                profile_window)

SUBPASTA = "13_curadoria"
ARQ_ADICIONADOS = "rnaseq_bams.json"
MAX_REGIAO = 2_000_000          # nt; larger regions are refused (the page asks for one cluster or segment)
MAX_BINS = 1200

_trava = threading.Lock()
_refs = {}                      # (bam, mtime, contigs) -> BamReference
_origens = {}                   # (outdir, species) -> calibrated segment origins
_cache = OrderedDict()          # request -> answer (last 64)


def _pastas(outdir: str) -> dict:
    from evetool.curation.curation import pastas_especies
    return pastas_especies(outdir)


def _arq_adicionados(outdir: str) -> str:
    return os.path.join(outdir, SUBPASTA, ARQ_ADICIONADOS)


def _ler_adicionados(outdir: str) -> dict:
    try:
        with open(_arq_adicionados(outdir), encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def amostras(outdir: str) -> dict:
    """{species: [{name, bam, index, library, source, added}]} for every species of the run."""
    adicionados = _ler_adicionados(outdir)
    saida = {}
    for sp, pasta in _pastas(outdir).items():
        lista, vistos = [], set()
        base = os.path.join(pasta, "10_rnaseq_expressao")
        registro = os.path.join(base, "samples.json")
        if os.path.exists(registro):
            try:
                with open(registro, encoding="utf-8") as fh:
                    for a in json.load(fh):
                        if a.get("bam") and os.path.exists(a["bam"]) and a["name"] not in vistos:
                            vistos.add(a["name"])
                            lista.append({**a, "added": False})
            except Exception:
                pass
        # Older runs: the BAMs left by the coverage step, without their library type
        for cam in sorted(glob.glob(os.path.join(base, "bam", "*.bam")) + glob.glob(os.path.join(base, "bam_index", "*.bam"))):
            nome = os.path.basename(cam)[:-4]
            if nome not in vistos:
                vistos.add(nome)
                lista.append({"name": nome, "bam": cam, "index": None, "library": "unknown",
                              "source": cam, "added": False})
        for a in adicionados.get(sp, []):
            if a["name"] not in vistos and os.path.exists(a["bam"]):
                vistos.add(a["name"])
                lista.append({**a, "added": True})
        if lista:
            saida[sp] = lista
    return saida


def resumo_amostras(outdir: str) -> dict:
    """What the page needs to know up front: species -> sample names and library types."""
    return {sp: [{"name": a["name"], "library": a.get("library") or "unknown", "added": a["added"]} for a in l]
            for sp, l in amostras(outdir).items()}


def adicionar(outdir: str, species: str, bam: str, library: str = "unknown") -> dict:
    """Add a BAM (aligned to the genome or to context windows) to one species of the run."""
    bam = os.path.abspath(os.path.expanduser(str(bam or "").strip()))
    if not os.path.exists(bam):
        raise RuntimeError(f"file not found: {bam}")
    if not is_bam(bam):
        raise RuntimeError(f"not a BAM file: {bam}")
    if species not in _pastas(outdir):
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    if library not in ("fr", "rf", "unstranded", "unknown"):
        raise RuntimeError(f"unknown library type: {library}")
    try:
        caminho, indice = prepare_bam(bam, os.path.join(outdir, SUBPASTA, "rnaseq_bam_index"))
    except Exception as erro:
        raise RuntimeError(f"could not index the BAM; it must be coordinate-sorted (samtools sort): {erro}")
    nome = os.path.basename(bam)
    nome = nome[:-4] if nome.lower().endswith(".bam") else nome
    dados = _ler_adicionados(outdir)
    lista = [a for a in dados.get(species, []) if a["name"] != nome]
    lista.append({"name": nome, "bam": os.path.abspath(caminho), "index": indice, "library": library, "source": bam})
    dados[species] = lista
    os.makedirs(os.path.dirname(_arq_adicionados(outdir)), exist_ok=True)
    with open(_arq_adicionados(outdir), "w", encoding="utf-8") as fh:
        json.dump(dados, fh, indent=1)
    _cache.clear()
    return {"species": species, "name": nome, "samples": resumo_amostras(outdir)}


TIPOS_LEITURA = {"sr": "Illumina short reads", "ont": "Oxford Nanopore", "pb": "PacBio CLR", "hifi": "PacBio HiFi"}


def mapear(outdir: str, species: str, r1: str, r2: str = "", read_type: str = "sr", library: str = "unknown",
           reference: str = "genome", flank: int = 20000, nome: str = "", threads: int = 4, passo=None) -> dict:
    """
    Map long-RNA (or any other RNA) reads of one species and add the BAM to its tracks, as
    ViralQuest maps its reads: `minimap2 -ax <preset> | samtools sort`, preset by read technology.
    `reference` = "genome" (the whole genome: no read is forced onto an EVE) or "windows" (the
    EVE loci +/- `flank`, much faster, but reads of similar sequences elsewhere land here too).
    """
    from evetool.annotation.rnaseq_coverage import (Sample, align_reads, build_windows, contig_sizes_for,
                                                    load_loci, sample_name, write_window_fasta, MM2_PRESET)
    from evetool.curation.curation import genomas_para_sintenia
    passo = passo or (lambda t: print(f"[RNAseq] {t}", flush=True))
    pastas = _pastas(outdir)
    if species not in pastas:
        raise RuntimeError(f"unknown species: {species or '(single genome)'}")
    if read_type not in MM2_PRESET:
        raise RuntimeError(f"unknown read type: {read_type}")
    if library not in ("fr", "rf", "unstranded", "unknown"):
        raise RuntimeError(f"unknown library type: {library}")
    arquivos = [os.path.abspath(os.path.expanduser(x.strip())) for x in (r1, r2) if x and x.strip()]
    if not arquivos:
        raise RuntimeError("give the path of the reads (FASTQ)")
    for a in arquivos:
        if not os.path.exists(a):
            raise RuntimeError(f"file not found: {a}")
    if len(arquivos) == 2 and read_type != "sr":
        raise RuntimeError("paired files are for short reads only")
    g = genomas_para_sintenia(outdir).get(species) or {}
    destino = os.path.join(outdir, SUBPASTA, "rnaseq_map", species or "genome")
    os.makedirs(destino, exist_ok=True)
    nome = nome or sample_name(arquivos[0])
    if reference == "windows":
        genoma = g.get("own")
        if not genoma or not os.path.exists(genoma):
            raise RuntimeError(f"genome of '{species or 'this run'}' not found")
        loci = load_loci(pastas[species], species or None)
        if not loci:
            raise RuntimeError("no EVE locus to build the windows from")
        janelas = build_windows(loci, int(flank), contig_sizes_for(genoma, pastas[species]))
        passo(f"writing {len(janelas)} EVE window(s) (± {int(flank):,} bp)")
        ref = write_window_fasta(genoma, janelas, os.path.join(destino, f"eve_windows_flank{int(flank)}.fasta"))
        nome += f".windows{int(flank)}"
    else:
        ref = g.get("genome")
        if not ref or not os.path.exists(ref):
            raise RuntimeError(f"genome of '{species or 'this run'}' not found")
    passo(f"minimap2 -ax {MM2_PRESET[read_type]} on {'the EVE windows' if reference == 'windows' else 'the whole genome'}")
    bam = align_reads(Sample(nome, r1=arquivos[0], r2=arquivos[1] if len(arquivos) > 1 else None),
                      ref, destino, threads, read_type, force=True)
    passo("adding the BAM to the tracks")
    r = adicionar(outdir, species, bam, library)
    return {**r, "bam": bam, "read_type": read_type, "reference": reference}


def remover(outdir: str, species: str, nome: str) -> dict:
    """Forget a BAM added on the page (the file itself is not touched)."""
    dados = _ler_adicionados(outdir)
    dados[species] = [a for a in dados.get(species, []) if a["name"] != nome]
    with open(_arq_adicionados(outdir), "w", encoding="utf-8") as fh:
        json.dump(dados, fh, indent=1)
    _cache.clear()
    return {"samples": resumo_amostras(outdir)}


def _referencia(a: dict, contigs: set) -> BamReference:
    chave = (a["bam"], os.path.getmtime(a["bam"]), tuple(sorted(contigs)))
    with _trava:
        if chave not in _refs:
            bam, indice = a["bam"], a.get("index")
            if not indice or not os.path.exists(indice):
                # never written next to the user's BAM: a missing index goes to the user cache
                bam, indice = prepare_bam(bam, os.path.join(os.path.expanduser("~"), ".cache", "eve_tool", "bam_index"))
            _refs[chave] = BamReference(bam, contigs, indice)
        return _refs[chave]


def _origem(outdir: str, species: str) -> dict:
    """Calibrated origins of the segments of a species run on extracted segments (else {})."""
    chave = (outdir, species)
    if chave not in _origens:
        from evetool.curation.curation import genomas_para_sintenia
        from evetool.comparative.flank_synteny import origens_de_segmentos, calibrar_origens
        g = genomas_para_sintenia(outdir).get(species) or {}
        _origens[chave] = calibrar_origens(origens_de_segmentos(g["own"]), g["own"], g["genome"]) \
            if g.get("segments") and g.get("own") and g.get("genome") else {}
    return _origens[chave]


def _baixa_cobertura(arr: np.ndarray, media: float, desloc: int, fracao: float = 0.20) -> list:
    """
    Internal runs where depth falls below 20% of the mean, as ViralQuest flags them: runs that
    touch either end or are shorter than max(20 nt, 1% of the region) are left out.
    [[start, end], ...], 0-based half-open, in the coordinates asked for.
    """
    if media <= 0:
        return []
    abaixo = arr < media * fracao
    if not abaixo.any():
        return []
    minimo = max(20, int(arr.size * 0.01))
    bordas = np.flatnonzero(np.diff(np.concatenate(([0], abaixo.view(np.int8), [0])))).reshape(-1, 2)
    return [[int(a) + desloc, int(b) + desloc] for a, b in bordas
            if a != 0 and b != arr.size and b - a >= minimo][:50]


def _binar(arr: np.ndarray, ini: np.ndarray, fn) -> list:
    return [round(float(x), 2) for x in fn.reduceat(arr, ini)]


def cobertura(outdir: str, species: str, contig: str, start: int, end: int, bins: int = 600,
              foco=None, min_mapq: int = 0) -> dict:
    """
    Binned read depth of every sample of `species` over contig:[start, end) (0-based, half-open),
    in the coordinates asked for (also when the reads sit on the original contig of a segment).
    `foco` = [[start, end], ...] spans (the EVEs) whose reads are classed as exonic (aligned bases
    inside) or intron-spanning (only a splice gap crosses it).
    """
    start, end = int(start), int(end)
    if end <= start:
        raise RuntimeError("empty region")
    if end - start > MAX_REGIAO:
        raise RuntimeError(f"region too large ({end - start:,} nt; at most {MAX_REGIAO:,})")
    bins = max(10, min(int(bins or 600), MAX_BINS, end - start))
    foco = [(int(a), int(b)) for a, b in (foco or []) if int(b) > int(a)]
    chave = (outdir, species, contig, start, end, bins, tuple(foco), min_mapq)
    if chave in _cache:
        return _cache[chave]

    lista = amostras(outdir).get(species, [])
    resposta = {"species": species, "contig": contig, "start": start, "end": end, "bins": bins,
                "samples": [], "placed_on": None}
    if not lista:
        resposta["note"] = "no RNA-seq BAM for this species"
        return resposta

    # Region of a segment run -> the same span on the original contig
    from evetool.comparative.flank_synteny import mapear_para_genoma
    g_contig, g0, g1, invertido = mapear_para_genoma(contig, start, end, _origem(outdir, species))
    if g_contig != contig:
        resposta["placed_on"] = {"contig": g_contig, "start": g0, "end": g1, "flipped": invertido}
    foco_g = [mapear_para_genoma(contig, a, b, _origem(outdir, species))[1:3] for a, b in foco]

    cortes = np.array_split(np.arange(g1 - g0), bins)
    ini = np.array([c[0] for c in cortes])
    for a in lista:
        item = {"name": a["name"], "library": a.get("library") or "unknown", "added": a.get("added", False)}
        try:
            ref = _referencia(a, {contig, g_contig})
            janela = Window("region", g_contig, g0, g1,
                            loci=[EveLocus(f"f{k}", g_contig, x, y, ".", f"f{k}") for k, (x, y) in enumerate(foco_g)])
            prof = profile_window(ref, janela, min_mapq)
        except Exception as erro:
            item["error"] = str(erro)
            resposta["samples"].append(item)
            continue
        mais, menos = prof.plus, prof.minus
        if item["library"] == "rf":
            mais, menos = menos, mais
        prof_d, in_ref = prof.depth, prof.in_ref.astype(float)
        if invertido:                         # back to the orientation of the segment asked for
            prof_d, in_ref = prof_d[::-1], in_ref[::-1]
            mais, menos = menos[::-1], mais[::-1]
        tam = np.diff(np.append(ini, len(prof_d)))
        # ViralQuest's track: bar height = depth, bar colour = mean base quality of the bin (over the
        # covered positions only), plus the uniformity figures of its caption
        qsum = prof.qsum[::-1] if invertido else prof.qsum
        soma_d = np.add.reduceat(prof_d, ini)
        qual = np.divide(np.add.reduceat(qsum, ini), soma_d, out=np.zeros_like(soma_d), where=soma_d > 0)
        media = float(prof_d.mean())
        item.update({
            "quality": [round(float(x), 1) for x in qual],
            "mean_quality": round(float(qsum.sum() / prof_d.sum()), 1) if prof_d.sum() > 0 else 0.0,
            "breadth_1x": round(float((prof_d >= 1).mean()), 4),
            "cv": round(float(prof_d.std() / media), 3) if media > 0 else 0.0,
            "mean_plus": round(float(mais.mean()), 2), "mean_minus": round(float(menos.mean()), 2),
            "low_cov": _baixa_cobertura(prof_d, media, start),
            "mean": [round(x, 2) for x in (np.add.reduceat(prof_d, ini) / tam).tolist()],
            "max": _binar(prof_d, ini, np.maximum),
            "plus": [round(x, 2) for x in (np.add.reduceat(mais, ini) / tam).tolist()],
            "minus": [round(x, 2) for x in (np.add.reduceat(menos, ini) / tam).tolist()],
            "in_ref": [round(x, 2) for x in (np.add.reduceat(in_ref, ini) / tam).tolist()],
            "region_mean": round(float(prof_d.mean()), 2),
            "in_reference_pct": round(100 * float(in_ref.mean()), 1),
            "focus": [{"start": x, "end": y, **prof.reads[f"f{k}"],
                       "mean_depth": round(float(prof.depth[max(0, x - g0):max(0, y - g0)].mean()), 2) if y > x else 0}
                      for k, (x, y) in enumerate(foco_g)],
        })
        # spans back in the coordinates asked for (segment or contig)
        item["focus"] = [{**f, "start": a0, "end": b0} for f, (a0, b0) in zip(item["focus"], foco)]
        resposta["samples"].append(item)
    resposta["max_depth"] = max([max(s["max"]) for s in resposta["samples"] if s.get("max")] or [0])
    with _trava:
        _cache[chave] = resposta
        while len(_cache) > 64:
            _cache.popitem(last=False)
    return resposta


_cache_matriz = OrderedDict()   # (outdir, species, data stamp, BAMs) -> answer (last 8)


def matriz(outdir: str, species: str) -> dict:
    """
    Transcription of every approved EVE locus of `species` in every RNA-seq sample, for the
    heatmap of the page: reads with aligned bases inside the locus (exonic), reads whose splice
    gap only spans it, mean depth, the share of the locus covered and the mean depth per strand.
    Cached until the curation, the data on disk or the list of BAMs changes.
    """
    from evetool.curation.curation import estado_dados, loci_aprovados
    lista = amostras(outdir).get(species, [])
    marca = tuple((a["name"], a["bam"], os.path.getmtime(a["bam"]) if os.path.exists(a["bam"]) else 0) for a in lista)
    chave = (outdir, species, estado_dados(outdir), marca)
    if chave in _cache_matriz:
        return _cache_matriz[chave]
    resposta = {"species": species, "loci": [],
                "samples": [{"name": a["name"], "library": a.get("library") or "unknown"} for a in lista]}
    if not lista:
        resposta["note"] = "no RNA-seq BAM for this species"
        return resposta
    for l in loci_aprovados(outdir, species):
        celulas = {}
        try:
            r = cobertura(outdir, species, l["contig"], l["start"], l["end"], 10, [[l["start"], l["end"]]])
            for s in r["samples"]:
                if s.get("error"):
                    celulas[s["name"]] = {"error": s["error"]}
                    continue
                f = (s.get("focus") or [{}])[0]
                celulas[s["name"]] = {"reads": f.get("exonic", 0), "spanning": f.get("spanning", 0),
                                      "mean_depth": f.get("mean_depth", 0), "breadth": s.get("breadth_1x"),
                                      "mean_plus": s.get("mean_plus"), "mean_minus": s.get("mean_minus")}
        except Exception as erro:
            l = {**l, "error": str(erro)}
        resposta["loci"].append({**l, "cells": celulas})
    with _trava:
        _cache_matriz[chave] = resposta
        while len(_cache_matriz) > 8:
            _cache_matriz.popitem(last=False)
    return resposta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--region", help="contig:start-end (1-based, inclusive)")
    ap.add_argument("--bins", type=int, default=50)
    a = ap.parse_args()
    print(json.dumps(resumo_amostras(a.outdir), indent=1))
    if a.region:
        contig, faixa = a.region.rsplit(":", 1)
        s, e = (int(x.replace(",", "")) for x in faixa.split("-"))
        r = cobertura(a.outdir, a.species, contig, s - 1, e, a.bins)
        for am in r["samples"]:
            print(am["name"], am.get("error") or f"mean {am['region_mean']} · in BAM {am['in_reference_pct']}% · max {max(am['max'])}")


if __name__ == "__main__":
    main()
