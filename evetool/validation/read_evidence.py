#!/usr/bin/env python3
"""
What the mapped reads say about each locus, for the verdict (validation/verdict.py):

  DNA reads of the assembly   integrated (covered like its flanks, reads across both junctions),
                              or a caveat: junction no read crosses, depth far above the flanks
                              and the genome (exogenous virus, contaminant, collapsed repeat),
                              not covered at all
  small RNA                   the locus produces piRNA or siRNA by the sRNAProfiler rules: the
                              host treats it as an insertion, not as one of its genes

Both are optional: a species without such reads simply gets no line. The thresholds are the
defaults of the Genomic DNA tab; population samples are not used here (an EVE absent from
another population is still an EVE).

The result is kept in 13_curadoria/read_evidence.json and recomputed only when a library is
added or removed or the loci change, so listing the candidates stays fast.
"""

import hashlib
import json
import os
import threading

ARQUIVO = os.path.join("13_curadoria", "read_evidence.json")
LARGURA_MIN = 0.8        # breadth of the EVE covered by assembly reads to call it integrated
JUNCAO_MIN = 2           # reads across each junction
EXCESSO = 3.0            # EVE depth over flank (and genome) depth
NAO_COBERTA = 0.2        # breadth at or below which the assembly reads do not support the locus
PROF_MIN = 3.0           # flank depth below which the library says nothing about the locus
SRNA_MIN = 10            # 1/NH-weighted reads for a small RNA class to count

_trava = threading.Lock()
_local = threading.local()
_memoria = {}


def _mtime(caminho: str) -> float:
    try:
        return os.path.getmtime(caminho)
    except OSError:
        return 0.0


def _loci(itens: list) -> dict:
    """(species, locus) -> span of the locus and its representative item."""
    grupos = {}
    for i in itens:
        if not i.get("contig") or not i.get("end"):
            continue
        g = grupos.setdefault((i["species"], i.get("locus") or i["query"]),
                              {"rep": i, "start": i["start"], "end": i["end"]})
        g["start"], g["end"] = min(g["start"], i["start"]), max(g["end"], i["end"])
        if i.get("locus_representative") and not g["rep"].get("locus_representative"):
            g["rep"] = i
    return grupos


def _dna(outdir: str, itens: list) -> dict:
    from evetool.annotation import dna_reads
    libs = {sp: [b for b in l if b.get("role") == "assembly"] for sp, l in dna_reads.bibliotecas(outdir).items()}
    libs = {sp: l for sp, l in libs.items() if l}
    saida = {}
    if not libs:
        return saida
    medidas = dna_reads.medir(outdir, itens=[i for i in itens if i["species"] in libs])
    for l in medidas["loci"]:
        melhor = None
        for b in libs.get(l["species"], []):
            m = l["dna"].get(b["name"])
            if m and (melhor is None or (m.get("flank_depth") or 0) > (melhor[1].get("flank_depth") or 0)):
                melhor = (b, m)
        if not melhor:
            continue
        b, m = melhor
        flanco, gen = m.get("flank_depth"), b.get("genome_depth")
        razao = m["eve_depth"] / flanco if flanco else None
        razao_g = m["eve_depth"] / gen if gen else None
        if m.get("flank_bp") and (flanco or 0) < PROF_MIN and m["eve_depth"] < PROF_MIN:
            continue                                   # no data
        r = {"library": b["name"], "eve_over_flank": None if razao is None else round(razao, 2),
             "breadth": m["breadth"], "span_left": m["span_left"], "span_right": m["span_right"]}
        abertas = [lado for lado, bp, n in (("left", m["left_bp"], m["span_left"]), ("right", m["right_bp"], m["span_right"]))
                   if bp and n < JUNCAO_MIN]
        if m["breadth"] <= NAO_COBERTA:
            r.update(state="not covered", text=f"the assembly reads cover only {m['breadth']:.0%} of it ({b['name']})")
        elif razao is not None and razao > EXCESSO and (razao_g is None or razao_g > EXCESSO):
            r.update(state="excess", text=f"depth {razao:.1f}x its flanks in the assembly reads: exogenous virus, "
                                          f"contaminant or collapsed repeat ({b['name']})")
        elif abertas:
            r.update(state="junction open", text=f"no assembly read across the {' and '.join(abertas)} junction ({b['name']})")
        elif not m["left_bp"] and not m["right_bp"]:
            continue                                   # the contig ends at the EVE: nothing to test
        elif m["breadth"] >= LARGURA_MIN:
            r.update(state="integrated",
                     text=f"integrated: assembly reads cover it like its flanks (EVE/flank "
                          f"{'n/a' if razao is None else f'{razao:.2f}'}) and cross the junctions "
                          f"({m['span_left']} / {m['span_right']} reads, {b['name']})")
        else:
            r.update(state="partial", text=f"assembly reads cover {m['breadth']:.0%} of it ({b['name']})")
        saida[(l["species"], l["locus"])] = r
    return saida


def _srna(outdir: str, loci: dict) -> dict:
    from evetool.annotation import srna
    libs = srna.bibliotecas(outdir)
    saida = {}
    for (sp, locus), g in loci.items():
        if not libs.get(sp):
            continue
        i = g["rep"]
        try:
            r = srna.perfil(outdir, sp, i["contig"], g["start"], g["end"], 10, [[g["start"], g["end"], i.get("strand", ".")]])
        except Exception:
            continue
        melhor = None
        for lib in r.get("libraries", []):
            f = (lib.get("focus") or [{}])[0]
            if (f.get("reads") or 0) >= SRNA_MIN and f.get("cls") in ("piRNA", "siRNA") and not f.get("review"):
                if melhor is None or f["reads"] > melhor[1]["reads"]:
                    melhor = (lib["name"], f)
        if melhor:
            nome, f = melhor
            anti = None if f.get("sense_frac") is None else 1 - f["sense_frac"]
            saida[(sp, locus)] = {
                "library": nome, "cls": f["cls"], "reads": round(f["reads"], 1), "antisense": anti,
                "text": f"{f['cls']}: {f['reads']:,.0f} small RNA reads"
                        + ("" if anti is None else f", {anti:.0%} antisense") + f" ({nome})"}
    return saida


def evidencia(outdir: str, itens: list) -> dict:
    """{(species, locus): {"dna": {...}?, "srna": {...}?}}; empty when no reads are mapped."""
    if getattr(_local, "ocupado", False):              # called back from inside its own computation
        return {}
    from evetool.annotation import dna_reads, srna
    loci = _loci(itens)
    h = hashlib.md5()
    for (sp, locus), g in sorted(loci.items()):
        h.update(f"{sp}|{locus}|{g['start']}|{g['end']}|{g['rep'].get('strand')};".encode())
    assinatura = [_mtime(dna_reads._arq_registro(outdir)), _mtime(srna._arq_registro(outdir)), h.hexdigest()]
    if not assinatura[0] and not assinatura[1]:
        return {}
    with _trava:
        guardado = _memoria.get(outdir)
        if guardado and guardado[0] == assinatura:
            return guardado[1]
        caminho = os.path.join(outdir, ARQUIVO)
        dados = None
        try:
            with open(caminho, encoding="utf-8") as fh:
                d = json.load(fh)
            if d.get("signature") == assinatura:
                dados = {tuple(k.split("\t", 1)): v for k, v in d["loci"].items()}
        except Exception:
            pass
        if dados is None:
            _local.ocupado = True
            try:
                dados = {}
                for fonte, valores in (("dna", _dna(outdir, itens)), ("srna", _srna(outdir, loci))):
                    for chave, v in valores.items():
                        dados.setdefault(chave, {})[fonte] = v
            finally:
                _local.ocupado = False
            try:
                os.makedirs(os.path.dirname(caminho), exist_ok=True)
                with open(caminho + ".tmp", "w", encoding="utf-8") as fh:
                    json.dump({"signature": assinatura, "loci": {"\t".join(k): v for k, v in dados.items()}}, fh)
                os.replace(caminho + ".tmp", caminho)
            except OSError:
                pass
        _memoria[outdir] = (assinatura, dados)
        return dados
