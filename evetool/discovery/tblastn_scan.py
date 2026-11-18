#!/usr/bin/env python3
"""
tblastn scan of the viral references over the host genome (the last analysis).

The ORF screening only sees what orfipy cut into reading frames, so a viral insertion broken
by stops, frameshifts or decay shows up as a few scattered candidates, or not at all. Here the
best viral hits of the first DIAMOND blastx (the viral screening) go back to the genome with
tblastn, which aligns each protein over the six frames of the whole genome:

  1. references   best viral hit(s) of every ORF in the viral screening
                  (02_alinhamento_diamond/triagem_viral_diamond_refseq_filtrado.tsv), with the
                  sequences read from the DIAMOND database of the run;
  2. tblastn      every reference against the genome; every HSP is kept, with the aligned
                  protein and genome sequences;
  3. segments     HSPs on the same contig closer than --gap (1000 nt by default, any strand)
                  are concatenated into one segment;
  4. ORFs         every orfipy ORF that overlaps a segment, already analysed (a candidate EVE)
                  or not. For the new ones the earlier searches are reported (best viral hit of
                  the screening and why it was dropped, HMM hit), and their peptides are
                  searched with CD-Search (local CDD) and Pfam, as the candidates were.

The page redraws the segments with any other distance (the HSPs travel with the data), and can
use the HSPs to bridge the genomic-context clusters.

Outputs, under 05_ancoragem_tblastn/ of each species:
  tblastn_viral_refs.faa              the reference proteins searched
  tblastn_viral_refs_vs_genome.tsv    every HSP (with qseq/sseq)
  tblastn_scan_segments.tsv           one segment per line, at --gap
  tblastn_scan_orfs.tsv               the ORFs of the segments, with every annotation
  tblastn_scan_orfs_pep.faa           peptides of the ORFs that were not candidates
  scan_cdsearch_*.tsv, scan_pfam_*.tsv   CD-Search and Pfam of those peptides (cached)
  tblastn_scan.json                   everything above, read by the curation page

Usage:
  python -m evetool.discovery.tblastn_scan -o results_eve -t 8
  python -m evetool.discovery.tblastn_scan -o results_eve --hits-per-orf 5 --evalue 1e-5 --gap 2000
  python -m evetool.discovery.tblastn_scan -o results_eve --sem-cdd --sem-pfam      # tblastn and ORFs only
"""
import os
import re
import csv
import glob
import json
import bisect
import shutil
import argparse
import hashlib
import subprocess
from collections import OrderedDict, defaultdict

from evetool import config
from evetool.annotation.cdsearch import pastas_especies, _gravar, ESCOPOS_ALVO

SUBPASTA = "05_ancoragem_tblastn"
NOME_JSON = "tblastn_scan.json"
COLS_TBLASTN = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "sframe", "evalue",
                "bitscore", "pident", "length", "qlen", "slen", "qseq", "sseq"]
RX_ORF = re.compile(r"\[(\d+)-(\d+)\]\(([+-])\)")
MAX_ALINHAMENTOS = 20000    # beyond this many HSPs the page gets them without the sequences


def log(msg):
    print(f"[tblastn-scan] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Inputs of a species folder
# ---------------------------------------------------------------------------
def _info(pasta: str) -> dict:
    try:
        with open(os.path.join(pasta, "pipeline_info.json"), encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def _orfs_fasta(pasta: str, tipo: str) -> str:
    """Full orfipy output (before the MMseqs2 deduplication): *_orfs_pep.fasta / *_orfs_nt.fasta."""
    achados = [f for f in glob.glob(os.path.join(pasta, "01_orfs_cds", f"*_orfs_{tipo}.fasta"))]
    return sorted(achados)[0] if achados else ""


def referencias_virais(pasta: str, hits_por_orf: int = 1) -> "OrderedDict[str, dict]":
    """
    {protein: {product, virus, genome_type, family, length, orfs, best_bitscore}} from the viral
    screening: the best `hits_por_orf` subjects of every ORF kept by the screening.
    """
    from evetool.validation.protein_map import _titulo
    caminho = next((c for c in (
        os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond_refseq_filtrado.tsv"),
        os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond.tsv"))
        if os.path.exists(c) and os.path.getsize(c) > 0), None)
    if not caminho:
        return OrderedDict()
    por_orf = defaultdict(dict)
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < 8 or linha.startswith("#"):
                continue
            try:
                bits = float(p[7])
            except ValueError:
                continue
            if bits > por_orf[p[0]].get(p[2], (0,))[0]:
                por_orf[p[0]][p[2]] = (bits, p[4], int(p[3]) if p[3].isdigit() else 0)
    refs = {}
    for orf, hits in por_orf.items():
        for sid, (bits, titulo, tam) in sorted(hits.items(), key=lambda x: -x[1][0])[:hits_por_orf]:
            r = refs.setdefault(sid, {"id": sid, **_titulo(titulo), "length": tam,
                                      "orfs": [], "best_bitscore": 0.0})
            r["orfs"].append(orf)
            r["best_bitscore"] = max(r["best_bitscore"], bits)
    return OrderedDict(sorted(refs.items(), key=lambda x: -x[1]["best_bitscore"]))


# ---------------------------------------------------------------------------
# tblastn
# ---------------------------------------------------------------------------
def _banco_genoma(genoma: str, dir_saida: str) -> str:
    """BLAST nucleotide database of the genome, reused while it is newer than the FASTA."""
    db = os.path.join(dir_saida, os.path.splitext(os.path.basename(genoma))[0] + "_blastdb")
    feitos = glob.glob(db + ".*")
    if feitos and min(os.path.getmtime(f) for f in feitos) >= os.path.getmtime(genoma):
        return db
    log(f"makeblastdb on {os.path.basename(genoma)}...")
    r = subprocess.run(["makeblastdb", "-in", genoma, "-dbtype", "nucl", "-out", db, "-parse_seqids"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        # -parse_seqids refuses some headers; the plain database is enough for tblastn
        r = subprocess.run(["makeblastdb", "-in", genoma, "-dbtype", "nucl", "-out", db],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"makeblastdb failed: {r.stderr.strip()[:400]}")
    return db


def rodar_tblastn(refs_faa: str, db: str, saida: str, evalue: float, threads: int, forcar: bool) -> str:
    """tblastn of the references on the genome; cached by the references + E-value."""
    with open(refs_faa, "rb") as fh:
        assinatura = hashlib.md5(fh.read() + f"|{evalue}|{db}".encode()).hexdigest()
    meta = saida + ".meta.json"
    if not forcar and os.path.exists(saida) and os.path.exists(meta):
        try:
            with open(meta) as fh:
                if json.load(fh).get("signature") == assinatura:
                    log("tblastn: same references and E-value, reusing the table")
                    return saida
        except Exception:
            pass
    if not shutil.which("tblastn"):
        raise RuntimeError("tblastn not found on the PATH (BLAST+ package)")
    log(f"tblastn of the references on the genome (E-value {evalue}, {threads} threads)...")
    tmp = saida + ".tmp"
    r = subprocess.run(["tblastn", "-query", refs_faa, "-db", db, "-evalue", str(evalue),
                        "-outfmt", "6 " + " ".join(COLS_TBLASTN), "-num_threads", str(threads),
                        "-max_target_seqs", "100000", "-out", tmp], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"tblastn failed: {r.stderr.strip()[:500]}")
    with open(tmp, encoding="utf-8") as fi, open(saida, "w", encoding="utf-8") as fo:
        fo.write("\t".join(COLS_TBLASTN) + "\n")
        shutil.copyfileobj(fi, fo)
    os.remove(tmp)
    with open(meta, "w") as fh:
        json.dump({"signature": assinatura, "evalue": evalue, "db": db}, fh)
    return saida


def _midline(q: str, s: str) -> str:
    """BLAST-style middle line: the residue when identical, '+' when positive, ' ' otherwise."""
    try:
        from Bio.Align import substitution_matrices
        m = _midline.matriz = getattr(_midline, "matriz", None) or substitution_matrices.load("BLOSUM62")
    except Exception:
        m = None
    out = []
    for a, b in zip(q, s):
        if a == b and a not in "-*":
            out.append(a)
        elif m is not None and a not in "-*" and b not in "-*" and a in m.alphabet and b in m.alphabet \
                and m[a][b] > 0:
            out.append("+")
        else:
            out.append(" ")
    return "".join(out)


def ler_hsps(tabela: str, refs: dict) -> list:
    """HSPs in genome coordinates (0-based start, end exclusive, like the ORF headers)."""
    hsps = []
    with open(tabela, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                ss, se, fr = int(r["sstart"]), int(r["send"]), int(r["sframe"])
                h = {"ref": r["qseqid"], "contig": r["sseqid"],
                     "start": min(ss, se) - 1, "end": max(ss, se), "strand": "-" if fr < 0 else "+",
                     "frame": fr, "qstart": int(r["qstart"]), "qend": int(r["qend"]),
                     "qlen": int(r["qlen"]), "evalue": float(r["evalue"]), "bitscore": float(r["bitscore"]),
                     "pident": float(r["pident"]), "length": int(r["length"]),
                     "qseq": r["qseq"], "sseq": r["sseq"]}
            except (TypeError, ValueError, KeyError):
                continue
            # BLAST names a contig "ref|NC_1.1|" when the database was built with -parse_seqids
            if "|" in h["contig"]:
                partes = [x for x in h["contig"].split("|") if x]
                h["contig"] = partes[1] if len(partes) > 1 else partes[0]
            ref = refs.get(h["ref"], {})
            h["family"] = ref.get("family", "")
            h["midline"] = _midline(h["qseq"], h["sseq"])
            hsps.append(h)
    hsps.sort(key=lambda h: (h["contig"], h["start"]))
    for k, h in enumerate(hsps):
        h["id"] = k
    return hsps


def segmentos(hsps: list, gap: int = 1000) -> list:
    """HSPs of one contig closer than `gap` (any strand) concatenated into segments."""
    segs = []
    for h in hsps:                                    # already sorted by contig, start
        s = segs[-1] if segs else None
        if s and s["contig"] == h["contig"] and h["start"] - s["end"] <= gap:
            s["end"] = max(s["end"], h["end"])
            s["hsps"].append(h["id"])
        else:
            segs.append({"contig": h["contig"], "start": h["start"], "end": h["end"], "hsps": [h["id"]]})
    for k, s in enumerate(segs, 1):
        s["name"] = f"Seg_{k}"
    return segs


# ---------------------------------------------------------------------------
# ORFs of the segments
# ---------------------------------------------------------------------------
def orfs_nos_segmentos(pasta: str, segs: list, hsps: list) -> "OrderedDict[str, dict]":
    """Every orfipy ORF overlapping a segment, with its overlap with the HSPs."""
    fasta = _orfs_fasta(pasta, "pep")
    if not fasta:
        raise RuntimeError("orfipy output (01_orfs_cds/*_orfs_pep.fasta) not found")
    por_contig = defaultdict(list)
    for s in segs:
        por_contig[s["contig"]].append(s)
    inicios = {c: [s["start"] for s in l] for c, l in por_contig.items()}
    hsps_contig = defaultdict(list)
    for h in hsps:
        hsps_contig[h["contig"]].append(h)

    orfs, atual = OrderedDict(), None
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                atual = None
                nome = linha[1:].split()[0]
                m = RX_ORF.search(linha)
                if not m or "_ORF." not in nome:
                    continue
                contig = nome.rsplit("_ORF.", 1)[0]
                if contig not in por_contig:
                    continue
                a, b, fita = int(m.group(1)), int(m.group(2)), m.group(3)
                l = por_contig[contig]
                # segments are sorted and disjoint: only the last one starting before the
                # ORF ends can overlap it (an ORF spanning two goes to the second)
                k = bisect.bisect_left(inicios[contig], b) - 1
                if k < 0 or l[k]["end"] <= a:
                    continue
                seg = l[k]
                orfs[nome] = {"orf": nome, "contig": contig, "start": a, "end": b, "strand": fita,
                              "segment": seg["name"], "pep": ""}
                atual = orfs[nome]
            elif atual is not None:
                atual["pep"] += linha.strip().replace("*", "")

    # overlap with the HSPs: largest overlap, and whether one HSP shares strand and frame
    for o in orfs.values():
        maior, quadro, refs = 0, False, set()
        for h in hsps_contig[o["contig"]]:
            if h["start"] >= o["end"]:
                break
            sob = min(o["end"], h["end"]) - max(o["start"], h["start"])
            if sob <= 0:
                continue
            maior = max(maior, sob)
            refs.add(h["ref"])
            if h["strand"] == o["strand"]:
                fase = (h["start"] - o["start"]) % 3 if o["strand"] == "+" else (o["end"] - h["end"]) % 3
                quadro = quadro or fase == 0
        o.update({"length_nt": o["end"] - o["start"], "length_aa": len(o["pep"]),
                  "hsp_overlap_nt": maior, "hsp_in_frame": quadro, "hsp_refs": sorted(refs)})
    return orfs


def _representantes(pasta: str, orfs: set) -> dict:
    """ORF -> MMseqs2 representative (the screening only searched the representatives)."""
    caminhos = glob.glob(os.path.join(pasta, "01_orfs_cds", "*_nr_cluster.tsv"))
    rep = {o: o for o in orfs}
    if not caminhos:
        return rep
    with open(caminhos[0], encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] in orfs:
                rep[p[1]] = p[0]
    return rep


def anotacoes_anteriores(pasta: str, orfs: dict, candidatas: dict) -> None:
    """Earlier results of each ORF: candidate status, screening decision, best viral hit, HMM."""
    from evetool.validation.protein_map import _titulo
    rep = _representantes(pasta, set(orfs))
    alvo = set(rep.values())

    decisoes = {}
    c = os.path.join(pasta, "02_alinhamento_diamond", "viral_screening_decision_per_orf.tsv")
    if os.path.exists(c):
        with open(c, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf") in alvo:
                    decisoes[r["orf"]] = r

    melhor = {}
    c = os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond.tsv")
    if os.path.exists(c):
        with open(c, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                q = linha.split("\t", 1)[0]
                if q not in alvo:
                    continue
                p = linha.rstrip("\n").split("\t")
                try:
                    bits = float(p[7])
                except (IndexError, ValueError):
                    continue
                if bits > melhor.get(q, {}).get("bitscore", -1):
                    t = _titulo(p[4])
                    melhor[q] = {"id": p[2], "bitscore": bits, "pident": float(p[5]),
                                 "evalue": p[6], "product": t["product"], "virus": t["virus"],
                                 "family": t["family"], "genome_type": t["genome_type"]}

    hmm = {}
    c = os.path.join(pasta, "03_triagem_hmm", "triagem_profunda_hmm.tsv")
    if os.path.exists(c):
        with open(c, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith("#"):
                    continue
                p = linha.split()
                if len(p) > 13 and p[0] in alvo:
                    try:
                        sc = float(p[13])
                    except ValueError:
                        continue
                    if sc > hmm.get(p[0], (None, -1))[1]:
                        hmm[p[0]] = (p[3], sc)

    for nome, o in orfs.items():
        r = rep.get(nome, nome)
        cand = candidatas.get(nome)
        d, b, hm = decisoes.get(r), melhor.get(r), hmm.get(r)
        if cand:
            status = "candidate"
        elif d:
            status = "screening: " + (d.get("decision") or "?")
        elif b:
            status = "viral hit, not kept"
        else:
            status = "no viral hit"
        o.update({
            "representative": r if r != nome else "",
            "status": status, "candidate": bool(cand),
            "query": cand["query"] if cand else "",
            "group": cand.get("group", "") if cand else "",
            "screening_reason": (d or {}).get("reason", ""),
            "screening_best_class": (d or {}).get("best_hit_class", ""),
            "viral_best_hit": (b or {}).get("id", ""), "viral_bitscore": (b or {}).get("bitscore", ""),
            "viral_pident": (b or {}).get("pident", ""), "viral_evalue": (b or {}).get("evalue", ""),
            "viral_product": (b or {}).get("product", ""), "viral_virus": (b or {}).get("virus", ""),
            "viral_family": (b or {}).get("family", ""), "viral_genome_type": (b or {}).get("genome_type", ""),
            "hmm_model": hm[0] if hm else "", "hmm_score": hm[1] if hm else "",
        })


def em_escopo(o: dict) -> bool:
    """RNA- or ssDNA-virus domain in CD-Search or Pfam."""
    return any(o.get(k) in ESCOPOS_ALVO for k in ("cdd_scope", "pfam_scope"))


def _candidatas(outdir: str, especie: str) -> dict:
    """{orf: row of eves_classified.tsv} of one species."""
    c = os.path.join(outdir, "eves_classified.tsv")
    saida = {}
    if os.path.exists(c):
        with open(c, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if not especie or r.get("species") == especie:
                    saida[r.get("orf", "")] = r
    return saida


def _cdd_pfam(pasta: str, dir_saida: str, faa: str, threads: int, usar_cdd: bool, usar_pfam: bool) -> tuple:
    """CD-Search and Pfam rows of the new ORFs ({orf: summary}, {orf: [domains]}), cached."""
    resumo, dominios = defaultdict(dict), defaultdict(list)
    if usar_cdd:
        try:
            from evetool.annotation import cdsearch as cd
            desc = cd.carregar_dominios_locais(cd.DOMINIOS_PADRAO)
            r = cd.anotar_pasta(pasta, desc, threads=threads, fasta=faa, dir_saida=dir_saida, prefixo="scan_")
            for x in r.get("summary_rows", []):
                resumo[x["orf"]].update({k: x[k] for k in ("cdd_scope", "cdd_best", "cdd_best_scope",
                                         "cdd_best_evalue", "cdd_best_bitscore", "cdd_best_description",
                                         "cdd_in_scope_domains")})
            for h in r.get("hit_rows", []):
                dominios[h["orf"]].append({"src": "CDD", "name": h["short_name"], "scope": h["scope"],
                                           "from": int(h["from"]), "to": int(h["to"]),
                                           "bits": float(h["bitscore"] or 0), "desc": h["description"][:120]})
        except Exception as erro:
            log(f"[warning] CD-Search skipped: {erro}")
    if usar_pfam:
        try:
            from evetool.annotation import pfam as pf
            r = pf.anotar_pasta(pasta, threads=threads, fasta=faa, dir_saida=dir_saida, prefixo="scan_")
            for x in r.get("summary_rows", []):
                resumo[x["orf"]].update({k: x[k] for k in ("pfam_scope", "pfam_best", "pfam_best_scope",
                                         "pfam_best_evalue", "pfam_best_bitscore", "pfam_best_description",
                                         "pfam_in_scope_domains")})
            for h in r.get("hit_rows", []):
                dominios[h["orf"]].append({"src": "Pfam", "name": h["short_name"], "scope": h["scope"],
                                           "from": int(h["from"]), "to": int(h["to"]),
                                           "bits": float(h["bitscore"] or 0), "desc": h["description"][:120]})
        except Exception as erro:
            log(f"[warning] Pfam skipped: {erro}")
    return resumo, dominios


# ---------------------------------------------------------------------------
# One species folder
# ---------------------------------------------------------------------------
def varrer_pasta(outdir: str, especie: str, pasta: str, hits_por_orf: int = 1, evalue: float = 1e-5,
                 gap: int = 1000, threads: int = 4, forcar: bool = False, usar_cdd: bool = True,
                 usar_pfam: bool = True, dmnd_db: str = None) -> dict:
    from evetool.validation.protein_map import banco_da_corrida, buscar_sequencias
    rotulo = especie or os.path.basename(os.path.abspath(pasta))
    genoma = _info(pasta).get("genome")
    if not genoma or not os.path.exists(genoma):
        raise RuntimeError(f"{rotulo}: genome not found (pipeline_info.json)")
    dir_s = os.path.join(pasta, SUBPASTA)
    os.makedirs(dir_s, exist_ok=True)

    # 1. references
    refs = referencias_virais(pasta, hits_por_orf)
    if not refs:
        raise RuntimeError(f"{rotulo}: no viral hit in the screening (02_alinhamento_diamond)")
    db_viral = banco_da_corrida(outdir, dmnd_db)
    seqs = buscar_sequencias(db_viral, set(refs))
    faltam = [r for r in refs if r not in seqs]
    if faltam:
        log(f"{rotulo}: {len(faltam)} reference(s) not in {db_viral} (e.g. {', '.join(faltam[:3])})")
    if not seqs:
        raise RuntimeError(f"{rotulo}: no reference sequence found in {db_viral}")
    refs_faa = os.path.join(dir_s, "tblastn_viral_refs.faa")
    with open(refs_faa, "w") as fh:
        for r in refs:
            if r in seqs:
                fh.write(f">{r} {refs[r]['product']} [{refs[r]['virus']}]\n{seqs[r]}\n")
                refs[r]["length"] = len(seqs[r])
    log(f"{rotulo}: {len(seqs)} viral references (best {hits_por_orf} hit(s) of "
        f"{len({o for r in refs.values() for o in r['orfs']})} screened ORFs)")

    # 2. tblastn
    db = _banco_genoma(genoma, dir_s)
    tabela = rodar_tblastn(refs_faa, db, os.path.join(dir_s, "tblastn_viral_refs_vs_genome.tsv"),
                           evalue, threads, forcar)
    hsps = ler_hsps(tabela, refs)
    log(f"{rotulo}: {len(hsps)} HSPs")

    # 3. segments and 4. their ORFs
    segs = segmentos(hsps, gap)
    orfs = orfs_nos_segmentos(pasta, segs, hsps)
    anotacoes_anteriores(pasta, orfs, _candidatas(outdir, especie))
    novas = [o for o in orfs.values() if not o["candidate"] and o["pep"]]
    log(f"{rotulo}: {len(segs)} segments at {gap} nt, {len(orfs)} ORFs inside them "
        f"({len(orfs) - len(novas)} candidates, {len(novas)} not analysed before)")

    faa = os.path.join(dir_s, "tblastn_scan_orfs_pep.faa")
    with open(faa, "w") as fh:
        fh.writelines(f">{o['orf']}\n{o['pep']}\n" for o in novas)
    resumo, dominios = (_cdd_pfam(pasta, dir_s, faa, threads, usar_cdd, usar_pfam) if novas
                        else ({}, {}))
    for o in orfs.values():
        o.update(resumo.get(o["orf"], {}))
        o["domains"] = sorted(dominios.get(o["orf"], []), key=lambda d: d["from"])

    # Tables
    cols_orf = ["orf", "segment", "contig", "start", "end", "strand", "length_nt", "length_aa", "status",
                "query", "group", "representative", "hsp_overlap_nt", "hsp_in_frame", "hsp_refs",
                "screening_reason", "screening_best_class", "viral_best_hit", "viral_product", "viral_virus",
                "viral_family", "viral_genome_type", "viral_bitscore", "viral_pident", "viral_evalue",
                "hmm_model", "hmm_score",
                "cdd_scope", "cdd_best", "cdd_best_scope", "cdd_best_evalue", "cdd_best_bitscore",
                "cdd_best_description", "cdd_in_scope_domains",
                "pfam_scope", "pfam_best", "pfam_best_scope", "pfam_best_evalue", "pfam_best_bitscore",
                "pfam_best_description", "pfam_in_scope_domains"]
    _gravar(os.path.join(dir_s, "tblastn_scan_orfs.tsv"), cols_orf,
            [{**o, "start": o["start"] + 1, "hsp_refs": ",".join(o["hsp_refs"]),
              "hsp_in_frame": "yes" if o["hsp_in_frame"] else "no"} for o in orfs.values()])
    por_seg = defaultdict(list)
    for o in orfs.values():
        por_seg[o["segment"]].append(o)
    linhas_seg = []
    for s in segs:
        hs = [hsps[i] for i in s["hsps"]]
        os_ = por_seg[s["name"]]
        linhas_seg.append({
            "segment": s["name"], "contig": s["contig"], "start": s["start"] + 1, "end": s["end"],
            "length": s["end"] - s["start"], "n_hsps": len(hs),
            "strands": "/".join(sorted({h["strand"] for h in hs})),
            "best_bitscore": max(h["bitscore"] for h in hs),
            "references": ",".join(OrderedDict.fromkeys(h["ref"] for h in sorted(hs, key=lambda h: -h["bitscore"]))),
            "families": ",".join(sorted({h["family"] for h in hs if h["family"]})),
            "n_orfs": len(os_), "n_candidates": sum(o["candidate"] for o in os_),
            "n_new": sum(not o["candidate"] for o in os_),
            "new_in_scope": sum(1 for o in os_ if not o["candidate"] and em_escopo(o)),
            "candidates": ",".join(o["orf"] for o in os_ if o["candidate"]),
        })
    _gravar(os.path.join(dir_s, "tblastn_scan_segments.tsv"), list(linhas_seg[0]) if linhas_seg else ["segment"],
            linhas_seg)

    # JSON for the page (the ORF peptides stay out: they are in the FASTA)
    com_alinhamento = len(hsps) <= MAX_ALINHAMENTOS
    dados = {
        "species": especie, "genome": genoma, "gap": gap, "evalue": evalue, "hits_per_orf": hits_por_orf,
        "db": db_viral, "references": {r: {k: v for k, v in refs[r].items() if k != "orfs"} | {"n_orfs": len(refs[r]["orfs"])}
                                       for r in refs if r in seqs},
        "hsps": [{k: v for k, v in h.items() if com_alinhamento or k not in ("qseq", "sseq", "midline")}
                 for h in hsps],
        "orfs": [{k: v for k, v in o.items() if k != "pep"} for o in orfs.values()],
        "cdd": usar_cdd, "pfam": usar_pfam,
    }
    with open(os.path.join(dir_s, NOME_JSON) + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(dados, fh, ensure_ascii=False)
    os.replace(os.path.join(dir_s, NOME_JSON) + ".tmp", os.path.join(dir_s, NOME_JSON))
    return {"species": especie, "references": len(seqs), "hsps": len(hsps), "segments": len(segs),
            "orfs": len(orfs), "new_orfs": len(novas),
            "segments_without_candidate": sum(1 for l in linhas_seg if not l["n_candidates"])}


def varrer(outdir: str, **opcoes) -> dict:
    total = {"species": [], "references": 0, "hsps": 0, "segments": 0, "orfs": 0, "new_orfs": 0,
             "segments_without_candidate": 0, "errors": []}
    for especie, pasta in pastas_especies(outdir).items():
        try:
            r = varrer_pasta(outdir, especie, pasta, **opcoes)
        except Exception as erro:
            log(f"[warning] {especie or pasta}: {erro}")
            total["errors"].append(f"{especie or os.path.basename(pasta)}: {erro}")
            continue
        total["species"].append(r)
        for k in ("references", "hsps", "segments", "orfs", "new_orfs", "segments_without_candidate"):
            total[k] += r[k]
    if not total["species"] and total["errors"]:
        raise RuntimeError("; ".join(total["errors"]))
    return total


def carregar(outdir: str) -> list:
    """The saved scans of every species (for the curation page); [] when never run."""
    saida = []
    for especie, pasta in pastas_especies(outdir).items():
        c = os.path.join(pasta, SUBPASTA, NOME_JSON)
        if os.path.exists(c):
            try:
                with open(c, encoding="utf-8") as fh:
                    d = json.load(fh)
                d["species"] = especie
                saida.append(d)
            except Exception as erro:
                log(f"[warning] {c}: {erro}")
    return saida


def arquivos(outdir: str) -> list:
    return [os.path.join(p, SUBPASTA, NOME_JSON) for p in pastas_especies(outdir).values()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--hits-per-orf", type=int, default=1,
                    help="Viral hits of each screened ORF used as references (default 1: the best hit)")
    ap.add_argument("--evalue", type=float, default=1e-5, help="tblastn E-value (default 1e-5)")
    ap.add_argument("--gap", type=int, default=1000,
                    help="HSPs closer than this (nt) form one segment (default 1000; the page can change it)")
    ap.add_argument("--dmnd-db", help="DIAMOND database of the screening (default: the run's)")
    ap.add_argument("--sem-cdd", action="store_true", help="Do not run CD-Search on the new ORFs")
    ap.add_argument("--sem-pfam", action="store_true", help="Do not run Pfam on the new ORFs")
    ap.add_argument("--forcar", action="store_true", help="Run tblastn again even with the same references")
    ap.add_argument("-t", "--threads", type=int, default=4)
    a = ap.parse_args()
    r = varrer(a.outdir, hits_por_orf=a.hits_per_orf, evalue=a.evalue, gap=a.gap, threads=a.threads,
               forcar=a.forcar, usar_cdd=not a.sem_cdd, usar_pfam=not a.sem_pfam, dmnd_db=a.dmnd_db)
    print(f"\n{r['references']} references | {r['hsps']} HSPs | {r['segments']} segments "
          f"({r['segments_without_candidate']} without any candidate) | {r['orfs']} ORFs in the segments, "
          f"{r['new_orfs']} not analysed before")
    for e in r["errors"]:
        print(f"  warning: {e}")


if __name__ == "__main__":
    main()
