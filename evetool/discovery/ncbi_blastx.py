#!/usr/bin/env python3
"""
Online BLASTX (NCBI) of curated segments, in small batches that can always be recovered.

The NCBI queue can take from minutes to hours, and a single large request is easily lost
(timeout, closed browser, server restarted). Here every request is small and its state is
written to disk before and after each step, so a run can stop at any point and be resumed:

  * batches of at most 20 sequences and 10,000 nt (sequences longer than that are split into
    overlapping windows, mapped back to the segment coordinates afterwards);
  * each batch keeps its RID in blastx_jobs.json; resuming polls the same RID instead of
    submitting it again (an RID lives about 24 h at NCBI);
  * a failed or expired batch is resubmitted, split in two after a failure;
  * results are cached by sequence (md5 + program/database/E-value), so a segment that was
    already searched is never sent again;
  * NCBI usage rules: one submission per 10 s, each RID polled at most once a minute, and
    at most 3 batches in the queue at a time.

Outputs, under <results>/13_curadoria/blastx_ncbi/:
  blastx_jobs.json         batches, RIDs and their state (the resume point)
  blastx_queries.tsv       every segment searched: id, origin, coordinates, length
  blastx_hits.tsv          one HSP per line, in segment coordinates, with title and taxon
  blastx_summary.tsv       one segment per line: best hit, best viral and best non-viral hit
  blastx_progress.json     what the curation page shows while it runs
  cache/<md5>.json         hits of one sequence for one program/database/E-value

Usage:
  python -m evetool.discovery.ncbi_blastx -o results --fasta segments.fasta
  python -m evetool.discovery.ncbi_blastx -o results --region NC_074709.1:9517461-9517811 --species Lutzomyia_longipalpis
  python -m evetool.discovery.ncbi_blastx -o results --resume            # only polls the batches already sent
"""
import os
import re
import csv
import json
import time
import hashlib
import argparse
import subprocess
import urllib.parse
import urllib.request

URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
MAX_SEQS = 20
MAX_NT = 10000
SOBREPOSICAO = 300          # overlap between windows of a long segment
MAX_NA_FILA = 3             # batches waiting at NCBI at the same time
INTERVALO_ENVIO = 10        # s between submissions
INTERVALO_CONSULTA = 60     # s between status checks of one RID
LIMITE_RID = 6 * 3600       # an RID still waiting after max(this, 2 x NCBI estimate) is resubmitted
TETO_RID = 30 * 3600        # ... but never later than this (NCBI keeps an RID about 36 h)
TENTATIVAS = 4
SUBPASTA = os.path.join("13_curadoria", "blastx_ncbi")

RX_VIRAL = re.compile(r"(?<!anti)vir(us|al|idae|ales|inae|icetes|icota)|phage|riboviria", re.I)


def log(msg):
    print(f"[blastx-ncbi] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Segments: ids, sequences and the windows that are actually sent
# ---------------------------------------------------------------------------
def _md5(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def ler_fasta(caminho: str) -> list:
    seqs, nome, partes = [], None, []
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if nome:
                    seqs.append((nome, "".join(partes)))
                nome, partes = linha[1:].strip(), []
            else:
                partes.append(linha.strip())
    if nome:
        seqs.append((nome, "".join(partes)))
    return [(n, s.replace("*", "").upper()) for n, s in seqs if s]


def extrair_regiao(genoma: str, contig: str, ini: int, fim: int) -> str:
    """1-based inclusive region of the genome (samtools faidx)."""
    p = subprocess.run(["samtools", "faidx", genoma, f"{contig}:{ini}-{fim}"],
                       capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        raise RuntimeError(f"samtools faidx {contig}:{ini}-{fim}: {p.stderr.strip()[:300]}")
    return "".join(l.strip() for l in p.stdout.splitlines() if not l.startswith(">")).upper()


def janelas(segmento: dict) -> list:
    """Pieces of at most MAX_NT nt, each with its offset in the segment."""
    seq, n = segmento["seq"], len(segmento["seq"])
    if n <= MAX_NT:
        return [{"wid": segmento["id"], "seg": segmento["id"], "offset": 0, "seq": seq}]
    passo, pecas, ini = MAX_NT - SOBREPOSICAO, [], 0
    while True:
        pecas.append({"wid": f"{segmento['id']}__w{len(pecas) + 1}", "seg": segmento["id"],
                      "offset": ini, "seq": seq[ini:ini + MAX_NT]})
        if ini + MAX_NT >= n:
            return pecas
        ini += passo


def montar_lotes(pecas: list) -> list:
    lotes, atual, nt = [], [], 0
    for p in pecas:
        if atual and (len(atual) >= MAX_SEQS or nt + len(p["seq"]) > MAX_NT):
            lotes.append(atual)
            atual, nt = [], 0
        atual.append(p)
        nt += len(p["seq"])
    if atual:
        lotes.append(atual)
    return lotes


# ---------------------------------------------------------------------------
# NCBI Common URL API
# ---------------------------------------------------------------------------
def _http(params: dict, metodo: str = "GET", timeout: int = 120) -> bytes:
    dados = urllib.parse.urlencode(params).encode()
    ultimo = None
    for k in range(TENTATIVAS):
        try:
            if metodo == "POST":
                req = urllib.request.Request(URL, dados, headers={"User-Agent": "EVE-pipeline/1.0 (blastx_ncbi)"})
            else:
                req = urllib.request.Request(f"{URL}?{dados.decode()}",
                                             headers={"User-Agent": "EVE-pipeline/1.0 (blastx_ncbi)"})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as erro:
            ultimo = erro
            time.sleep(15 * (k + 1))
    raise RuntimeError(f"NCBI BLAST unreachable: {ultimo}")


def enviar(fasta: str, cfg: dict) -> tuple:
    params = {"CMD": "Put", "PROGRAM": "blastx", "DATABASE": cfg["database"], "QUERY": fasta,
              "EXPECT": str(cfg["evalue"]), "HITLIST_SIZE": str(cfg["max_hits"]),
              "TOOL": "eve_pipeline_curation"}
    if cfg.get("email"):
        params["EMAIL"] = cfg["email"]
    texto = _http(params, "POST", timeout=180).decode("utf-8", "replace")
    rid = re.search(r"RID = (\S+)", texto)
    if not rid:
        raise RuntimeError("NCBI did not return an RID: " + re.sub(r"<[^>]+>", " ", texto)[:300])
    rtoe = re.search(r"RTOE = (\d+)", texto)
    return rid.group(1), int(rtoe.group(1)) if rtoe else 60


def situacao(rid: str) -> str:
    texto = _http({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}).decode("utf-8", "replace")
    m = re.search(r"Status=(\w+)", texto)
    return m.group(1) if m else "UNKNOWN"


def baixar(rid: str) -> dict:
    corpo = _http({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2_S"}, timeout=300)
    return json.loads(corpo.decode("utf-8", "replace"))


def hits_do_json(dados: dict) -> dict:
    """{query window id: [hsp dict, ...]} from a JSON2_S report."""
    por_query = {}
    for bloco in dados.get("BlastOutput2", []):
        busca = ((bloco.get("report") or {}).get("results") or {}).get("search") or {}
        titulo = busca.get("query_title") or busca.get("query_id") or ""
        wid = titulo.split()[0] if titulo else ""
        lista = por_query.setdefault(wid, [])
        for h in busca.get("hits", []):
            desc = (h.get("description") or [{}])[0]
            for hsp in h.get("hsps", []):
                lista.append({
                    "accession": desc.get("accession", ""), "title": desc.get("title", ""),
                    "sciname": desc.get("sciname", ""), "taxid": desc.get("taxid", ""),
                    "subject_len": h.get("len", ""), "pident": round(100.0 * hsp.get("identity", 0) /
                                                                   max(1, hsp.get("align_len", 1)), 1),
                    "align_len": hsp.get("align_len", 0), "qstart": hsp.get("query_from", 0),
                    "qend": hsp.get("query_to", 0), "frame": hsp.get("query_frame", ""),
                    "sstart": hsp.get("hit_from", 0), "send": hsp.get("hit_to", 0),
                    "evalue": hsp.get("evalue", ""), "bitscore": hsp.get("bit_score", 0)})
    return por_query


# ---------------------------------------------------------------------------
# State on disk
# ---------------------------------------------------------------------------
def _gravar_json(caminho: str, dados):
    with open(caminho + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(dados, fh, ensure_ascii=False, indent=1)
    os.replace(caminho + ".tmp", caminho)


def _ler_json(caminho: str, padrao):
    try:
        with open(caminho, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return padrao


def _chave_cache(seq: str, cfg: dict) -> str:
    return _md5(f"{seq}|blastx|{cfg['database']}|{cfg['evalue']}|{cfg['max_hits']}")


class Execucao:
    def __init__(self, outdir: str):
        self.dir = os.path.join(outdir, SUBPASTA)
        self.cache = os.path.join(self.dir, "cache")
        os.makedirs(self.cache, exist_ok=True)
        self.cam_jobs = os.path.join(self.dir, "blastx_jobs.json")
        self.cam_prog = os.path.join(self.dir, "blastx_progress.json")
        self.estado = _ler_json(self.cam_jobs, {"config": {}, "segments": {}, "batches": []})

    def salvar(self, msg: str = None):
        _gravar_json(self.cam_jobs, self.estado)
        lotes = self.estado["batches"]
        cont = {}
        for b in lotes:
            cont[b["status"]] = cont.get(b["status"], 0) + 1
        prog = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "message": msg or "",
                "segments": len(self.estado["segments"]), "batches": len(lotes), "per_status": cont,
                "done": all(b["status"] in ("done", "failed") for b in lotes),
                "waiting": [{"batch": b["n"], "rid": b.get("rid"), "since": b.get("submitted"),
                             "estimate_min": round(b.get("rtoe", 0) / 60),
                             "sequences": len(b["windows"])} for b in lotes if b["status"] == "waiting"]}
        _gravar_json(self.cam_prog, prog)
        if msg:
            log(msg)

    def em_cache(self, seq: str, cfg: dict):
        return _ler_json(os.path.join(self.cache, _chave_cache(seq, cfg) + ".json"), None)

    def guardar_cache(self, seq: str, cfg: dict, hits: list):
        _gravar_json(os.path.join(self.cache, _chave_cache(seq, cfg) + ".json"), hits)


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def executar(outdir: str, segmentos: list = None, database: str = "ClusteredNR", evalue: float = 1e-3,
             max_hits: int = 20, email: str = None, deve_parar=lambda: False) -> dict:
    """
    segmentos: [{"id", "seq", "origin", "species", "contig", "start", "end", "strand"}].
    Without segmentos the batches already on disk are resumed.
    """
    ex = Execucao(outdir)
    cfg = {"database": database, "evalue": evalue, "max_hits": max_hits, "email": email}

    if segmentos:
        # A new request replaces the plan, but keeps the RIDs of batches still useful
        antigos = {tuple(sorted(w["wid"] for w in b["windows"])): b for b in ex.estado["batches"]
                   if b["status"] == "waiting" and ex.estado.get("config", {}).get("database") == database
                   and ex.estado.get("config", {}).get("evalue") == evalue}
        ex.estado = {"config": cfg, "segments": {}, "batches": [],
                     "created": time.strftime("%Y-%m-%d %H:%M:%S")}
        pecas = []
        for s in segmentos:
            s = dict(s)
            s["seq"] = re.sub(r"[^A-Za-z]", "", s["seq"]).upper()
            if not s["seq"]:
                continue
            s["length"] = len(s["seq"])
            ex.estado["segments"][s["id"]] = {k: v for k, v in s.items() if k != "seq"}
            for p in janelas(s):
                if ex.em_cache(p["seq"], cfg) is None:
                    pecas.append(p)
                else:
                    ex.estado["segments"][s["id"]].setdefault("cached_windows", 0)
                    ex.estado["segments"][s["id"]]["cached_windows"] += 1
                ex.estado["segments"][s["id"]].setdefault("windows", []).append(
                    {"wid": p["wid"], "offset": p["offset"], "len": len(p["seq"]),
                     "cache": _chave_cache(p["seq"], cfg)})
        for i, lote in enumerate(montar_lotes(pecas), 1):
            chave = tuple(sorted(p["wid"] for p in lote))
            velho = antigos.get(chave)
            ex.estado["batches"].append({
                "n": i, "status": "waiting" if velho else "pending", "rid": velho and velho["rid"],
                "submitted": velho and velho["submitted"], "last_check": 0, "tries": 0,
                "windows": [{"wid": p["wid"], "seq": p["seq"]} for p in lote]})
        n_cache = sum(len(v.get("windows", [])) for v in ex.estado["segments"].values()) - len(pecas)
        ex.salvar(f"{len(segmentos)} segment(s): {len(pecas)} piece(s) to search in "
                  f"{len(ex.estado['batches'])} batch(es), {n_cache} already in the cache")
    else:
        cfg = ex.estado.get("config") or cfg
        if not ex.estado["batches"] and not ex.estado["segments"]:
            raise RuntimeError("Nothing to resume: no blastx_jobs.json in " + ex.dir)
        ex.salvar("resuming the batches on disk")

    ultimo_envio = 0
    while not deve_parar():
        lotes = ex.estado["batches"]
        if all(b["status"] in ("done", "failed") for b in lotes):
            break
        agora = time.time()
        na_fila = [b for b in lotes if b["status"] == "waiting"]

        # Submit while there is room in the queue
        pendentes = [b for b in lotes if b["status"] == "pending"]
        if pendentes and len(na_fila) < MAX_NA_FILA and agora - ultimo_envio >= INTERVALO_ENVIO:
            b = pendentes[0]
            fasta = "".join(f">{w['wid']}\n{w['seq']}\n" for w in b["windows"])
            try:
                b["rid"], rtoe = enviar(fasta, cfg)
                b.update(status="waiting", submitted=time.time(), rtoe=rtoe,
                         last_check=time.time() + min(rtoe, 600) - INTERVALO_CONSULTA)
                ex.salvar(f"batch {b['n']} sent ({len(b['windows'])} seq, "
                          f"{sum(len(w['seq']) for w in b['windows'])} nt): RID {b['rid']}, "
                          f"NCBI estimate {rtoe / 60:.0f} min")
            except Exception as erro:
                b["tries"] += 1
                if b["tries"] >= TENTATIVAS:
                    b["status"], b["error"] = "failed", str(erro)
                ex.salvar(f"batch {b['n']}: submission failed ({erro})")
            ultimo_envio = time.time()
            continue

        # Check the batches in the queue, each at most once a minute
        for b in na_fila:
            if deve_parar() or time.time() - b["last_check"] < INTERVALO_CONSULTA:
                continue
            b["last_check"] = time.time()
            try:
                st = situacao(b["rid"])
            except Exception as erro:
                ex.salvar(f"batch {b['n']}: status check failed ({erro}); will try again")
                continue
            espera = time.time() - (b["submitted"] or time.time())
            if st == "READY":
                try:
                    por_query = hits_do_json(baixar(b["rid"]))
                except Exception as erro:
                    ex.salvar(f"batch {b['n']}: download failed ({erro}); will try again")
                    continue
                for w in b["windows"]:
                    ex.guardar_cache(w["seq"], cfg, por_query.get(w["wid"], []))
                b["status"], b["finished"] = "done", time.time()
                ex.salvar(f"batch {b['n']} done after {espera / 60:.0f} min")
            elif st == "WAITING" and espera < min(TETO_RID, max(LIMITE_RID, 2 * b.get("rtoe", 0))):
                continue
            else:
                motivo = {"WAITING": "too long in the queue", "FAILED": "failed at NCBI",
                          "UNKNOWN": "expired or unknown at NCBI"}.get(st, st)
                _reenviar(ex, b, motivo)
        ex.salvar()
        time.sleep(5)

    resumo = consolidar(outdir)
    ex.salvar("stopped; resume to continue" if deve_parar() and not resumo["done"] else
              f"finished: {resumo['with_hits']} of {resumo['segments']} segment(s) with hits, "
              f"{resumo['viral_best']} with a viral best hit")
    return resumo


def _reenviar(ex: Execucao, b: dict, motivo: str):
    """Put a lost batch back in line; after a failure it goes in two halves."""
    b["tries"] += 1
    if b["tries"] >= TENTATIVAS:
        b["status"], b["error"] = "failed", motivo
        ex.salvar(f"batch {b['n']}: {motivo}; giving up after {b['tries']} tries")
        return
    if motivo == "failed at NCBI" and len(b["windows"]) > 1:
        meio = len(b["windows"]) // 2
        novo = {"n": max(x["n"] for x in ex.estado["batches"]) + 1, "status": "pending", "rid": None,
                "submitted": None, "last_check": 0, "tries": b["tries"], "windows": b["windows"][meio:]}
        b["windows"] = b["windows"][:meio]
        ex.estado["batches"].append(novo)
    b.update(status="pending", rid=None, submitted=None)
    ex.salvar(f"batch {b['n']}: {motivo}; sent again")


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def consolidar(outdir: str) -> dict:
    ex = Execucao(outdir)
    segs = ex.estado.get("segments", {})
    linhas, resumo = [], []
    for sid, s in segs.items():
        hits, buscado = [], True
        for w in s.get("windows", []):
            h = _ler_json(os.path.join(ex.cache, w["cache"] + ".json"), None)
            if h is None:
                buscado = False
                continue
            for x in h:
                x = dict(x)
                x["qstart"] = int(x["qstart"]) + w["offset"]
                x["qend"] = int(x["qend"]) + w["offset"]
                hits.append(x)
        # Overlapping windows report the same HSP twice
        vistos, unicos = set(), []
        for x in sorted(hits, key=lambda x: -float(x["bitscore"] or 0)):
            k = (x["accession"], min(x["qstart"], x["qend"]) // 30, x["sstart"] // 10)
            if k not in vistos:
                vistos.add(k)
                unicos.append(x)
        for x in unicos:
            x["viral"] = "yes" if RX_VIRAL.search(f"{x['sciname']} {x['title']}") else "no"
            gc = s.get("start")
            if gc not in (None, "") and s.get("contig"):
                ini = int(gc)
                x["genome_start"], x["genome_end"] = ini + min(x["qstart"], x["qend"]) - 1, \
                    ini + max(x["qstart"], x["qend"]) - 1
            linhas.append({"segment": sid, "contig": s.get("contig", ""), **x})
        melhor = unicos[0] if unicos else None
        viral = next((x for x in unicos if x["viral"] == "yes"), None)
        celular = next((x for x in unicos if x["viral"] == "no"), None)
        resumo.append({
            "segment": sid, "origin": s.get("origin", ""), "species": s.get("species", ""),
            "contig": s.get("contig", ""), "start": s.get("start", ""), "end": s.get("end", ""),
            "length": s.get("length", ""), "searched": "yes" if buscado else "pending",
            "n_hits": len({x["accession"] for x in unicos}),
            "best_hit": melhor["accession"] if melhor else "", "best_title": melhor["title"] if melhor else "",
            "best_taxon": melhor["sciname"] if melhor else "",
            "best_evalue": melhor["evalue"] if melhor else "", "best_bitscore": melhor["bitscore"] if melhor else "",
            "best_is_viral": melhor["viral"] if melhor else "",
            "best_viral": viral["accession"] if viral else "", "best_viral_title": viral["title"] if viral else "",
            "best_viral_evalue": viral["evalue"] if viral else "",
            "best_viral_bitscore": viral["bitscore"] if viral else "",
            "best_nonviral": celular["accession"] if celular else "",
            "best_nonviral_title": celular["title"] if celular else "",
            "best_nonviral_bitscore": celular["bitscore"] if celular else "",
        })
    cols_h = ["segment", "contig", "accession", "title", "sciname", "taxid", "viral", "pident", "align_len", "qstart",
              "qend", "frame", "genome_start", "genome_end", "sstart", "send", "subject_len", "evalue", "bitscore"]
    _gravar_tsv(os.path.join(ex.dir, "blastx_hits.tsv"), cols_h, linhas)
    cols_r = list(resumo[0]) if resumo else ["segment"]
    _gravar_tsv(os.path.join(ex.dir, "blastx_summary.tsv"), cols_r, resumo)
    _gravar_tsv(os.path.join(ex.dir, "blastx_queries.tsv"),
                ["segment", "origin", "species", "contig", "start", "end", "length"], resumo)
    return {"segments": len(resumo), "with_hits": sum(1 for r in resumo if r["n_hits"]),
            "viral_best": sum(1 for r in resumo if r["best_is_viral"] == "yes"),
            "pending": sum(1 for r in resumo if r["searched"] != "yes"),
            "done": all(r["searched"] == "yes" for r in resumo),
            "summary": os.path.join(ex.dir, "blastx_summary.tsv"),
            "hits": os.path.join(ex.dir, "blastx_hits.tsv")}


def _gravar_tsv(caminho: str, colunas: list, linhas: list):
    with open(caminho + ".tmp", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=colunas, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(linhas)
    os.replace(caminho + ".tmp", caminho)


def progresso(outdir: str) -> dict:
    return _ler_json(os.path.join(outdir, SUBPASTA, "blastx_progress.json"), {})


def hits_do_segmento(outdir: str, contig: str, ini: int, fim: int, query: str = None) -> list:
    """BLASTX hits of segments that overlap a region (or that are this EVE)."""
    cam = os.path.join(outdir, SUBPASTA, "blastx_hits.tsv")
    if not os.path.exists(cam):
        return []
    saida = []
    with open(cam, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mesmo = query and r["segment"] == query
            sobrepoe = (r.get("genome_start") and r.get("contig") == contig and
                        int(r["genome_start"]) <= fim and int(r["genome_end"]) >= ini)
            if mesmo or sobrepoe:
                saida.append(r)
    saida.sort(key=lambda r: -float(r["bitscore"] or 0))
    return saida


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--fasta", help="Nucleotide FASTA with the segments")
    ap.add_argument("--region", action="append", default=[], help="contig:start-end (1-based), repeatable")
    ap.add_argument("--species", default="", help="Species folder of --region (multi-species runs)")
    ap.add_argument("--genome", help="Genome FASTA for --region (default: from pipeline_info.json)")
    ap.add_argument("--resume", action="store_true", help="Only resume the batches already on disk")
    ap.add_argument("--database", default="ClusteredNR", help="ClusteredNR (default), nr, refseq_protein, swissprot")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--max-hits", type=int, default=20)
    ap.add_argument("--email", default=None, help="Contact e-mail sent to NCBI with the requests")
    a = ap.parse_args()

    segmentos = []
    if a.fasta:
        segmentos += [{"id": n.split()[0], "seq": s, "origin": "fasta"} for n, s in ler_fasta(a.fasta)]
    if a.region:
        genoma = a.genome
        if not genoma:
            from evetool.curation.curation import pastas_especies, _info_pipeline
            pastas = pastas_especies(a.outdir)
            genoma = _info_pipeline(pastas.get(a.species) or next(iter(pastas.values()))).get("genome")
        for r in a.region:
            contig, _, faixa = r.rpartition(":")
            ini, fim = (int(x.replace(",", "")) for x in faixa.split("-"))
            segmentos.append({"id": f"{a.species + '|' if a.species else ''}{contig}:{ini}-{fim}",
                              "seq": extrair_regiao(genoma, contig, ini, fim), "origin": "region",
                              "species": a.species, "contig": contig, "start": ini, "end": fim})
    if not segmentos and not a.resume:
        ap.error("give --fasta, --region or --resume")
    r = executar(a.outdir, segmentos or None, database=a.database, evalue=a.evalue,
                 max_hits=a.max_hits, email=a.email)
    print(f"\n{r['segments']} segment(s), {r['with_hits']} with hits, {r['viral_best']} with a viral best hit"
          f"{', ' + str(r['pending']) + ' still pending' if r['pending'] else ''}\n  {r['summary']}")


if __name__ == "__main__":
    main()
