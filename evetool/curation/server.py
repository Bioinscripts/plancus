#!/usr/bin/env python3
"""
Local server for the EVE page: curation and overview charts on one page.

  GET  /, /index.html, /curadoria.html -> the page, always read from evetool/curation/web/curation.html,
                           so an edit to its code (colours, layout) shows on reload without
                           touching the data
  GET  /v2, /curation_v2.html -> the same tool with the report layout (web/curation_v2.html): per-EVE
                           cards and the read heatmaps; the page above stays as it is
  GET  /api/rnaseq-matrix -> ?sp=: every approved locus x RNA-seq sample (exonic reads, depth, strands)
  GET  /api/srna-matrix -> ?sp=: every approved locus x small RNA library (sRNAProfiler features and
                           the piRNA-sized reads per strand: 5' U, A10, ping-pong)
  GET  /api/eve-spacing -> gap from every approved locus to the nearest one of its genome, per species,
                           against random placement on the assembly (comparative/eve_spacing.py)
  GET  /api/candidates  -> classified EVEs, curation state, review flags and sequence metrics
  GET  /api/state       -> version stamp of the data on disk (the page polls it to refresh itself)
  GET  /api/detail      -> inspect one EVE: position, neighbours, GC and alignments
  GET  /api/protein-map -> where each EVE falls on its viral protein (?force=1 realigns)
  GET  /api/rnaseq      -> binned RNA-seq depth of every sample over ?sp=&c=&s=&e= (0-based) [&focus=a-b,...]
  GET  /api/rnaseq-samples -> RNA-seq BAMs known per species
  POST /api/rnaseq-bam  -> add ({species, bam, library}) or remove ({op: "remove", species, name}) a BAM
  POST /api/rnaseq-map  -> map long / other RNA reads with minimap2, as ViralQuest ({species, r1, r2, read_type,
                           library, reference, flank}); runs in the background
  GET  /api/srna        -> small RNA reads over ?sp=&c=&s=&e= [&focus=a-b:strand,...]: depth by size class and
                           the sRNAProfiler features of each span
  GET  /api/srna-libraries -> small RNA libraries per species
  POST /api/srna-library -> map trimmed small RNA reads with Bowtie1 ({species, reads, name, ncrna_index}),
                           register a BAM ({species, bam}) or forget one ({op: "remove", species, name})
  GET  /api/dna         -> genomic DNA reads over every locus: depth, breadth, junction and empty-site reads per
                           library (raw measurements; presence / absence is called on the page)
  POST /api/dna-library -> map DNA reads with minimap2 on the whole genome, keeping the EVE windows ({species, r1,
                           r2, read_type, role: assembly|population, sample_type: individual|pool, group, name}),
                           register a BAM (r1 ending in .bam, or {bam}), edit ({op: "edit", ...}) or forget one
  GET  /api/map-jobs    -> the read-mapping jobs of this session and the step each one is in
  GET  /api/primers     -> Primer3 pairs for an EVE (?q=<query>) or a region (?sp=&c=&s=&e=, 0-based):
                           internal and junction (EVE + host flank)
  POST /api/export      -> save a chart as SVG and, when possible, convert it to PDF
  GET  /api/circulating -> circulating-virus sets and the last composition analysis (codon_composition.py)
  POST /api/regions     -> multi-FASTA of several genome regions (copy / download FASTA)
  POST /api/curation    -> save the selection and run the requested action:
                           save | fasta_nt | fasta_pep | annotate_hmm | cluster |
                           multispecies | context | metrics | repeats | cdsearch | pfam |
                           tblastn_scan | circ_add | circ_remove | circ_run | circ_gff |
                           primers_tsv | srna_tsv | dna_tsv | host_genes

After every change the data snapshot (12_dashboard_web/dados/snapshot.js) is rewritten, so the
page also opens with its data straight from disk. The old review page (revisao.html) is gone:
its ORFs are flagged on the curation page and decided there like any other row.

Requests are served in parallel, so the page keeps answering while a long action runs.
Listens on 127.0.0.1 only. The results directory is --outdir (default: the served directory,
or its parent when serving 12_dashboard_web).
"""
import os
import sys
import json
import argparse
import http.server
import socketserver
import threading
import traceback

PAGINAS = ("/", "/index.html", "/curadoria.html")
PAGINAS_V2 = ("/v2", "/v2/", "/curation_v2.html")
_estado = {"outdir": ".", "rodando": False, "acao": None, "ultimo_resumo": None, "base_code": "", "threads": 4,
          "blastx": None}      # online BLASTX runs in its own thread: it can take hours at NCBI
_trava = threading.Lock()


def _situacao_blastx() -> dict:
    """Progress of the online BLASTX plus its summary table, for the page."""
    import csv
    from urllib.parse import quote
    bx = _importar("evetool.discovery.ncbi_blastx")
    pasta = os.path.join(_estado["outdir"], bx.SUBPASTA)
    t = _estado.get("blastx")
    res = {"running": bool(t and t.is_alive()), "error": _estado.get("blastx_erro"),
           "progress": bx.progresso(_estado["outdir"]), "rows": []}
    cam = os.path.join(pasta, "blastx_summary.tsv")
    if os.path.exists(cam):
        with open(cam, encoding="utf-8", errors="replace") as fh:
            res["rows"] = list(csv.DictReader(fh, delimiter="\t"))
        for nome in ("blastx_summary.tsv", "blastx_hits.tsv"):
            res[nome.split(".")[0].split("_")[1]] = "api/download?f=" + quote(os.path.abspath(os.path.join(pasta, nome)))
    return res


def _importar(nome: str):
    """Import a pipeline module (evetool.*), lazily, so the server starts fast."""
    import importlib
    return importlib.import_module(nome)


def _snapshot():
    """Rewrite the data snapshot; never lets a failure here break the request."""
    try:
        _importar("evetool.curation.curation").gravar_snapshot(_estado["outdir"])
    except Exception as erro:
        print(f"[server] snapshot not written: {erro}", file=sys.stderr)


def _urls_arquivos(caminhos):
    """Browser URLs: a relative path when inside the served folder, otherwise /api/download."""
    from urllib.parse import quote
    base = os.path.abspath(os.getcwd())
    urls = []
    for c in caminhos or []:
        a = os.path.abspath(c)
        if a.startswith(base + os.sep):
            urls.append(os.path.relpath(a, base).replace(os.sep, "/"))
        else:
            urls.append("api/download?f=" + quote(a))
    return urls


class Handler(http.server.SimpleHTTPRequestHandler):
    def _json(self, codigo: int, payload: dict):
        corpo = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(codigo)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(corpo)))
        self.end_headers()
        self.wfile.write(corpo)

    TIPOS_INLINE = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
                    ".svg": "image/svg+xml", ".pdf": "application/pdf"}
    SEM_INLINE = {".pdf"}   # PDFs go out as downloads, so we do not depend on the viewer

    def _arquivo(self, caminho: str, nome_download: str):
        """Serve a file from disk: images and PDFs inline, everything else as a download."""
        try:
            dados = open(caminho, "rb").read()
        except OSError as erro:
            return self._json(500, {"error": f"could not read {caminho}: {erro}"})
        ext = os.path.splitext(caminho)[1].lower()
        tipo = self.TIPOS_INLINE.get(ext)
        if ext in self.SEM_INLINE:
            tipo = None
        self.send_response(200)
        self.send_header("Content-Type", tipo or "application/octet-stream")
        if tipo:
            self.send_header("Content-Disposition", f'inline; filename="{nome_download}"')
        else:
            self.send_header("Content-Disposition", f'attachment; filename="{nome_download}"')
        self.send_header("Content-Length", str(len(dados)))
        self.end_headers()
        self.wfile.write(dados)

    def _pagina(self, arquivo: str = "curation.html"):
        caminho = os.path.join(_estado["base_code"], "web", arquivo)
        corpo = open(caminho, "rb").read()
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(corpo)))
        self.end_headers()
        self.wfile.write(corpo)

    def end_headers(self):
        # Data and the snapshot must never come from the browser cache
        if self.path.startswith("/api/") or "snapshot.js" in self.path:
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self):
        rota = self.path.split("?")[0]
        if rota in PAGINAS:
            return self._pagina()
        if rota in PAGINAS_V2:
            if rota != "/curation_v2.html":     # relative links (api/..., dados/...) must resolve at the root
                self.send_response(302)
                self.send_header("Location", "/curation_v2.html" + (self.path[len(rota):] or ""))
                self.end_headers()
                return
            return self._pagina("curation_v2.html")
        if rota == "/api/eve-spacing":
            from urllib.parse import parse_qs, urlparse
            juntar = (parse_qs(urlparse(self.path).query).get("merge") or ["0"])[0]
            try:
                return self._json(200, _importar("evetool.comparative.eve_spacing").distancias(
                    _estado["outdir"], juntar=int(juntar) if juntar.isdigit() else 0))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota in ("/api/rnaseq-matrix", "/api/srna-matrix"):
            from urllib.parse import parse_qs, urlparse
            sp = (parse_qs(urlparse(self.path).query).get("sp") or [""])[0]
            modulo = "evetool.annotation.rnaseq_tracks" if rota == "/api/rnaseq-matrix" else "evetool.annotation.srna"
            try:
                return self._json(200, _importar(modulo).matriz(_estado["outdir"], sp))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/revisao.html":
            self.send_response(302)
            self.send_header("Location", "curation.html#review")
            self.end_headers()
            return
        if rota == "/api/state":
            try:
                return self._json(200, {"stamp": _importar("evetool.curation.curation").estado_dados(_estado["outdir"]),
                                        "running": _estado["rodando"], "action": _estado["acao"]})
            except Exception as erro:
                return self._json(500, {"error": str(erro)})
        if rota == "/api/candidates":
            try:
                dados = _importar("evetool.curation.curation").listar_candidatas(_estado["outdir"])
                dados["running"] = _estado["rodando"]
                return self._json(200, dados)
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/detail":
            from urllib.parse import parse_qs, urlparse
            q = (parse_qs(urlparse(self.path).query).get("q") or [""])[0]
            if not q:
                return self._json(400, {"error": "provide ?q=<query>"})
            try:
                cur = _importar("evetool.curation.curation")
                detalhe = cur.detalhar_eve(_estado["outdir"], q)
                # verdict of this candidate (the inspection opens in its own tab, without the table)
                estado = cur.listar_candidatas(_estado["outdir"])
                todos = estado["items"]
                it = next((i for i in todos if i["query"] == q), None)
                if it and isinstance(detalhe, dict):
                    # the row of the candidate table, as saved, and the lists its two selects offer
                    detalhe["record"] = it
                    detalhe["lists"] = {k: estado.get(k) or [] for k in ("classes", "custom_classes", "families", "custom_families")}
                    detalhe["verdict"] = {k: it.get(k) for k in ("verdict", "detection", "viral_scope", "auto_reject", "read_support", "simple_seq",
                                                                 "verdict_for", "verdict_against")}
                    detalhe.update(query=q, locus=it["locus"], ctx_cluster=it["ctx_cluster"],
                                   saved_state={"removed": it["removed_flag"], "qualified": it["qualified"]})
                    # neighbours: one per locus (its representative frame), with what the page needs to
                    # show them as a cluster and open them; other frames of this EVE's own locus are left out
                    por_orf = {i["orf"]: i for i in todos if i["species"] == it["species"]}
                    vistos, lista = set(), []
                    for v in detalhe.get("neighbours") or []:
                        n = por_orf.get(v["orf"])
                        if not n or n["locus"] == it["locus"]:
                            continue
                        rep = next((x for x in por_orf.values() if x["locus"] == n["locus"] and x["locus_representative"]), n)
                        if rep["query"] in vistos:
                            continue
                        vistos.add(rep["query"])
                        ini = min(x["start"] for x in por_orf.values() if x["locus"] == n["locus"])
                        fim = max(x["end"] for x in por_orf.values() if x["locus"] == n["locus"])
                        dist = 0 if (ini < detalhe["end"] and fim > detalhe["start"]) else (
                            ini - detalhe["end"] if ini >= detalhe["end"] else detalhe["start"] - fim)
                        lista.append({"orf": rep["orf"], "query": rep["query"], "start": ini, "end": fim,
                                      "strand": rep["strand"], "distance": dist, "family": rep["family"],
                                      "protein_class": rep["protein_class"], "verdict": rep["verdict"],
                                      "ctx_cluster": rep["ctx_cluster"], "same_cluster": bool(rep["ctx_cluster"])
                                      and rep["ctx_cluster"] == it["ctx_cluster"],
                                      "removed": rep["removed_flag"], "qualified": rep["qualified"]})
                    detalhe["neighbours"] = sorted(lista, key=lambda v: v["start"])
                    detalhe["n_neighbours"] = len(lista)
                return self._json(200, detalhe)
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/protein-map":
            from urllib.parse import parse_qs, urlparse
            forcar = (parse_qs(urlparse(self.path).query).get("force") or [""])[0] == "1"
            try:
                cur = _importar("evetool.curation.curation")
                res = cur.mapa_de_proteinas(_estado["outdir"], forcar=forcar)
                try:
                    cur.gravar_snapshot(_estado["outdir"])   # the offline copy gets the map too
                except Exception:
                    pass
                return self._json(200, res)
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/region":
            from urllib.parse import parse_qs, urlparse
            q = {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}
            try:
                return self._json(200, _importar("evetool.curation.curation").sequencia_regiao(
                    _estado["outdir"], q.get("sp", ""), q["c"], int(q["s"]), int(q["e"]), q.get("n", "")))
            except Exception as erro:
                return self._json(500, {"error": str(erro)})
        if rota == "/api/rnaseq-samples":
            try:
                return self._json(200, _importar("evetool.annotation.rnaseq_tracks").resumo_amostras(_estado["outdir"]))
            except Exception as erro:
                return self._json(500, {"error": str(erro)})
        if rota == "/api/rnaseq":
            from urllib.parse import parse_qs, urlparse
            q = {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}
            try:
                foco = [tuple(int(x) for x in f.split("-")) for f in q.get("focus", "").split(",") if "-" in f]
                return self._json(200, _importar("evetool.annotation.rnaseq_tracks").cobertura(
                    _estado["outdir"], q.get("sp", ""), q["c"], int(q["s"]), int(q["e"]),
                    int(q.get("bins") or 600), foco, int(q.get("mapq") or 0)))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/srna-libraries":
            try:
                return self._json(200, _importar("evetool.annotation.srna").resumo(_estado["outdir"]))
            except Exception as erro:
                return self._json(500, {"error": str(erro)})
        if rota == "/api/srna":
            from urllib.parse import parse_qs, urlparse
            q = {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}
            try:
                foco = []
                for f in q.get("focus", "").split(","):
                    faixa, _, fita = f.partition(":")
                    if "-" in faixa:
                        a, b = faixa.split("-")
                        foco.append([int(a), int(b), fita or "."])
                return self._json(200, _importar("evetool.annotation.srna").perfil(
                    _estado["outdir"], q.get("sp", ""), q["c"], int(q["s"]), int(q["e"]),
                    int(q.get("bins") or 600), foco))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/dna":
            # genomic DNA reads over every locus: raw measurements, the call is made on the page
            try:
                return self._json(200, _importar("evetool.annotation.dna_reads").medir(_estado["outdir"]))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/map-jobs":
            return self._json(200, {"jobs": _importar("evetool.annotation.read_jobs").situacao()})
        if rota == "/api/primers":
            from urllib.parse import parse_qs, urlparse
            q = {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}
            try:
                pr = _importar("evetool.curation.primers")
                faixa = lambda ch: tuple(int(x) for x in q[ch].split("-")) if q.get(ch) else None
                opcoes = {"flank": q.get("flank"), "n": q.get("pairs"), "internal_size": faixa("isize"),
                          "junction_size": faixa("jsize"),
                          "avoid_masked": q.get("masked", "1") != "0"}
                if q.get("tm"):
                    t = float(q["tm"])
                    opcoes["tm"] = (t - 3.0, t, t + 3.0)
                if q.get("q"):
                    cur = _importar("evetool.curation.curation")
                    it = next((i for i in cur.listar_candidatas(_estado["outdir"])["items"] if i["query"] == q["q"]), None)
                    if not it or not it["end"]:
                        return self._json(404, {"error": f"{q['q']}: no coordinates for this candidate"})
                    alvo = (it["species"], it["contig"], it["start"], it["end"], it["orf"])
                else:
                    alvo = (q.get("sp", ""), q["c"], int(q["s"]), int(q["e"]), q.get("n", ""))
                return self._json(200, pr.para_regiao(_estado["outdir"], *alvo, opcoes=opcoes))
            except ImportError:
                return self._json(500, {"error": "primer3-py is not installed (pip install primer3-py)"})
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/blastx":
            return self._json(200, _situacao_blastx())
        if rota == "/api/circulating":
            try:
                return self._json(200, _importar("evetool.comparative.codon_composition").carregar(_estado["outdir"]))
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/download":
            from urllib.parse import parse_qs, urlparse
            alvo = (parse_qs(urlparse(self.path).query).get("f") or [""])[0]
            permitido = os.path.abspath(os.path.join(_estado["outdir"], "13_curadoria"))
            caminho = os.path.abspath(alvo)
            if not caminho.startswith(permitido + os.sep):
                return self._json(403, {"error": "path outside the curation folder"})
            return self._arquivo(caminho, os.path.basename(caminho))
        return super().do_GET()

    def do_POST(self):
        rota = self.path.split("?")[0]
        if rota == "/api/export":
            try:
                tamanho = int(self.headers.get("Content-Length") or 0)
                pedido = json.loads(self.rfile.read(tamanho) or b"{}")
                res = _importar("evetool.curation.curation").exportar_grafico(
                    _estado["outdir"], pedido.get("svg") or "", pedido.get("name") or "chart")
                from urllib.parse import quote
                return self._json(200, {
                    "pdf_error": res["error"],
                    "svg": "api/download?f=" + quote(os.path.abspath(res["svg"])),
                    "pdf": ("api/download?f=" + quote(os.path.abspath(res["pdf"]))) if res["pdf"] else None,
                    "path_svg": res["svg"], "path_pdf": res["pdf"]})
            except Exception as erro:
                return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        if rota == "/api/rnaseq-bam":
            try:
                tamanho = int(self.headers.get("Content-Length") or 0)
                pedido = json.loads(self.rfile.read(tamanho) or b"{}")
                rt = _importar("evetool.annotation.rnaseq_tracks")
                if pedido.get("op") == "remove":
                    return self._json(200, rt.remover(_estado["outdir"], pedido.get("species", ""), pedido.get("name", "")))
                return self._json(200, rt.adicionar(_estado["outdir"], pedido.get("species", ""), pedido.get("bam", ""),
                                                    pedido.get("library") or "unknown"))
            except Exception as erro:
                return self._json(400, {"error": str(erro)})
        if rota in ("/api/rnaseq-map", "/api/srna-library"):
            try:
                tamanho = int(self.headers.get("Content-Length") or 0)
                pedido = json.loads(self.rfile.read(tamanho) or b"{}")
                jobs = _importar("evetool.annotation.read_jobs")
                sp, threads, outdir = str(pedido.get("species") or ""), _estado["threads"], _estado["outdir"]
                if rota == "/api/rnaseq-map":
                    rt = _importar("evetool.annotation.rnaseq_tracks")
                    r1 = str(pedido.get("r1") or "").strip()
                    if not r1:
                        return self._json(400, {"error": "give the path of the reads (FASTQ)"})
                    job = jobs.iniciar("rna", sp, os.path.basename(r1), lambda passo: rt.mapear(
                        outdir, sp, r1, str(pedido.get("r2") or ""), str(pedido.get("read_type") or "sr"),
                        str(pedido.get("library") or "unknown"), str(pedido.get("reference") or "genome"),
                        int(pedido.get("flank") or 20000), str(pedido.get("name") or ""), threads, passo))
                    return self._json(200, {"job": job})
                sr = _importar("evetool.annotation.srna")
                if pedido.get("op") == "remove":
                    return self._json(200, sr.remover(outdir, sp, str(pedido.get("name") or "")))
                reads, bam = str(pedido.get("reads") or "").strip(), str(pedido.get("bam") or "").strip()
                if not reads and not bam:
                    return self._json(400, {"error": "give the path of the trimmed reads (FASTQ/FASTA) or of a BAM"})
                if reads.lower().endswith(".bam"):
                    reads, bam = "", reads
                if bam:
                    job = jobs.iniciar("srna", sp, os.path.basename(bam), lambda passo: sr.adicionar_bam(
                        outdir, sp, bam, str(pedido.get("name") or ""), passo))
                else:
                    job = jobs.iniciar("srna", sp, os.path.basename(reads), lambda passo: sr.mapear(
                        outdir, sp, reads, str(pedido.get("name") or ""), threads,
                        str(pedido.get("ncrna_index") or ""), passo))
                return self._json(200, {"job": job})
            except Exception as erro:
                return self._json(400, {"error": str(erro)})
        if rota == "/api/dna-library":
            try:
                tamanho = int(self.headers.get("Content-Length") or 0)
                pedido = json.loads(self.rfile.read(tamanho) or b"{}")
                jobs = _importar("evetool.annotation.read_jobs")
                dn = _importar("evetool.annotation.dna_reads")
                sp, threads, outdir = str(pedido.get("species") or ""), _estado["threads"], _estado["outdir"]
                nome = str(pedido.get("name") or "")
                if pedido.get("op") == "remove":
                    return self._json(200, dn.remover(outdir, sp, nome))
                if pedido.get("op") == "edit":
                    return self._json(200, dn.editar(outdir, sp, nome, role=pedido.get("role"),
                                                     sample_type=pedido.get("sample_type"), group=pedido.get("group")))
                tipo = str(pedido.get("sample_type") or "individual")
                papel = str(pedido.get("role") or "population")
                grupo = str(pedido.get("group") or "")
                dn._validar(tipo, papel)
                r1, bam = str(pedido.get("r1") or "").strip(), str(pedido.get("bam") or "").strip()
                if r1.lower().endswith(".bam"):
                    r1, bam = "", r1
                if not r1 and not bam:
                    return self._json(400, {"error": "give the path of the DNA reads (FASTQ) or of a BAM"})
                if bam:
                    job = jobs.iniciar("dna", sp, nome or os.path.basename(bam), lambda passo: dn.adicionar_bam(
                        outdir, sp, bam, nome, tipo, papel, grupo, str(pedido.get("read_type") or "sr"), passo))
                else:
                    job = jobs.iniciar("dna", sp, nome or os.path.basename(r1), lambda passo: dn.mapear(
                        outdir, sp, r1, str(pedido.get("r2") or ""), str(pedido.get("read_type") or "sr"), nome,
                        tipo, papel, grupo, threads, passo))
                return self._json(200, {"job": job})
            except Exception as erro:
                return self._json(400, {"error": str(erro)})
        if rota == "/api/regions":
            try:
                tamanho = int(self.headers.get("Content-Length") or 0)
                pedido = json.loads(self.rfile.read(tamanho) or b"{}")
                return self._json(200, _importar("evetool.curation.curation").sequencias_regioes(
                    _estado["outdir"], pedido.get("regions") or [], pedido.get("flank", 0)))
            except Exception as erro:
                return self._json(500, {"error": str(erro)})
        if rota == "/api/curation":
            return self._curadoria()
        return self._json(404, {"error": "unknown route"})

    def _blastx(self, acao, pedido, selecionados):
        """Start (or resume) the online BLASTX in a thread and answer at once."""
        with _trava:
            if _estado["blastx"] and _estado["blastx"].is_alive():
                return self._json(409, {"error": "an online BLASTX run is already in progress"})
        try:
            cur = _importar("evetool.curation.curation")
            bx = _importar("evetool.discovery.ncbi_blastx")
            segmentos = None
            if acao == "blastx":
                fontes = pedido.get("blastx_sources") or ["eves"]
                segmentos = cur.segmentos_para_blastx(
                    _estado["outdir"], queries=selecionados if "eves" in fontes else None,
                    clusters=pedido.get("clusters") if "clusters" in fontes else None,
                    regioes=pedido.get("regions") or "")
            opcoes = {"database": pedido.get("database") or "ClusteredNR",
                      "evalue": float(pedido.get("evalue") or 1e-3),
                      "max_hits": int(pedido.get("max_hits") or 20)}
        except Exception as erro:
            return self._json(400, {"error": str(erro)})

        def correr():
            try:
                bx.executar(_estado["outdir"], segmentos, **opcoes)
            except Exception as erro:
                print(f"[blastx-ncbi] stopped with an error: {erro}", flush=True)
                _estado["blastx_erro"] = str(erro)

        _estado["blastx_erro"] = None
        t = threading.Thread(target=correr, daemon=True, name="blastx_ncbi")
        _estado["blastx"] = t
        t.start()
        return self._json(200, {"action": acao, "started": True,
                                "segments": len(segmentos) if segmentos else None})

    def _curadoria(self):
        try:
            tamanho = int(self.headers.get("Content-Length") or 0)
            pedido = json.loads(self.rfile.read(tamanho) or b"{}")
        except Exception as erro:
            return self._json(400, {"error": f"invalid JSON: {erro}"})

        acao = str(pedido.get("action") or "save")
        removed = [str(q) for q in pedido.get("removed", [])]
        proteins = pedido.get("proteins")
        classes_novas = pedido.get("custom_classes") or []
        families = pedido.get("families")
        familias_novas = pedido.get("custom_families") or []
        selecionados = [str(q) for q in pedido.get("selected", [])] or None

        if acao in ("blastx", "blastx_resume"):
            return self._blastx(acao, pedido, selecionados)

        with _trava:
            if _estado["rodando"]:
                return self._json(409, {"error": "a run is already in progress"})
            _estado["rodando"] = True
            _estado["acao"] = acao
        try:
            cur = _importar("evetool.curation.curation")
            cur.salvar_selecao(_estado["outdir"], removed,
                               proteins=proteins if isinstance(proteins, dict) else None,
                               classes_novas=classes_novas,
                               families=families if isinstance(families, dict) else None,
                               familias_novas=familias_novas,
                               tblastn_removed=pedido.get("tblastn_removed")
                               if isinstance(pedido.get("tblastn_removed"), list) else None,
                               removed_via=pedido.get("removed_via")
                               if isinstance(pedido.get("removed_via"), dict) else None,
                               qualified=pedido.get("qualified")
                               if isinstance(pedido.get("qualified"), list) else None,
                               host_genes_removed=pedido.get("host_genes_removed")
                               if isinstance(pedido.get("host_genes_removed"), list) else None)
            resposta = {"action": acao, "removed": len(set(removed))}

            if acao in ("fasta_nt", "fasta_pep"):
                res = cur.exportar_fasta(_estado["outdir"], acao.split("_")[1], queries=selecionados)
                resposta.update(res)
                from urllib.parse import quote
                resposta["download"] = "api/download?f=" + quote(os.path.abspath(res["file"]))
            elif acao == "cluster":
                res = cur.rodar_clusterizacao(_estado["outdir"], threads=_estado["threads"],
                                              queries=selecionados)
                resposta.update({"sequences": res["sequences"], "clusters": res["clusters"],
                                 "images": _urls_arquivos(res["images"]),
                                 "images_abs": res["images"], "tables": res["tables"]})
            elif acao == "annotate_hmm":
                resposta.update(cur.rodar_anotacao_hmm(_estado["outdir"], threads=_estado["threads"]))
            elif acao == "multispecies":
                res = cur.rodar_multiespecies(
                    _estado["outdir"], threads=_estado["threads"], queries=selecionados,
                    min_identity=float(pedido.get("min_identity", 30.0)),
                    min_coverage=float(pedido.get("min_coverage", 30.0)),
                    sintenia=bool(pedido.get("synteny")),
                    proteoma_ref=pedido.get("proteoma_ref") or None,
                    share_identity=float(pedido.get("share_identity", 50.0)),
                    share_coverage=float(pedido.get("share_coverage", 50.0)),
                    modo=str(pedido.get("units_mode") or "orfs"), unidades=pedido.get("units"),
                    gap=pedido.get("units_gap"), flanco=int(pedido.get("synteny_flank") or 20000),
                    max_sitios=int(pedido.get("synteny_max_sites") or 3))
                resposta.update({k: res[k] for k in ("clusters", "shared", "supported", "species_count", "tables", "units", "synteny")})
                resposta.update({"images": _urls_arquivos(res["images"]),
                                 "images_abs": res["images"]})
            elif acao == "context":
                res = cur.rodar_contexto(_estado["outdir"], queries=selecionados)
                resposta.update({"warnings": res["warnings"], "images": _urls_arquivos(res["images"]),
                                 "images_abs": res["images"]})
            elif acao == "cdsearch":
                resposta.update(cur.rodar_cdsearch(_estado["outdir"],
                                                   evalue=float(pedido.get("evalue", 0.01)),
                                                   threads=_estado["threads"],
                                                   modo=pedido.get("mode", "local")))
            elif acao == "pfam":
                resposta.update(cur.rodar_pfam(_estado["outdir"], threads=_estado["threads"]))
            elif acao == "tblastn_scan":
                resposta.update(cur.rodar_varredura_tblastn(
                    _estado["outdir"], threads=_estado["threads"],
                    hits_por_orf=max(1, int(pedido.get("scan_hits_per_orf", 1))),
                    evalue=float(pedido.get("scan_evalue", 1e-5)),
                    gap=max(0, int(pedido.get("scan_gap", 1000))), forcar=bool(pedido.get("scan_force"))))
            elif acao == "repeats":
                resposta.update(cur.rodar_repeticoes(_estado["outdir"], threads=_estado["threads"],
                                                     flanco=int(pedido.get("repeats_flank") or 20000),
                                                     queries=selecionados))
            elif acao == "primers_tsv":
                from urllib.parse import quote
                o = pedido.get("primer_options") or {}
                res = _importar("evetool.curation.primers").para_aprovadas(
                    _estado["outdir"], queries=selecionados,
                    opcoes={k: (tuple(v) if isinstance(v, list) else v) for k, v in o.items()})
                resposta.update(res)
                resposta["download"] = "api/download?f=" + quote(os.path.abspath(res["file"]))
            elif acao == "srna_tsv":
                from urllib.parse import quote
                res = _importar("evetool.annotation.srna").tabela_aprovadas(_estado["outdir"], queries=selecionados)
                resposta.update(res)
                resposta["download"] = "api/download?f=" + quote(os.path.abspath(res["file"]))
            elif acao == "host_genes":
                prots = [str(p) for p in pedido.get("host_proteomes") or [] if p]
                resposta.update(cur.rodar_genes_hospedeiro(_estado["outdir"], prots,
                                                           float(pedido.get("host_evalue") or 1e-5), _estado["threads"]))
            elif acao == "dna_tsv":
                from urllib.parse import quote
                res = _importar("evetool.annotation.dna_reads").tabela(_estado["outdir"])
                resposta.update(res)
                resposta["download"] = "api/download?f=" + quote(os.path.abspath(res["file"]))
            elif acao == "metrics":
                resposta.update(cur.rodar_metricas(_estado["outdir"], threads=_estado["threads"]))
            elif acao in ("circ_add", "circ_remove", "circ_run", "circ_gff"):
                cc = _importar("evetool.comparative.codon_composition")
                if acao == "circ_add":
                    resposta["set"] = cc.adicionar_conjunto(
                        _estado["outdir"], pedido.get("circ_fasta") or "", nome=pedido.get("circ_name") or "",
                        tipo=pedido.get("circ_kind") or "circulating", familia_padrao=pedido.get("circ_family") or "",
                        familia_auto=bool(pedido.get("circ_auto_family", True)), threads=_estado["threads"])
                elif acao == "circ_remove":
                    resposta.update(cc.remover_conjunto(_estado["outdir"], str(pedido.get("circ_set") or "")))
                elif acao == "circ_gff":
                    resposta.update(cc.definir_gff(_estado["outdir"], str(pedido.get("circ_species") or ""),
                                                   str(pedido.get("circ_gff_path") or "")))
                else:
                    o = pedido.get("circ_options") or {}
                    r = cc.analisar(_estado["outdir"], top=int(o.get("top") or 6),
                                    discriminar=str(o.get("discriminate") or "family_circ"),
                                    conjunto_features=str(o.get("feature_set") or "top"), ajuste=str(o.get("fit") or "all"),
                                    min_codons=int(o.get("min_codons") or 0), max_por_especie=int(o.get("max_per_species") or 0),
                                    max_hospedeiro=int(o.get("max_host") or 0), classe=str(o.get("protein_class") or ""),
                                    incluir_orf=bool(o.get("include_orf_frame")), perplexidade=float(o.get("perplexity") or 30),
                                    hospedeiro_anotacao=bool(o.get("host_annotation", True)),
                                    flanco_vizinhos=int(o.get("neighbour_flank") or 50000))
                    resposta.update({k: r.get(k) for k in ("counts", "excluded", "messages", "used_features", "amelioration")})
            elif acao != "save":
                return self._json(400, {"error": f"unknown action: {acao}"})

            dados = cur.listar_candidatas(_estado["outdir"])
            resposta["summary"] = dados["summary"]
            resposta["stamp"] = dados["stamp"]
            _estado["ultimo_resumo"] = resposta
            _snapshot()
            return self._json(200, resposta)
        except Exception as erro:
            return self._json(500, {"error": str(erro), "trace": traceback.format_exc()[-1500:]})
        finally:
            _estado["rodando"] = False
            _estado["acao"] = None

    def log_message(self, formato, *args):
        # args[0] is the request line for normal logs, but an HTTPStatus code for send_error
        if "/api/" in (str(args[0]) if args else ""):
            super().log_message(formato, *args)


def main():
    parser = argparse.ArgumentParser(description="Local server for the EVE dashboard")
    parser.add_argument("--port", type=int, default=8080, help="Server port (default: 8080)")
    parser.add_argument("--dir", type=str, default=".", help="Directory to serve")
    parser.add_argument("--threads", type=int, default=4, help="Threads for clustering and context figures")
    parser.add_argument("--outdir", type=str, default=None,
                        help="Pipeline results folder used by the re-analysis "
                             "(default: the served directory, or its parent for 12_dashboard_web)")
    args = parser.parse_args()

    _estado["base_code"] = os.path.dirname(os.path.abspath(__file__))
    _estado["threads"] = args.threads
    dir_servido = os.path.abspath(args.dir)

    if args.outdir:
        _estado["outdir"] = os.path.abspath(args.outdir)
    elif os.path.basename(dir_servido) == "12_dashboard_web":
        _estado["outdir"] = os.path.dirname(dir_servido)
    else:
        _estado["outdir"] = dir_servido

    os.makedirs(dir_servido, exist_ok=True)
    # Refresh the copies of the page on disk too (12_dashboard_web), so opening them
    # without the server shows the current page and data
    try:
        _importar("evetool.curation.curation").instalar_pagina(_estado["outdir"])
    except Exception as erro:
        print(f"[server] page not installed: {erro}", file=sys.stderr)
        _snapshot()

    os.chdir(dir_servido)
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    socketserver.ThreadingTCPServer.daemon_threads = True
    try:
        with socketserver.ThreadingTCPServer(("127.0.0.1", args.port), Handler) as httpd:
            print(f"PLANCUS · EVE page (curation + overview): http://127.0.0.1:{args.port}/")
            print(f"Results used by the re-analysis: {_estado['outdir']}")
            print("Press Ctrl+C to stop.")
            httpd.serve_forever()
    except OSError as e:
        print(f"Could not start the server on port {args.port}: {e}", file=sys.stderr)
    except KeyboardInterrupt:
        print("\nServer stopped.")


if __name__ == "__main__":
    main()
