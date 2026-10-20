#!/usr/bin/env python3
"""
Local server for the EVE page: curation and overview charts on one page.

  GET  /, /index.html, /curadoria.html -> the page, always read from evetool/curation/web/curation.html,
                           so an edit to its code (colours, layout) shows on reload without
                           touching the data
  GET  /api/candidates  -> classified EVEs, curation state, review flags and sequence metrics
  GET  /api/state       -> version stamp of the data on disk (the page polls it to refresh itself)
  GET  /api/detail      -> inspect one EVE: position, neighbours, GC and alignments
  GET  /api/protein-map -> where each EVE falls on its viral protein (?force=1 realigns)
  GET  /api/rnaseq      -> binned RNA-seq depth of every sample over ?sp=&c=&s=&e= (0-based) [&focus=a-b,...]
  GET  /api/rnaseq-samples -> RNA-seq BAMs known per species
  POST /api/rnaseq-bam  -> add ({species, bam, library}) or remove ({op: "remove", species, name}) a BAM
  POST /api/export      -> save a chart as SVG and, when possible, convert it to PDF
  POST /api/regions     -> multi-FASTA of several genome regions (copy / download FASTA)
  POST /api/curation    -> save the selection and run the requested action:
                           save | fasta_nt | fasta_pep | annotate_hmm | cluster |
                           multispecies | context | metrics | repeats | cdsearch | pfam |
                           tblastn_scan

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

    def _pagina(self):
        caminho = os.path.join(_estado["base_code"], "web", "curation.html")
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
                return self._json(200, _importar("evetool.curation.curation").detalhar_eve(_estado["outdir"], q))
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
        if rota == "/api/blastx":
            return self._json(200, _situacao_blastx())
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
                               if isinstance(pedido.get("tblastn_removed"), list) else None)
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
            elif acao == "metrics":
                resposta.update(cur.rodar_metricas(_estado["outdir"], threads=_estado["threads"]))
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
            print(f"EVE page (curation + overview): http://127.0.0.1:{args.port}/")
            print(f"Results used by the re-analysis: {_estado['outdir']}")
            print("Press Ctrl+C to stop.")
            httpd.serve_forever()
    except OSError as e:
        print(f"Could not start the server on port {args.port}: {e}", file=sys.stderr)
    except KeyboardInterrupt:
        print("\nServer stopped.")


if __name__ == "__main__":
    main()
