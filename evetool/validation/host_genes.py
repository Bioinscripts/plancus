#!/usr/bin/env python3
"""
Host genes hit by the candidate EVEs.

A host gene is not always a non-EVE (an old insertion may be annotated as a gene of the host,
and some cellular genes are hard to tell from viral ones), so nothing is removed here. The
peptides of every candidate are searched against reference host proteomes (DIAMOND blastp) and
the hits are grouped by host gene (isoforms collapsed). The curation page lists those genes;
the researcher removes the ones they judge to be real host genes, which rejects every candidate
ORF that hit them, and exports the list of removed genes to report in the methods.

Each candidate keeps its best hit in each proteome. Flags that help the decision:

  uncharacterised   the gene has no described function ("uncharacterized protein LOC..."):
                    it may be the EVE itself, annotated as a gene
  same species      the hit is in the proteome of the candidate's own species at >= 95 %
                    identity: the annotation of that genome simply contains this ORF
  host >= viral     the host hit scores at least as much as the best viral hit

Files (13_curadoria/host_genes/): host_hits.tsv (one line per ORF x proteome),
host_genes.json (what the page reads).
"""

import argparse
import csv
import datetime
import json
import os
import re
import shutil
import subprocess

SUBPASTA = os.path.join("13_curadoria", "host_genes")
RX_NAO_CARACTERIZADA = re.compile(
    r"uncharacteri[sz]ed|hypothetical|unknown function|unnamed protein|predicted protein|"
    r"\bDUF\d+|\bUPF\d+|\bLOC\d+", re.I)
COLUNAS = ["query", "proteome", "protein", "gene", "description", "organism", "pident", "length",
           "evalue", "bitscore", "qcov", "scov"]


def log(msg):
    print(f"[host-genes] {msg}", flush=True)


def _pasta(outdir: str) -> str:
    return os.path.join(outdir, SUBPASTA)


def _peptideos(outdir: str, destino: str) -> int:
    """Peptides of every candidate of the run in one FASTA, named `species|orf` like the queries."""
    from evetool.curation.curation import pastas_especies
    multi = os.path.isdir(os.path.join(outdir, "especies"))
    n = 0
    with open(destino, "w") as saida:
        for species, pasta in pastas_especies(outdir).items():
            cam = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
            if not os.path.exists(cam):
                continue
            with open(cam, encoding="utf-8", errors="replace") as fh:
                for linha in fh:
                    if linha.startswith(">"):
                        orf = linha[1:].split()[0]
                        saida.write(f">{species}|{orf}\n" if (multi and species) else f">{orf}\n")
                        n += 1
                    else:
                        saida.write(linha.rstrip("*\n").replace("*", "X") + "\n")
    return n


def _organismo(descricao: str) -> str:
    m = re.search(r"\[([^\]]+)\]\s*$", descricao or "")
    return m.group(1) if m else ""


def rodar(outdir: str, proteomas=None, evalue: float = 1e-5, threads: int = 4, passo=None) -> dict:
    """Search the candidates against the host proteomes (all the installed ones by default)."""
    from evetool.comparative.flank_synteny import descrever, preparar_proteoma
    from evetool.curation.curation import proteomas_sintenia
    passo = passo or (lambda t: None)
    if not shutil.which("diamond"):
        raise RuntimeError("diamond not found in the PATH")
    instalados = proteomas_sintenia()
    if proteomas:
        quais = set(proteomas)
        escolhidos = [p for p in instalados if p["path"] in quais] + [
            {"path": p, "label": os.path.basename(p)} for p in proteomas
            if p not in {x["path"] for x in instalados} and os.path.exists(p)]
    else:
        escolhidos = instalados
    if not escolhidos:
        raise RuntimeError("no host proteome found (expected *.faa in <hmm dir>/proteomas_calibracao)")
    os.makedirs(_pasta(outdir), exist_ok=True)
    fasta = os.path.join(_pasta(outdir), "candidates_pep.fasta")
    n = _peptideos(outdir, fasta)
    if not n:
        raise RuntimeError("no candidate peptide found in this run")
    linhas, usados = [], []
    for k, prot in enumerate(escolhidos, 1):
        passo(f"host proteome {k} of {len(escolhidos)}: {prot['label']}")
        log(f"{prot['label']}: preparing")
        ref = preparar_proteoma(prot["path"], threads)
        tsv = os.path.join(_pasta(outdir), f"tmp_{k}.tsv")
        r = subprocess.run(
            ["diamond", "blastp", "-q", fasta, "-d", ref["dmnd"], "-o", tsv, "--evalue", str(evalue),
             "--max-target-seqs", "1", "--sensitive", "--threads", str(threads), "--quiet", "--outfmt", "6",
             "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovhsp", "scovhsp"],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"diamond failed on {prot['label']}: {(r.stderr or '')[-300:]}")
        organismo = ""
        achados = 0
        with open(tsv, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                p = linha.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                desc = ref["nome"].get(p[1], p[1])
                gene = ref["grupo"].get(p[1], p[1])
                organismo = organismo or _organismo(desc)
                linhas.append({"query": p[0], "proteome": os.path.basename(prot["path"]), "protein": p[1],
                               "gene": gene, "description": descrever(ref["nome"].get(gene, desc)),
                               "organism": _organismo(desc) or organismo, "pident": float(p[2]),
                               "length": int(p[3]), "evalue": p[4], "bitscore": float(p[5]),
                               "qcov": float(p[6]), "scov": float(p[7])})
                achados += 1
        os.remove(tsv)
        usados.append({"file": os.path.basename(prot["path"]), "label": prot["label"],
                       "organism": organismo or prot["label"], "proteins": len(ref["nome"]), "orfs_with_hit": achados})
        log(f"{prot['label']}: {achados} of {n} candidate ORFs hit a protein")
    with open(os.path.join(_pasta(outdir), "host_hits.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUNAS, delimiter="\t")
        w.writeheader()
        w.writerows(linhas)
    dados = {"ran_at": datetime.datetime.now().isoformat(timespec="seconds"), "evalue": evalue,
             "candidates": n, "proteomes": usados, "hits": linhas,
             "method": f"DIAMOND blastp --sensitive, e-value <= {evalue}, best hit per ORF in each proteome; "
                       f"isoforms collapsed with MMseqs2 (90% identity, 80% coverage)"}
    with open(os.path.join(_pasta(outdir), "host_genes.json"), "w", encoding="utf-8") as fh:
        json.dump(dados, fh)
    os.remove(fasta)
    return {"candidates": n, "orfs_with_hit": len({l["query"] for l in linhas}),
            "genes": len({(l["proteome"], l["gene"]) for l in linhas}), "proteomes": len(usados)}


def carregar(outdir: str, itens: list = None):
    """
    The genes for the page: [{id, proteome, organism, gene, description, uncharacterised, orfs:
    [{query, bitscore, pident, qcov, scov, protein}], best_bitscore, max_pident}], strongest first.
    With `itens`, the best cellular hit of the general database (Swiss-Prot) is added as one
    more source, so those genes can be removed the same way. None when nothing is known.
    """
    cam = os.path.join(_pasta(outdir), "host_genes.json")
    dados = None
    if os.path.exists(cam):
        try:
            with open(cam, encoding="utf-8") as fh:
                dados = json.load(fh)
        except Exception:
            dados = None
    hits = list((dados or {}).get("hits", []))
    for i in itens or []:
        if i.get("cellular_bitscore", 0) > 0 and i.get("best_cellular") not in ("", "n/a", None):
            titulo = i.get("cellular_title") or ""
            m = re.search(r"OS=(.+?)(?:\s+OX=|\s+GN=|$)", titulo)
            hits.append({"query": i["query"], "proteome": "general database", "protein": i["best_cellular"],
                         "gene": i["best_cellular"],
                         "description": re.sub(r"\s+OS=.*$", "", titulo.split(" ", 1)[1] if " " in titulo else titulo),
                         "organism": m.group(1) if m else "", "pident": None, "length": None, "evalue": "",
                         "bitscore": i["cellular_bitscore"], "qcov": None, "scov": None})
    if not hits:
        return None if dados is None else {**{k: v for k, v in dados.items() if k != "hits"}, "genes": []}
    genes = {}
    for h in hits:
        chave = f"{h['proteome']}|{h['gene']}"
        g = genes.setdefault(chave, {"id": chave, "proteome": h["proteome"], "organism": h["organism"],
                                     "gene": h["gene"], "description": h["description"],
                                     "uncharacterised": bool(RX_NAO_CARACTERIZADA.search(h["description"] or ""))
                                     or not h["description"], "orfs": []})
        g["orfs"].append({k: h[k] for k in ("query", "bitscore", "pident", "qcov", "scov", "protein")})
    for g in genes.values():
        g["orfs"].sort(key=lambda o: -o["bitscore"])
        g["best_bitscore"] = g["orfs"][0]["bitscore"]
        g["max_pident"] = max((o["pident"] for o in g["orfs"] if o["pident"] is not None), default=None)
    base = {k: v for k, v in (dados or {}).items() if k != "hits"}
    return {**base, "genes": sorted(genes.values(), key=lambda g: (-len(g["orfs"]), -g["best_bitscore"]))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="results folder of the run")
    ap.add_argument("--proteome", action="append", help="host proteome (protein FASTA); repeat for several. "
                                                        "Default: every installed proteome")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("-t", "--threads", type=int, default=4)
    a = ap.parse_args()
    print(json.dumps(rodar(a.outdir, a.proteome, a.evalue, a.threads, log), indent=1))


if __name__ == "__main__":
    main()
