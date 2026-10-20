#!/usr/bin/env python3
"""
Pfam annotation of the candidate EVEs (an independent layer next to the viral HMM and CDD).

Every candidate peptide is scanned with hmmscan against Pfam-A, using the gathering
thresholds curated by Pfam (--cut_ga), so each family keeps its own cutoff instead of one
E-value for all. CDD already carries Pfam models as rpsblast PSSMs; the HMM search here is
more sensitive for divergent sequences and gives cleaner domain boundaries.

Each domain gets the same viral scope as the CD-Search layer (cdsearch.py), read
from the Pfam family name and description:

  RNA virus / ssDNA virus   in scope
  Retro (excluded)          reverse transcriptase, integrase, gag, retrotransposons...
  Other virus               viral, but neither of the above
  Non-viral                 cellular domains

Outputs, under 04_candidatas_eves/ of each species:
  pfam_raw.tsv       raw domain hits (cache: sequences already scanned are not scanned again)
  pfam_hits.tsv      one domain per line, with its scope, positions and description
  pfam_summary.tsv   one ORF per line: scope, best in-scope domain, all in-scope domains

Usage:
  python -m evetool.annotation.pfam -o results_eve
  python -m evetool.annotation.pfam -o results_eve --pfam-db /path/Pfam-A.hmm -t 8 --forcar
"""
import os
import csv
import shutil
import argparse
import tempfile
import subprocess
from collections import OrderedDict

from evetool.annotation.cdsearch import (escopo_do_dominio, ler_fasta, _assinatura, _gravar, _rival,
                             pastas_especies, ESCOPOS_ALVO)
from evetool import config

PFAM_PADRAO = config.PFAM_HMM
ORDEM_ESCOPO = ["RNA virus", "ssDNA virus", "Retro (excluded)", "Other virus", "Non-viral", "No Pfam hit"]
COLS_RAW = ["orf", "seq_md5", "db", "accession", "short_name", "model_length", "hmm_from", "hmm_to",
            "from", "to", "env_from", "env_to", "ievalue", "bitscore", "description"]


def log(msg):
    print(f"[pfam] {msg}", flush=True)


def _versao_db(pfam_db: str) -> str:
    """Name + size of the database: a new Pfam release invalidates the cache."""
    return f"{os.path.basename(pfam_db)}:{os.path.getsize(pfam_db)}"


def buscar(seqs: "OrderedDict[str, str]", pfam_db: str, threads: int = 4) -> list:
    """Raw domain rows (dicts without orf md5/db) from hmmscan --cut_ga on Pfam-A."""
    if not shutil.which("hmmscan"):
        raise RuntimeError("hmmscan not found on the PATH (HMMER package)")
    if not os.path.exists(pfam_db):
        raise RuntimeError(f"Pfam database not found: {pfam_db}")
    if not os.path.exists(pfam_db + ".h3m"):
        raise RuntimeError(f"Pfam database is not pressed: run  hmmpress {pfam_db}")
    with tempfile.TemporaryDirectory(prefix="pfam_") as tmp:
        query = os.path.join(tmp, "query.fasta")
        tabela = os.path.join(tmp, "domtbl.txt")
        with open(query, "w") as fh:
            fh.writelines(f">{q}\n{s}\n" for q, s in seqs.items())
        cmd = ["hmmscan", "--cut_ga", "--noali", "--cpu", str(threads),
               "--domtblout", tabela, pfam_db, query]
        r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"hmmscan failed: {r.stderr.strip()[:500]}")
        linhas = []
        with open(tabela, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith("#"):
                    continue
                p = linha.rstrip("\n").split(None, 22)
                if len(p) < 22:
                    continue
                linhas.append({"orf": p[3], "accession": p[1], "short_name": p[0], "model_length": p[2],
                               "hmm_from": p[15], "hmm_to": p[16], "from": p[17], "to": p[18],
                               "env_from": p[19], "env_to": p[20], "ievalue": p[12], "bitscore": p[13],
                               "description": p[22].strip() if len(p) > 22 else ""})
    return linhas


# ---------------------------------------------------------------------------
# One species folder
# ---------------------------------------------------------------------------
def anotar_pasta(pasta: str, pfam_db: str = PFAM_PADRAO, forcar: bool = False, threads: int = 4,
                 fasta: str = None, dir_saida: str = None, prefixo: str = "") -> dict:
    """fasta/dir_saida/prefixo scan another peptide set and write <prefixo>pfam_*.tsv to dir_saida."""
    dir_c = dir_saida or os.path.join(pasta, "04_candidatas_eves")
    fasta = fasta or os.path.join(dir_c, "eve_candidatas_unificadas_pep.fasta")
    if not os.path.exists(fasta):
        return {"orfs": 0, "searched": 0}
    seqs = ler_fasta(fasta)
    versao = _versao_db(pfam_db)
    raw_path = os.path.join(dir_c, prefixo + "pfam_raw.tsv")

    # Cache: an ORF already scanned with the same sequence and the same Pfam is not scanned again
    feitos, brutos = set(), []
    if os.path.exists(raw_path) and not forcar:
        with open(raw_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["orf"] in seqs and r["seq_md5"] == _assinatura(seqs[r["orf"]]) and r["db"] == versao:
                    feitos.add(r["orf"])
                    brutos.append(r)      # no-hit rows too: they mark the ORF as scanned
    faltam = OrderedDict((q, s) for q, s in seqs.items() if q not in feitos)
    log(f"{os.path.basename(pasta.rstrip('/')) or pasta}: {len(seqs)} ORFs, "
        f"{len(feitos)} cached, {len(faltam)} to scan")

    if faltam:
        achados = buscar(faltam, pfam_db, threads=threads)
        com_hit = set()
        for h in achados:
            if h["orf"] not in faltam:
                continue
            com_hit.add(h["orf"])
            brutos.append({**h, "seq_md5": _assinatura(faltam[h["orf"]]), "db": versao})
        for q in faltam:
            if q not in com_hit:
                brutos.append({"orf": q, "seq_md5": _assinatura(faltam[q]), "db": versao})
    _gravar(raw_path, COLS_RAW, brutos)

    hits = [{"orf": r["orf"], "accession": r["accession"], "short_name": r["short_name"],
             "from": r["from"], "to": r["to"], "env_from": r["env_from"], "env_to": r["env_to"],
             "hmm_from": r["hmm_from"], "hmm_to": r["hmm_to"], "model_length": r["model_length"],
             "model_coverage": round(100 * (int(r["hmm_to"]) - int(r["hmm_from"]) + 1)
                                     / max(1, int(r["model_length"])), 1),
             "ievalue": r["ievalue"], "bitscore": r["bitscore"],
             "scope": escopo_do_dominio(r["short_name"], r["description"]),
             "description": r["description"][:200]}
            for r in brutos if r.get("accession")]
    cols_hits = ["orf", "scope", "accession", "short_name", "from", "to", "env_from", "env_to",
                 "hmm_from", "hmm_to", "model_length", "model_coverage", "ievalue", "bitscore", "description"]
    _gravar(os.path.join(dir_c, prefixo + "pfam_hits.tsv"), cols_hits, hits)

    # One row per ORF, with the same rule as CD-Search: a viral domain beaten on the same
    # region by a stronger non-viral one does not make the ORF viral
    por_orf = {}
    for h in hits:
        por_orf.setdefault(h["orf"], []).append(h)
    resumo = []
    for orf in seqs:
        lista = por_orf.get(orf, [])
        vencidos = {id(h): _rival(h, lista) for h in lista if h["scope"] in ESCOPOS_ALVO}
        vencidos = {k: v for k, v in vencidos.items() if v}
        validos = [h for h in lista if id(h) not in vencidos]
        escopos = {h["scope"] for h in validos}
        escopo = next((e for e in ORDEM_ESCOPO if e in escopos), "No Pfam hit")
        alvo = [h for h in validos if h["scope"] in ESCOPOS_ALVO]
        ordem = sorted(alvo or validos, key=lambda h: -float(h["bitscore"] or 0))
        melhor = ordem[0] if ordem else None
        nomes_alvo = list(OrderedDict.fromkeys(h["short_name"] for h in sorted(
            alvo, key=lambda h: -float(h["bitscore"] or 0))))
        resumo.append({
            "orf": orf, "pfam_scope": escopo, "pfam_in_scope": "yes" if alvo else "no",
            "pfam_n_domains": len(validos), "pfam_n_in_scope": len(alvo),
            "pfam_best": melhor["short_name"] if melhor else "",
            "pfam_best_accession": melhor["accession"] if melhor else "",
            "pfam_best_scope": melhor["scope"] if melhor else "",
            "pfam_best_evalue": melhor["ievalue"] if melhor else "",
            "pfam_best_bitscore": melhor["bitscore"] if melhor else "",
            "pfam_best_description": melhor["description"] if melhor else "",
            "pfam_in_scope_domains": "; ".join(nomes_alvo),
            "pfam_outcompeted": "; ".join(OrderedDict.fromkeys(
                f"{h['short_name']} < {vencidos[id(h)]['short_name']}" for h in lista if id(h) in vencidos)),
        })
    cols_res = list(resumo[0]) if resumo else ["orf"]
    _gravar(os.path.join(dir_c, prefixo + "pfam_summary.tsv"), cols_res, resumo)
    n_alvo = sum(1 for r in resumo if r["pfam_in_scope"] == "yes")
    log(f"  {n_alvo}/{len(resumo)} ORFs with an RNA- or ssDNA-virus Pfam domain")
    return {"orfs": len(resumo), "searched": len(faltam), "in_scope": n_alvo,
            "summary_rows": resumo, "hit_rows": hits,
            "with_domain": sum(1 for r in resumo if r["pfam_n_domains"]),
            "per_scope": {e: sum(1 for r in resumo if r["pfam_scope"] == e) for e in ORDEM_ESCOPO}}


def anotar(outdir: str, pfam_db: str = PFAM_PADRAO, forcar: bool = False, threads: int = 4) -> dict:
    log(f"Pfam: {pfam_db}")
    total = {"orfs": 0, "searched": 0, "in_scope": 0, "with_domain": 0,
             "per_scope": {e: 0 for e in ORDEM_ESCOPO}}
    for _sp, pasta in pastas_especies(outdir).items():
        r = anotar_pasta(pasta, pfam_db=pfam_db, forcar=forcar, threads=threads)
        for k in ("orfs", "searched", "in_scope", "with_domain"):
            total[k] += r.get(k, 0)
        for e, n in (r.get("per_scope") or {}).items():
            total["per_scope"][e] += n
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--pfam-db", default=PFAM_PADRAO, help="Pfam-A.hmm, pressed with hmmpress")
    ap.add_argument("--forcar", action="store_true", help="Scan every ORF again, ignoring the cache")
    ap.add_argument("-t", "--threads", type=int, default=4, help="hmmscan threads")
    a = ap.parse_args()
    r = anotar(a.outdir, pfam_db=a.pfam_db, forcar=a.forcar, threads=a.threads)
    print(f"\n{r['orfs']} ORFs ({r['searched']} scanned now) | with a Pfam domain: {r['with_domain']}"
          f" | in scope (RNA/ssDNA virus): {r['in_scope']}")
    for e, n in r["per_scope"].items():
        print(f"  {n:5}  {e}")


if __name__ == "__main__":
    main()
