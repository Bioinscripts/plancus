#!/usr/bin/env python3
"""
HMM annotation of the EVEs already identified.

During screening the HMM models only see the ORFs that did NOT align in DIAMOND: they
exist to rescue what the alignment misses. Here the calibrated database is applied to ALL
candidates, including those confirmed by DIAMOND, as an extra layer of information:
which viral domains each segment carries, at what score, and which family and genome
type those models belong to.

It corroborates (or contradicts) the assignment made by the alignment and helps
characterise segments of uncertain origin.

Outputs, under 04_candidatas_eves/:
  hmm_candidates_domtbl.tsv  raw hmmsearch output (domtblout)
  hmm_annotation.tsv         one domain per line, with family and genome type
  hmm_summary.tsv            one ORF per line, with the best domain and the consensus

Usage:
  python -m evetool.annotation.hmm_domains -o results_eve --hmm-db calibrated_database.hmm -t 5
"""
import os
import csv
import glob
import argparse
from collections import defaultdict, Counter
from evetool import config

BASE = config.HMM_DIR
HMM_PADRAO = f"{BASE}/COMBINED_DATABASE/combined_viral_calibrado.hmm"
ANOTACOES_PADRAO = f"{BASE}/COMBINED_DATABASE/combined_viral_annotations.tsv"


def log(msg):
    print(f"[hmm-anot] {msg}", flush=True)


def carregar_tc(hmm_db: str) -> dict:
    """
    TC cutoff of each model in the calibrated database, cached on disk.
    It lets a weak hit be reported together with how far below the cutoff it fell.
    """
    cache = os.path.splitext(hmm_db)[0] + "_tc.json"
    if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(hmm_db):
        try:
            import json
            return json.load(open(cache, encoding="utf-8"))
        except Exception:
            pass
    import json
    import pyhmmer
    tcs = {}
    log("reading the TC cutoffs from the database (once; then cached)...")
    with pyhmmer.plan7.HMMFile(hmm_db) as hf:
        for h in hf:
            nome = h.name.decode() if isinstance(h.name, bytes) else str(h.name)
            t = getattr(h.cutoffs, "trusted", None)
            if t:
                tcs[nome] = float(t[0])
    try:
        json.dump(tcs, open(cache, "w", encoding="utf-8"))
    except Exception:
        pass
    return tcs


def carregar_anotacoes(caminho: str) -> dict:
    """{HMM_ID: {family, genome, source, keywords}} from the combined database."""
    mapa = {}
    if not caminho or not os.path.exists(caminho):
        return mapa
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mapa[r.get("HMM_ID", "")] = {
                "family": r.get("Family", "") or "",
                "genoma": r.get("Genome_Type", "") or "",
                "source": r.get("Source_DB", "") or "",
                "palavras": (r.get("Protein_Keywords", "") or "")[:70],
            }
    return mapa


def ler_domtbl(caminho: str) -> dict:
    """
    HMMER domtblout: space-separated columns. Uses target(0), model(3), model length(5),
    i-Evalue(12), domain score(13) and the model (15,16) and target (17,18) coordinates.
    """
    por_orf = defaultdict(list)
    if not os.path.exists(caminho):
        return por_orf
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            p = linha.split()
            if len(p) < 19:
                continue
            try:
                por_orf[p[0]].append({
                    "model": p[3], "tlen_modelo": int(p[5]),
                    "ievalue": float(p[12]), "score": float(p[13]),
                    "hmm_de": int(p[15]), "hmm_ate": int(p[16]),
                    "orf_from": int(p[17]), "orf_to": int(p[18]),
                })
            except ValueError:
                continue
    return por_orf


def anotar(outdir: str, hmm_db: str = HMM_PADRAO, anotacoes: str = ANOTACOES_PADRAO,
           threads: int = 4, pasta_especie: str = None, forcar: bool = False,
           evalue_permissivo: float = 1.0) -> dict:
    import sys
    from evetool.discovery.hmm_scanner import escanear_hmm

    pasta = pasta_especie or outdir
    dir_cand = os.path.join(pasta, "04_candidatas_eves")
    pep = os.path.join(dir_cand, "eve_candidatas_unificadas_pep.fasta")
    if not os.path.exists(pep) or os.path.getsize(pep) == 0:
        raise RuntimeError(f"Candidate FASTA not found: {pep}")
    n_cand = sum(1 for l in open(pep) if l.startswith(">"))

    domtbl = os.path.join(dir_cand, "hmm_candidates_domtbl.tsv")
    if forcar or not os.path.exists(domtbl) or os.path.getsize(domtbl) == 0:
        log(f"scanning {n_cand} candidates against the HMM database (permissive mode, "
            f"e-value <= {evalue_permissivo})...")
        # Permissive scan: every candidate gets a score, even below the TC. The TC is then
        # applied as a flag instead of erasing the row.
        gerado = escanear_hmm(pep, dir_cand, hmm_db=hmm_db, threads=threads,
                              usar_tc=False, evalue=evalue_permissivo)
        if os.path.abspath(gerado) != os.path.abspath(domtbl):
            os.replace(gerado, domtbl)
    else:
        log(f"reusing {domtbl} (pass --forcar to redo it)")

    mapa = carregar_anotacoes(anotacoes)
    tcs = carregar_tc(hmm_db)
    por_orf = ler_domtbl(domtbl)

    linhas, resumo = [], []
    for orf in (l[1:].split()[0] for l in open(pep) if l.startswith(">")):
        dominios = sorted(por_orf.get(orf, []), key=lambda d: -d["score"])
        for d in dominios:
            a = mapa.get(d["model"], {})
            tc = tcs.get(d["model"])
            linhas.append({
                "orf": orf, "model": d["model"], "score": round(d["score"], 1),
                "model_tc": round(tc, 1) if tc else "",
                "above_tc": bool(tc and d["score"] >= tc),
                "ievalue": f"{d['ievalue']:.2g}",
                "model_coverage": round(100.0 * (d["hmm_ate"] - d["hmm_de"] + 1) / max(1, d["tlen_modelo"]), 1),
                "orf_from": d["orf_from"], "orf_to": d["orf_to"],
                "hmm_family": a.get("family", ""), "hmm_genome": a.get("genoma", ""),
                "source": a.get("source", ""), "keywords": a.get("palavras", ""),
            })
        acima = [d for d in dominios if tcs.get(d["model"]) and d["score"] >= tcs[d["model"]]]
        if dominios:
            # The best hit above the TC wins; if none passed, still report the best one
            melhor = (acima or dominios)[0]
            a = mapa.get(melhor["model"], {})
            base = acima or dominios
            families = Counter(mapa.get(d["model"], {}).get("family", "") for d in base)
            families.pop("", None)
            consenso, n_cons = (families.most_common(1)[0] if families else ("", 0))
            tc = tcs.get(melhor["model"])
            resumo.append({
                "orf": orf, "n_domains": len(acima), "n_domains_total": len(dominios),
                "best_model": melhor["model"], "best_score": round(melhor["score"], 1),
                "model_tc": round(tc, 1) if tc else "",
                "above_tc": bool(tc and melhor["score"] >= tc),
                "tc_margin": round(melhor["score"] - tc, 1) if tc else "",
                "best_ievalue": f"{melhor['ievalue']:.2g}",
                "best_model_family": a.get("family", ""), "best_model_genome": a.get("genoma", ""),
                "hmm_consensus_family": consenso,
                "consensus_support": f"{n_cons}/{len(base)}",
            })
        else:
            resumo.append({"orf": orf, "n_domains": 0, "n_domains_total": 0, "best_model": "",
                           "best_score": 0.0, "model_tc": "", "above_tc": False,
                           "tc_margin": "", "best_ievalue": "", "best_model_family": "",
                           "best_model_genome": "", "hmm_consensus_family": "", "consensus_support": ""})

    def grava(nome, dados):
        caminho = os.path.join(dir_cand, nome)
        if dados:
            with open(caminho, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(dados[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(dados)
        return caminho

    com = sum(1 for r in resumo if r["n_domains"])
    com_qualquer = sum(1 for r in resumo if r["n_domains_total"])
    log(f"{com}/{n_cand} above the TC cutoff | {com_qualquer}/{n_cand} with any hit "
        f"({len(linhas)} domains recorded)")
    return {"domtbl": domtbl, "anotacao": grava("hmm_annotation.tsv", linhas),
            "resumo": grava("hmm_summary.tsv", resumo),
            "n_candidatas": n_cand, "n_com_dominio": com, "n_domains": len(linhas)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--hmm-db", default=HMM_PADRAO)
    ap.add_argument("--anotacoes", default=ANOTACOES_PADRAO)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--forcar", action="store_true", help="Rerun the scan even when the output already exists")
    ap.add_argument("--evalue-permissivo", type=float, default=1.0,
                    help="Scan e-value; the TC cutoff becomes a flag, not a filter (default: 1.0)")
    cfg = ap.parse_args()

    pastas = [cfg.outdir]
    especies = os.path.join(cfg.outdir, "especies")
    if os.path.isdir(especies):
        pastas = sorted(p for p in glob.glob(os.path.join(especies, "*")) if os.path.isdir(p))
    for pasta in pastas:
        try:
            r = anotar(cfg.outdir, cfg.hmm_db, cfg.anotacoes, cfg.threads,
                       pasta_especie=pasta, forcar=cfg.forcar,
                       evalue_permissivo=cfg.evalue_permissivo)
            print(f"{os.path.basename(pasta)}: {r['n_com_dominio']}/{r['n_candidatas']} with a domain "
                  f"-> {r['resumo']}")
        except Exception as erro:
            print(f"[warning] {os.path.basename(pasta)}: {erro}")


if __name__ == "__main__":
    main()
