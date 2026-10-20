#!/usr/bin/env python3
"""
Audit of the viral database: find reference proteins that are, at origin, host genes
captured by the virus (horizontal transfer from host to virus).

An ORF aligning to such a protein looks like an EVE but is a cellular gene. Instead of
filtering ORF by ORF afterwards, the audit flags the problem in the reference database.

Method: align the viral proteins of interest (non-retroviral RNA + ssDNA) against a
cellular set (Swiss-Prot without viruses + dipteran proteomes) and record, for each viral
protein, the best cellular similarity and the TAXONOMIC BREADTH of the hits.

Breadth is what separates the two scenarios; identity alone does not:
  - cellular hits spread over many phyla   -> an ancient cellular gene captured by the virus
  - hits restricted to a few close taxa    -> possibly a domesticated EVE; do NOT discard

Nothing is removed automatically: the output is a table for curation plus a file of
suggested accessions, which is only used if the researcher approves it.

Usage:
  python -m evetool.databases.audit_viral_db --threads 6
  python -m evetool.databases.audit_viral_db --min-pident 40 --min-filos 4 --threads 6
"""
import os
import re
import csv
import sys
import tarfile
import argparse
import subprocess
from collections import defaultdict, Counter
from evetool import config

BASE = config.HMM_DIR
PADROES = {
    "viral": f"{BASE}/COMBINED_DATABASE/calibracao/3_positivos_refseq_RNA_ssDNA.faa",
    "celular": f"{BASE}/COMBINED_DATABASE/calibracao/2_negativos_nr70.faa",
    "taxdump": f"{BASE}/RVDB/taxdump.tar.gz",
    "saida": f"{BASE}/COMBINED_DATABASE/auditoria_banco",
}
# Protein families known to be cellular that viruses often capture
REGEX_FAMILIA_CELULAR = re.compile(
    r"ham1|itpase|inosine triphosphate|dITP|pyrophosphatase|phosphatase|"
    r"chaperonin|groel|calcineurin|cdc48|thymidylate|dUTPase|ribonucleotide reductase|"
    r"thioredoxin|glutaredoxin|ubiquitin|kinase|actin|tubulin|histone|"
    r"ankyrin|elongation factor|ribosomal", re.IGNORECASE)


def log(msg):
    print(f"[auditoria] {msg}", flush=True)


def linhagens(taxdump: str) -> tuple:
    """Return (parent, names, rank) so each taxid lineage can be rebuilt."""
    pai, nomes, rank = {}, {}, {}
    with tarfile.open(taxdump) as tar:
        for l in tar.extractfile("nodes.dmp"):
            p = l.decode().split("\t|\t")
            pai[p[0]], rank[p[0]] = p[1], p[2]
        for l in tar.extractfile("names.dmp"):
            p = l.decode().split("\t|\t")
            if p[3].startswith("scientific name"):
                nomes[p[0]] = p[1]
    return pai, nomes, rank


def filo_de(taxid: str, pai: dict, nomes: dict, rank: dict, cache: dict) -> str:
    """Phylum (or the highest rank available) of a taxid, to measure taxonomic breadth."""
    if taxid in cache:
        return cache[taxid]
    t, visto, achado = taxid, 0, None
    while t and t in pai and t != "1" and visto < 60:
        if rank.get(t) in ("phylum", "kingdom", "superkingdom"):
            achado = nomes.get(t)
            if rank.get(t) == "phylum":
                break
        t, visto = pai[t], visto + 1
    cache[taxid] = achado or "?"
    return cache[taxid]


def organismo_e_taxid(stitle: str) -> tuple:
    """Swiss-Prot carries 'OS=... OX=...'; the dipteran proteomes carry '[Organism]'."""
    ox = re.search(r"OX=(\d+)", stitle)
    os_ = re.search(r"OS=(.+?)(?: OX=| GN=| PE=|$)", stitle)
    col = re.search(r"\[([^\]]+)\]", stitle)
    return (os_.group(1).strip() if os_ else (col.group(1) if col else "?"),
            ox.group(1) if ox else None)


def preparar_viral(cfg) -> str:
    """
    With --filtrar-genoma-alvo, cut from the raw viral FASTA only the viruses of interest
    (non-retroviral RNA + ssDNA, minus the excluded families), using the same genome-type
    criterion as the pipeline. When new proteins are added to the database, just point at
    the full FASTA and rerun the audit.
    """
    if not cfg.filtrar_genoma_alvo:
        return cfg.viral
    from evetool.core.utils import classe_genoma_viral
    from evetool.databases.calibrate_hmm import FAMILIAS_EXCLUIDAS

    os.makedirs(cfg.saida, exist_ok=True)
    destino = os.path.join(cfg.saida, "viral_alvo.faa")
    n, gravar = 0, False
    with open(cfg.viral, encoding="utf-8", errors="replace") as fin, open(destino, "w") as fout:
        for linha in fin:
            if linha.startswith(">"):
                campos = linha.split("|")
                genoma = campos[2].strip() if len(campos) > 2 else ""
                family = campos[5].strip() if len(campos) > 5 else ""
                gravar = (classe_genoma_viral(genoma) in ("RNA", "ssDNA")
                          and family not in FAMILIAS_EXCLUIDAS)
                n += gravar
            if gravar:
                fout.write(linha)
    log(f"viruses of interest kept: {n} proteins -> {destino}")
    return destino


def accessions_ja_auditados(cfg) -> set:
    """Accessions already present in a previous audit, so only new ones are audited."""
    anterior = os.path.join(cfg.saida, "viral_database_audit.tsv")
    if not cfg.apenas_novas or not os.path.exists(anterior):
        return set()
    with open(anterior, encoding="utf-8") as fh:
        return {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}


def alinhar(cfg) -> str:
    """DIAMOND blastp of the viral proteins against the cellular set."""
    os.makedirs(cfg.saida, exist_ok=True)
    tsv = os.path.join(cfg.saida, "viral_vs_celular.tsv")
    if os.path.exists(tsv) and os.path.getsize(tsv) > 0 and not cfg.refazer:
        log(f"alignment already exists (pass --refazer to redo it): {tsv}")
        return tsv
    db = os.path.join(cfg.saida, "celular_db")
    if not os.path.exists(db + ".dmnd"):
        log("formatting the cellular database...")
        subprocess.run(["diamond", "makedb", "--in", cfg.celular, "-d", db],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    log(f"aligning the viral proteins against the cellular set (threads={cfg.threads})...")
    subprocess.run(
        ["diamond", "blastp", "-q", cfg.viral, "-d", db, "-o", tsv,
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
         "qlen", "slen", "qcovhsp", "stitle",
         "--max-target-seqs", str(cfg.hits_por_proteina), "--evalue", str(cfg.evalue),
         "--more-sensitive", "--threads", str(cfg.threads)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return tsv


def auditar(cfg):
    cfg.viral = preparar_viral(cfg)
    ja_vistos = accessions_ja_auditados(cfg)
    if ja_vistos:
        log(f"incremental mode: {len(ja_vistos)} already audited accessions will be skipped")
    tsv = alinhar(cfg)
    log("loading taxonomy...")
    pai, nomes, rank = linhagens(cfg.taxdump)
    cache = {}

    # Viral protein headers: ACC |product [organism]|genome|species|genus|family|len
    info_viral = {}
    with open(cfg.viral, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                acc = linha[1:].split()[0]
                campos = linha.split("|")
                info_viral[acc] = {
                    "product": campos[1].strip() if len(campos) > 1 else "",
                    "genoma": campos[2].strip() if len(campos) > 2 else "",
                    "family": campos[5].strip() if len(campos) > 5 else "",
                }

    hits = defaultdict(list)
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            org, taxid = organismo_e_taxid(p[9])
            hits[p[0]].append({"pident": float(p[2]), "bits": float(p[5]),
                               "qcov": float(p[8]), "org": org, "taxid": taxid,
                               "titulo": p[9][:90]})

    linhas = []
    for acc, lista in hits.items():
        if acc in ja_vistos:
            continue
        lista.sort(key=lambda h: -h["bits"])
        melhor = lista[0]
        organismos = {h["org"] for h in lista}
        phyla = {filo_de(h["taxid"], pai, nomes, rank, cache) for h in lista if h["taxid"]}
        phyla.discard("?")
        iv = info_viral.get(acc, {})

        amplo = len(phyla) >= cfg.min_filos
        forte = melhor["pident"] >= cfg.min_pident and melhor["qcov"] >= cfg.min_qcov
        familia_conhecida = bool(REGEX_FAMILIA_CELULAR.search(iv.get("product", "")))

        if forte and amplo:
            verdict = "SUSPECT_cellular_HGT"
        elif forte and familia_conhecida:
            verdict = "SUSPECT_cellular_family"
        elif forte:
            verdict = "REVIEW_restricted_hits"   # possibly a domesticated EVE; do not discard
        else:
            verdict = "ok_weak_similarity"

        linhas.append({
            "accession": acc, "verdict": verdict,
            "viral_product": iv.get("product", "")[:90], "viral_family": iv.get("family", ""),
            "viral_genome": iv.get("genoma", ""),
            "cellular_pident": round(melhor["pident"], 1),
            "cellular_qcov": round(melhor["qcov"], 1),
            "cellular_bits": round(melhor["bits"], 1),
            "best_organism": melhor["org"][:40],
            "best_cellular_title": melhor["titulo"],
            "n_organisms": len(organismos), "n_phyla": len(phyla),
            "phyla": ";".join(sorted(phyla)[:6]),
        })

    linhas.sort(key=lambda d: (d["verdict"], -d["cellular_pident"]))
    destino = os.path.join(cfg.saida, "viral_database_audit.tsv")
    with open(destino, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(linhas)

    sugeridos = [d for d in linhas if d["verdict"].startswith("SUSPEITO")]
    lista_sug = os.path.join(cfg.saida, "accessions_suggested_removal.txt")
    with open(lista_sug, "w") as fh:
        fh.write("# Suggestions from the viral database audit - NOT applied automatically.\n")
        fh.write("# Review them before passing the file to --lista-remocao-manual in main.py.\n")
        for d in sugeridos:
            fh.write(f"{d['accession']}\t# {d['verdict']} | {d['cellular_pident']}% com "
                     f"{d['best_organism']} | {d['viral_product'][:50]}\n")

    total_viral = len(info_viral)
    log(f"{total_viral} viral proteins | {len(hits)} with some cellular hit")
    for v, n in Counter(d["verdict"] for d in linhas).most_common():
        log(f"  {n:6}  {v}")
    log(f"table: {destino}")
    log(f"removal suggestions ({len(sugeridos)}): {lista_sug}")
    return linhas


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k, v in PADROES.items():
        ap.add_argument(f"--{k}", dest=k, default=v)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--hits-por-protein", type=int, default=50)
    ap.add_argument("--min-pident", type=float, default=35.0,
                    help="Minimum identity to a cellular protein for the hit to count as strong")
    ap.add_argument("--min-qcov", type=float, default=50.0,
                    help="Minimum coverage of the viral protein in the cellular alignment")
    ap.add_argument("--min-phyla", type=int, default=3,
                    help="Number of distinct cellular phyla that marks an ancient cellular gene")
    ap.add_argument("--filtrar-genoma-alvo", action="store_true",
                    help="Keep only the viruses of interest (RNA without RT + ssDNA) from the viral FASTA "
                         "before auditing; use it when pointing at the full viral database")
    ap.add_argument("--apenas-novas", action="store_true",
                    help="Audit only the accessions missing from the previous audit")
    ap.add_argument("--refazer", action="store_true", help="Redo the alignment even if it already exists")
    auditar(ap.parse_args())


if __name__ == "__main__":
    main()
