#!/usr/bin/env python3
"""
Where each candidate EVE falls on the viral protein it resembles (the "protein map").

The DIAMOND screening keeps the best viral hits of every ORF, but not the coordinates on
the subject protein. Here every candidate peptide is realigned (local Smith-Waterman,
BLOSUM62) to a reference viral protein, which gives the stretch of the protein it covers.

Two ways of choosing the reference:

  shared     inside each viral family, the reference proteins are chosen greedily: the
             protein found among the DIAMOND hits of the most candidates comes first, and
             takes all of them; then the next one, for the candidates left. EVEs of one
             family therefore share an axis, as in the figures of EVE papers, and the
             fragments of one broken gene line up along it.
  best       every candidate on its own best DIAMOND hit.

Candidates supported only by the viral HMM have no DIAMOND hit and are listed apart.

The reference sequences come from the DIAMOND database of the run (diamond getseq), read
once and cached with the alignments in 13_curadoria/protein_map_cache.json.

Usage:
  python -m evetool.validation.protein_map -o results_eve                      # summary on the terminal
  python -m evetool.validation.protein_map -o results_eve --tsv map.tsv        # one alignment per line
  python -m evetool.validation.protein_map -o results_eve --dmnd-db mydb.dmnd --forcar
"""
import os
import csv
import json
import glob
import shutil
import argparse
import subprocess
from collections import OrderedDict, defaultdict

from evetool.annotation.cdsearch import ler_fasta, _assinatura
from evetool import config

DMND_PADRAO = config.diamond_path("mydb.dmnd")
CACHE = os.path.join("13_curadoria", "protein_map_cache.json")
MAX_HITS = 25          # DIAMOND hits per ORF considered when choosing a shared reference
VERSAO_CACHE = 2


def log(msg):
    print(f"[protein-map] {msg}", flush=True)


def _pastas(outdir: str) -> "OrderedDict[str, str]":
    pastas = OrderedDict()
    especies_dir = os.path.join(outdir, "especies")
    if os.path.isdir(especies_dir):
        for p in sorted(glob.glob(os.path.join(especies_dir, "*"))):
            if os.path.isdir(p):
                pastas[os.path.basename(p)] = p
    else:
        pastas[""] = outdir
    return pastas


def _titulo(stitle: str) -> dict:
    """'acc |product [virus]|genome|species|genus|family|length' -> its parts."""
    p = [x.strip() for x in stitle.split("|")]
    prod_virus = p[1] if len(p) > 1 else stitle
    produto, virus = prod_virus, ""
    if prod_virus.endswith("]") and "[" in prod_virus:
        k = prod_virus.rfind("[")
        produto, virus = prod_virus[:k].strip(), prod_virus[k + 1:-1].strip()
    return {"product": produto, "virus": virus,
            "genome_type": p[2] if len(p) > 2 else "",
            # acc|product [virus]|genome|species|genus|family|length
            "family": p[5] if len(p) > 5 else ""}


# ---------------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------------
def banco_da_corrida(outdir: str, dmnd_db: str = None) -> str:
    """DIAMOND database used by the run: argument, pipeline_info.json, pipeline.log, default."""
    if dmnd_db:
        return dmnd_db
    for pasta in [outdir, *_pastas(outdir).values()]:
        caminho = os.path.join(pasta, "pipeline_info.json")
        if os.path.exists(caminho):
            try:
                with open(caminho, encoding="utf-8") as fh:
                    db = json.load(fh).get("diamond_db")
                if db and os.path.exists(db):
                    return db
            except Exception:
                pass
    caminho_log = os.path.join(outdir, "pipeline.log")
    if os.path.exists(caminho_log):
        with open(caminho_log, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith("#") and "--dmnd-db" in linha:
                    p = linha.split()
                    db = p[p.index("--dmnd-db") + 1]
                    if os.path.exists(db):
                        return db
                if not linha.startswith("#"):
                    break
    return DMND_PADRAO


def buscar_sequencias(db: str, ids: set) -> dict:
    """Protein sequences of `ids` from a .dmnd (diamond getseq) or from FASTA file(s)."""
    falta, seqs = set(ids), {}
    if not falta:
        return seqs
    if os.path.isdir(db):
        arquivos = sorted(f for f in glob.glob(os.path.join(db, "*"))
                          if f.lower().endswith((".fa", ".faa", ".fasta", ".fas")))
    elif db.endswith(".dmnd"):
        arquivos = None
    else:
        arquivos = [db]

    def ler(linhas):
        nome = None
        for linha in linhas:
            if linha.startswith(">"):
                if len(seqs) == len(ids):
                    break
                nome = linha[1:].split()[0] if linha[1:].split() else ""
                nome = nome if nome in falta else None
                if nome:
                    seqs[nome] = []
            elif nome:
                seqs[nome].append(linha.strip())

    if arquivos is None:
        if not shutil.which("diamond"):
            raise RuntimeError("diamond not found on the PATH")
        if not os.path.exists(db):
            raise RuntimeError(f"DIAMOND database not found: {db}")
        proc = subprocess.Popen(["diamond", "getseq", "--db", db], stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
        try:
            ler(proc.stdout)
        finally:
            proc.stdout.close()
            proc.kill()
            proc.wait()
    else:
        for f in arquivos:
            with open(f, encoding="utf-8", errors="replace") as fh:
                ler(fh)
    return {k: "".join(v).replace("*", "").upper() for k, v in seqs.items()}


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------
_ALINHADOR = None


def _alinhador():
    global _ALINHADOR
    if _ALINHADOR is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        a = PairwiseAligner()
        a.mode = "local"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -11, -1
        _ALINHADOR = a
    return _ALINHADOR


def _seis_quadros(nt: str) -> list:
    """(frame, translation) of the six reading frames of a nucleotide stretch."""
    from Bio.Seq import Seq
    nt = "".join(c for c in nt.upper() if c.isalpha())
    rc = str(Seq(nt).reverse_complement())
    saida = []
    for sinal, s in (("+", nt), ("-", rc)):
        for f in range(3):
            trecho = s[f:len(s) - ((len(s) - f) % 3)]
            if len(trecho) >= 3:
                saida.append((f"{sinal}{f + 1}", str(Seq(trecho).translate())))
    return saida


def alinhar_trecho(qseq_nt: str, pep_orf: str, prot: str) -> dict:
    """Aligns the six translations of the stretch DIAMOND aligned (blastx qseq) and keeps
    the best one. DIAMOND translates the ORF in every frame, so the viral signal can sit on
    the opposite strand or in another frame than the ORF itself; `orf_frame` says whether
    the aligned frame is the ORF's own reading frame (positions on the ORF only then)."""
    melhor = None
    for quadro, trad in _seis_quadros(qseq_nt):
        a = alinhar(trad, prot)
        if a and (melhor is None or a["score"] > melhor["score"]):
            melhor = {**a, "frame": quadro, "_trad": trad}
    if not melhor:
        return {}
    trad = melhor.pop("_trad")
    # the translated stretch inside the ORF peptide: same frame as the ORF
    nucleo = trad[melhor["qfrom"] - 1:melhor["qto"]].replace("*", "")
    k = pep_orf.find(nucleo) if nucleo else -1
    melhor["orf_frame"] = k >= 0
    if k >= 0:
        melhor["qfrom"], melhor["qto"] = k + 1, k + len(nucleo)
    else:
        melhor.pop("qfrom"), melhor.pop("qto")
    melhor["stops"] = trad[melhor.get("qfrom", 1) - 1:].count("*") if k >= 0 else nucleo.count("*")
    return melhor


def alinhar(pep: str, prot: str) -> dict:
    """Best local alignment of the peptide on the protein: positions (1-based) and identity."""
    pep = "".join(c for c in pep.upper() if c.isalpha())
    aln = _alinhador().align(pep, prot)[0]
    blocos_q, blocos_s = aln.aligned
    if not len(blocos_q):
        return {}
    iguais = colunas = 0
    for (q0, q1), (s0, s1) in zip(blocos_q, blocos_s):
        iguais += sum(1 for x, y in zip(pep[q0:q1], prot[s0:s1]) if x == y)
        colunas += q1 - q0
    # gaps count as columns too, as in BLAST identity
    q_ini, q_fim = int(blocos_q[0][0]), int(blocos_q[-1][1])
    s_ini, s_fim = int(blocos_s[0][0]), int(blocos_s[-1][1])
    colunas += (q_fim - q_ini - sum(b - a for a, b in blocos_q)) + (s_fim - s_ini - sum(b - a for a, b in blocos_s))
    return {"from": s_ini + 1, "to": s_fim, "qfrom": q_ini + 1, "qto": q_fim,
            "pident": round(100 * iguais / max(1, colunas), 1), "score": float(aln.score)}


# ---------------------------------------------------------------------------
# Map
# ---------------------------------------------------------------------------
def _hits_virais(pasta: str, orfs: set) -> dict:
    """ORF -> its DIAMOND viral hits, best first (one per subject)."""
    hits = defaultdict(OrderedDict)
    for nome in ("triagem_viral_diamond_refseq_filtrado.tsv", "triagem_viral_diamond.tsv"):
        caminho = os.path.join(pasta, "02_alinhamento_diamond", nome)
        if not os.path.exists(caminho):
            continue
        with open(caminho, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.split("\t", 1)[0] not in orfs:
                    continue
                p = linha.rstrip("\n").split("\t")
                if len(p) < 8 or p[2] in hits[p[0]]:
                    continue
                try:
                    hits[p[0]][p[2]] = {"id": p[2], "length": int(p[3]), "bitscore": float(p[7]),
                                        "pident": float(p[5]), **_titulo(p[4]),
                                        # aligned stretch of the ORF, in nucleotides (blastx qseq)
                                        "qseq": p[10] if len(p) > 10 else ""}
                except ValueError:
                    continue
    return {o: sorted(h.values(), key=lambda x: -x["bitscore"])[:MAX_HITS] for o, h in hits.items()}


def _ler_cache(outdir: str) -> dict:
    caminho = os.path.join(outdir, CACHE)
    try:
        with open(caminho, encoding="utf-8") as fh:
            c = json.load(fh)
        if c.get("version") == VERSAO_CACHE:
            return c
    except Exception:
        pass
    return {"version": VERSAO_CACHE, "lengths": {}, "alignments": {}}


def _gravar_cache(outdir: str, cache: dict):
    caminho = os.path.join(outdir, CACHE)
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    tmp = caminho + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(cache, fh)
    os.replace(tmp, caminho)


def _chave_resultado(outdir: str, candidatos: list) -> str:
    """Changes when the candidates, their families or the DIAMOND tables change."""
    import hashlib
    partes = sorted(f"{c['query']}|{c['species']}|{c['orf']}|{c.get('family', '')}" for c in candidatos)
    for pasta in _pastas(outdir).values():
        for nome in ("triagem_viral_diamond_refseq_filtrado.tsv", "triagem_viral_diamond.tsv"):
            f = os.path.join(pasta, "02_alinhamento_diamond", nome)
            if os.path.exists(f):
                partes.append(f"{f}:{os.path.getsize(f)}:{int(os.path.getmtime(f))}")
    return hashlib.md5("\n".join(partes).encode()).hexdigest()


def mapa(outdir: str, dmnd_db: str = None, forcar: bool = False, candidatos: list = None,
         so_cache: bool = False) -> dict:
    """Alignment of every candidate on a shared and on its best reference protein.

    `candidatos`: dicts with query/species/orf/family (the curation rows); when omitted,
    read from eves_classified.tsv. `so_cache`: return the last map if still valid, else
    only {"pending": True} (the curation page asks for the full map separately)."""
    if candidatos is None:
        candidatos = []
        with open(os.path.join(outdir, "eves_classified.tsv"), encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                candidatos.append({"query": r.get("query") or r.get("orf", ""), "species": r.get("species", ""),
                                   "orf": r.get("orf", ""), "family": r.get("viral_family", "")})
    chave_res = _chave_resultado(outdir, candidatos)
    cache = {"version": VERSAO_CACHE, "lengths": {}, "alignments": {}} if forcar else _ler_cache(outdir)
    if not forcar and cache.get("result_key") == chave_res and cache.get("result"):
        return cache["result"]
    if so_cache:
        return {"pending": True, "items": [], "references": {}, "no_alignment": []}
    pastas = _pastas(outdir)
    pasta_de = (lambda sp: pastas.get(sp)) if len(pastas) > 1 else (lambda sp: next(iter(pastas.values())))

    # Peptides and DIAMOND hits, per species folder
    peps, hits = {}, {}
    por_pasta = defaultdict(set)
    for c in candidatos:
        if pasta_de(c["species"]):
            por_pasta[pasta_de(c["species"])].add(c["orf"])
    for pasta, orfs in por_pasta.items():
        fasta = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
        seqs = ler_fasta(fasta) if os.path.exists(fasta) else {}
        h = _hits_virais(pasta, orfs)
        for o in orfs:
            peps[(pasta, o)] = seqs.get(o, "")
            hits[(pasta, o)] = h.get(o, [])

    # Shared references: greedy, inside each family (the curated family when edited by hand)
    def familia(c, hs):
        f = c.get("family") or ""
        return f if f and f != "n/a" else (hs[0]["family"] if hs and hs[0]["family"] else "Unassigned")

    grupos = defaultdict(list)
    for c in candidatos:
        k = (pasta_de(c["species"]), c["orf"])
        if hits.get(k):
            grupos[familia(c, hits[k])].append(c)
    compartilhada = {}
    for fam, membros in grupos.items():
        restantes = list(membros)
        while restantes:
            conta, peso, info = defaultdict(int), defaultdict(float), {}
            for c in restantes:
                for h in hits[(pasta_de(c["species"]), c["orf"])]:
                    conta[h["id"]] += 1
                    peso[h["id"]] += h["bitscore"]
                    info[h["id"]] = h
            ref = max(conta, key=lambda i: (conta[i], peso[i], info[i]["length"]))
            ficam = []
            for c in restantes:
                if any(h["id"] == ref for h in hits[(pasta_de(c["species"]), c["orf"])]):
                    compartilhada[c["query"]] = info[ref]
                else:
                    ficam.append(c)
            restantes = ficam

    # Pairs to align: (peptide, reference) for both modes
    pares = {}
    for c in candidatos:
        k = (pasta_de(c["species"]), c["orf"])
        if not hits.get(k) or not peps.get(k):
            continue
        for ref in {hits[k][0]["id"]: hits[k][0], compartilhada[c["query"]]["id"]: compartilhada[c["query"]]}.values():
            h = next(x for x in hits[k] if x["id"] == ref["id"])
            pares[(c["query"], ref["id"])] = (peps[k], h.get("qseq", ""))

    chave = lambda pep, qseq, rid: f"{_assinatura(pep + '|' + qseq)}|{rid}"
    faltam = {rid for (_q, rid), (pep, qs) in pares.items() if chave(pep, qs, rid) not in cache["alignments"]}
    erro_db = None
    if faltam:
        db = banco_da_corrida(outdir, dmnd_db)
        log(f"{len(faltam)} reference protein(s) to fetch from {db}")
        try:
            seqs_ref = buscar_sequencias(db, faltam)
        except Exception as e:
            seqs_ref, erro_db = {}, str(e)
        if len(seqs_ref) < len(faltam):
            erro_db = erro_db or (f"{len(faltam) - len(seqs_ref)} reference protein(s) not found in {db}")
        novos = 0
        for (q, rid), (pep, qs) in pares.items():
            ck = chave(pep, qs, rid)
            if ck not in cache["alignments"] and rid in seqs_ref:
                if qs:
                    cache["alignments"][ck] = alinhar_trecho(qs, pep, seqs_ref[rid])
                else:   # older outputs without qseq: the ORF peptide itself
                    cache["alignments"][ck] = {**alinhar(pep, seqs_ref[rid]), "frame": "ORF", "orf_frame": True}
                cache["lengths"][rid] = len(seqs_ref[rid])
                novos += 1
        if novos:
            log(f"{novos} alignment(s) computed")

    itens, sem_hit = [], []
    for c in candidatos:
        k = (pasta_de(c["species"]), c["orf"])
        hs = hits.get(k) or []
        if not hs:
            sem_hit.append(c["query"])
            continue
        pep = peps.get(k, "")
        lin = {"query": c["query"], "family": familia(c, hs), "peptide_length": len(pep),
               "best_ref": hs[0]["id"]}
        for modo, ref in (("shared", compartilhada[c["query"]]), ("best", hs[0])):
            dh = next((h for h in hs if h["id"] == ref["id"]), {})
            aln = cache["alignments"].get(chave(pep, dh.get("qseq", ""), ref["id"])) if pep else None
            lin[modo] = {"ref": ref["id"], "diamond_pident": dh.get("pident"),
                         "diamond_bitscore": dh.get("bitscore"), **(aln or {})} if aln else {"ref": ref["id"]}
        itens.append(lin)

    refs = {}
    for c in candidatos:
        for h in hits.get((pasta_de(c["species"]), c["orf"])) or []:
            if h["id"] in cache["lengths"] and h["id"] not in refs:
                refs[h["id"]] = {k: v for k, v in h.items() if k not in ("qseq", "bitscore", "pident")}
                refs[h["id"]]["length"] = cache["lengths"][h["id"]]
    resultado = {"items": itens, "references": refs, "no_alignment": sem_hit,
                 "error": erro_db, "db": banco_da_corrida(outdir, dmnd_db)}
    if not erro_db:   # an incomplete map is not kept: the next call tries the database again
        cache["result_key"], cache["result"] = chave_res, resultado
    _gravar_cache(outdir, cache)
    return resultado


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--dmnd-db", help="DIAMOND database of the run (default: pipeline_info.json, then pipeline.log)")
    ap.add_argument("--modo", choices=["shared", "best"], default="shared", help="Reference choice")
    ap.add_argument("--tsv", help="Write one alignment per line to this file")
    ap.add_argument("--forcar", action="store_true", help="Realign everything, ignoring the cache")
    a = ap.parse_args()
    r = mapa(a.outdir, dmnd_db=a.dmnd_db, forcar=a.forcar)
    if r["error"]:
        log(f"warning: {r['error']}")
    linhas = []
    for i in r["items"]:
        m = i[a.modo]
        ref = r["references"].get(m["ref"], {})
        linhas.append({"query": i["query"], "family": i["family"], "reference": m["ref"],
                       "product": ref.get("product", ""), "virus": ref.get("virus", ""),
                       "ref_length": ref.get("length", ""), "from": m.get("from", ""), "to": m.get("to", ""),
                       "orf_from": m.get("qfrom", ""), "orf_to": m.get("qto", ""),
                       "pident": m.get("pident", ""), "frame": m.get("frame", ""),
                       "orf_frame": m.get("orf_frame", ""), "is_best_hit": m["ref"] == i["best_ref"]})
    if a.tsv:
        with open(a.tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(linhas[0]) if linhas else ["query"], delimiter="\t")
            w.writeheader()
            w.writerows(linhas)
        log(f"written: {a.tsv}")
    for l in linhas:
        print(f"{l['family'][:18]:18} {l['query'][:30]:30} {l['reference']:16} {str(l['from']):>5}-{str(l['to']):<5}"
              f" / {l['ref_length']:<5} {l['pident']:>5}%  {l['frame']:3} {'' if l['orf_frame'] else 'OTHER-FRAME'}"
              f"  {l['product'][:30]}")
    if r["no_alignment"]:
        print(f"\n{len(r['no_alignment'])} candidate(s) without a DIAMOND hit (HMM only): "
              + ", ".join(r["no_alignment"]))


if __name__ == "__main__":
    main()
