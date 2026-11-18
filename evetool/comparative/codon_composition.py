#!/usr/bin/env python3
"""
Composition of the EVEs next to circulating viruses: GC, codon usage and dinucleotides.

The researcher adds coding sequences (CDS) of circulating viruses, and optionally host genes
as a reference. Every approved EVE is then measured on the same features as those sequences:

  composition    %A %C %G %T, %GC, GC1, GC2, GC3, GC12 and GC3s (GC at the synonymous third
                 positions: Met, Trp and stops left out)
  codon usage    RSCU of the 59 codons of the amino acids with synonyms, and the effective
                 number of codons (ENC, Wright 1990, codonW conventions)
  dinucleotides  relative abundance rho = f(XY) / (f(X) f(Y)) of the 16 dinucleotides (Karlin &
                 Mrazek 1997): over the whole coding sequence, and at each codon position,
                 1-2 and 2-3 inside the codon, 3-1 across the boundary of two codons

An EVE is read in the frame of the viral protein, not in the frame of its ORF. The stretch that
DIAMOND aligned to the viral protein (blastx qseq) is translated in the six frames, and the frame
that aligns best to the reference protein (protein_map.py) gives the codons. In many EVEs the
longest ORF is another frame of the insertion, and its codons say nothing about the virus.
Overlapping ORFs of one locus are one insertion: the locus is measured once.

The features that best separate the groups (Kruskal-Wallis, a feature too correlated with one
already chosen is skipped) go into a PCA and a t-SNE. The page colours the points by viral
family, draws EVEs, circulating viruses and host genes with different symbols, and fills the
points that carry a viral HMM domain (the same calibrated HMMs for EVEs and circulating CDS).

Files, under 13_curadoria/circulantes/:
  sets.json                 the sequence sets added
  set_<id>.fasta            the sequences of a set, as uploaded
  set_<id>.json             per sequence: family, product, HMM domain, checks and features
  analysis.json             the last analysis (points, ranking, PCA, t-SNE), read by the page
  composition_features.tsv  every feature of every sequence of the last analysis
  composition_pca.tsv       PCA and t-SNE coordinates of the last analysis

Usage:
  python -m evetool.comparative.codon_composition -o results --add rhabdo_cds.fasta --name Rhabdoviridae
  python -m evetool.comparative.codon_composition -o results --add host_cds.fna --kind host --name L_longipalpis
  python -m evetool.comparative.codon_composition -o results --run --top 6
  python -m evetool.comparative.codon_composition -o results --list
"""
import os
import re
import csv
import json
import math
import time
import random
import shutil
import argparse
import datetime
import subprocess
from collections import defaultdict, Counter

SUBPASTA = os.path.join("13_curadoria", "circulantes")
VERSAO_FEATURES = 1
ORIGENS = ("EVE", "circulating", "host")

# ---------------------------------------------------------------------------
# Genetic code and feature names
# ---------------------------------------------------------------------------
_B = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODIGO = {a + b + c: _AAS[16 * i + 4 * j + k]
          for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
STOPS = {c for c, a in CODIGO.items() if a == "*"}
SINONIMOS = defaultdict(list)
for _c, _a in sorted(CODIGO.items()):
    if _a != "*":
        SINONIMOS[_a].append(_c)
# amino acids with synonymous codons, in one-letter order; Met and Trp have one codon each
AA_RSCU = [a for a in "ACDEFGHIKLNPQRSTVY" if len(SINONIMOS[a]) > 1]
CODONS_RSCU = [c for a in AA_RSCU for c in SINONIMOS[a]]           # 59 codons
DINUCS = [x + y for x in "ACGT" for y in "ACGT"]
POSICOES = {"": "all", "12": "1-2", "23": "2-3", "31": "3-1"}
ACGT = set("ACGT")

COMPOSICAO = ["GC", "GC1", "GC2", "GC3", "GC12", "GC3s", "A", "C", "G", "T"]
FEATURES = (COMPOSICAO + ["ENC"] + ["RSCU_" + c for c in CODONS_RSCU]
            + [f"rho{p}_{d}" for p in POSICOES for d in DINUCS])
GRUPOS_FEATURES = {
    "composition": COMPOSICAO,
    "codon": ["ENC"] + ["RSCU_" + c for c in CODONS_RSCU],
    "dinucleotide": [f"rho{p}_{d}" for p in POSICOES for d in DINUCS],
}


def log(msg):
    print(f"[composition] {msg}", flush=True)


def pasta(outdir: str) -> str:
    return os.path.join(outdir, SUBPASTA)


def limpar(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper().replace("U", "T")


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def traduzir(nt: str) -> str:
    s = limpar(nt)
    return "".join(CODIGO.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


# ---------------------------------------------------------------------------
# Features of one coding sequence
# ---------------------------------------------------------------------------
def enc_wright(cont: Counter):
    """Effective number of codons (Wright 1990): 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6, at most 61.
    F of an amino acid is its codon homozygosity, (n sum p^2 - 1)/(n - 1), needing n >= 2; a
    missing F3 (Ile) is the mean of F2 and F4, any other missing class leaves ENC undefined."""
    F = defaultdict(list)
    for a, cods in SINONIMOS.items():
        k = len(cods)
        if k < 2:
            continue
        n = sum(cont[c] for c in cods)
        if n < 2:
            continue
        soma = sum((cont[c] / n) ** 2 for c in cods)
        F[k].append((n * soma - 1) / (n - 1))
    media = {k: sum(v) / len(v) for k, v in F.items() if v}
    if 3 not in media and 2 in media and 4 in media:
        media[3] = (media[2] + media[4]) / 2
    if any(media.get(k, 0) <= 0 for k in (2, 3, 4, 6)):
        return None
    return min(61.0, 2 + 9 / media[2] + 1 / media[3] + 5 / media[4] + 3 / media[6])


def _rho(pares) -> dict:
    """rho of the 16 dinucleotides from a list of (X, Y) pairs; X and Y frequencies are the
    margins of the same pairs (for 1-2 and 2-3, the base frequencies at those codon positions)."""
    pares = [p for p in pares if p[0] in ACGT and p[1] in ACGT]
    n = len(pares)
    if n < 10:
        return {d: None for d in DINUCS}
    xy = Counter(a + b for a, b in pares)
    fx = Counter(a for a, _ in pares)
    fy = Counter(b for _, b in pares)
    out = {}
    for d in DINUCS:
        den = (fx[d[0]] / n) * (fy[d[1]] / n)
        out[d] = (xy[d] / n) / den if den else None
    return out


def caracteristicas(nt: str):
    """(features, info) of a coding sequence read from its first base, or (None, info).
    The terminal stop is dropped; internal stops (EVEs) are counted and left out of the codon
    usage, but stay in the base composition and the dinucleotides; codons with N are skipped."""
    s = limpar(nt)
    s = s[:len(s) - len(s) % 3]
    codons = [s[i:i + 3] for i in range(0, len(s), 3)]
    if codons and codons[-1] in STOPS:
        codons = codons[:-1]
    validos = [c if set(c) <= ACGT else None for c in codons]
    sentido = [c for c in validos if c and c not in STOPS]
    info = {"codons": len(sentido), "stops": sum(1 for c in validos if c in STOPS),
            "ambiguous": sum(1 for c in validos if c is None)}
    if len(sentido) < 10:
        return None, info
    out = {}
    nts = "".join(c for c in validos if c)
    for b in "ACGT":
        out[b] = 100 * nts.count(b) / len(nts)
    out["GC"] = out["G"] + out["C"]
    for k in range(3):
        out[f"GC{k + 1}"] = 100 * sum(c[k] in "GC" for c in sentido) / len(sentido)
    out["GC12"] = (out["GC1"] + out["GC2"]) / 2
    terc = [c[2] for c in sentido if len(SINONIMOS[CODIGO[c]]) > 1]
    out["GC3s"] = 100 * sum(x in "GC" for x in terc) / len(terc) if terc else None

    cont = Counter(sentido)
    for a in AA_RSCU:
        cods = SINONIMOS[a]
        n = sum(cont[c] for c in cods)
        for c in cods:   # an amino acid absent from the sequence has no RSCU (None)
            out["RSCU_" + c] = cont[c] * len(cods) / n if n else None
    out["ENC"] = enc_wright(cont)

    seq_frame = "".join(c if c else "NNN" for c in validos)
    for pos, pares in (("", zip(seq_frame, seq_frame[1:])),
                       ("12", ((c[0], c[1]) for c in validos if c)),
                       ("23", ((c[1], c[2]) for c in validos if c)),
                       ("31", ((a[2], b[0]) for a, b in zip(validos, validos[1:]) if a and b))):
        for d, v in _rho(list(pares)).items():
            out[f"rho{pos}_{d}"] = v
    return {k: (round(v, 4) if v is not None else None) for k, v in out.items()}, info


# ---------------------------------------------------------------------------
# Sequence sets (circulating viruses, host genes)
# ---------------------------------------------------------------------------
RX_FAMILIA = re.compile(r"\b([A-Z][a-z]+viridae)\b")
RX_TAG = re.compile(r"\[?\b(family|protein|product|gene|organism|species|virus)\s*=\s*([^\]|;]+)\]?", re.I)


def ler_fasta_texto(texto: str):
    cab, partes = None, []
    for linha in (texto or "").splitlines():
        if linha.startswith(">"):
            if cab is not None:
                yield cab, "".join(partes)
            cab, partes = linha[1:].strip(), []
        elif cab is not None:
            partes.append(linha.strip())
    if cab is not None:
        yield cab, "".join(partes)


def ler_cabecalho(cab: str) -> dict:
    """Family, species and product from a FASTA header. Understands the NCBI Virus custom
    definition line (acc |title|molecule|species|genus|family|length, the same as the
    RefSeq viral proteins of the tool), NCBI CDS tags ([protein=...] [gene=...]), key=value
    tags (family=..., species=...) and, as a last resort, any word ending in -viridae."""
    info = {"family": "", "species": "", "product": "", "accession": cab.split()[0] if cab.split() else ""}
    campos = [x.strip() for x in cab.split("|")]
    if len(campos) >= 6 and (campos[5].endswith("viridae") or re.search(r"(ss|ds)(RNA|DNA)", campos[2] or "")):
        titulo = campos[1]
        if titulo.endswith("]") and "[" in titulo:
            titulo = titulo[:titulo.rfind("[")].strip()
        info.update(product=titulo, species=campos[3], family=campos[5])
    tags = {k.lower(): v.strip() for k, v in RX_TAG.findall(cab)}
    info["family"] = info["family"] or tags.get("family", "")
    info["species"] = info["species"] or tags.get("species") or tags.get("organism") or tags.get("virus") or ""
    info["product"] = info["product"] or tags.get("protein") or tags.get("product") or tags.get("gene") or ""
    if not info["family"]:
        m = RX_FAMILIA.search(cab)
        info["family"] = m.group(1) if m else ""
    if not info["species"]:
        m = re.search(r"\[([^\[\]=]+)\]\s*$", cab)
        info["species"] = m.group(1).strip() if m else ""
    if not info["product"] and len(cab.split(None, 1)) > 1 and "|" not in cab and "[" not in cab:
        info["product"] = cab.split(None, 1)[1].strip()
    return info


def _ler_json(caminho, padrao):
    try:
        with open(caminho, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return padrao


def _gravar_json(caminho, dados):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(dados, fh, ensure_ascii=False)
    os.replace(caminho + ".tmp", caminho)


def listar_conjuntos(outdir: str) -> list:
    return _ler_json(os.path.join(pasta(outdir), "sets.json"), [])


def _bancos(outdir: str) -> dict:
    """HMM and DIAMOND databases of the run (pipeline_info.json), else the installed ones."""
    from evetool import config
    from evetool.curation.curation import pastas_especies, _info_pipeline
    hmm, dmnd = None, None
    for p in [outdir, *pastas_especies(outdir).values()]:
        info = _info_pipeline(p)
        hmm = hmm or (info.get("hmm_db") if info.get("hmm_db") and os.path.exists(info["hmm_db"]) else None)
        dmnd = dmnd or (info.get("diamond_db") if info.get("diamond_db") and os.path.exists(info["diamond_db"]) else None)
    return {"hmm": hmm or config.HMM_DB, "dmnd": dmnd or config.VIRAL_DMND}


def _familias_diamond(proteinas: dict, db: str, threads: int, trabalho: str) -> dict:
    """{id: hit} from DIAMOND blastp of the translated CDS on the viral database of the run:
    family, species and product come from the RefSeq title of the best hit."""
    from evetool.validation.protein_map import _titulo
    if not proteinas:
        return {}
    if not shutil.which("diamond"):
        raise RuntimeError("diamond not found on the PATH (needed to assign the family automatically)")
    if not db or not os.path.exists(db):
        raise RuntimeError(f"viral DIAMOND database not found: {db}")
    faa, tsv = os.path.join(trabalho, "familia_query.faa"), os.path.join(trabalho, "familia_hits.tsv")
    with open(faa, "w") as fh:
        for k, p in proteinas.items():
            fh.write(f">{k}\n{p}\n")
    subprocess.run(["diamond", "blastp", "--db", db, "--query", faa, "--out", tsv, "--threads", str(max(1, threads)),
                    "--outfmt", "6", "qseqid", "sseqid", "stitle", "pident", "bitscore", "evalue",
                    "--max-target-seqs", "1", "--evalue", "1e-5", "--quiet"], check=True)
    hits = {}
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < 6 or p[0] in hits:
                continue
            t = _titulo(p[2])
            hits[p[0]] = {"ref": p[1], "family": t["family"], "species": t["virus"], "product": t["product"],
                          "pident": float(p[3]), "bitscore": float(p[4])}
    for f in (faa, tsv):
        if os.path.exists(f):
            os.remove(f)
    return hits


def _dominios_hmm(proteinas: dict, hmm_db: str, threads: int, trabalho: str) -> dict:
    """{id: {"n": domains, "model": best, "score": bits, "label": family · genome}} with the
    calibrated viral HMMs (TC cut-offs), the same search that annotates the EVEs."""
    from evetool.visualization import dashboard
    if not proteinas:
        return {}
    faa = os.path.join(trabalho, "dominio_query.faa")
    domtbl = os.path.join(trabalho, "dominio.domtbl")
    with open(faa, "w") as fh:
        for k, p in proteinas.items():
            fh.write(f">{k}\n{p}\n")
    if os.path.exists(domtbl):
        os.remove(domtbl)
    saida = dashboard.hmmsearch_todas([faa], domtbl, hmm_db, threads)
    hits = defaultdict(dict)
    dashboard.ler_domtbl(saida, hits)
    rotulos = dashboard.anotacao_modelos(os.path.join(os.path.dirname(hmm_db), "combined_viral_annotations.tsv"))
    out = {}
    for k, modelos in hits.items():
        melhor = max(modelos, key=modelos.get)
        out[k] = {"n": len(modelos), "model": melhor, "score": round(modelos[melhor], 1),
                  "label": rotulos.get(melhor, "")}
    for f in (faa, domtbl):
        if os.path.exists(f):
            os.remove(f)
    return out if saida else {"_error": "HMM search unavailable (pyhmmer or the HMM database missing)"}


def adicionar_conjunto(outdir: str, texto_fasta: str, nome: str = "", tipo: str = "circulating",
                       familia_padrao: str = "", familia_auto: bool = True, threads: int = 4) -> dict:
    """Store a FASTA of coding sequences as a new set and annotate each sequence."""
    from evetool.validation.classify_eves import classe_da_proteina
    if tipo not in ("circulating", "host"):
        raise ValueError("kind must be 'circulating' or 'host'")
    registros = list(ler_fasta_texto(texto_fasta))
    if not registros:
        raise ValueError("no FASTA record found (headers must start with '>')")
    base = pasta(outdir)
    os.makedirs(base, exist_ok=True)
    conjuntos = listar_conjuntos(outdir)
    nome = (nome or "").strip() or f"set {len(conjuntos) + 1}"
    raiz = re.sub(r"[^A-Za-z0-9]+", "_", nome).strip("_")[:40] or "set"
    ids = {c["id"] for c in conjuntos}
    sid, k = raiz, 2
    while sid in ids:
        sid, k = f"{raiz}_{k}", k + 1
    log(f"set {sid}: {len(registros)} sequence(s), kind {tipo}")

    seqs, avisos = [], Counter()
    for n, (cab, seq) in enumerate(registros):
        s = limpar(seq)
        cab_info = ler_cabecalho(cab) if tipo == "circulating" else {"family": "", "species": nome, "product": "",
                                                                       "accession": cab.split()[0] if cab.split() else ""}
        r = {"key": f"q{n:06d}", "id": cab_info["accession"] or f"seq{n + 1}", "header": cab[:300],
             "length": len(s), "family": cab_info["family"], "family_source": "header" if cab_info["family"] else "",
             "species": cab_info["species"], "product": cab_info["product"], "problems": []}
        if len(s) % 3:
            r["problems"].append("length not a multiple of 3 (trailing bases ignored)")
        if s[:3] != "ATG":
            avisos["no ATG start"] += 1
        feats, info = caracteristicas(s)
        r.update(info)
        if info["stops"]:
            r["problems"].append(f"{info['stops']} internal stop codon(s): wrong frame or not a CDS?")
        if feats is None:
            r["problems"].append("fewer than 10 sense codons")
        r["features"] = feats
        r["protein"] = traduzir(s).rstrip("*")
        seqs.append(r)
        for p in r["problems"]:
            avisos[p.split(":")[0]] += 1

    trabalho = os.path.join(base, "tmp_" + sid)
    os.makedirs(trabalho, exist_ok=True)
    bancos = _bancos(outdir)
    msgs = []
    try:
        if tipo == "circulating":
            for r in seqs:
                if not r["family"] and familia_padrao:
                    r["family"], r["family_source"] = familia_padrao.strip(), "form"
            faltam = {r["key"]: r["protein"].replace("*", "") for r in seqs
                      if r["protein"] and (familia_auto and not r["family"] or not r["product"])}
            if faltam and familia_auto:
                try:
                    log(f"DIAMOND blastp of {len(faltam)} CDS without family or product on {bancos['dmnd']}")
                    hits = _familias_diamond(faltam, bancos["dmnd"], threads, trabalho)
                    for r in seqs:
                        h = hits.get(r["key"])
                        if not h:
                            continue
                        r["diamond"] = h
                        if not r["family"] and h["family"]:
                            r["family"], r["family_source"] = h["family"], "DIAMOND"
                        r["species"] = r["species"] or h["species"]
                        r["product"] = r["product"] or h["product"]
                except Exception as erro:
                    msgs.append(f"family by DIAMOND failed: {erro}")
            log(f"viral HMM domains of {len(seqs)} CDS ({bancos['hmm']})")
            dom = _dominios_hmm({r["key"]: r["protein"].replace("*", "") for r in seqs if r["protein"]},
                                bancos["hmm"], threads, trabalho)
            if "_error" in dom:
                msgs.append(dom.pop("_error"))
            for r in seqs:
                d = dom.get(r["key"])
                r["domain"] = d if d else None
                r["has_domain"] = bool(d) if "_error" not in dom else None
        for r in seqs:
            r["family"] = r["family"] or ("" if tipo == "host" else "Unassigned")
            r["protein_class"] = classe_da_proteina(r["product"]) if r["product"] else "n/a"
            r.pop("protein", None)
    finally:
        shutil.rmtree(trabalho, ignore_errors=True)

    with open(os.path.join(base, f"set_{sid}.fasta"), "w", encoding="utf-8") as fh:
        for cab, seq in registros:
            fh.write(f">{cab}\n{limpar(seq)}\n")
    _gravar_json(os.path.join(base, f"set_{sid}.json"),
                 {"version": VERSAO_FEATURES, "id": sid, "sequences": seqs})
    validas = [r for r in seqs if r["features"]]
    resumo = {"id": sid, "name": nome, "kind": tipo, "n": len(seqs), "valid": len(validas),
              "families": dict(Counter(r["family"] for r in validas if r["family"]).most_common()),
              "with_domain": sum(1 for r in validas if r.get("has_domain")),
              "family_sources": dict(Counter(r["family_source"] or "none" for r in validas)),
              "warnings": dict(avisos), "messages": msgs,
              "added": datetime.datetime.now().isoformat(timespec="seconds")}
    conjuntos.append(resumo)
    _gravar_json(os.path.join(base, "sets.json"), conjuntos)
    log(f"set {sid}: {len(validas)}/{len(seqs)} usable, {resumo['with_domain']} with a viral HMM domain")
    return resumo


def remover_conjunto(outdir: str, sid: str) -> dict:
    base = pasta(outdir)
    conjuntos = [c for c in listar_conjuntos(outdir) if c["id"] != sid]
    for ext in (".fasta", ".json"):
        caminho = os.path.join(base, f"set_{sid}{ext}")
        if os.path.exists(caminho):
            os.remove(caminho)
    _gravar_json(os.path.join(base, "sets.json"), conjuntos)
    return {"removed": sid, "sets": len(conjuntos)}


def _sequencias_do_conjunto(outdir: str, c: dict) -> list:
    """Records of a set, recomputing the features when the stored ones are from an older version."""
    base = pasta(outdir)
    dados = _ler_json(os.path.join(base, f"set_{c['id']}.json"), {})
    seqs = dados.get("sequences", [])
    if dados.get("version") != VERSAO_FEATURES:
        nts = {}
        caminho = os.path.join(base, f"set_{c['id']}.fasta")
        if os.path.exists(caminho):
            with open(caminho, encoding="utf-8") as fh:
                nts = {k: s for k, (_cab, s) in enumerate(ler_fasta_texto(fh.read()))}
        for n, r in enumerate(seqs):
            r["features"], info = caracteristicas(nts.get(n, ""))
            r.update(info)
        dados.update(version=VERSAO_FEATURES, sequences=seqs)
        _gravar_json(os.path.join(base, f"set_{c['id']}.json"), dados)
    return seqs


# ---------------------------------------------------------------------------
# Host genes from the genome annotation (GFF3 or GTF)
# ---------------------------------------------------------------------------
NOME_ANOT = "host_annotation.json"


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s or "genome").strip("_")[:60] or "genome"


def _gff_ao_lado(genoma: str):
    """A GFF/GTF next to the genome FASTA: same file stem, or same NCBI accession (GCF_/GCA_)."""
    if not genoma or not os.path.exists(genoma):
        return None
    d, base = os.path.split(os.path.abspath(genoma))
    stem = re.sub(r"\.(fna|fa|fasta|fas)(\.gz)?$", "", base)
    acc = re.match(r"(GC[AF]_\d+\.\d+)", base)
    exts = (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")
    try:
        arquivos = sorted(os.listdir(d))
    except OSError:
        return None
    for f in arquivos:
        if f.endswith(exts) and (f.startswith(stem) or (acc and f.startswith(acc.group(1)))):
            return os.path.join(d, f)
    return None


def _cache_hospedeiro(outdir: str, sp: str) -> str:
    return os.path.join(pasta(outdir), f"host_genes_{_slug(sp)}.json")


def anotacoes_hospedeiro(outdir: str) -> dict:
    """{species: genome, GFF and where the GFF came from, and the genes already extracted}.

    The GFF of each species is, in order: the one set on the page, the third column of
    synteny_genomes.tsv, the GFF of the pipeline run (whole-genome runs only: a run on segments
    has no annotation of its own), or a GFF next to the genome FASTA with the same name or accession."""
    from evetool.curation import curation as cur
    manual = _ler_json(os.path.join(pasta(outdir), NOME_ANOT), {})
    terceira = {}
    cam = os.path.join(outdir, "synteny_genomes.tsv")
    if os.path.exists(cam):
        with open(cam, encoding="utf-8") as fh:
            for linha in fh:
                p = linha.split("#", 1)[0].rstrip("\n").split("\t")
                if len(p) >= 3 and p[0].strip() and p[2].strip():
                    terceira[p[0].strip()] = p[2].strip()
    pastas = cur.pastas_especies(outdir)
    saida = {}
    for sp, g in cur.genomas_para_sintenia(outdir).items():
        info = cur._info_pipeline(pastas.get(sp, outdir))
        candidatos = [(manual.get(sp), "set on the page"), (terceira.get(sp), "synteny_genomes.tsv"),
                      (None if g["segments"] else info.get("gff"), "pipeline run"),
                      (_gff_ao_lado(g["genome"]), "found next to the genome")]
        gff, fonte = next(((c, f) for c, f in candidatos if c and os.path.exists(c)), (None, ""))
        meta = _ler_json(_cache_hospedeiro(outdir, sp)[:-5] + ".meta.json", {})
        atual = meta.get("gff") == gff and meta.get("genome") == g["genome"] and meta.get("version") == VERSAO_FEATURES
        saida[sp] = {"genome": g["genome"], "own": g["own"], "segments": g["segments"], "gff": gff, "source": fonte,
                     "manual": manual.get(sp) or "", "genes": meta.get("stats") if atual else None}
    return saida


def definir_gff(outdir: str, especie: str, caminho: str) -> dict:
    """Set (or clear, with an empty path) the GFF/GTF of one species."""
    caminho = os.path.abspath(os.path.expanduser(caminho.strip())) if caminho and caminho.strip() else ""
    if caminho and not os.path.exists(caminho):
        raise ValueError(f"file not found: {caminho}")
    if caminho and not re.search(r"\.(gff3?|gtf)(\.gz)?$", caminho):
        raise ValueError("the annotation must be a .gff, .gff3 or .gtf file (optionally .gz)")
    arq = os.path.join(pasta(outdir), NOME_ANOT)
    manual = _ler_json(arq, {})
    if caminho:
        manual[especie] = caminho
    else:
        manual.pop(especie, None)
    _gravar_json(arq, manual)
    return {"species": especie, "gff": caminho}


def _atributos(campo: str, gtf: bool) -> dict:
    if gtf:
        return dict(re.findall(r'(\S+) "([^"]*)"', campo))
    from urllib.parse import unquote
    out = {}
    for par in campo.strip().split(";"):
        if "=" in par:
            k, v = par.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def _ler_gff_cds(gff: str):
    """(models, stats): the CDS of every transcript, keeping the longest isoform of each gene.
    Transcripts on two contigs or strands (trans-splicing) and CDS under another genetic code
    (transl_table, e.g. mitochondrial) are left out."""
    import gzip
    gtf = bool(re.search(r"\.gtf(\.gz)?$", gff))
    abrir = gzip.open if gff.endswith(".gz") else open
    pai, nome_gene, cds = {}, {}, {}
    stats = Counter()
    with abrir(gff, "rt", encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#"):
                if linha.startswith("##FASTA"):
                    break
                continue
            p = linha.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("CDS", "mRNA", "transcript", "gene"):
                continue
            at = _atributos(p[8], gtf)
            if gtf:
                tid = at.get("transcript_id", "")
                if p[2] != "CDS" or not tid:
                    continue
                pai[tid] = at.get("gene_id") or tid
                nome_gene.setdefault(pai[tid], at.get("gene_name") or pai[tid])
            else:
                ident, par = at.get("ID", ""), at.get("Parent", "").split(",")[0]
                if p[2] in ("mRNA", "transcript"):
                    if ident:
                        pai[ident] = par or ident
                    continue
                if p[2] == "gene":
                    if ident:
                        nome_gene[ident] = at.get("Name") or at.get("gene") or ident
                    continue
                tid = par or ident
            if not tid:
                continue
            e = cds.setdefault(tid, {"tid": tid, "contig": p[0], "strand": p[6], "parts": [], "ok": True,
                                     "product": at.get("product", ""), "protein": at.get("protein_id", ""),
                                     "gene": at.get("gene") or at.get("gene_name") or ""})
            if e["contig"] != p[0] or e["strand"] != p[6] or at.get("transl_table", "1") != "1":
                e["ok"] = False
            try:
                e["parts"].append((int(p[3]) - 1, int(p[4]), int(p[7]) if p[7].isdigit() else 0))
            except ValueError:
                e["ok"] = False
    melhor = {}
    for tid, e in cds.items():
        if not e["ok"] or not e["parts"]:
            stats["transcripts left out (trans-spliced or another genetic code)"] += 1
            continue
        gid = pai.get(tid, tid)
        e["gid"], e["gene"] = gid, e["gene"] or nome_gene.get(gid, "")
        e["len"] = sum(b - a for a, b, _ in e["parts"])
        if gid not in melhor or e["len"] > melhor[gid]["len"]:
            melhor[gid] = e
    stats["transcripts with CDS"] = len(cds)
    stats["genes"] = len(melhor)
    return list(melhor.values()), stats


def _ler_contigs(genoma: str, nomes: set):
    """(name, sequence) of the contigs in `nomes`: random access through the .fai when the FASTA
    is plain (samtools faidx builds it when missing), a single pass over the file otherwise."""
    import gzip
    fai = genoma + ".fai"
    if not genoma.endswith(".gz") and not os.path.exists(fai) and shutil.which("samtools"):
        subprocess.run(["samtools", "faidx", genoma], capture_output=True)
    if not genoma.endswith(".gz") and os.path.exists(fai):
        idx = {}
        with open(fai) as fh:
            for linha in fh:
                p = linha.split("\t")
                if len(p) >= 5:
                    idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
        with open(genoma, "rb") as fh:
            for n in nomes:
                if n not in idx:
                    continue
                tam, ini, lb, lw = idx[n]
                fh.seek(ini)
                bruto = fh.read(tam + ((tam + lb - 1) // lb) * (lw - lb))
                yield n, bruto.replace(b"\n", b"").replace(b"\r", b"")[:tam].decode("ascii", "replace").upper()
        return
    abrir = gzip.open if genoma.endswith(".gz") else open
    nome, partes = None, []
    with abrir(genoma, "rt", encoding="ascii", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if nome in nomes:
                    yield nome, "".join(partes).upper()
                nome, partes = linha[1:].split()[0] if linha[1:].split() else "", []
            elif nome in nomes:
                partes.append(linha.strip())
    if nome in nomes:
        yield nome, "".join(partes).upper()


def genes_do_hospedeiro(outdir: str, especie: str, genoma: str, gff: str) -> tuple:
    """(genes, stats) of one species: the CDS of the longest isoform of every annotated gene, read
    from the genome, with its features. Cached per species until the genome or the GFF change."""
    chave = (f"{os.path.abspath(genoma)}|{os.path.abspath(gff)}|{os.path.getmtime(genoma):.0f}|"
             f"{os.path.getmtime(gff):.0f}|{VERSAO_FEATURES}")
    cache = _cache_hospedeiro(outdir, especie)
    d = _ler_json(cache, {})
    if d.get("key") == chave:
        return d["genes"], d["stats"]
    log(f"{especie or 'genome'}: reading the CDS of {os.path.basename(gff)}")
    modelos, stats = _ler_gff_cds(gff)
    por_contig = defaultdict(list)
    for m in modelos:
        por_contig[m["contig"]].append(m)
    genes = []
    for nome, s in _ler_contigs(genoma, set(por_contig)):
        for m in por_contig.pop(nome):
            partes = sorted(m["parts"])
            nt = "".join(s[a:b] for a, b, _ in partes)
            fase = partes[-1][2] if m["strand"] == "-" else partes[0][2]
            nt = (revcomp(nt) if m["strand"] == "-" else nt)[fase:]
            feats, info = caracteristicas(nt)
            if feats is None:
                stats["genes under 10 codons"] += 1
                continue
            if info["stops"]:
                stats["genes with internal stops (left out)"] += 1
                continue
            genes.append({"id": m["protein"] or m["tid"], "gene": m["gene"], "product": m["product"][:120],
                          "contig": m["contig"], "start": partes[0][0], "end": partes[-1][1], "strand": m["strand"],
                          "cds": [[a, b] for a, b, _ in partes], "length": len(nt), "codons": info["codons"],
                          "features": feats})
    stats["genes on contigs missing from the genome"] = sum(len(v) for v in por_contig.values())
    stats["genes measured"] = len(genes)
    stats = {k: v for k, v in stats.items() if v}
    _gravar_json(cache, {"key": chave, "genes": genes, "stats": stats})
    _gravar_json(cache[:-5] + ".meta.json", {"gff": gff, "genome": genoma, "version": VERSAO_FEATURES, "stats": stats})
    log(f"{especie or 'genome'}: {len(genes)} host genes measured")
    return genes, stats


def _genes_vizinhos(outdir: str, eves: list, flanco: int, max_fundo: int, min_codons: int, msgs: list, excl: Counter):
    """Host gene points: those within `flanco` bp of an EVE locus ("host_near", all of them) and a
    sample of the rest of the genome ("host"). Each EVE gets the list of its neighbours (nb).
    A gene model whose CDS covers the EVE is probably the EVE itself, annotated as a host gene: it is
    kept out of both groups and named on the EVE (annotated_as)."""
    from evetool.comparative.flank_synteny import origens_de_segmentos, calibrar_origens, mapear_para_genoma
    anot = anotacoes_hospedeiro(outdir)
    por_sp = defaultdict(list)
    for r in eves:
        por_sp[r["species"]].append(r)
    pontos, sem_gff = [], []
    for sp, lista in sorted(por_sp.items()):
        a = anot.get(sp) or {}
        if not a.get("gff") or not a.get("genome"):
            sem_gff.append(sp.replace("_", " ") or "genome")
            continue
        try:
            genes, _stats = genes_do_hospedeiro(outdir, sp, a["genome"], a["gff"])
        except Exception as erro:
            msgs.append(f"{sp or 'genome'}: host genes not read ({erro})")
            continue
        origens = {}
        if a["segments"] and a["own"]:
            origens = calibrar_origens(origens_de_segmentos(a["own"]), a["own"], a["genome"])
        por_contig = defaultdict(list)
        for r in lista:
            if not r.get("contig") or (a["segments"] and r["contig"] not in origens):
                continue
            c, x, y, _inv = mapear_para_genoma(r["contig"], r["start"], r["end"], origens)
            r["g_loc"] = f"{c}:{x + 1}-{y}"
            por_contig[c].append((r, x, y))
        perto, anotados = {}, set()
        for g in genes:
            for r, x, y in por_contig.get(g["contig"], ()):
                if sum(max(0, min(b, y) - max(a_, x)) for a_, b in g["cds"]) >= 30:
                    anotados.add(g["id"])
                    r.setdefault("annotated_as", []).append(g["gene"] or g["id"])
                    continue
                d = max(0, g["start"] - y, x - g["end"])
                if d <= flanco:
                    perto.setdefault(g["id"], []).append((r, d))
        if anotados:
            msgs.append(f"{sp.replace('_', ' ') or 'genome'}: {len(anotados)} host gene model(s) cover an EVE and were "
                        "left out (probably the EVE itself, annotated as a host gene)")
        rotulo = sp.replace("_", " ") or "genome"
        fundo = []
        for g in genes:
            if g["id"] in anotados:
                continue
            if g["codons"] < min_codons:
                excl[f"host genes < {min_codons} codons"] += 1
                continue
            p = {"id": f"{sp}|{g['id']}", "origin": "host_near" if g["id"] in perto else "host", "family": "Host",
                 "taxon": "", "species": sp, "gene": g["gene"], "product": g["product"], "protein_class": "n/a",
                 "has_domain": None, "domain": None, "length": g["length"], "codons": g["codons"], "stops": 0,
                 "g_loc": f"{g['contig']}:{g['start'] + 1}-{g['end']}", "features": g["features"]}
            if g["id"] in perto:
                viz = sorted(perto[g["id"]], key=lambda t: t[1])
                p.update(set=f"{rotulo} · near EVEs", near=[r["id"] for r, _d in viz], dist=viz[0][1])
                for r, d in viz:
                    r.setdefault("nb", []).append([p["id"], d])
                pontos.append(p)
            else:
                p["set"] = f"{rotulo} · genome"
                fundo.append(p)
        amostra, fora = _amostrar(fundo, lambda p: p["species"], max_fundo)
        if fora:
            excl[f"{rotulo}: genome genes over the sampling cap"] += fora
        pontos += amostra
    if sem_gff:
        msgs.append("no host annotation (GFF) for " + ", ".join(sem_gff[:6]) + (" and others" if len(sem_gff) > 6 else "")
                    + ": set one in the host annotation table to compare those EVEs with their neighbouring genes")
    return pontos


# ---------------------------------------------------------------------------
# EVEs in the frame of the viral protein
# ---------------------------------------------------------------------------
def _no_quadro(qseq_nt: str, quadro: str) -> str:
    """Nucleotides of the stretch read in a frame of protein_map ('+1'..'-3')."""
    s = limpar(qseq_nt)
    if quadro.startswith("-"):
        s = revcomp(s)
    s = s[int(quadro[1]) - 1:]
    return s[:len(s) - len(s) % 3]


def eves_no_quadro_viral(outdir: str, incluir_orf: bool = False) -> tuple:
    """(records, counts) of the approved EVE loci, each read in the frame of its viral protein.

    One record per locus (overlapping ORFs of one insertion): the member with the longest
    stretch in the viral frame. A locus counts as carrying a viral HMM domain when any of its
    ORFs does. EVEs supported by the HMM alone have no aligned viral protein: they are left
    out, or read in the ORF frame with incluir_orf (flagged)."""
    from evetool.curation import curation as cur
    from evetool.validation import protein_map as pm
    from evetool.annotation.cdsearch import ler_fasta
    dados = cur.listar_candidatas(outdir)
    if dados.get("error"):
        raise RuntimeError(dados["error"])
    todos = dados["items"]
    aprovados = [i for i in todos if not i["removed_flag"]]
    mapa = cur.mapa_de_proteinas(outdir, todos)
    alin = {m["query"]: m.get("best") or {} for m in (mapa.get("items") or [])}
    pastas = cur.pastas_especies(outdir)
    pasta_de = (lambda sp: pastas.get(sp)) if len(pastas) > 1 else (lambda sp: next(iter(pastas.values())))

    hits, orfs_nt = {}, {}
    por_pasta = defaultdict(set)
    for i in aprovados:
        if pasta_de(i["species"]):
            por_pasta[pasta_de(i["species"])].add(i["orf"])
    for p, orfs in por_pasta.items():
        for o, hs in pm._hits_virais(p, orfs).items():
            hits[(p, o)] = hs
        fa = os.path.join(p, "04_candidatas_eves", "eve_candidatas_unificadas.fasta")
        if incluir_orf and os.path.exists(fa):
            for o, s in ler_fasta(fa).items():
                if o in orfs:
                    orfs_nt[(p, o)] = s

    contagem = Counter()
    por_locus = defaultdict(list)
    for i in aprovados:
        por_locus[(i["species"], i["locus"] or i["query"])].append(i)
    registros = []
    for (sp, locus), membros in por_locus.items():
        melhor = None
        for i in membros:
            a = alin.get(i["query"]) or {}
            k = (pasta_de(i["species"]), i["orf"])
            seq, fonte = "", ""
            if a.get("frame"):
                h = next((x for x in hits.get(k, []) if x["id"] == a.get("ref")), None)
                if a["frame"] == "ORF" and orfs_nt.get(k):
                    seq, fonte = orfs_nt[k], "ORF (in the viral frame)"
                elif h and h.get("qseq"):
                    seq, fonte = _no_quadro(h["qseq"], a["frame"]), "viral frame"
            if not seq and incluir_orf and orfs_nt.get(k):
                seq, fonte = orfs_nt[k], "ORF frame (no viral alignment)"
            if seq and (melhor is None or len(seq) > len(melhor[1])):
                melhor = (i, seq, fonte, a)
        if not melhor:
            contagem["no viral frame"] += 1
            continue
        i, seq, fonte, a = melhor
        # an EVE classified only to a higher rank (Riboviria, an order...) has no family to compare
        taxon = i["family"] if i["family"] not in ("", "n/a") else ""
        fam = taxon if taxon.endswith("viridae") else "Unassigned"
        dominio = next((m for m in membros if m["hmm_n_domains"] > 0), None)
        # span of the locus (every ORF of the insertion) on the analysed genome, 0-based
        com_coord = [m for m in membros if m.get("contig") and m.get("end")]
        span = {"contig": com_coord[0]["contig"], "start": min(m["start"] for m in com_coord),
                "end": max(m["end"] for m in com_coord)} if com_coord else {"contig": "", "start": 0, "end": 0}
        registros.append({**span,
            "id": i["query"], "orf": i["orf"], "locus": locus, "species": sp, "family": fam, "taxon": taxon,
            "product": i["product"], "protein_class": i["protein_class"] or "n/a",
            "has_domain": dominio is not None,
            "domain": {"n": dominio["hmm_n_domains"], "model": dominio["hmm_best_model"],
                       "score": dominio["hmm_best_score"], "label": dominio["hmm_consensus_family"]} if dominio else None,
            "frame_source": fonte, "frame": a.get("frame", ""), "orf_frame": a.get("orf_frame"),
            "ref": a.get("ref", ""), "pident": a.get("pident"), "length": len(seq), "seq": seq,
            "n_orfs": len(membros)})
    contagem["loci"] = len(por_locus)
    return registros, contagem, mapa.get("error")


# ---------------------------------------------------------------------------
# Analysis: feature ranking, PCA and t-SNE
# ---------------------------------------------------------------------------
def _amostrar(lista, chave, maximo, semente=0):
    """At most `maximo` items per value of chave(item), drawn reproducibly."""
    if not maximo or maximo <= 0:
        return lista, 0
    por = defaultdict(list)
    for x in lista:
        por[chave(x)].append(x)
    rng = random.Random(semente)
    saida, fora = [], 0
    for k in sorted(por, key=str):
        grupo = sorted(por[k], key=lambda x: x["id"])
        if len(grupo) > maximo:
            rng.shuffle(grupo)
            fora += len(grupo) - maximo
            grupo = grupo[:maximo]
        saida.extend(grupo)
    return saida, fora


def _escolher(np, X, nomes, grupos, top, limite_r=0.9):
    """Features ranked by the Kruskal-Wallis H across the groups (groups with >= 3 members);
    a feature correlated (|Pearson r| > limite_r) with one already chosen is skipped."""
    from scipy.stats import kruskal
    cont = Counter(grupos)
    validos = [g for g in cont if cont[g] >= 3]
    if len(validos) < 2:
        return None, [], []
    idx = {g: [k for k, x in enumerate(grupos) if x == g] for g in validos}
    n = sum(len(v) for v in idx.values())
    ranking = []
    for j, f in enumerate(nomes):
        col = X[:, j]
        amostras = [col[idx[g]] for g in validos]
        if np.nanstd(np.concatenate(amostras)) == 0:
            continue
        try:
            H, p = kruskal(*amostras)
        except ValueError:
            continue
        if math.isnan(H):
            continue
        ranking.append({"feature": f, "H": round(float(H), 3), "p": float(p),
                        "eps2": round(float(H) / max(1, n - 1), 4)})
    ranking.sort(key=lambda r: -r["H"])
    linhas = sorted(i for v in idx.values() for i in v)
    escolhidas = []
    for r in ranking:
        j = nomes.index(r["feature"])
        if any(abs(np.corrcoef(X[linhas, j], X[linhas, nomes.index(e["feature"])])[0, 1]) > limite_r
               for e in escolhidas):
            r["skipped"] = "redundant"
            continue
        escolhidas.append(r)
        if len(escolhidas) >= top:
            break
    return {"groups": {g: len(v) for g, v in idx.items()}, "n": n}, ranking, escolhidas


def analisar(outdir: str, top: int = 6, discriminar: str = "family_circ", conjunto_features: str = "top",
             ajuste: str = "all", min_codons: int = 50, max_por_especie: int = 20, max_hospedeiro: int = 300,
             classe: str = "", incluir_orf: bool = False, perplexidade: float = 30.0,
             hospedeiro_anotacao: bool = True, flanco_vizinhos: int = 50000) -> dict:
    """Features of the EVEs, of their neighbouring host genes and of every added set, the ranking,
    PCA, t-SNE and the amelioration index; writes analysis.json."""
    import numpy as np
    base = pasta(outdir)
    os.makedirs(base, exist_ok=True)
    msgs = []
    excl = Counter()

    eves, cont_eves, erro_mapa = eves_no_quadro_viral(outdir, incluir_orf=incluir_orf)
    if erro_mapa:
        msgs.append(f"protein map: {erro_mapa}")
    excl["EVE loci without a viral frame"] = cont_eves.get("no viral frame", 0)
    pontos = []
    for r in eves:
        feats, info = caracteristicas(r.pop("seq"))
        r.update(info)
        if feats is None or info["codons"] < min_codons:
            excl[f"EVE loci < {min_codons} codons"] += 1
            continue
        pontos.append({**r, "origin": "EVE", "set": "EVEs", "features": feats})

    # Host genes of the genome annotation: neighbours of each EVE and a sample of the genome
    if hospedeiro_anotacao:
        pontos += _genes_vizinhos(outdir, [p for p in pontos if p["origin"] == "EVE"], max(0, int(flanco_vizinhos)),
                                  max_hospedeiro, min_codons, msgs, excl)

    for c in listar_conjuntos(outdir):
        seqs = [r for r in _sequencias_do_conjunto(outdir, c) if r.get("features")]
        curtas = [r for r in seqs if r["codons"] < min_codons]
        excl[f"{c['name']}: < {min_codons} codons"] += len(curtas)
        seqs = [r for r in seqs if r["codons"] >= min_codons]
        if c["kind"] == "host":
            seqs, fora = _amostrar(seqs, lambda r: c["id"], max_hospedeiro)
        else:
            seqs, fora = _amostrar(seqs, lambda r: r.get("species") or r["id"], max_por_especie)
        if fora:
            excl[f"{c['name']}: over the sampling cap"] += fora
        for r in seqs:
            pontos.append({"id": r["id"], "origin": c["kind"], "set": c["name"], "family": r["family"] or "Host",
                           "taxon": r["family"],
                           "species": r.get("species", ""), "product": r.get("product", ""),
                           "protein_class": r.get("protein_class", "n/a"),
                           "has_domain": r.get("has_domain") if c["kind"] == "circulating" else None,
                           "domain": r.get("domain"), "family_source": r.get("family_source", ""),
                           "length": r["length"], "codons": r["codons"], "stops": r["stops"],
                           "problems": r.get("problems", []), "features": r["features"]})

    if classe:
        antes = len(pontos)
        pontos = [p for p in pontos if p["origin"] in ("host", "host_near") or p.get("protein_class") == classe]
        excl[f"other protein classes than {classe}"] = antes - len(pontos)

    n_origem = Counter(p["origin"] for p in pontos)
    resultado = {"created": datetime.datetime.now().isoformat(timespec="seconds"),
                 "params": {"top": top, "discriminate": discriminar, "feature_set": conjunto_features,
                            "fit": ajuste, "min_codons": min_codons, "max_per_species": max_por_especie,
                            "max_host": max_hospedeiro, "protein_class": classe, "include_orf_frame": incluir_orf,
                            "perplexity": perplexidade, "host_annotation": hospedeiro_anotacao,
                            "neighbour_flank": flanco_vizinhos},
                 "features": FEATURES, "feature_groups": GRUPOS_FEATURES,
                 "counts": dict(n_origem), "eve_loci": cont_eves.get("loci", 0),
                 "excluded": {k: v for k, v in excl.items() if v}, "messages": msgs,
                 "points": [], "ranking": [], "selected": [], "pca": None, "tsne": None}
    if len(pontos) < 4:
        msgs.append("fewer than 4 sequences pass the filters: nothing to project")
        _gravar_resultado(outdir, resultado, pontos)
        return resultado

    # Matrix: an amino acid absent from a sequence has RSCU 1 (no preference); any other gap
    # takes the median of the feature; features constant over every sequence are dropped
    X = np.array([[p["features"].get(f) if p["features"].get(f) is not None else np.nan for f in FEATURES]
                  for p in pontos], dtype=float)
    for j, f in enumerate(FEATURES):
        col = X[:, j]
        if np.isnan(col).any():
            col[np.isnan(col)] = 1.0 if f.startswith("RSCU_") else (np.nanmedian(col) if (~np.isnan(col)).any() else 0.0)
    variaveis = [j for j in range(len(FEATURES)) if np.std(X[:, j]) > 1e-9]

    # Which sequences and labels define the "best discriminators"
    origem = np.array([p["origin"] for p in pontos])
    familia = [p["family"] for p in pontos]
    modos = {
        "family_circ": ("viral family of the circulating viruses",
                        [k for k in range(len(pontos)) if origem[k] == "circulating"], familia),
        "family_all": ("viral family (EVEs and circulating viruses)",
                       [k for k in range(len(pontos)) if origem[k] in ("EVE", "circulating")], familia),
        "origin": ("origin (EVE / circulating / host near EVEs / host genome)", list(range(len(pontos))), list(origem)),
        "family_origin": ("viral family x origin", [k for k in range(len(pontos)) if origem[k] in ("EVE", "circulating")],
                          [f"{familia[k]} · {origem[k]}" for k in range(len(pontos))]),
    }
    ordem_modos = [discriminar] + [m for m in ("family_circ", "family_all", "family_origin", "origin") if m != discriminar]
    nomes_var = [FEATURES[j] for j in variaveis]
    selecao = None
    for m in ordem_modos:
        rot, linhas, rotulos = modos.get(m, modos["family_circ"])
        if len(linhas) < 6:
            continue
        info, ranking, escolhidas = _escolher(np, X[np.ix_(linhas, variaveis)], nomes_var,
                                              [rotulos[k] for k in linhas], max(2, top))
        if info:
            selecao = (m, rot, info, ranking, escolhidas)
            break
    if selecao:
        m, rot, info, ranking, escolhidas = selecao
        if m != discriminar:
            msgs.append(f"not enough groups to rank by '{modos[discriminar][0]}' (at least two groups with "
                        f"3 sequences): ranked by {rot} instead")
        resultado.update(discriminate_used=m, discriminate_label=rot, rank_groups=info["groups"],
                         ranking=ranking[:40], selected=escolhidas)
    else:
        msgs.append("no grouping has two groups with 3 sequences: the PCA uses every feature")

    if conjunto_features == "top" and selecao:
        usar = [e["feature"] for e in selecao[4]]
    elif conjunto_features in GRUPOS_FEATURES:
        usar = [f for f in GRUPOS_FEATURES[conjunto_features] if f in nomes_var]
    else:
        usar = nomes_var
    resultado["used_features"] = usar
    J = [FEATURES.index(f) for f in usar]
    Xs = X[:, J]

    # PCA on z-scores (mean and sd of the fitted set); the rest is projected on its axes
    ajuste_idx = [k for k in range(len(pontos)) if origem[k] == "circulating"] if ajuste == "circulating" else list(range(len(pontos)))
    if len(ajuste_idx) < 3:
        msgs.append("fewer than 3 circulating sequences: the PCA was fitted on every sequence")
        ajuste_idx, ajuste = list(range(len(pontos))), "all"
    media, dp = Xs[ajuste_idx].mean(axis=0), Xs[ajuste_idx].std(axis=0, ddof=1)
    dp[dp == 0] = 1.0
    Z = (Xs - media) / dp
    _u, s, vt = np.linalg.svd(Z[ajuste_idx], full_matrices=False)
    autoval = s ** 2 / max(1, len(ajuste_idx) - 1)
    ncomp = min(4, len(usar))
    scores = Z @ vt[:ncomp].T
    loadings = (vt[:ncomp].T * np.sqrt(autoval[:ncomp]))      # correlation of each feature with each PC
    resultado["pca"] = {"fit": ajuste, "explained": [round(float(v), 4) for v in (autoval / autoval.sum())[:ncomp]],
                        "loadings": {f: [round(float(x), 4) for x in loadings[k]] for k, f in enumerate(usar)}}

    resultado["amelioration"] = _amelioracao(np, X, variaveis, J, pontos)

    # t-SNE on the same z-scores (Barnes-Hut, fixed seed)
    tsne = None
    if len(pontos) >= 10:
        try:
            from sklearn.manifold import TSNE
            perp = float(max(2.0, min(perplexidade, (len(pontos) - 1) / 3.0)))
            emb = TSNE(n_components=2, perplexity=perp, init="pca", learning_rate="auto",
                       random_state=0).fit_transform(Z)
            tsne = emb
            resultado["tsne"] = {"perplexity": round(perp, 1)}
            if perp < perplexidade:
                msgs.append(f"t-SNE perplexity lowered to {perp:.0f} (at most a third of the points)")
        except ImportError:
            msgs.append("t-SNE needs scikit-learn (conda install scikit-learn)")
        except Exception as erro:
            msgs.append(f"t-SNE failed: {erro}")
    else:
        msgs.append("t-SNE needs at least 10 sequences")

    for k, p in enumerate(pontos):
        p["pc"] = [round(float(v), 4) for v in scores[k]]
        p["ts"] = [round(float(v), 3) for v in tsne[k]] if tsne is not None else None
    _gravar_resultado(outdir, resultado, pontos)
    log(f"{n_origem.get('EVE', 0)} EVE loci, {n_origem.get('circulating', 0)} circulating and "
        f"{n_origem.get('host', 0)} host CDS; PCA on {len(usar)} feature(s)")
    return resultado


def _amelioracao(np, X, variaveis, usadas, pontos) -> dict:
    """Where each EVE sits on the axis from the circulating viruses of its family to the host genes:
    A = (E - V)·(H - V) / |H - V|², with V the centroid of the circulating viruses of the same family
    (at least 3) and H the centroid of its neighbouring host genes ("local") or of the host genome
    genes ("genome": the genes of its species, else the host sets added). 0 = at the virus centroid,
    1 = at the host centroid; below 0 or above 1 is beyond them. A projection, not a ratio of
    distances: the noise of a short EVE off the axis does not pull it to the middle. In z-scores
    over every sequence, on every feature ("all") or on the features of the PCA ("pca")."""
    def z(M):
        dp = M.std(axis=0, ddof=1)
        dp[dp == 0] = 1.0
        return (M - M.mean(axis=0)) / dp
    espacos = {"all": z(X[:, variaveis]), "pca": z(X[:, usadas])}
    idx = {p["id"]: k for k, p in enumerate(pontos)}
    circ, fundo, enviados = defaultdict(list), defaultdict(list), []
    for k, p in enumerate(pontos):
        if p["origin"] == "circulating" and p["family"] != "Unassigned":
            circ[p["family"]].append(k)
        elif p["origin"] == "host":
            (fundo[p["species"]] if p.get("g_loc") else enviados).append(k)
    cont = Counter()
    for k, p in enumerate(pontos):
        if p["origin"] != "EVE":
            continue
        V = circ.get(p["family"], [])
        L = [idx[i] for i, _d in p.get("nb", []) if i in idx]
        G = fundo.get(p["species"]) if len(fundo.get(p["species"], [])) >= 3 else (enviados if len(enviados) >= 3 else [])
        p["am"] = {}
        for nome, Z in espacos.items():
            if len(V) < 3:
                p["am"][nome] = None
                continue
            cv = Z[V].mean(axis=0)
            par = []
            for H in (L, G):
                if (H is L and not H) or (H is G and len(H) < 3):
                    par.append(None)
                    continue
                eixo = Z[H].mean(axis=0) - cv
                n2 = float(eixo @ eixo)
                par.append(round(float((Z[k] - cv) @ eixo) / n2, 4) if n2 > 0 else None)
            p["am"][nome] = par
        a = p["am"].get("all")
        cont["with a family reference"] += a is not None
        cont["local"] += bool(a and a[0] is not None)
        cont["genome"] += bool(a and a[1] is not None)
    return dict(cont)


def _gravar_resultado(outdir: str, resultado: dict, pontos: list):
    """analysis.json for the page (feature values as one array per point) and the TSVs."""
    base = pasta(outdir)
    campos = ["id", "origin", "set", "family", "taxon", "species", "protein_class", "product", "has_domain", "length",
              "codons", "stops", "frame_source", "pc", "ts", "gene", "near", "dist", "nb", "am", "annotated_as", "g_loc"]
    resultado["points"] = []
    for p in pontos:
        d = {k: p.get(k) for k in campos}
        d["dom"] = (p["domain"] or {}).get("model") if p.get("domain") else None
        d["x"] = [p["features"].get(f) for f in FEATURES]
        resultado["points"].append(d)
    _gravar_json(os.path.join(base, "analysis.json"), resultado)
    meta = ["id", "origin", "set", "family", "taxon", "species", "protein_class", "product", "has_domain", "domain_model",
            "length", "codons", "internal_stops", "frame_source", "genome_location", "host_gene", "near_eves",
            "neighbour_genes", "nearest_neighbour_bp", "amelioration_local", "amelioration_genome", "annotated_as_host_gene"]
    with open(os.path.join(base, "composition_features.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(meta + FEATURES)
        for p in pontos:
            w.writerow([p.get("id"), p.get("origin"), p.get("set"), p.get("family"), p.get("taxon", ""), p.get("species"),
                        p.get("protein_class"), p.get("product"), p.get("has_domain"),
                        (p.get("domain") or {}).get("model", ""), p.get("length"), p.get("codons"), p.get("stops"),
                        p.get("frame_source", ""), p.get("g_loc", ""), p.get("gene", ""), ",".join(p.get("near") or []),
                        len(p.get("nb") or []) if p.get("origin") == "EVE" else "",
                        min((d for _i, d in p.get("nb") or []), default=""),
                        *[("" if not (p.get("am") or {}).get("all") or p["am"]["all"][j] is None else p["am"]["all"][j]) for j in (0, 1)],
                        "; ".join(p.get("annotated_as") or [])]
                       + ["" if p["features"].get(f) is None else p["features"][f] for f in FEATURES])
    if resultado.get("pca"):
        with open(os.path.join(base, "composition_pca.tsv"), "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            npc = len(resultado["pca"]["explained"])
            w.writerow(["id", "origin", "set", "family", "has_domain"] + [f"PC{k + 1}" for k in range(npc)]
                       + ["tSNE1", "tSNE2"])
            for p in pontos:
                w.writerow([p["id"], p["origin"], p["set"], p["family"], p.get("has_domain")] + list(p.get("pc") or [])
                           + list(p.get("ts") or ["", ""]))


def carregar(outdir: str) -> dict:
    """What the page reads: the sets added and the last analysis (or None)."""
    base = pasta(outdir)
    an = _ler_json(os.path.join(base, "analysis.json"), None)
    if an is not None:
        an["files"] = {k: os.path.join(base, n) for k, n in
                       (("features", "composition_features.tsv"), ("pca", "composition_pca.tsv"))
                       if os.path.exists(os.path.join(base, n))}
    try:
        anot = anotacoes_hospedeiro(outdir)
    except Exception as erro:
        anot = {"_error": str(erro)}
    return {"sets": listar_conjuntos(outdir), "analysis": an, "host_annotation": anot}


def arquivos(outdir: str) -> list:
    """Files whose change must refresh the page (version stamp)."""
    base = pasta(outdir)
    return [os.path.join(base, "sets.json"), os.path.join(base, "analysis.json"), os.path.join(base, NOME_ANOT)]


def main():
    ap = argparse.ArgumentParser(description="Codon usage and dinucleotide composition of EVEs vs circulating viruses")
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--add", help="FASTA of coding sequences (CDS) to add as a set")
    ap.add_argument("--name", default="", help="Name of the set added")
    ap.add_argument("--kind", choices=("circulating", "host"), default="circulating")
    ap.add_argument("--family", default="", help="Family of the sequences whose header has none")
    ap.add_argument("--no-auto-family", action="store_true", help="Do not assign the missing families with DIAMOND")
    ap.add_argument("--remove", help="Id of a set to remove")
    ap.add_argument("--list", action="store_true", help="List the sets")
    ap.add_argument("--run", action="store_true", help="Run the analysis")
    ap.add_argument("--top", type=int, default=6, help="Number of discriminating features in the PCA (default 6)")
    ap.add_argument("--discriminate", default="family_circ",
                    choices=("family_circ", "family_all", "origin", "family_origin"))
    ap.add_argument("--features", default="top", choices=("top", "all", "composition", "codon", "dinucleotide"))
    ap.add_argument("--fit", default="all", choices=("all", "circulating"))
    ap.add_argument("--min-codons", type=int, default=50)
    ap.add_argument("--max-per-species", type=int, default=20)
    ap.add_argument("--max-host", type=int, default=300)
    ap.add_argument("--protein-class", default="")
    ap.add_argument("--include-orf-frame", action="store_true",
                    help="Also measure the EVEs without a viral alignment, in the frame of their ORF")
    ap.add_argument("--perplexity", type=float, default=30.0)
    ap.add_argument("--no-host-annotation", action="store_true",
                    help="Do not read the host genes of the genome annotation (neighbours of the EVEs and genome sample)")
    ap.add_argument("--neighbour-flank", type=int, default=50000,
                    help="Host genes within this distance of an EVE are its neighbours (default 50000 bp)")
    ap.add_argument("--gff", nargs=2, metavar=("SPECIES", "GFF"),
                    help="Set the GFF/GTF of one species (species folder name, or '' in a single-species run)")
    ap.add_argument("-t", "--threads", type=int, default=4)
    a = ap.parse_args()
    if a.add:
        with open(a.add, encoding="utf-8", errors="replace") as fh:
            print(json.dumps(adicionar_conjunto(a.outdir, fh.read(), a.name or os.path.basename(a.add), a.kind,
                                                a.family, not a.no_auto_family, a.threads), indent=1))
    if a.gff:
        print(json.dumps(definir_gff(a.outdir, a.gff[0], a.gff[1])))
    if a.remove:
        print(json.dumps(remover_conjunto(a.outdir, a.remove)))
    if a.list:
        for c in listar_conjuntos(a.outdir):
            print(f"{c['id']}\t{c['kind']}\t{c['valid']}/{c['n']}\t{c['name']}\t"
                  + ", ".join(f"{k} {v}" for k, v in list(c["families"].items())[:6]))
    if a.run:
        r = analisar(a.outdir, top=a.top, discriminar=a.discriminate, conjunto_features=a.features, ajuste=a.fit,
                     min_codons=a.min_codons, max_por_especie=a.max_per_species, max_hospedeiro=a.max_host,
                     classe=a.protein_class, incluir_orf=a.include_orf_frame, perplexidade=a.perplexity,
                     hospedeiro_anotacao=not a.no_host_annotation, flanco_vizinhos=a.neighbour_flank)
        print(json.dumps({k: r.get(k) for k in ("counts", "excluded", "messages", "selected", "pca", "tsne", "amelioration")},
                         indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
