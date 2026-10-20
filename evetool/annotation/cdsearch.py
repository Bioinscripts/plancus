#!/usr/bin/env python3
"""
CD-Search annotation of the candidate EVEs (an independent layer next to the HMM).

Every candidate peptide is searched against CDD, by default locally with rpsblast+ on the
CDD database in CDD/rpsdb (offline, no NCBI queue), or with --modo online through the
NCBI Batch CD-Search service (the same search the NCBI web page runs). Each
conserved-domain hit is placed in a viral scope from the domain's own CDD description;
in local mode the source taxon of the domain (Cdd_domains.tsv) is used as well:

  RNA virus          RdRp, nucleocapsid, glycoprotein... of RNA viruses (in scope)
  ssDNA virus        Rep, NS1, capsid... of non-retroviral ssDNA viruses (in scope)
  Retro (excluded)   reverse transcriptase, integrase, gag, RNase H, retrotransposons,
                     retroviruses and pararetroviruses (hepadna-, caulimovirus)
  Other virus        viral, but neither of the above (dsDNA viruses, phages...)
  Non-viral          cellular domains

CD-Search has no taxonomic filter, so the scope comes from the description (plus, in local
mode, the domain's source taxon: a domain built from viral proteins is never "Non-viral",
and one built only from cellular proteins is not called RNA/ssDNA virus unless its title
says it is viral). The domain names, descriptions and taxa travel with the results so any
call can be checked by hand.

Local mode does not classify hits as specific / non-specific / superfamily (that needs
NCBI's rpsbproc); its hits are typed "rpsblast".

Outputs, under 04_candidatas_eves/ of each species:
  cdsearch_raw.tsv      raw hits from NCBI (cache: sequences already searched are not resent)
  cdsearch_raw_local.tsv  the same for the local rpsblast+ search
  cdsearch_hits.tsv     one domain per line, with its scope and description
  cdsearch_summary.tsv  one ORF per line: scope, best in-scope domain, all in-scope domains

Usage:
  python -m evetool.annotation.cdsearch -o results_eve                    # local rpsblast+ (default)
  python -m evetool.annotation.cdsearch -o results_eve --evalue 0.001 --forcar -t 8
  python -m evetool.annotation.cdsearch -o results_eve --modo online      # NCBI Batch CD-Search
"""
import os
import re
import csv
import glob
import gzip
import time
import shutil
import hashlib
import argparse
import tempfile
import subprocess
import urllib.parse
import urllib.request
from collections import OrderedDict
from evetool import config

URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
CDDID_URL = "https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz"
CDDID_PADRAO = config.hmm_path("CDD", "cddid_all.tbl.gz")
RPSDB_PADRAO = config.hmm_path("CDD", "rpsdb", "Cdd")
DOMINIOS_PADRAO = os.path.join(os.path.dirname(RPSDB_PADRAO), "Cdd_domains.tsv")
LOTE = 1000            # the service takes up to 4000 queries per request; smaller is gentler
ESPERA = 5             # seconds between status checks

ESCOPOS_ALVO = ("RNA virus", "ssDNA virus")
ORDEM_ESCOPO = ["RNA virus", "ssDNA virus", "Retro (excluded)", "Other virus", "Non-viral", "No CDD hit"]

# Checked on the domain name + the first sentence of its description (its title), where
# the defining words are; the rest of a description often cites other groups in passing.
RX_RETRO = re.compile(
    r"retrovir|reverse transcriptase|\brt_|\brvt|rna-dependent dna pol|integrase core|\brve\b|"
    r"retroviral|retrotranspos|retroelement|\bltr\b|non-ltr|lentivir|spumavir|\bhiv\b|\bsiv\b|"
    r"\bhtlv|\bherv|\bgag\b|gag[_-]|retrotran|retropepsin|hepadnavir|hepatitis b\b|caulimovir|"
    r"pararetrovir|metavir|pseudovir|belpao|\bty[13]\b|copia|gypsy", re.I)
RX_SSDNA = re.compile(
    r"ssdna vir|single[- ]strand(ed)? dna vir|parvovir|parvo_|\bpv_nsp|densovir|denso_|"
    r"erythrovir|bocavir|dependovir|circovir|circo_|geminivir|gemini_|nanovir|anellovir|"
    r"bidnavir|smacovir|genomovir|bacilladnavir|redondovir|microvirus|microviridae|microvir_|"
    r"\binovir|spiravir", re.I)
RX_CRESS = re.compile(r"\bCRESS\b")            # upper case only: never the plant "cress"
RX_RNA = re.compile(
    r"rna vir|rna-containing vir|rna containing vir|viral rna|rna-dependent rna pol|\brdrp|"
    r"rna_dep_rnap|[+-]ssrna|ssrna vir|ssrna phage|dsrna vir|negative[- ]s(trand|ense)|positive[- ]s(trand|ense)|"
    r"riboviria|orthornavir|rhabdovir|rhabdo_|mononeg|paramyxo|pneumovir|filovir|ebola|marburg ?vir|bornavir|"
    r"nyamivir|chuvir|bunya|phlebovir|phlebo_|phenui|hantavir|hanta_|nairovir|peribunya|arenavir|"
    r"orthomyxo|(?<!haemophilus )influenza(?!e)|\bflu_|tenuivir|emaravir|ophiovir|flavivir|flavi_|pestivir|hepacivir|\bhcv|"
    r"togavir|alphavir|alpha_e|rubivir|coronavir|\bcov\b|\bcov_|-cov\b|nidovir|arterivir|picorna|"
    r"calicivir|astrovir|nodavir|tombusvir|luteovir|potyvir|bromovir|tymovir|tobamovir|"
    r"virgavir|closterovir|hepevir|reovir|rotavir|orbivir|birnavir|totivir|partitivir|narnavir|"
    r"mitovir|ourmiavir|dicistrovir|iflavir|negevir|levivir|cystovir|tolivir|martellivir|"
    r"amarillovir|ghabrivir|durnavir|corona_|corona nucleocap|mosaic vir|ringspot|mottle vir|\bcmv_|\brabies|measles|vesicular stomatitis|dengue|"
    r"hepatitis [ace] vir|poliovir|\bsars|\bmers\b", re.I)
# dsDNA viruses and phages: viral, but outside the scope
RX_DSDNA = re.compile(
    r"herpes|poxvir|\bpox_|vaccinia|baculovir|adenovir|polyomavir|papillomavir|iridovir|asfarvir|"
    r"african swine|mimivir|phycodnavir|nudivir|polydnavir|ascovir|hytrosavir|nimavir|"
    r"whispovir|epstein-barr|\bebv|cytomegalo|phage|caudovir", re.I)
RX_CELULAR = re.compile(r"arabidopsis|eukaryot|\brdr\d|rdrp1|qde-?1|ego-1|\brrf-|cellular|"
                        r"plant rna-dependent|fungal", re.I)
RX_VIRAL = re.compile(r"(?<!anti)vir(us|al|idae|ales|ion)|phage|capsid|virion", re.I)
# phrases that only borrow a virus name ("BTB/POZ: poxvirus and zinc finger", COG categories)
RX_RUIDO = re.compile(r"poxvirus and zinc finger|\[mobilome:[^\]]*\]", re.I)


def log(msg):
    print(f"[cdsearch] {msg}", flush=True)


# ---------------------------------------------------------------------------
# CDD descriptions and the viral scope of each domain
# ---------------------------------------------------------------------------
def carregar_descricoes(caminho: str = CDDID_PADRAO) -> dict:
    """{accession and short name: (short name, description)} from cddid_all.tbl(.gz)."""
    if not os.path.exists(caminho):
        log(f"downloading the CDD descriptions ({CDDID_URL})...")
        os.makedirs(os.path.dirname(caminho), exist_ok=True)
        urllib.request.urlretrieve(CDDID_URL, caminho + ".tmp")
        os.replace(caminho + ".tmp", caminho)
    abrir = gzip.open if caminho.endswith(".gz") else open
    desc = {}
    with abrir(caminho, "rt", encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) >= 4:
                desc[p[1]] = (p[2], p[3])
    return desc


def titulo(nome: str, descricao: str) -> str:
    """Domain name + first sentence of its description."""
    primeira = re.split(r"(?<=[a-z0-9)\]])\. ", descricao or "", maxsplit=1)[0]
    return f"{nome} {primeira}"


def escopo_do_dominio(nome: str, descricao: str, taxon: str = "", origem: str = "") -> str:
    """Viral scope of a domain; taxon/origin (local mode) are the CDD source taxon of the model."""
    if origem == "viral":
        # The taxon names the group the title often leaves out (Coronaviridae, Poxviridae...)
        escopo = escopo_do_dominio(f"{nome} {taxon}", descricao)
        return "Other virus" if escopo == "Non-viral" else escopo
    escopo = _escopo_pela_descricao(nome, descricao)
    if origem == "cellular" and escopo in ESCOPOS_ALVO and not RX_VIRAL.search(titulo(nome, descricao)):
        return "Non-viral"
    return escopo


def _escopo_pela_descricao(nome: str, descricao: str) -> str:
    t = RX_RUIDO.sub(" ", titulo(nome, descricao))
    d = RX_RUIDO.sub(" ", f"{nome} {descricao}")
    rna, ss = RX_RNA.search(t), RX_SSDNA.search(t) or RX_CRESS.search(t)
    # An RT title often names retroviruses, but never the RdRp of RNA viruses
    if RX_RETRO.search(t) and not re.search(r"rna-dependent rna pol|\brdrp", t, re.I):
        return "Retro (excluded)"
    if RX_DSDNA.search(t) and not (rna or ss):
        return "Other virus"
    # Host RNA-dependent RNA polymerases (plant/fungal RDR, QDE-1, EGO-1) are cellular
    if rna and not RX_VIRAL.search(t) and RX_CELULAR.search(d):
        return "Non-viral"
    if rna:
        return "RNA virus"
    if ss:
        return "ssDNA virus"
    # A viral title that names no group: the full description decides
    if RX_VIRAL.search(t):
        if RX_RETRO.search(d):
            return "Retro (excluded)"
        if RX_DSDNA.search(d):
            return "Other virus"
        if RX_RNA.search(d):
            return "RNA virus"
        if RX_SSDNA.search(d) or RX_CRESS.search(d):
            return "ssDNA virus"
        return "Other virus"
    return "Non-viral"


# ---------------------------------------------------------------------------
# NCBI Batch CD-Search
# ---------------------------------------------------------------------------
def _pedido(dados: dict, tentativas: int = 4) -> str:
    corpo = urllib.parse.urlencode(dados).encode()
    for k in range(tentativas):
        try:
            req = urllib.request.Request(URL, corpo, headers={"User-Agent": "EVE-pipeline/1.0 (anotar_cdsearch)"})
            return urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
        except Exception as erro:
            if k == tentativas - 1:
                raise RuntimeError(f"NCBI CD-Search unreachable: {erro}")
            time.sleep(10 * (k + 1))


def buscar_lote(seqs: "OrderedDict[str, str]", evalue: float, limite_s: int = 1800) -> list:
    """Raw hit rows ([query, type, pssm, from, to, evalue, bits, acc, name, incomplete, superfamily])."""
    consulta = "".join(f">{q}\n{s}\n" for q, s in seqs.items())
    r = _pedido({"queries": consulta, "db": "cdd", "smode": "auto", "useid1": "true",
                 "compbasedadj": "T", "filter": "true", "evalue": str(evalue), "maxhit": "500",
                 "dmode": "full", "tdata": "hits"})
    m = re.search(r"#cdsid\s+(\S+)", r)
    if not m:
        raise RuntimeError("CD-Search did not accept the request: " + r[:300])
    cdsid, inicio = m.group(1), time.time()
    while True:
        time.sleep(ESPERA)
        r = _pedido({"cdsid": cdsid, "tdata": "hits", "dmode": "full"})
        st = re.search(r"#status\s+(\d+)", r)
        st = st.group(1) if st else "?"
        if st == "0":
            break
        if st != "3":                      # 3 = still running; anything else is an error
            raise RuntimeError(f"CD-Search job {cdsid} failed (status {st}): {r[:300]}")
        if time.time() - inicio > limite_s:
            raise RuntimeError(f"CD-Search job {cdsid} still running after {limite_s} s")
    linhas = []
    for linha in r.splitlines():
        if not linha.startswith("Q#"):
            continue
        p = linha.split("\t")
        # "Q#3 - >ORF_ID" -> ORF_ID
        p[0] = p[0].split(">", 1)[-1].strip().split()[0]
        linhas.append(p[:11])
    return linhas


# ---------------------------------------------------------------------------
# Local search: rpsblast+ on the CDD database
# ---------------------------------------------------------------------------
def carregar_dominios_locais(caminho: str = DOMINIOS_PADRAO) -> dict:
    """{PSSM tag and accession: (accession, name, description, taxon, origin)} from Cdd_domains.tsv."""
    if not os.path.exists(caminho):
        raise RuntimeError(f"CDD domain table not found: {caminho}")
    dom = {}
    with open(caminho, encoding="utf-8", errors="replace", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k.strip(): (v or "").strip() for k, v in r.items() if k}
            info = (r["accession"], r["name"], r["description"], r.get("taxname", ""), r.get("origin", ""))
            dom[r["tag"]] = dom[r["accession"]] = info
    return dom


def _rpsblast() -> str:
    for exe in ("rpsblast", "rpsblast+"):
        if shutil.which(exe):
            return exe
    raise RuntimeError("rpsblast+ not found on the PATH (BLAST+ package)")


def buscar_local(seqs: "OrderedDict[str, str]", evalue: float, dominios: dict,
                 rpsdb: str = RPSDB_PADRAO, threads: int = 4) -> list:
    """Raw hit rows in the same layout as buscar_lote, from rpsblast+ on the local CDD."""
    if not os.path.exists(rpsdb + ".pal") and not os.path.exists(rpsdb + ".rps"):
        raise RuntimeError(f"Local CDD database not found: {rpsdb}")
    with tempfile.TemporaryDirectory(prefix="cdsearch_") as tmp:
        query = os.path.join(tmp, "query.fasta")
        with open(query, "w") as fh:
            fh.writelines(f">{q}\n{s}\n" for q, s in seqs.items())
        cmd = [_rpsblast(), "-query", query, "-db", rpsdb, "-evalue", str(evalue),
               "-outfmt", "6 qseqid sseqid qstart qend evalue bitscore",
               "-max_target_seqs", "500", "-num_threads", str(threads), "-seg", "no",
               "-comp_based_stats", "1"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"rpsblast+ failed: {r.stderr.strip()[:500]}")
    linhas = []
    for linha in r.stdout.splitlines():
        p = linha.split("\t")
        if len(p) < 6:
            continue
        tag = p[1].rsplit("|", 1)[-1].split(":")[-1]          # gnl|CDD|395756 -> 395756
        acc, nome = dominios.get(tag, (tag, "", "", "", ""))[:2]
        linhas.append([p[0], "rpsblast", tag, p[2], p[3], p[4], p[5], acc, nome, "", ""])
    return linhas


def ler_fasta(caminho: str) -> "OrderedDict[str, str]":
    seqs, nome = OrderedDict(), None
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                nome = linha[1:].split()[0]
                seqs[nome] = ""
            elif nome:
                seqs[nome] += linha.strip().replace("*", "")
    return OrderedDict((k, v) for k, v in seqs.items() if v)


def _assinatura(seq: str) -> str:
    return hashlib.md5(seq.encode()).hexdigest()[:16]


COLS_RAW = ["orf", "seq_md5", "hit_type", "pssm_id", "from", "to", "evalue", "bitscore",
            "accession", "short_name", "incomplete", "superfamily"]


# ---------------------------------------------------------------------------
# One species folder
# ---------------------------------------------------------------------------
def anotar_pasta(pasta: str, descricoes: dict, evalue: float = 0.01, forcar: bool = False,
                 local: bool = True, rpsdb: str = RPSDB_PADRAO, threads: int = 4,
                 fasta: str = None, dir_saida: str = None, prefixo: str = "") -> dict:
    """
    descricoes: {accession: (name, description[, taxon, origin])}.
    fasta/dir_saida/prefixo search another peptide set (e.g. the tblastn scan ORFs) and write
    <prefixo>cdsearch_*.tsv to dir_saida; by default the candidates of 04_candidatas_eves.
    """
    dir_c = dir_saida or os.path.join(pasta, "04_candidatas_eves")
    fasta = fasta or os.path.join(dir_c, "eve_candidatas_unificadas_pep.fasta")
    if not os.path.exists(fasta):
        return {"orfs": 0, "searched": 0}
    seqs = ler_fasta(fasta)
    # One cache per mode, so local and NCBI hits are never mixed
    raw_path = os.path.join(dir_c, prefixo + ("cdsearch_raw_local.tsv" if local else "cdsearch_raw.tsv"))

    # Cache: an ORF whose sequence was already searched is not sent again
    feitos, brutos = set(), []
    if os.path.exists(raw_path) and not forcar:
        with open(raw_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["orf"] in seqs and r["seq_md5"] == _assinatura(seqs[r["orf"]]):
                    feitos.add(r["orf"])
                    brutos.append(r)      # no-hit rows too: they mark the ORF as searched
    faltam = OrderedDict((q, s) for q, s in seqs.items() if q not in feitos)
    log(f"{os.path.basename(pasta.rstrip('/')) or pasta}: {len(seqs)} ORFs, "
        f"{len(feitos)} cached, {len(faltam)} to search")

    nomes = list(faltam)
    for i in range(0, len(nomes), LOTE):
        lote = OrderedDict((q, faltam[q]) for q in nomes[i:i + LOTE])
        if local:
            log(f"  batch {i // LOTE + 1}: {len(lote)} sequences, rpsblast+ on the local CDD...")
            achados = buscar_local(lote, evalue, descricoes, rpsdb=rpsdb, threads=threads)
        else:
            log(f"  batch {i // LOTE + 1}: {len(lote)} sequences sent to NCBI...")
            achados = buscar_lote(lote, evalue)
        com_hit = set()
        for p in achados:
            if p[0] not in lote:
                continue
            com_hit.add(p[0])
            brutos.append(dict(zip(COLS_RAW, [p[0], _assinatura(lote[p[0]])] + p[1:11])))
        # ORFs without any hit are recorded too, so they count as searched
        for q in lote:
            if q not in com_hit:
                brutos.append(dict(zip(COLS_RAW, [q, _assinatura(lote[q])] + [""] * 10)))
        _gravar(raw_path, COLS_RAW, brutos)      # saved after every batch: a crash keeps them
    if not faltam:
        _gravar(raw_path, COLS_RAW, brutos)

    # Domains with scope and description
    hits = []
    for r in brutos:
        if not r.get("hit_type"):
            continue
        info = descricoes.get(r["accession"])
        if not info:
            nome, desc, taxon, origem = r["short_name"], "", "", ""
        elif local:
            _acc, nome, desc, taxon, origem = info
        else:
            (nome, desc), taxon, origem = info, "", ""
        primeira = re.split(r"(?<=[a-z0-9)\]])\. ", desc, maxsplit=1)[0]
        hits.append({"orf": r["orf"], "hit_type": r["hit_type"], "accession": r["accession"],
                     "short_name": r["short_name"] or nome, "from": r["from"], "to": r["to"],
                     "evalue": r["evalue"], "bitscore": r["bitscore"],
                     "incomplete": (r["incomplete"] or "").strip(" -"),
                     "superfamily": (r["superfamily"] or "").strip(" -"),
                     "scope": escopo_do_dominio(r["short_name"] or nome, desc, taxon, origem),
                     "description": primeira[:200], "taxon": taxon})
    cols_hits = ["orf", "scope", "hit_type", "accession", "short_name", "from", "to", "evalue",
                 "bitscore", "incomplete", "superfamily", "description", "taxon"]
    _gravar(os.path.join(dir_c, prefixo + "cdsearch_hits.tsv"), cols_hits, hits)

    # One row per ORF. Superfamily rows only repeat a specific hit, so they count only
    # when the ORF has nothing more specific.
    peso_tipo = {"specific": 0, "non-specific": 1, "rpsblast": 1, "superfamily": 2}
    por_orf = {}
    for h in hits:
        por_orf.setdefault(h["orf"], []).append(h)
    resumo = []
    for orf in seqs:
        lista = por_orf.get(orf, [])
        especificos = [h for h in lista if h["hit_type"] != "superfamily"] or lista
        # A viral domain beaten on the same region by a stronger non-viral one is a cellular
        # protein that only resembles the viral family (e.g. DNA2/Helz helicases hitting the
        # coronavirus Nsp13 helicase): it does not make the ORF viral
        vencidos = {id(h): _rival(h, especificos) for h in especificos if h["scope"] in ESCOPOS_ALVO}
        vencidos = {k: v for k, v in vencidos.items() if v}
        superados = [h for h in especificos if id(h) in vencidos]
        especificos = [h for h in especificos if id(h) not in vencidos]
        escopos = {h["scope"] for h in especificos}
        escopo = next((e for e in ORDEM_ESCOPO if e in escopos), "No CDD hit")
        alvo = [h for h in especificos if h["scope"] in ESCOPOS_ALVO]
        ordem = sorted(alvo or especificos,
                       key=lambda h: (peso_tipo.get(h["hit_type"], 3), -float(h["bitscore"] or 0)))
        melhor = ordem[0] if ordem else None
        nomes_alvo = list(OrderedDict.fromkeys(h["short_name"] for h in sorted(
            alvo, key=lambda h: -float(h["bitscore"] or 0))))
        resumo.append({
            "orf": orf, "cdd_scope": escopo, "cdd_in_scope": "yes" if alvo else "no",
            "cdd_n_domains": len(especificos), "cdd_n_in_scope": len(alvo),
            "cdd_best": melhor["short_name"] if melhor else "",
            "cdd_best_accession": melhor["accession"] if melhor else "",
            "cdd_best_scope": melhor["scope"] if melhor else "",
            "cdd_best_evalue": melhor["evalue"] if melhor else "",
            "cdd_best_bitscore": melhor["bitscore"] if melhor else "",
            "cdd_best_description": melhor["description"] if melhor else "",
            "cdd_in_scope_domains": "; ".join(nomes_alvo),
            "cdd_outcompeted": "; ".join(OrderedDict.fromkeys(
                f"{h['short_name']} < {vencidos[id(h)]['short_name']}" for h in superados)),
        })
    cols_res = list(resumo[0]) if resumo else ["orf"]
    _gravar(os.path.join(dir_c, prefixo + "cdsearch_summary.tsv"), cols_res, resumo)
    n_alvo = sum(1 for r in resumo if r["cdd_in_scope"] == "yes")
    log(f"  {n_alvo}/{len(resumo)} ORFs with an RNA- or ssDNA-virus domain")
    return {"orfs": len(resumo), "searched": len(faltam), "in_scope": n_alvo,
            "summary_rows": resumo, "hit_rows": hits,
            "per_scope": {e: sum(1 for r in resumo if r["cdd_scope"] == e) for e in ORDEM_ESCOPO}}


def _rival(h: dict, hits: list, sobreposicao: float = 0.5):
    """Strongest non-viral hit covering >= half of h (or of itself) with a higher bitscore."""
    a, b, bits = int(h["from"]), int(h["to"]), float(h["bitscore"] or 0)
    melhor = None
    for o in hits:
        if o["scope"] in ESCOPOS_ALVO or o["scope"].startswith("Other") or float(o["bitscore"] or 0) <= bits:
            continue
        c, d = int(o["from"]), int(o["to"])
        if min(b, d) - max(a, c) + 1 >= sobreposicao * min(b - a + 1, d - c + 1):
            if not melhor or float(o["bitscore"]) > float(melhor["bitscore"]):
                melhor = o
    return melhor


def _gravar(caminho: str, colunas: list, linhas: list):
    with open(caminho + ".tmp", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=colunas, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(linhas)
    os.replace(caminho + ".tmp", caminho)


def pastas_especies(outdir: str) -> "OrderedDict[str, str]":
    especies = os.path.join(outdir, "especies")
    if os.path.isdir(especies):
        return OrderedDict((os.path.basename(p), p) for p in sorted(glob.glob(os.path.join(especies, "*")))
                           if os.path.isdir(p))
    return OrderedDict([("", outdir)])


def anotar(outdir: str, evalue: float = 0.01, forcar: bool = False, cddid: str = CDDID_PADRAO,
           modo: str = "local", rpsdb: str = RPSDB_PADRAO, dominios: str = DOMINIOS_PADRAO,
           threads: int = 4) -> dict:
    """modo: "local" (rpsblast+ on rpsdb, the default) or "online" (NCBI Batch CD-Search)."""
    local = modo == "local"
    if local:
        log(f"local CDD: {rpsdb}")
        descricoes = carregar_dominios_locais(dominios)
    else:
        descricoes = carregar_descricoes(cddid)
    total = {"orfs": 0, "searched": 0, "in_scope": 0, "mode": modo,
             "per_scope": {e: 0 for e in ORDEM_ESCOPO}}
    for _sp, pasta in pastas_especies(outdir).items():
        r = anotar_pasta(pasta, descricoes, evalue=evalue, forcar=forcar, local=local,
                         rpsdb=rpsdb, threads=threads)
        for k in ("orfs", "searched", "in_scope"):
            total[k] += r.get(k, 0)
        for e, n in (r.get("per_scope") or {}).items():
            total["per_scope"][e] += n
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--evalue", type=float, default=0.01, help="CD-Search E-value cutoff (default 0.01, as the web page)")
    ap.add_argument("--forcar", action="store_true", help="Search every ORF again, ignoring the cache")
    ap.add_argument("--modo", choices=("local", "online"), default="local",
                    help="local: rpsblast+ on the local CDD (default); online: NCBI Batch CD-Search")
    ap.add_argument("--rpsdb", default=RPSDB_PADRAO, help="Local CDD database for rpsblast+ (local mode)")
    ap.add_argument("--dominios", default=DOMINIOS_PADRAO,
                    help="Cdd_domains.tsv with the descriptions and taxa of the local CDD")
    ap.add_argument("-t", "--threads", type=int, default=4, help="rpsblast+ threads (local mode)")
    ap.add_argument("--cddid", default=CDDID_PADRAO, help="cddid_all.tbl(.gz) with the domain descriptions (online mode)")
    a = ap.parse_args()
    r = anotar(a.outdir, evalue=a.evalue, forcar=a.forcar, cddid=a.cddid, modo=a.modo,
               rpsdb=a.rpsdb, dominios=a.dominios, threads=a.threads)
    print(f"\n{r['orfs']} ORFs ({r['searched']} searched now) | in scope (RNA/ssDNA virus): {r['in_scope']}")
    for e, n in r["per_scope"].items():
        print(f"  {n:5}  {e}")


if __name__ == "__main__":
    main()
