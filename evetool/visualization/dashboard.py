#!/usr/bin/env python3
"""
EVE dashboard (the "Overview" page), built from scratch out of the pipeline outputs.

For each candidate EVE it computes:
  * length (nt) and %GC of the ORF;
  * delta GC between the EVE and the host genome;
  * delta GC between the flank and the EVE for 10, 25 and 50 kb windows on each side
    (the EVE itself stays out of the window);
  * viral protein of the best viral hit (general blastx; failing that, the primary RefSeq
    viral screen);
  * presence of a viral domain (hmmsearch over ALL candidates; no hit = "Non_alignment");
  * taxonomy of the best viral hit in the general blastx (nr by default), by ICTV rank.

Output: <outdir>/12_dashboard_web/dados/eves_dashboard.tsv (the table behind the charts),
        read by the single curation page, which is installed as 12_dashboard_web/index.html
        together with a data snapshot so it also opens straight from disk.

Usage (reusing an existing run):
  python -m evetool.visualization.dashboard -o Species_phlebotominae --species-tsv lista_especies.tsv \\
         --genomes-dir . --nomes lista_especies_corrigida.tsv -t 6
  python -m evetool.visualization.dashboard -o resultados_eve -i genoma.fna --species "Lutzomyia longipalpis"
"""
import os
import re
import csv
import glob
import json
import gzip
import shutil
import tarfile
import argparse
from collections import defaultdict, Counter
from evetool import config

BASE_HMM = config.HMM_DIR
HMM_PADRAO = f"{BASE_HMM}/COMBINED_DATABASE/combined_viral_curated.hmm"
ANOTACOES_HMM = f"{BASE_HMM}/COMBINED_DATABASE/combined_viral_annotations.tsv"
TAXDUMP = f"{BASE_HMM}/RVDB/taxdump.tar.gz"
PLOTLY_JS = os.path.join(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "curation", "web"), "assets", "plotly-2.24.1.min.js")

JANELAS_KB = (10, 25, 50)
RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "genus")
COLUNAS_GERAL = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore", "qlen", "slen", "positive", "ppos",
                 "qcovhsp", "stitle", "sseq"]
IG = {c: i for i, c in enumerate(COLUNAS_GERAL)}
REGEX_VIRAL = re.compile(r"vir(us|al|idae|ales|icota|icetes|ia)\b|phage|viroid", re.IGNORECASE)
REGEX_COORD = re.compile(r"\[(\d+)-(\d+)\]\(([+-])\)")

GRUPOS = {
    "1_VIRAL_DOMINANT": "Viral dominant",
    "2_VIRAL_PREFERRED": "Viral preferred",
    "3_AMBIGUOUS": "Ambiguous",
    "4_CELLULAR_DOMINANT": "Cellular dominant",
    "5_HMM_ONLY": "HMM only",
    "6_NO_EVIDENCE": "No evidence",
    "n/a": "Unclassified",
}
GRUPOS_FORA_POR_PADRAO = ["4_CELLULAR_DOMINANT", "6_NO_EVIDENCE"]

# Functional category of the viral protein, tested in order (first pattern that matches wins)
CATEGORIAS_PROTEINA = [
    ("RdRp / polymerase", r"rna-dependent rna polymerase|\brdrp\b|polymerase|replicase|\bl protein\b|large protein|\bpb[12]\b|\bns5\b"),
    ("Glycoprotein", r"^g$|glycoprotein|envelope|spike|\bg protein\b|\bgp\d+\b|hemagglutinin"),
    ("Nucleoprotein", r"nucleoprotein|nucleocapsid|\bn protein\b"),
    ("Capsid / structural protein", r"capsid|coat|\bvp\d\b|structural protein(?!.*non)"),
    ("Rep / NS1", r"\brep\b|replication.*protein|replication initiator|\bns1\b|non-?structural protein 1\b"),
    ("Helicase", r"helicase"),
    ("Polyprotein", r"polyprotein"),
    ("Matrix", r"matrix"),
    ("Movement", r"movement"),
    ("Non-structural", r"non-?structural|\bnsp?\d"),
    ("Hypothetical", r"hypothetical|uncharacterized|unknown|putative protein$|^protein$|^orf\d*"),
]
CATEGORIAS_PROTEINA = [(nome, re.compile(p, re.IGNORECASE)) for nome, p in CATEGORIAS_PROTEINA]
SEM_HIT_BLASTX = "No viral hit in the blastx (HMM only)"


def log(msg):
    print(f"[dashboard] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Taxonomy (NCBI taxdump, restricted to the Viruses branch)
# ---------------------------------------------------------------------------
class TaxonomiaViral:
    def __init__(self, taxdump=TAXDUMP, cache=None):
        self.pai, self.rank, self.nome, self.por_nome = {}, {}, {}, {}
        if cache and os.path.exists(cache):
            with open(cache) as fh:
                d = json.load(fh)
            self.pai, self.rank, self.nome, self.por_nome = d["pai"], d["rank"], d["nome"], d["por_nome"]
            return
        if not os.path.exists(taxdump):
            log(f"taxdump missing ({taxdump}); taxonomy from the title text only")
            return
        log("Reading the taxdump (Viruses branch)...")
        filhos, pai, rank = defaultdict(list), {}, {}
        with tarfile.open(taxdump) as tar:
            for linha in tar.extractfile("nodes.dmp"):
                p = linha.decode().split("\t|\t", 3)
                filhos[p[1]].append(p[0])
                pai[p[0]], rank[p[0]] = p[1], p[2]
            virais, pilha = set(), ["10239"]
            while pilha:
                t = pilha.pop()
                virais.add(t)
                pilha.extend(filhos.get(t, []))
            self.pai = {t: pai[t] for t in virais}
            self.rank = {t: rank[t] for t in virais}
            for linha in tar.extractfile("names.dmp"):
                p = linha.decode().split("\t|\t")
                if p[0] not in virais:
                    continue
                nome, classe = p[1], p[3].rstrip("\t|\n")
                if classe == "scientific name":
                    self.nome[p[0]] = nome
                self.por_nome.setdefault(nome.lower(), p[0])
        if cache:
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            with open(cache, "w") as fh:
                json.dump({"pai": self.pai, "rank": self.rank, "nome": self.nome, "por_nome": self.por_nome}, fh)

    def eh_viral(self, taxid):
        return taxid in self.pai

    def linhagem(self, taxid) -> dict:
        out, t, passos = {}, taxid, 0
        while t in self.pai and passos < 60:
            if self.rank.get(t) in RANKS:
                out[self.rank[t]] = self.nome.get(t, t)
            if t == "10239":
                break
            t, passos = self.pai[t], passos + 1
        return out


def interpretar_titulo(sseqid: str, stitle: str, tax: TaxonomiaViral):
    """(is_viral, taxid|None, organism, product) for nr, Swiss-Prot or pipeline RefSeq viral titles."""
    campos = stitle.split("|")
    # Pipeline RefSeq viral: "ACC |product [organism]|genome|species|genus|family|len"
    if len(campos) >= 6 and re.fullmatch(r"\s*(ss|ds)(RNA|DNA)[()+\-/ ]*\s*", campos[2] or ""):
        m = re.match(r"\s*(.*?)\s*\[([^\]]+)\]\s*$", campos[1])
        produto, organismo = (m.group(1), m.group(2)) if m else (campos[1].strip(), "")
        taxid = tax.por_nome.get(campos[3].strip().lower()) or tax.por_nome.get(organismo.lower())
        return True, taxid, organismo, produto
    # Swiss-Prot: "sp|ACC|ID product OS=org OX=taxid ..."
    m = re.search(r"OX=(\d+)", stitle)
    if m:
        os_ = re.search(r"OS=(.+?)(?: OX=|$)", stitle)
        corpo = stitle.split(" ", 1)[1] if " " in stitle else stitle
        produto = corpo.split(" OS=")[0]
        taxid = m.group(1)
        viral = tax.eh_viral(taxid) if tax.pai else bool(REGEX_VIRAL.search(os_.group(1) if os_ else ""))
        return viral, (taxid if viral else None), (os_.group(1) if os_ else ""), produto
    # nr: "ACC product [organism]"
    corpo = stitle.split(" ", 1)[1] if stitle.startswith(sseqid) and " " in stitle else stitle
    orgs = re.findall(r"\[([^\[\]]+)\]", corpo)
    organismo = orgs[-1] if orgs else ""
    produto = re.sub(r"\s*\[[^\[\]]+\]\s*$", "", corpo).strip()
    taxid = tax.por_nome.get(organismo.lower()) if organismo else None
    if taxid:
        return True, taxid, organismo, produto
    return bool(REGEX_VIRAL.search(organismo)), None, organismo, produto


def categoria_proteina(produto: str) -> str:
    p = re.sub(r"\b(putative|probable|partial|predicted|truncated)\b", "", produto or "", flags=re.I).strip(" ,;")
    for nome, rx in CATEGORIAS_PROTEINA:
        if rx.search(p):
            return nome
    return "Other"


# ---------------------------------------------------------------------------
# Pipeline inputs
# ---------------------------------------------------------------------------
def ler_fasta(caminho):
    """Yields (id, header, sequence)."""
    abrir = gzip.open if caminho.endswith(".gz") else open
    cab, partes = None, []
    with abrir(caminho, "rt") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if cab is not None:
                    yield cab.split()[0], cab, "".join(partes)
                cab, partes = linha[1:].strip(), []
            else:
                partes.append(linha.strip())
    if cab is not None:
        yield cab.split()[0], cab, "".join(partes)


def gc_de(seq: str):
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    acgt = gc + s.count("A") + s.count("T")
    return gc, acgt


def ler_tabela_especies(caminho):
    """{species column name: genome (file or accession)} from a two-column TSV in either order."""
    mapa = {}
    if not caminho:
        return mapa
    with open(caminho) as fh:
        for linha in fh:
            p = [t.strip() for t in linha.rstrip("\n").split("\t")]
            if len(p) < 2 or p[0].lower() in ("species", "especie", "espécie"):  # old headers still accepted
                continue
            a, b = p[0], p[1]
            genoma_primeiro = re.match(r"GC[AF]_", a) or a.endswith((".fna", ".fa", ".fasta", ".gz"))
            nome, genoma = (b, a) if genoma_primeiro else (a, b)
            mapa[nome] = genoma
    return mapa


def resolver_genoma(genoma: str, pastas):
    if genoma and os.path.exists(genoma):
        return genoma
    acc = os.path.basename(genoma or "")
    for pasta in pastas:
        if not pasta:
            continue
        achados = sorted(glob.glob(os.path.join(pasta, f"{acc}*.fna*")) +
                         glob.glob(os.path.join(pasta, f"{acc}*.fa*")) +
                         glob.glob(os.path.join(pasta, "**", f"{acc}*_genomic.fna*"), recursive=True))
        if achados:
            return achados[0]
    return None


def acc_de(genoma: str):
    m = re.search(r"(GC[AF]_\d+\.\d+)", os.path.basename(genoma or ""))
    return m.group(1) if m else None


def carregar_candidatas(pasta_especie):
    nt = os.path.join(pasta_especie, "04_candidatas_eves", "eve_candidatas_unificadas.fasta")
    pep = os.path.join(pasta_especie, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
    eves = []
    if not os.path.exists(nt):
        return eves, pep
    for orf, cab, seq in ler_fasta(nt):
        m = REGEX_COORD.search(cab)
        contig = orf.rsplit("_ORF.", 1)[0]
        ini, fim, fita = (int(m.group(1)), int(m.group(2)), m.group(3)) if m else (None, None, ".")
        gc, acgt = gc_de(seq)
        eves.append({"orf": orf, "contig": contig, "start": ini, "end": fim, "strand": fita,
                     "length": len(seq), "gc": round(100 * gc / acgt, 2) if acgt else None})
    return eves, pep


def metricas_genoma(caminho_genoma, eves):
    """Fills gc_genome, delta_gc_genome and delta_gc_flank_{N}kb on each EVE (in place)."""
    por_contig = defaultdict(list)
    for e in eves:
        if e["start"] is not None:
            por_contig[e["contig"]].append(e)
    gc_tot = acgt_tot = 0
    for contig, _cab, seq in ler_fasta(caminho_genoma):
        gc, acgt = gc_de(seq)
        gc_tot += gc
        acgt_tot += acgt
        for e in por_contig.get(contig, []):
            ini, fim = e["start"], e["end"]
            for kb in JANELAS_KB:
                n = kb * 1000
                g, a = gc_de(seq[max(0, ini - n):ini] + seq[fim:fim + n])
                e[f"gc_flank_{kb}kb"] = round(100 * g / a, 2) if a else None
    gc_genoma = round(100 * gc_tot / acgt_tot, 2) if acgt_tot else None
    for e in eves:
        e["gc_genome"] = gc_genoma
        e["delta_gc_genome"] = round(e["gc"] - gc_genoma, 2) if (e["gc"] is not None and gc_genoma) else None
        for kb in JANELAS_KB:
            f = e.get(f"gc_flank_{kb}kb")
            e[f"delta_gc_flank_{kb}kb"] = round(f - e["gc"], 2) if (f is not None and e["gc"] is not None) else None
    return gc_genoma


def melhores_hits_virais(caminho, tax):
    """{query: (bits, sseqid, taxid, organism, product, stitle)} - best viral hit per query in the general blastx."""
    melhor = {}
    if not caminho or not os.path.exists(caminho):
        return melhor
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) <= IG["stitle"]:
                continue
            try:
                bits = float(p[IG["bitscore"]])
            except ValueError:
                continue
            q = p[IG["qseqid"]]
            if q in melhor and melhor[q][0] >= bits:
                continue
            viral, taxid, org, prod = interpretar_titulo(p[IG["sseqid"]], p[IG["stitle"]], tax)
            if viral:
                melhor[q] = (bits, p[IG["sseqid"]], taxid, org, prod, p[IG["stitle"]])
    return melhor


def melhores_hits_triagem(caminho, tax):
    """Best hit of the primary screen (RefSeq viral): {orf: (bits, taxid, organism, product)}."""
    melhor = {}
    if not caminho or not os.path.exists(caminho):
        return melhor
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            try:
                bits = float(p[7])
            except ValueError:
                continue
            if p[0] in melhor and melhor[p[0]][0] >= bits:
                continue
            _v, taxid, org, prod = interpretar_titulo(p[2], p[4], tax)
            melhor[p[0]] = (bits, taxid, org, prod)
    return melhor


def ler_domtbl(caminho, destino):
    """Adds the hits of a domtblout table (target = ORF) into destino as {orf: {model: score}}."""
    if not caminho or not os.path.exists(caminho):
        return
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            p = linha.split()
            if len(p) < 14:
                continue
            try:
                score = float(p[7])
            except ValueError:
                continue
            atual = destino[p[0]].get(p[3], 0.0)
            destino[p[0]][p[3]] = max(atual, score)


def hmmsearch_todas(fastas_pep, saida, hmm_db, threads):
    """hmmsearch over every candidate (same rule as hmm_scanner: TC when calibrated, else E<=1e-3)."""
    if os.path.exists(saida) and os.path.getsize(saida) > 0:
        log(f"Reusing the HMM scan of every candidate: {saida}")
        return saida
    try:
        import pyhmmer
    except ImportError:
        log("pyhmmer unavailable; using the pipeline HMM tables only")
        return None
    if not hmm_db or not os.path.exists(hmm_db):
        log(f"HMM database not found ({hmm_db}); using the pipeline HMM tables only")
        return None
    alfabeto = pyhmmer.easel.Alphabet.amino()
    alvos = []
    for fasta in fastas_pep:
        if not os.path.exists(fasta):
            continue
        for orf, _cab, seq in ler_fasta(fasta):
            seq = seq.replace("*", "").upper()
            if seq:
                alvos.append(pyhmmer.easel.TextSequence(name=orf.encode(), sequence=seq).digitize(alfabeto))
    if not alvos:
        return None
    alvos = pyhmmer.easel.DigitalSequenceBlock(alfabeto, alvos)
    with pyhmmer.plan7.HMMFile(hmm_db) as hf:
        hmms = list(hf.optimized_profiles()) if hf.is_pressed() else list(hf)
    tem_tc = all(getattr(h.cutoffs, "trusted", None) is not None for h in hmms)
    criterio = {"bit_cutoffs": "trusted"} if tem_tc else {"E": 1e-3, "domE": 1e-3}
    log(f"hmmsearch: {len(hmms)} models x {len(alvos)} candidates ({'TC' if tem_tc else 'E<=1e-3'})...")
    os.makedirs(os.path.dirname(saida), exist_ok=True)
    primeiro = True
    with open(saida, "wb") as fh:
        for hits in pyhmmer.hmmer.hmmsearch(hmms, alvos, cpus=max(1, threads), **criterio):
            hits.write(fh, format="domains", header=primeiro)
            primeiro = False
    return saida


def anotacao_modelos(caminho=ANOTACOES_HMM):
    rotulos = {}
    if os.path.exists(caminho):
        with open(caminho, encoding="utf-8", errors="replace") as fh:
            for linha in csv.DictReader(fh, delimiter="\t"):
                rotulos[linha["HMM_ID"]] = f'{linha.get("Family") or linha.get("Tax_Name") or "?"} · {linha.get("Genome_Type") or "?"}'
    return rotulos


def ler_classificacao(outdir):
    caminho = os.path.join(outdir, "eves_classified.tsv")
    grupos = {}
    if os.path.exists(caminho):
        with open(caminho) as fh:
            for linha in csv.DictReader(fh, delimiter="\t"):
                grupos[linha["query"]] = linha["group"]
    return grupos


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def gerar_dashboard(outdir, genomas=None, nomes=None, hmm_db=HMM_PADRAO, hmm_dir=None,
                    threads=4, pastas_genomas=(), titulo="EVEs - Overview", recalcular=False):
    """
    genomas: {species folder (or "" in single-species mode): genome path or accession}
    nomes:   {accession: name to display} (fixes labels from older folders)
    """
    genomas, nomes = genomas or {}, nomes or {}
    destino = os.path.join(outdir, "12_dashboard_web")
    dados = os.path.join(destino, "dados")
    os.makedirs(dados, exist_ok=True)

    multi = os.path.isdir(os.path.join(outdir, "especies"))
    pastas = sorted(glob.glob(os.path.join(outdir, "especies", "*"))) if multi else [outdir]
    if hmm_dir and not (hmm_db and os.path.exists(hmm_db)):  # combined database built by hmm_scanner
        combinado = os.path.join(pastas[0], "03_triagem_hmm", "banco_viral_customizado.hmm")
        hmm_db = combinado if os.path.exists(combinado) else hmm_db

    tax = TaxonomiaViral(cache=os.path.join(dados, ".cache_taxonomia_viral.json"))
    geral = next((c for c in (os.path.join(outdir, "04_candidatas_eves", "diamond_blastx_referencia_geral.tsv"),
                              os.path.join(outdir, "02_alinhamento_diamond", "diamond_blastx_referencia_geral.tsv"))
                  if os.path.exists(c)), None)
    log(f"general blastx: {geral or 'not found'}")
    hits_gerais = melhores_hits_virais(geral, tax)
    grupos = ler_classificacao(outdir)

    # Candidates and GC metrics (cached per species)
    registros, fastas_pep, tabelas_hmm = [], [], []
    for pasta in pastas:
        chave = os.path.basename(pasta) if multi else ""
        eves, pep = carregar_candidatas(pasta)
        if not eves:
            continue
        fastas_pep.append(pep)
        tabelas_hmm.append(os.path.join(pasta, "03_triagem_hmm", "triagem_profunda_hmm.tsv"))
        genoma_cfg = genomas.get(chave) or (next(iter(genomas.values())) if not multi and genomas else None)
        acc = acc_de(genoma_cfg) if genoma_cfg else None
        rotulo = nomes.get(acc) or nomes.get(chave) or chave or "Genome"
        rotulo = rotulo.replace("_", " ")

        cache = os.path.join(dados, f"gc_{chave or 'genome'}.json")
        em_cache = {}
        if os.path.exists(cache) and not recalcular:
            with open(cache) as fh:
                em_cache = json.load(fh)
        # A cache written before the columns were renamed holds the old key names, so it
        # has to be recomputed rather than merged into records that expect the new ones.
        esperadas = {"gc_genome", f"gc_flank_{JANELAS_KB[0]}kb"}
        if em_cache and not esperadas.issubset(set(next(iter(em_cache.values()), {}))):
            em_cache = {}
        if em_cache and all(e["orf"] in em_cache for e in eves):
            for e in eves:
                e.update(em_cache[e["orf"]])
        else:
            caminho = resolver_genoma(genoma_cfg, [*pastas_genomas, os.path.join(outdir, "genomas_ncbi")])
            if caminho:
                log(f"{rotulo}: genome and flank GC ({os.path.basename(caminho)})...")
                metricas_genoma(caminho, eves)
                with open(cache, "w") as fh:
                    json.dump({e["orf"]: {k: v for k, v in e.items() if "gc_" in k} for e in eves}, fh)
            else:
                log(f"{rotulo}: genome not located ({genoma_cfg}); delta GC unavailable")

        triagem = melhores_hits_triagem(next((c for c in (
            os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond_refseq_filtrado.tsv"),
            os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond.tsv")) if os.path.exists(c)), None), tax)
        for e in eves:
            query = f"{chave}|{e['orf']}" if multi else e["orf"]
            e["query"], e["species"], e["group"] = query, rotulo, grupos.get(query, "n/a")
            hg = hits_gerais.get(query)
            e["tax"] = tax.linhagem(hg[2]) if hg and hg[2] else ({} if hg else None)
            if hg and not e["tax"]:
                fam = re.search(r"\b([A-Z][a-z]+viridae)\b", hg[5])
                e["tax"] = {"family": fam.group(1)} if fam else {}
            e["general_hit"] = f"{hg[1]} {hg[4]} [{hg[3]}]" if hg else ""
            fonte = hg[4] if hg else (triagem[e["orf"]][3] if e["orf"] in triagem else None)
            e["product"] = fonte or ""
            e["protein"] = categoria_proteina(fonte) if fonte else SEM_HIT_BLASTX
        registros.extend(eves)
        log(f"{rotulo}: {len(eves)} candidates")

    if not registros:
        raise SystemExit("No candidate found in " + outdir)

    # Viral domains: hmmsearch over every candidate + the pipeline tables
    hmm_hits = defaultdict(dict)
    ler_domtbl(hmmsearch_todas(fastas_pep, os.path.join(dados, "hmm_todas_candidatas.domtbl"), hmm_db, threads), hmm_hits)
    for t in tabelas_hmm:
        ler_domtbl(t, hmm_hits)
    rotulos_hmm = anotacao_modelos()
    for e in registros:
        modelos = sorted(hmm_hits.get(e["orf"], {}).items(), key=lambda kv: -kv[1])
        e["domain"] = "Viral domain" if modelos else "Non_alignment"
        e["hmm_models"] = [f"{m} ({rotulos_hmm.get(m, '?')}) {s:.1f} bits" for m, s in modelos[:3]]

    # One row per EVE
    colunas = ["query", "species", "group", "contig", "start", "end", "strand", "length", "gc", "gc_genome",
               "delta_gc_genome"] + [f"gc_flank_{kb}kb" for kb in JANELAS_KB] + \
              [f"delta_gc_flank_{kb}kb" for kb in JANELAS_KB] + \
              ["protein", "product", "domain", "hmm_models", "general_hit"] + [f"tax_{r}" for r in RANKS]
    with open(os.path.join(dados, "eves_dashboard.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(colunas)
        for e in registros:
            linha = dict(e, hmm_models="; ".join(e["hmm_models"]),
                         **{f"tax_{r}": (e["tax"] or {}).get(r, "") for r in RANKS})
            w.writerow(["" if linha.get(c) is None else linha.get(c) for c in colunas])

    # The overview charts now live on the single curation page, which reads these metrics
    # (eves_dashboard.tsv) through the server, or from the snapshot when opened from disk.
    import sys
    from evetool.curation.curation import instalar_pagina
    saida = instalar_pagina(outdir)
    log(f"Dashboard: {saida}")
    return saida


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("-i", "--input", help="Host genome (single-species mode)")
    ap.add_argument("--species", help="Species name to display (single-species mode)")
    ap.add_argument("--species-tsv", help="TSV used in the multi-species run (folder -> genome/accession)")
    ap.add_argument("--genomes-dir", action="append", default=[], help="Folder(s) to search for the genomes by accession")
    ap.add_argument("--nomes", help="TSV accession -> correct species name (e.g. lista_especies_corrigida.tsv)")
    ap.add_argument("--hmm-db", default=HMM_PADRAO, help="HMM database to scan every candidate against")
    ap.add_argument("--hmm-dir", help="Folder of .hmm files (uses the combined database built by the pipeline)")
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--recalcular", action="store_true", help="Ignore the GC cache")
    a = ap.parse_args()

    genomas = ler_tabela_especies(a.species_tsv)
    if a.input:
        genomas = {"": a.input}
    nomes = {}
    if a.nomes:
        nomes = {acc_de(g) or g: n for n, g in ler_tabela_especies(a.nomes).items()}
    if a.species:
        nomes[acc_de(a.input or "") or ""] = a.species
    pastas = a.genomes_dir or ([os.path.dirname(os.path.abspath(a.species_tsv))] if a.species_tsv else [])
    gerar_dashboard(a.outdir, genomas=genomas, nomes=nomes, hmm_db=a.hmm_db, hmm_dir=a.hmm_dir,
                    threads=a.threads, pastas_genomas=pastas, recalcular=a.recalcular)


if __name__ == "__main__":
    main()
