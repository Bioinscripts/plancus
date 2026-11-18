#!/usr/bin/env python3
"""
Score-based classification of the candidate EVEs, from the DIAMOND run against the
general database.

For each candidate it compares, on the same query, the best viral and the best cellular hit:

    razao_v_c = bitscore_viral / bitscore_celular

Viral evidence is the higher bitscore between the primary screening (DIAMOND vs RefSeq
viral) and the best viral hit in the general database. The cellular hit always comes from
the general database.

    VIRAL_DOMINANT     ratio >= 1.5, or no cellular hit at all
    VIRAL_PREFERRED    1.1 <= ratio < 1.5
    AMBIGUOUS          0.9 <= ratio < 1.1
    CELLULAR_DOMINANT  ratio < 0.9                     -> most likely a host gene
    HMM_ONLY           no viral hit in the general database, but an HMM hit above the calibrated TC
    NO_EVIDENCE        neither a viral hit nor an HMM hit

Within the viral groups, the conservation tier comes from the best viral hit:

    CONSERVED      pident >= 50 and qcovhsp >= 50
    INTERMEDIATE   pident >= 30
    DIVERGENT      pident <  30

The rank of the first viral hit is kept as an auxiliary column: how many cellular
proteins come before the first virus.

Uso:
  python -m evetool.validation.classify_eves -o results_eve                    # single species
  python -m evetool.validation.classify_eves -o results_multispecies --multi   # multi-species mode
"""
import os
import re
import csv
import glob
import tarfile
import argparse
from collections import defaultdict, Counter
from evetool import config

TAXDUMP = config.hmm_path("RVDB", "taxdump.tar.gz")

# Columns of the general blastx, as written by sequence_searcher.buscar_proteinas_referencia_diamond
COLUNAS_GERAL = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore", "qlen", "slen", "positive", "ppos",
                 "qcovhsp", "stitle", "sseq"]
I = {nome: idx for idx, nome in enumerate(COLUNAS_GERAL)}

# RefSeq viral header used in the database: ACC |product [organism]|genome|species|genus|family|len
REGEX_ORGANISMO_VIRAL = re.compile(
    r"vir(us|al|idae|ales|icota|icetes|ia)\b|phage|viroid|satellite nucleic acid", re.IGNORECASE)


def log(msg):
    print(f"[classificar] {msg}", flush=True)


def taxids_virais(taxdump=TAXDUMP) -> set:
    """Every taxid descending from Viruses (10239)."""
    if not os.path.exists(taxdump):
        return set()
    filhos = defaultdict(list)
    with tarfile.open(taxdump) as tar:
        for linha in tar.extractfile("nodes.dmp"):
            p = linha.decode().split("\t|\t")
            filhos[p[1]].append(p[0])
    virais, pilha = set(), ["10239"]
    while pilha:
        t = pilha.pop()
        virais.add(t)
        pilha.extend(filhos.get(t, []))
    return virais


def eh_viral(sseqid: str, stitle: str, virais: set) -> bool:
    """Decide whether a hit is viral: Swiss-Prot taxonomy (OX=) when present, otherwise the title text."""
    m = re.search(r"OX=(\d+)", stitle)
    if m and virais:
        return m.group(1) in virais
    if sseqid.startswith("sp|"):  # Swiss-Prot without a readable OX=
        return bool(REGEX_ORGANISMO_VIRAL.search(stitle))
    # Pipeline RefSeq viral: field 3 of the title is the genome type (ssRNA(+), dsRNA, ...)
    campos = stitle.split("|")
    if len(campos) > 2 and re.fullmatch(r"\s*(ss|ds)(RNA|DNA)[()+\-/ ]*\s*", campos[2] or ""):
        return True
    return bool(REGEX_ORGANISMO_VIRAL.search(stitle))


def carregar_hits_gerais(caminho: str, virais: set, especies_ignoradas=()) -> dict:
    """{query: {'viral': melhor hit viral, 'celular': melhor hit celular, 'rank_viral': n, 'n_hits': n}}"""
    bruto = defaultdict(list)
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < len(COLUNAS_GERAL) - 1:
                continue
            sseqid = p[I["sseqid"]]
            # Skip the species' own proteome: the ORF would find itself at ~100% identity
            if any(sseqid.startswith(pref) for pref in especies_ignoradas):
                continue
            try:
                bits = float(p[I["bitscore"]])
            except ValueError:
                continue
            bruto[p[I["qseqid"]]].append((bits, p))

    resultado = {}
    for query, lista in bruto.items():
        lista.sort(key=lambda t: -t[0])
        melhor_viral = best_cellular = None
        rank_viral = 0
        for pos, (_bits, p) in enumerate(lista, start=1):
            viral = eh_viral(p[I["sseqid"]], p[I["stitle"]], virais)
            if viral and melhor_viral is None:
                melhor_viral, rank_viral = p, pos
            elif not viral and best_cellular is None:
                best_cellular = p
            if melhor_viral is not None and best_cellular is not None:
                break
        resultado[query] = {"viral": melhor_viral, "celular": best_cellular,
                            "rank_viral": rank_viral, "n_hits": len(lista)}
    return resultado


def carregar_resumo_hmm(caminho: str) -> dict:
    """
    hmm_resumo.tsv (hmm_domains.py): varredura do banco calibrado sobre TODAS as candidatas.
    This is the preferred source; the screening table only covers ORFs without a DIAMOND hit.
    """
    import csv as _csv
    resumo = {}
    if caminho and os.path.exists(caminho):
        with open(caminho, encoding="utf-8", errors="replace") as fh:
            for r in _csv.DictReader(fh, delimiter="\t"):
                if int(r.get("n_domains") or 0):
                    resumo[r.get("orf", "")] = r
    return resumo


def carregar_hmm(caminho: str) -> dict:
    """{orf: (best_model, best_score, n_models)} read from the domtblout."""
    melhor, modelos = {}, defaultdict(set)
    if not caminho or not os.path.exists(caminho):
        return {}
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            p = linha.split()
            if len(p) < 14:
                continue
            orf, model = p[0], p[3]
            try:
                score = float(p[13])
            except ValueError:
                continue
            modelos[orf].add(model)
            if orf not in melhor or score > melhor[orf][1]:
                melhor[orf] = (model, score)
    return {orf: (m, s, len(modelos[orf])) for orf, (m, s) in melhor.items()}


_TAXONOMIA = None
RANKS_TAXON = ("family", "order", "class", "phylum", "kingdom", "realm")


def _taxonomia():
    """NCBI taxonomy of the Viruses branch (taxdump), loaded once and cached on disk."""
    global _TAXONOMIA
    if _TAXONOMIA is None:
        try:
            from evetool.visualization.dashboard import TaxonomiaViral
            cache = os.path.join(os.path.expanduser("~"), ".cache", "eve_tool", "taxonomia_viral.json")
            _TAXONOMIA = TaxonomiaViral(cache=cache)
        except Exception as erro:
            log(f"[warning] viral taxonomy unavailable ({erro}); families only from the RefSeq header")
            _TAXONOMIA = False
    return _TAXONOMIA


def viral_family(stitle: str) -> str:
    """
    Family from the RefSeq viral header (field 6). Viruses that RefSeq leaves without a family
    (e.g. 'Wuhan insect virus 23') get the lowest rank NCBI assigns them (order, class... down
    to the realm, often 'Riboviria'), never the virus name, which would read as a family of its own.
    """
    campos = stitle.split("|")
    if len(campos) >= 6 and campos[5].strip().endswith(("viridae", "virales", "viricetes")):
        return campos[5].strip()
    nomes = [campos[3].strip()] if len(campos) >= 6 and campos[3].strip() else []
    m = re.search(r"\[([^\]]+)\]", stitle) or re.search(r"OS=(.+?)(?: OX=| GN=| PE=|$)", stitle)
    if m:
        nomes.append(m.group(1).strip())
    if not nomes:
        return "n/a"
    tax = _taxonomia()
    if not tax or not tax.por_nome:
        return nomes[-1]                       # no taxdump: the organism, as before
    for nome in nomes:
        taxid = tax.por_nome.get(nome.lower())
        if taxid:
            lin = tax.linhagem(taxid)
            for rank in RANKS_TAXON:
                if lin.get(rank):
                    return lin[rank]
    return "unclassified"


# Functional classes of viral protein, in priority order: the first match wins.
# They group the RefSeq products, whose wording varies a lot, into categories that are
# comparable across viral families.
CLASSES_PROTEINA = [
    ("RdRp / replicase", r"rna[- ]dependent rna[- ]?polymerase|\brdrp\b|replicase|rna[- ]replicase"),
    ("polymerase",       r"\bpolymerase\b|\bl protein\b|\bpol\b"),
    ("glycoprotein",    r"glycoprotein|\benvelope\b|\bspike\b|\bg protein\b|\bgn\b|\bgc\b"),
    ("nucleoprotein",   r"nucleoprotein|nucleocapsid|\bn protein\b"),
    ("capsid / coat",  r"capsid|\bcoat\b|\bvp\d|virion protein"),
    ("matrix",           r"\bmatrix\b|\bm protein\b"),
    ("phosphoprotein",    r"phosphoprotein|\bp protein\b"),
    ("protease",         r"protease|peptidase|proteinase"),
    ("helicase",         r"helicase"),
    ("methyltransferase", r"methyltransferase"),
    ("polyprotein",     r"polyprotein"),
    ("non-structural",   r"non[- ]?structural|\bns\d|\bnss\b|\bnsm\b"),
    ("movement",        r"movement protein"),
    ("hypothetical",       r"hypothetical|uncharacteri|unknown|putative protein|\borf\s*\w*\b"),
]
_CLASSES = [(nome, re.compile(padrao, re.IGNORECASE)) for nome, padrao in CLASSES_PROTEINA]


def classe_da_proteina(product: str) -> str:
    """Functional category of the viral product; anything unmatched becomes 'other'."""
    texto = re.sub(r"\[.*?\]", " ", product or "")          # drop the organism in square brackets
    texto = re.sub(r"^\s*MAG:\s*", "", texto, flags=re.I)   # prefix used by metagenome-assembled genomes
    for nome, rx in _CLASSES:
        if rx.search(texto):
            return nome
    return "other" if texto.strip() else "n/a"


def produto_de(stitle: str) -> str:
    """Product name in the RefSeq viral header: 'ACC |product [organism]|genome|...'."""
    campos = stitle.split("|")
    bruto = campos[1].strip() if len(campos) > 1 else stitle
    return re.sub(r"\s*\[.*?\]\s*$", "", bruto).strip()


# Colunas de triagem_viral_diamond.tsv (sequence_searcher.preparar_e_buscar)
I_VIRAL = {"qseqid": 0, "qlen": 1, "sseqid": 2, "slen": 3, "stitle": 4,
           "pident": 5, "evalue": 6, "bitscore": 7, "qcovhsp": 8}


def carregar_triagem_viral(caminho: str) -> dict:
    """Best hit per ORF in the primary RefSeq viral screening: {orf: (bits, pident, qcov, family, sseqid)}."""
    melhor = {}
    if not caminho or not os.path.exists(caminho):
        return melhor
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) <= I_VIRAL["qcovhsp"]:
                continue
            try:
                bits = float(p[I_VIRAL["bitscore"]])
                pid = float(p[I_VIRAL["pident"]])
                cov = float(p[I_VIRAL["qcovhsp"]])
            except ValueError:
                continue
            orf = p[I_VIRAL["qseqid"]]
            if orf not in melhor or bits > melhor[orf][0]:
                melhor[orf] = (bits, pid, cov, viral_family(p[I_VIRAL["stitle"]]),
                               p[I_VIRAL["sseqid"]], produto_de(p[I_VIRAL["stitle"]]))
    return melhor


def classificar(hits_gerais: dict, triagem_viral: dict, hmm: dict, candidatas: list, cortes,
                resumo_hmm: dict = None) -> list:
    resumo_hmm = resumo_hmm or {}
    linhas = []
    for query in candidatas:
        orf = query.split("|", 1)[1] if "|" in query else query
        species = query.split("|", 1)[0] if "|" in query else ""
        info = hits_gerais.get(query) or hits_gerais.get(orf) or {}
        h_viral, h_celular = info.get("viral"), info.get("celular")
        hmm_hit = hmm.get(orf)
        res_hmm = resumo_hmm.get(orf)
        primario = triagem_viral.get(orf)

        # Viral evidence: the better of the primary RefSeq viral screening and the general-database viral hit
        bits_geral = float(h_viral[I["bitscore"]]) if h_viral else 0.0
        bits_v = max(bits_geral, primario[0] if primario else 0.0)
        bits_c = float(h_celular[I["bitscore"]]) if h_celular else 0.0
        razao = (bits_v / bits_c) if (bits_v and bits_c) else (float("inf") if bits_v else 0.0)

        if bits_v:
            if razao >= cortes.razao_dominante or not bits_c:
                group = "1_VIRAL_DOMINANT"
            elif razao >= cortes.razao_preferido:
                group = "2_VIRAL_PREFERRED"
            elif razao >= cortes.razao_ambiguo:
                group = "3_AMBIGUOUS"
            else:
                group = "4_CELLULAR_DOMINANT"
        elif hmm_hit:
            group = "5_HMM_ONLY"
        else:
            group = "6_NO_EVIDENCE"

        # Conservation metrics come from whichever source gave the best viral bitscore
        if primario and primario[0] >= bits_geral:
            pid, cov = primario[1], primario[2]
            family, protein, product = primario[3], primario[4], primario[5]
        elif h_viral:
            pid = float(h_viral[I["pident"]])
            cov = float(h_viral[I["qcovhsp"]]) if h_viral[I["qcovhsp"]] else 0.0
            family, protein = viral_family(h_viral[I["stitle"]]), h_viral[I["sseqid"]]
            product = produto_de(h_viral[I["stitle"]])
        else:
            pid = cov = 0.0
            family = protein = product = "n/a"

        if bits_v:
            if pid >= cortes.pident_conservada and cov >= cortes.qcov_conservada:
                conservation = "CONSERVED"
            elif pid >= cortes.pident_intermediaria:
                conservation = "INTERMEDIATE"
            else:
                conservation = "DIVERGENT"
        else:
            conservation = "n/a"

        linhas.append({
            "query": query, "species": species, "orf": orf,
            "group": group, "conservation": conservation,
            "viral_cellular_ratio": "inf" if razao == float("inf") else round(razao, 3),
            "viral_bitscore": round(bits_v, 1), "cellular_bitscore": round(bits_c, 1),
            "bitscore_delta": round(bits_v - bits_c, 1),
            "viral_pident": round(pid, 1), "viral_qcov": round(cov, 1),
            "viral_family": family, "viral_protein": protein,
            "viral_product": product[:90], "protein_class": classe_da_proteina(product),
            "viral_evidence_source": ("triagem_primaria" if primario and primario[0] >= bits_geral
                                      else "banco_geral" if h_viral else "n/a"),
            "first_viral_rank": info.get("rank_viral", 0),
            "n_hits_general": info.get("n_hits", 0),
            "best_cellular": h_celular[I["sseqid"]] if h_celular else "n/a",
            "cellular_title": (h_celular[I["stitle"]][:90] if h_celular else "n/a"),
            # Single HMM source: the full annotation when it exists, otherwise whatever the
            # screening found (which only covers ORFs without a DIAMOND hit).
            "hmm_n_domains": int(res_hmm["n_domains"]) if res_hmm else (hmm_hit[2] if hmm_hit else 0),
            "hmm_best_model": (res_hmm["best_model"] if res_hmm
                                  else (hmm_hit[0] if hmm_hit else "n/a")),
            "hmm_best_score": float(res_hmm["best_score"]) if res_hmm
                                else (round(hmm_hit[1], 1) if hmm_hit else 0.0),
            "hmm_consensus_family": res_hmm.get("hmm_consensus_family", "") if res_hmm else "",
            "hmm_model_genome": res_hmm.get("best_model_genome", "") if res_hmm else "",
            "hmm_consensus_support": res_hmm.get("consensus_support", "") if res_hmm else "",
            "hmm_source": "full_annotation" if res_hmm else ("screening" if hmm_hit else "n/a"),
        })
    return linhas


def agrupar_por_locus(linhas: list, fasta_pep: str, sobreposicao_min: float = 0.5) -> int:
    """
    Collapse ORFs that describe the SAME genomic segment into one locus.

    orfipy reports every reading frame, so one insertion shows up as several overlapping
    ORFs (the same protein, different frames, lower bitscores). Two ORFs join the same
    locus when they overlap by at least `sobreposicao_min` of the shorter one; merely
    touching at the edges is not enough, otherwise two neighbouring genes would chain
    into a single locus through a few shared bases. Strand is ignored on purpose: the
    reverse frame of an insertion is the same segment, not a second one.

    Adds `locus` and `locus_representative` to each row and returns the number of loci.
    The rows are kept: only the representative is counted, nothing is discarded.
    """
    coords = {}
    for cab in (l[1:].strip() for l in open(fasta_pep) if l.startswith(">")) if fasta_pep else []:
        orf = cab.split()[0]
        m = re.search(r"\[(\d+)-(\d+)\]", cab)
        if m:
            coords[orf] = (orf.rsplit("_ORF.", 1)[0], int(m.group(1)), int(m.group(2)))

    por_contig = defaultdict(list)
    for linha in linhas:
        c = coords.get(linha["orf"])
        if c:
            por_contig[(linha.get("species", ""), c[0])].append((c[1], c[2], linha))
        else:                                   # no coordinates: a locus of its own
            por_contig[("", linha["orf"])].append((0, 0, linha))

    n_loci = 0
    for (species, contig), itens in sorted(por_contig.items()):
        grupos = []
        for ini, fim, linha in sorted(itens, key=lambda t: (t[0], -t[1])):
            destino = None
            for g in grupos:
                for ini2, fim2, _ in g:
                    sobre = min(fim, fim2) - max(ini, ini2)
                    if sobre > 0 and sobre >= sobreposicao_min * max(1, min(fim - ini, fim2 - ini2)):
                        destino = g
                        break
                if destino:
                    break
            (destino if destino is not None else grupos.append([]) or grupos[-1]).append((ini, fim, linha))

        for g in grupos:
            n_loci += 1
            nome = f"{contig}:{min(i for i, _, _ in g)}-{max(f for _, f, _ in g)}"
            # The representative is the strongest evidence: viral bitscore first, then the
            # HMM score (which is what carries an HMM-only segment), then the longest ORF.
            melhor = max(g, key=lambda t: (float(t[2].get("viral_bitscore") or 0),
                                           float(t[2].get("hmm_best_score") or 0),
                                           t[1] - t[0]))
            for ini, fim, linha in g:
                linha["locus"] = nome
                linha["locus_representative"] = (linha is melhor[2])
    return n_loci


def ler_ids_fasta(caminho: str) -> list:
    if not caminho or not os.path.exists(caminho):
        return []
    with open(caminho) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def primeiro_existente(*caminhos):
    for c in caminhos:
        if c and os.path.exists(c) and os.path.getsize(c) > 0:
            return c
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True, help="Pipeline results folder")
    ap.add_argument("--multi", action="store_true", help="Multi-species mode (reads outdir/especies/*)")
    ap.add_argument("--geral-tsv", help="DIAMOND TSV against the general database (auto-detected when omitted)")
    ap.add_argument("--ignorar-proteoma", nargs="*", default=[],
                    help="sseqid prefixes to ignore (e.g. GCF_024334085.1 for the species' own proteome)")
    ap.add_argument("--razao-dominante", type=float, default=1.5)
    ap.add_argument("--razao-preferido", type=float, default=1.1)
    ap.add_argument("--razao-ambiguo", type=float, default=0.9)
    ap.add_argument("--pident-conservada", type=float, default=50.0)
    ap.add_argument("--qcov-conservada", type=float, default=50.0)
    ap.add_argument("--pident-intermediaria", type=float, default=30.0)
    ap.add_argument("--sobreposicao-locus", type=float, default=0.5,
                    help="Overlap fraction of the shorter ORF for two ORFs to count as the "
                         "same locus (default: 0.5; 0 disables the grouping)")
    cfg = ap.parse_args()

    geral = cfg.geral_tsv or primeiro_existente(
        os.path.join(cfg.outdir, "04_candidatas_eves", "diamond_blastx_referencia_geral.tsv"),
        os.path.join(cfg.outdir, "02_alinhamento_diamond", "diamond_blastx_referencia_geral.tsv"),
    )
    if not geral:
        raise SystemExit("General blastx TSV not found; pass --geral-tsv.")

    log(f"Reading viral taxonomy from {os.path.basename(TAXDUMP)}...")
    virais = taxids_virais()
    log(f"Reading general hits from {geral}...")
    hits_gerais = carregar_hits_gerais(geral, virais, tuple(cfg.ignorar_proteoma))

    # Candidates and HMM tables: one folder per species in multi mode, the root otherwise
    if cfg.multi:
        pastas = sorted(glob.glob(os.path.join(cfg.outdir, "especies", "*")))
    else:
        pastas = [cfg.outdir]

    todas = []
    for pasta in pastas:
        species = os.path.basename(pasta) if cfg.multi else ""
        cand = ler_ids_fasta(primeiro_existente(
            os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")))
        if not cand:
            continue
        hmm = carregar_hmm(primeiro_existente(
            os.path.join(pasta, "03_triagem_hmm", "triagem_profunda_hmm.tsv")))
        triagem_viral = carregar_triagem_viral(primeiro_existente(
            os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond_refseq_filtrado.tsv"),
            os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond.tsv")))
        queries = [f"{species}|{c}" if species else c for c in cand]
        resumo_hmm = carregar_resumo_hmm(primeiro_existente(
            os.path.join(pasta, "04_candidatas_eves", "hmm_summary.tsv")))
        linhas = classificar(hits_gerais, triagem_viral, hmm, queries, cfg, resumo_hmm)
        for linha in linhas:
            linha["species"] = species or linha["species"]
        # Overlapping reading frames of one insertion become a single locus; the ORFs of
        # the species are grouped together so a locus never spans two species.
        n_loci = agrupar_por_locus(linhas, primeiro_existente(
            *sorted(glob.glob(os.path.join(pasta, "01_orfs_cds", "*_orfs_pep.fasta"))),
            os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")),
            cfg.sobreposicao_locus)
        log(f"{species or 'single species'}: {len(linhas)} ORFs in {n_loci} distinct loci")
        todas.extend(linhas)

    if not todas:
        raise SystemExit("No candidate found to classify.")

    destino = os.path.join(cfg.outdir, "eves_classified.tsv")
    with open(destino, "w", newline="") as fh:
        escritor = csv.DictWriter(fh, fieldnames=list(todas[0].keys()), delimiter="\t")
        escritor.writeheader()
        escritor.writerows(todas)

    log(f"\nTable: {destino}")
    contagem = Counter(l["group"] for l in todas)
    n_loci = len({l.get("locus") or l["orf"] for l in todas})
    n_rep = sum(1 for l in todas if l.get("locus_representative"))
    log(f"Total: {len(todas)} candidate ORFs in {n_loci} distinct loci "
        f"({len(todas) - n_rep} collapsed as overlapping reading frames)")
    for group in sorted(contagem):
        conserv = Counter(l["conservation"] for l in todas if l["group"] == group and l["conservation"] != "n/a")
        extra = f"  ({', '.join(f'{k}={v}' for k, v in sorted(conserv.items()))})" if conserv else ""
        log(f"  {group:20} {contagem[group]:5}{extra}")


if __name__ == "__main__":
    main()
