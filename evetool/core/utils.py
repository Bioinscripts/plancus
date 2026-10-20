import os
import shutil

# Blacklist of RefSeq viral accessions with known false positives: bacteriophages, host
# genes captured by viruses (HAM1/GroEL/calcineurin), contaminants and spurious vertebrate viruses.
FALSOS_POSITIVOS_REFSEQ_ACCESSION = {
    "NP_040937.1", "NP_620053.1", "NP_776266.1", "YP_001876447.1", "YP_003359177.1",
    "YP_003773387.1", "YP_004063681.1", "YP_004063983.1", "YP_004598981.2", "YP_004767305.1",
    "YP_005271224.1", "YP_005352845.1", "YP_005352853.1", "YP_007027011.1", "YP_007032446.1",
    "YP_007697629.1", "YP_007697636.1", "YP_007697642.1", "YP_009026378.1",
    "YP_009111319.1", "YP_009333345.1", "YP_009333466.1", "YP_009337442.1", "YP_009342322.1",
    "YP_009345058.1", "YP_009448185.1", "YP_009505587.1", "YP_009505589.1",
    "YP_009513020.1", "YP_009551696.1", "YP_009551698.1", "YP_010085116.1", "YP_010086809.1",
    "YP_010087319.1", "YP_010356086.1", "YP_010356144.1", "YP_010356448.1", "YP_010356503.1",
    "YP_010769290.1", "YP_010800375.1", "YP_010800379.1", "YP_010800854.1", "YP_010840026.1",
    "YP_460036.1"
}
# Versions without the decimal suffix, for robust matching
FALSOS_POSITIVOS_REFSEQ_ACCESSION.update({a.split('.')[0] for a in list(FALSOS_POSITIVOS_REFSEQ_ACCESSION)})

# Disallowed HMM / conserved domains (RVDB FAM / eggNOG, curated by the researcher)
DOMINIOS_EXCLUSAO_HMM = {
    "FAM006507", "FAM010466", "FAM005753", "4R0ZU", "FAM003478", "FAM003906", "FAM007500",
    "FAM020061", "FAM021751", "FAM033562", "4R0XZ", "FAM008720", "4QP6Y", "4R0UT", "4R13J",
    "FAM001232", "FAM001338", "FAM002630", "FAM024525", "FAM024750", "FAM023132", "FAM024758",
    "FAM025425", "FAM025723", "FAM030806", "FAM020950", "FAM022720", "FAM023030", "FAM010785",
    "FAM010808", "FAM019848", "FAM009125", "FAM010068", "FAM008389", "FAM008332", "FAM003101",
    "FAM005852", "FAM005924", "FAM006212", "FAM006283", "FAM006950", "FAM007012", "FAM007251",
    "FAM007270", "FAM008135", "4QP5T", "4R0V2", "4R10Q", "FAM003092", "FAM006504", "FAM006619", "FAM011606", "FAM008567"
}
DOMINIOS_EXCLUSAO_HMM.update({d.split('.')[0].upper() for d in list(DOMINIOS_EXCLUSAO_HMM)})



# Viral products that are, at origin, cellular genes captured by horizontal transfer.
# They apply even when the hit is to an RNA virus: HAM1/ITPase in RNA viruses (cassava)
# is the case that motivates this list.
import re as _re

REGEX_PRODUTO_HGT = _re.compile(
    r"ham1|itpase|dITP|XTP.*pyrophosphatase|inosine triphosphate|"
    r"calcineurin|chaperonin|groel|cdc48|thymidylate kinase|"
    r"ribonucleotide reductase|thymidine kinase|dUTPase|"
    r"glutaredoxin|thioredoxin|ubiquitin[- ]conjugating",
    _re.IGNORECASE,
)

# Organisms/lineages that cannot be the source of an arthropod EVE
REGEX_ORGANISMO_EXCLUIDO = _re.compile(
    r"phage|prophage|archaeal virus|medusavirus|kaumoebavirus|mollivirus|"
    r"pacmanvirus|pandoravirus|mimiviridae|pithovirus",
    _re.IGNORECASE,
)

# Column indices of triagem_viral_diamond.tsv (sequence_searcher.preparar_e_buscar)
_I_VIRAL = {"qseqid": 0, "qlen": 1, "sseqid": 2, "slen": 3, "stitle": 4,
            "pident": 5, "evalue": 6, "bitscore": 7, "qcovhsp": 8}


def classe_genoma_viral(tipo_genoma: str) -> str:
    """Map the genome type from the RefSeq viral header onto the classes used in the decision."""
    g = (tipo_genoma or "").strip()
    if not g or g.lower() in ("unknown", "n/a", "-"):
        return "unknown"
    if g in ("ssRNA-RT", "dsDNA-RT") or g.endswith("-RT"):
        return "RT"
    if g.startswith("dsDNA"):
        return "dsDNA"
    if g.startswith("ssDNA"):
        return "ssDNA"
    if "RNA" in g:
        return "RNA"
    return "unknown"


def _partes_stitle(stitle: str):
    """Header layout: 'ACC |product [organism]|genome|species|genus|family|len'."""
    campos = stitle.split("|")
    product = campos[1].strip() if len(campos) > 1 else ""
    genoma = campos[2].strip() if len(campos) > 2 else ""
    family = campos[5].strip() if len(campos) > 5 else ""
    return product, genoma, family


def carregar_lista_remocao_manual(caminhos) -> dict:
    """
    Accessions or terms (one per line, '#' starts a comment) marked for removal.
    Accepts one path or several: the lists ADD to the built-in filter (accession
    blacklist plus the genome-type rule); they never replace it.
    Returns {term: source} so the recorded reason can name the file it came from.
    """
    if not caminhos:
        return {}
    if isinstance(caminhos, str):
        caminhos = [caminhos]
    termos = {}
    for caminho in caminhos:
        if not caminho or not os.path.exists(caminho):
            continue
        origem = os.path.splitext(os.path.basename(caminho))[0]
        with open(caminho, encoding="utf-8") as fh:
            for linha in fh:
                linha = linha.split("#", 1)[0].strip()
                if linha:
                    # Allows 'accession<TAB>comment', as in the audit file
                    termos[linha.split("\t")[0].strip().lower()] = origem
    return termos


def filtrar_taxonomia_refseq(tsv_entrada: str, diretorio_saida: str, is_custom_db: bool = False,
                             manual_removal_list=None, margem_ambiguo: float = 0.9) -> str:
    """
    Decide per ORF (not per line) whether the viral evidence is genuine or a host gene.

    The class of the BEST viral hit (highest bitscore) decides the fate of the whole ORF:
      - dsDNA or RT  -> removed. Large dsDNA viruses and retroviruses carry host genes they
                        captured (kinases, trypsin, P450, actin; ABL/SRC/v-kit for
                        retroviruses), so such a hit points to a host gene.
      - RNA or ssDNA -> kept, unless the product matches REGEX_PRODUTO_HGT (e.g. HAM1).
      - unknown      -> kept and flagged for manual review, unless it appears in one of the
                        removal lists given to --lista-remocao-manual (one or more; they add
                        to the built-in filter, for instance the database audit list written
                        by audit_viral_db.py).

    When the best dsDNA/RT hit nearly ties with the best RNA/ssDNA hit (ratio >=
    margem_ambiguo) the ORF is removed as well: the tie points to a cellular gene
    present in both viral lineages through horizontal transfer.

    Writes two helper files in the output directory:
      - triagem_viral_decisao_por_orf.tsv      (one line per ORF, with the reason)
      - triagem_viral_desconhecidos_revisar.tsv (candidates for manual curation)
    """
    if is_custom_db:
        return tsv_entrada

    os.makedirs(diretorio_saida, exist_ok=True)
    tsv_filtrado = os.path.join(diretorio_saida, "triagem_viral_diamond_refseq_filtrado.tsv")
    tsv_decisao = os.path.join(diretorio_saida, "viral_screening_decision_per_orf.tsv")
    tsv_revisar = os.path.join(diretorio_saida, "viral_screening_unknown_to_review.tsv")
    termos_manuais = carregar_lista_remocao_manual(manual_removal_list)

    # 1. Group the hits by ORF
    hits_por_orf = {}
    with open(tsv_entrada, encoding="utf-8", errors="replace") as f_in:
        for linha in f_in:
            partes = linha.rstrip("\n").split("\t")
            if len(partes) <= _I_VIRAL["qcovhsp"]:
                continue
            try:
                bitscore = float(partes[_I_VIRAL["bitscore"]])
            except ValueError:
                continue
            product, genoma, family = _partes_stitle(partes[_I_VIRAL["stitle"]])
            hits_por_orf.setdefault(partes[_I_VIRAL["qseqid"]], []).append({
                "linha": linha, "bitscore": bitscore, "classe": classe_genoma_viral(genoma),
                "product": product, "genoma": genoma, "family": family,
                "sseqid": partes[_I_VIRAL["sseqid"]], "pident": partes[_I_VIRAL["pident"]],
            })

    # 2. Decide per ORF
    decisoes, mantidas, revisar = [], set(), []
    for orf, hits in hits_por_orf.items():
        hits.sort(key=lambda h: -h["bitscore"])
        melhor = hits[0]
        contagem = {}
        for h in hits:
            contagem[h["classe"]] = contagem.get(h["classe"], 0) + 1

        melhor_alvo = next((h for h in hits if h["classe"] in ("RNA", "ssDNA")), None)
        melhor_mau = next((h for h in hits if h["classe"] in ("dsDNA", "RT")), None)

        acc = melhor["sseqid"]
        acc_base = acc.split(".")[0]
        texto_manual = f"{acc} {melhor['product']} {melhor['family']}".lower()

        if acc in FALSOS_POSITIVOS_REFSEQ_ACCESSION or acc_base in FALSOS_POSITIVOS_REFSEQ_ACCESSION:
            decision, reason = "removed", "accession_blacklist"
        elif termos_manuais and any(t in texto_manual for t in termos_manuais):
            origem = next(o for t, o in termos_manuais.items() if t in texto_manual)
            decision, reason = "removed", f"lista:{origem}"
        elif melhor["classe"] == "dsDNA":
            decision, reason = "removed", "HGT_best_hit_dsDNA"
        elif melhor["classe"] == "RT":
            decision, reason = "removed", "HGT_best_hit_retro"
        elif REGEX_ORGANISMO_EXCLUIDO.search(melhor["product"] + " " + melhor["family"]):
            decision, reason = "removed", "excluded_organism"
        elif melhor["classe"] == "unknown":
            decision, reason = "kept", "UNKNOWN_review"
            revisar.append((orf, melhor))
        elif REGEX_PRODUTO_HGT.search(melhor["product"]):
            decision, reason = "removed", "HGT_cellular_product"
        elif melhor_mau and melhor_alvo and melhor_mau["bitscore"] >= margem_ambiguo * melhor_alvo["bitscore"]:
            # A dsDNA/RT bitscore nearly tied with the RNA one: the segment is likely a
            # cellular gene present in both viral lineages through horizontal transfer.
            decision, reason = "removed", "HGT_dsDNA_RT_close_to_RNA"
        else:
            decision, reason = "kept", "viral_RNA_ssDNA"

        if decision == "kept":
            mantidas.add(orf)

        decisoes.append({
            "orf": orf, "decision": decision, "reason": reason,
            "best_hit_class": melhor["classe"], "best_bitscore": round(melhor["bitscore"], 1),
            "best_pident": melhor["pident"], "best_product": melhor["product"][:80],
            "best_family": melhor["family"], "best_acc": acc,
            "best_bitscore_RNA_ssDNA": round(melhor_alvo["bitscore"], 1) if melhor_alvo else 0.0,
            "best_bitscore_dsDNA_RT": round(melhor_mau["bitscore"], 1) if melhor_mau else 0.0,
            "n_hits": len(hits),
            "class_composition": ";".join(f"{k}={v}" for k, v in sorted(contagem.items())),
        })

    # 3. Write the filtered TSV (every line of the kept ORFs) and the reports
    with open(tsv_entrada, encoding="utf-8", errors="replace") as f_in, open(tsv_filtrado, "w") as f_out:
        for linha in f_in:
            qseqid = linha.split("\t", 1)[0]
            if qseqid in mantidas:
                f_out.write(linha)

    import csv as _csv
    if decisoes:
        with open(tsv_decisao, "w", newline="") as fh:
            escritor = _csv.DictWriter(fh, fieldnames=list(decisoes[0].keys()), delimiter="\t")
            escritor.writeheader()
            escritor.writerows(sorted(decisoes, key=lambda d: (d["decision"], d["reason"])))

    with open(tsv_revisar, "w") as fh:
        fh.write("orf\tacc\tproduct\tfamily\traw_genome\tbitscore\n")
        for orf, h in sorted(revisar, key=lambda x: -x[1]["bitscore"]):
            fh.write(f"{orf}\t{h['sseqid']}\t{h['product']}\t{h['family']}\t{h['genoma']}\t{h['bitscore']:.1f}\n")

    contagem_motivos = {}
    for d in decisoes:
        chave = f"{d['decision']}:{d['reason']}"
        contagem_motivos[chave] = contagem_motivos.get(chave, 0) + 1
    print(f"[Viral filter] {len(hits_por_orf)} ORFs evaluated -> {len(mantidas)} kept")
    for chave in sorted(contagem_motivos, key=lambda k: -contagem_motivos[k]):
        print(f"               {contagem_motivos[chave]:6}  {chave}")
    if revisar:
        print(f"[Viral filter] {len(revisar)} ORFs with unknown genome type listed in: {tsv_revisar}")

    return tsv_filtrado

    return tsv_filtrado


def expandir_membros_do_cluster(ids_representantes: set, tsv_clusters: str) -> tuple:
    """
    Return the original IDs of every ORF whose MMseqs2 representative had a hit.

    Without this, only the representative of a cluster of copies of the same EVE would
    move on, and the copy number per genome - central to the analysis - would be
    underestimated.
    """
    if not tsv_clusters or not os.path.exists(tsv_clusters):
        return set(ids_representantes), 0
    membros_por_rep = {}
    with open(tsv_clusters, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            partes = linha.rstrip("\n").split("\t")
            if len(partes) >= 2:
                membros_por_rep.setdefault(partes[0], []).append(partes[1])
    expandido = set(ids_representantes)
    for rep in ids_representantes:
        expandido.update(membros_por_rep.get(rep, []))
    return expandido, len(expandido) - len(ids_representantes)


def unir_e_extrair_orfs(tsv_diamond: str, tsv_hmm: str, nt_fasta: str, diretorio_saida: str,
                        pep_fasta: str = None, tsv_clusters: str = None,
                        nt_completo: str = None, pep_completo: str = None):
    """
    Merge the ORF IDs validated by DIAMOND (primary screening) and by hmmsearch (deep
    screening), drop duplicates and extract the matching sequences from the non-redundant
    nucleotide and peptide FASTAs. Returns {"nt": path, "pep": path} when pep_fasta is
    given, or just the nucleotide path in the legacy single-argument mode.
    """
    ids_unificados = set()
    orfs_excluidas_dominios = set()

    # 1. IDs from the DIAMOND table (custom outfmt 6: qseqid is the first column)
    if os.path.exists(tsv_diamond) and os.path.getsize(tsv_diamond) > 0:
        with open(tsv_diamond, 'r') as f_diamond:
            for linha in f_diamond:
                linha = linha.strip()
                if not linha or linha.startswith("#"):
                    continue
                qseqid = linha.split("\t")[0]
                ids_unificados.add(qseqid)

    # 2. IDs from the hmmsearch --domtblout table (space-separated, not tab-separated;
    # column 1 is the target name, i.e. the ORF ID, and column 4 is the domain name)
    if os.path.exists(tsv_hmm) and os.path.getsize(tsv_hmm) > 0:
        with open(tsv_hmm, 'r') as f_hmm:
            for linha in f_hmm:
                if not linha.strip() or linha.startswith("#"):
                    continue
                parts = linha.split()
                if len(parts) >= 4:
                    target_name = parts[0]
                    dom_id = parts[3].strip().split('.')[0].upper()
                    if dom_id in DOMINIOS_EXCLUSAO_HMM or parts[3].strip().upper() in DOMINIOS_EXCLUSAO_HMM:
                        orfs_excluidas_dominios.add(target_name)
                    else:
                        ids_unificados.add(target_name)

    # Drop outright every ORF carrying any of the disallowed domains
    ids_unificados.difference_update(orfs_excluidas_dominios)

    if not ids_unificados:
        raise RuntimeError("No candidate ORF was validated by DIAMOND or hmmsearch after the filters.")

    # Recover the remaining copies of each MMseqs2 cluster whose representative was
    # validated. Member sequences only exist in the original, pre-deduplication FASTAs.
    n_representantes = len(ids_unificados)
    ids_unificados, n_extras = expandir_membros_do_cluster(ids_unificados, tsv_clusters)
    if n_extras:
        print(f"[run] Cluster members recovered: {n_representantes} representatives "
              f"-> {len(ids_unificados)} ORFs (+{n_extras} copies).")
        nt_fasta = nt_completo or nt_fasta
        pep_fasta = pep_completo or pep_fasta

    # 3. Scan the non-redundant nucleotide FASTA and keep only the unified ORFs
    fasta_final_nt = os.path.join(diretorio_saida, "eve_candidatas_unificadas.fasta")
    with open(nt_fasta, 'r') as f_in, open(fasta_final_nt, 'w') as f_out:
        gravar_sequencia = False
        for linha in f_in:
            if linha.startswith(">"):
                id_seq = linha.strip().lstrip(">").split()[0]
                if id_seq in ids_unificados:
                    gravar_sequencia = True
                    f_out.write(linha)
                else:
                    gravar_sequencia = False
            elif gravar_sequencia:
                f_out.write(linha)

    # 4. When the non-redundant peptide FASTA is given, extract the amino-acid versions
    # too: tblastn needs a protein query to anchor the candidates back on the genome.
    fasta_final_pep = None
    if pep_fasta and os.path.exists(pep_fasta):
        fasta_final_pep = os.path.join(diretorio_saida, "eve_candidatas_unificadas_pep.fasta")
        with open(pep_fasta, 'r') as f_in, open(fasta_final_pep, 'w') as f_out:
            gravar_sequencia = False
            for linha in f_in:
                if linha.startswith(">"):
                    id_seq = linha.strip().lstrip(">").split()[0]
                    if id_seq in ids_unificados:
                        gravar_sequencia = True
                        f_out.write(linha)
                    else:
                        gravar_sequencia = False
                elif gravar_sequencia:
                    f_out.write(linha)

    if pep_fasta:
        return {"nt": fasta_final_nt, "pep": fasta_final_pep}
    return fasta_final_nt


def verificar_dependencias(ferramentas: list):
    """Make sure every external program is installed and on the PATH."""
    faltantes = [f for f in ferramentas if shutil.which(f) is None]
    if faltantes:
        raise RuntimeError(f"[Fatal] Tools missing from the PATH: {', '.join(faltantes)}\nInstall them before continuing.")
    print(f"[run] All {len(ferramentas)} dependencies found.")


# --------------------------------------------------------------------------
# Best-hit classification, shared by the funnel report and the dashboard
# --------------------------------------------------------------------------
ORGANISMOS_CELULARES_SWISSPROT = {
    "MOUSE", "HUMAN", "DROME", "XENLA", "XENTR", "DANRE", "RAT", "BOVIN", "CHICK",
    "ARATH", "YEAST", "CAEEL", "MACFA", "PIG", "RABIT", "ANOGA", "AEDAE",
}


def is_hit_celular(top_virus_id) -> bool:
    """True when the best hit is a Swiss-Prot protein from a cellular organism (e.g. sp|Q811F1|ZBT41_MOUSE)."""
    top = str(top_virus_id or "")
    if not top.startswith("sp|"):
        return False
    return top.split("|")[-1].rsplit("_", 1)[-1] in ORGANISMOS_CELULARES_SWISSPROT


def classificar_familia(top_virus_id, familia_diamond=None) -> str:
    """
    Family label for figures and the dashboard:
    1) the taxonomic family from the DIAMOND/RefSeq viral stitle, when present;
    2) 'Cellular hit (Swiss-Prot)' when the best hit is a cellular protein;
    3) 'No viral hit' / 'Viral (family n/a)' otherwise.
    """
    if familia_diamond:
        return familia_diamond
    top = str(top_virus_id or "")
    if top in ("", "Nenhum", "nan"):
        return "No viral hit"
    if is_hit_celular(top):
        return "Cellular hit (Swiss-Prot)"
    return "Viral (family n/a)"
