#!/usr/bin/env python3
"""
Interactive curation of the candidate EVEs, served by the local dashboard.

Reads `eves_classified.tsv`, lets ORFs be marked as false positives and, from the
selection that remains, runs:

  - FASTA export (nucleotide or peptide) of the approved EVEs only;
  - hierarchical clustering (all-vs-all tblastx + hierarchical_clustering);
  - genomic context figures for the selected loci.

The decisions live in `curation_selection.json`, at the root of the results, and never
touch the pipeline's own tables.

Usage:
  python -m evetool.curation.curation -o results_eve --listar
  python -m evetool.curation.curation -o results_eve --remover ORF_1 ORF_2
  python -m evetool.curation.curation -o results_eve --exportar-fasta nt
  python -m evetool.curation.curation -o results_eve --clusterizar
  python -m evetool.curation.curation -o results_eve --html-offline   # single file, opens with a double click
"""
import re
import os
import csv
import glob
import json
import shutil
import argparse
import subprocess
from collections import Counter, OrderedDict

NOME_SELECAO = "curation_selection.json"
TABELA = "eves_classified.tsv"
SUBPASTA = "13_curadoria"


def log(msg):
    print(f"[curadoria] {msg}", flush=True)


def pastas_especies(outdir: str) -> "OrderedDict[str, str]":
    pastas = OrderedDict()
    especies_dir = os.path.join(outdir, "especies")
    if os.path.isdir(especies_dir):
        for p in sorted(glob.glob(os.path.join(especies_dir, "*"))):
            if os.path.isdir(p):
                pastas[os.path.basename(p)] = p
    else:
        pastas[""] = outdir
    return pastas


def caminho_selecao(outdir: str) -> str:
    return os.path.join(outdir, NOME_SELECAO)


def _ler_selecao(outdir: str) -> dict:
    caminho = caminho_selecao(outdir)
    if os.path.exists(caminho):
        try:
            with open(caminho, encoding="utf-8") as fh:
                d = json.load(fh)
                return d if isinstance(d, dict) else {}
        except Exception:
            pass
    return {}


def carregar_selecao(outdir: str) -> set:
    return set(_ler_selecao(outdir).get("removed", []))


def classes_padrao() -> list:
    """Functional classes recognised automatically, in the order they are tested."""
    import sys
    from evetool.validation.classify_eves import CLASSES_PROTEINA
    return [nome for nome, _ in CLASSES_PROTEINA] + ["other"]


def salvar_selecao(outdir: str, removed, proteins: dict = None, classes_novas=None,
                   families: dict = None, familias_novas=None, tblastn_removed=None) -> str:
    """
    Store the curation state: rejected ORFs and, when the researcher overrides it by hand,
    the protein class assigned to each segment. Classes they create are recorded so they
    show up in the list next time.
    """
    atual = _ler_selecao(outdir)
    atual["removed"] = sorted(set(removed))
    if proteins is not None:
        # an empty value drops the manual annotation and returns the segment to the automatic class
        atual["proteins"] = {k: v for k, v in proteins.items() if v}
    if families is not None:
        atual["families"] = {k: v for k, v in families.items() if v}
    if tblastn_removed is not None:
        # tblastn HSPs discarded by hand ("species|contig|ref|start|end"): the segments are rebuilt
        # from the HSPs at any distance, so the removal survives a change of the merge distance
        atual["tblastn_removed"] = sorted(set(map(str, tblastn_removed)))
    registradas = set(atual.get("custom_classes", [])) | set(classes_novas or [])
    registradas |= {v for v in atual.get("proteins", {}).values() if v not in classes_padrao()}
    atual["custom_classes"] = sorted(registradas)
    atual["custom_families"] = sorted(set(atual.get("custom_families", []))
                                              | set(familias_novas or [])
                                              | set(atual.get("families", {}).values()))
    with open(caminho_selecao(outdir), "w", encoding="utf-8") as fh:
        json.dump(atual, fh, ensure_ascii=False, indent=1)
    log(f"{len(atual['removed'])} ORFs rejected, "
        f"{len(atual.get('proteins', {}))} proteins and "
        f"{len(atual.get('families', {}))} families annotated by hand")
    _exportar_anotadas(outdir)
    return caminho_selecao(outdir)


def _exportar_anotadas(outdir: str) -> str:
    """TSV of the approved EVEs with the final class (automatic or hand-corrected)."""
    dados = listar_candidatas(outdir)
    if dados.get("error"):
        return ""
    destino_dir = os.path.join(outdir, SUBPASTA)
    os.makedirs(destino_dir, exist_ok=True)
    destino = os.path.join(destino_dir, "eves_curated_annotated.tsv")
    cols = ["query", "species", "orf", "locus", "locus_representative", "contig", "start", "end",
            "group", "conservation", "family", "family_manual",
            "protein_class", "class_manual", "product", "protein",
            "viral_bitscore", "cellular_bitscore", "pident", "qcov",
            "hmm_n_domains", "hmm_best_model", "hmm_best_score",
            "hmm_consensus_family", "hmm_source"]
    with open(destino, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for i in dados["items"]:
            if not i["removed_flag"]:
                w.writerow(i)
    return destino


def _float(valor, padrao=0.0):
    try:
        return float(valor)
    except (TypeError, ValueError):
        return padrao


def listar_candidatas(outdir: str) -> dict:
    """Linhas de eves_classificadas.tsv enriquecidas com o estado da curadoria."""
    caminho = os.path.join(outdir, TABELA)
    if not os.path.exists(caminho):
        return {"error": f"{TABELA} not found in {outdir}. Run the classification first.",
                "items": [], "summary": {}}

    # HMM annotation summary (the extra layer over all candidates), per species
    hmm_por_orf = {}
    for _especie, pasta in pastas_especies(outdir).items():
        caminho_hmm = os.path.join(pasta, "04_candidatas_eves", "hmm_summary.tsv")
        if os.path.exists(caminho_hmm):
            with open(caminho_hmm, encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    hmm_por_orf[r.get("orf", "")] = r

    # Coordinates of each candidate, so the browser can measure the distances between EVEs
    coords = {}
    for _especie, pasta in pastas_especies(outdir).items():
        for contig, ini, end, strand, nome in _todos_os_loci(pasta):
            coords[nome] = {"contig": contig, "start": ini, "end": end, "strand": strand}
    # Contig lengths: needed to draw the whole chromosome
    tamanhos = {}
    for _especie, pasta in pastas_especies(outdir).items():
        info = _info_pipeline(pasta)
        cache = os.path.join(outdir, SUBPASTA, "cache_gc.json")
        if os.path.exists(cache):
            try:
                with open(cache, encoding="utf-8") as fh:
                    for c, v in (json.load(fh).get("contigs") or {}).items():
                        tamanhos.setdefault(c, v.get("tamanho", 0))
            except Exception:
                pass
        fai = (info.get("genome") or "") + ".fai"
        if os.path.exists(fai):
            with open(fai, encoding="utf-8", errors="replace") as fh:
                for linha in fh:
                    p_ = linha.split("\t")
                    if len(p_) >= 2:
                        tamanhos.setdefault(p_[0], int(p_[1]))

    # Genomic-context cluster of each candidate (bedtools merge of the anchored ORFs)
    ctx_cluster = {}
    for species, pasta in pastas_especies(outdir).items():
        regioes = _clusters_contexto(pasta)
        for contig, ini, end, _strand, orf in _todos_os_loci(pasta):
            for c in regioes.get(contig, ()):
                if c["start"] < end and c["end"] > ini:
                    ctx_cluster[(species, orf)] = c["name"]
                    break

    # NCBI CD-Search layer (cdsearch.py), per species
    cdd_por_orf = {}
    for species, pasta in pastas_especies(outdir).items():
        caminho_cdd = os.path.join(pasta, "04_candidatas_eves", "cdsearch_summary.tsv")
        if os.path.exists(caminho_cdd):
            with open(caminho_cdd, encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cdd_por_orf[(species, r.get("orf", ""))] = r
    tem_cdd = bool(cdd_por_orf)

    # Pfam layer (pfam.py), per species
    pfam_por_orf = {}
    for species, pasta in pastas_especies(outdir).items():
        caminho_pfam = os.path.join(pasta, "04_candidatas_eves", "pfam_summary.tsv")
        if os.path.exists(caminho_pfam):
            with open(caminho_pfam, encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    pfam_por_orf[(species, r.get("orf", ""))] = r
    tem_pfam = bool(pfam_por_orf)

    revisao = _orfs_em_revisao(outdir)
    metricas = metricas_sequencia(outdir)

    estado = _ler_selecao(outdir)
    removed = set(estado.get("removed", []))
    manuais = estado.get("proteins", {})
    manuais_fam = estado.get("families", {})
    itens = []
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            query = r.get("query") or r.get("orf", "")
            itens.append({
                "query": query,
                "species": r.get("species", ""),
                "orf": r.get("orf", ""),
                # Overlapping reading frames of one insertion share a locus; the page can
                # show one row per locus instead of one per ORF.
                "locus": r.get("locus", ""),
                "locus_representative": str(r.get("locus_representative", "True")) == "True",
                "group": r.get("group", ""),
                "conservation": r.get("conservation", ""),
                "family": manuais_fam.get(query) or r.get("viral_family", ""),
                "family_auto": r.get("viral_family", ""),
                "family_manual": bool(manuais_fam.get(query)),
                "protein": r.get("viral_protein", ""),
                "product": r.get("viral_product", ""),
                "protein_class": manuais.get(query) or r.get("protein_class", ""),
                "class_auto": r.get("protein_class", ""),
                "class_manual": bool(manuais.get(query)),
                "viral_bitscore": _float(r.get("viral_bitscore")),
                "cellular_bitscore": _float(r.get("cellular_bitscore")),
                "ratio": r.get("viral_cellular_ratio", ""),
                "pident": _float(r.get("viral_pident")),
                "qcov": _float(r.get("viral_qcov")),
                "hmm_model": r.get("hmm_best_model", r.get("hmm_model", "")),
                "hmm_score": _float(r.get("hmm_best_score", r.get("hmm_score"))),
                "best_cellular": r.get("best_cellular", ""),
                "cellular_title": r.get("cellular_title", ""),
                "removed_flag": query in removed,
                "ctx_cluster": ctx_cluster.get((r.get("species", ""), r.get("orf", "")), ""),
                **_campos_cdd(cdd_por_orf.get((r.get("species", ""), r.get("orf", ""))), tem_cdd),
                **_campos_pfam(pfam_por_orf.get((r.get("species", ""), r.get("orf", ""))), tem_pfam),
                "contig": coords.get(r.get("orf", ""), {}).get("contig", ""),
                "start": coords.get(r.get("orf", ""), {}).get("start", 0),
                "end": coords.get(r.get("orf", ""), {}).get("end", 0),
                "strand": coords.get(r.get("orf", ""), {}).get("strand", "."),
                "hmm_n_domains": int(r.get("hmm_n_domains")
                                      or hmm_por_orf.get(r.get("orf", ""), {}).get("n_domains") or 0),
                "hmm_consensus_family": (r.get("hmm_consensus_family")
                                         or hmm_por_orf.get(r.get("orf", ""), {}).get("hmm_consensus_family", "")),
                "hmm_best_model": (r.get("hmm_best_model")
                                      or hmm_por_orf.get(r.get("orf", ""), {}).get("best_model", "")),
                "hmm_best_score": _float(r.get("hmm_best_score")
                                           or hmm_por_orf.get(r.get("orf", ""), {}).get("best_score")),
                "hmm_source": r.get("hmm_source", ""),
                # Target with genome type "unknown": formerly decided on a separate review
                # page, now just flagged here and approved or rejected like any other row.
                "review": revisao.get((r.get("species", ""), r.get("orf", ""))),
                # Sequence metrics of the overview charts (length, GC, flanks, domains, taxonomy)
                "seq": metricas["by_query"].get(query),
            })

    # With one row per locus only the representative ORF shows; a flagged ORF hidden under it
    # must still be visible, so the flag is carried up to the representative of its locus.
    por_locus = {}
    for i in itens:
        if i["review"] and i["locus"]:
            por_locus.setdefault((i["species"], i["locus"]), dict(i["review"], via=i["orf"]))
    for i in itens:
        if not i["review"] and i["locus_representative"] and (i["species"], i["locus"]) in por_locus:
            i["review"] = por_locus[(i["species"], i["locus"])]

    aprovadas = [i for i in itens if not i["removed_flag"]]
    resumo = {
        "total": len(itens),
        "approved": len(aprovadas),
        "loci": len({i["locus"] or i["query"] for i in aprovadas}),
        "rejected": len(itens) - len(aprovadas),
        "per_group": dict(Counter(i["group"] for i in aprovadas)),
        "per_species": dict(Counter(i["species"] for i in aprovadas)),
        "per_conservation": dict(Counter(i["conservation"] for i in aprovadas if i["conservation"] != "n/a")),
        "per_evidence": dict(Counter(
            ("Alignment + HMM" if (i["viral_bitscore"] > 0 and i["hmm_n_domains"] > 0) else
             "Alignment only" if i["viral_bitscore"] > 0 else
             "HMM domain only" if i["hmm_n_domains"] > 0 else "No evidence")
            for i in aprovadas)),
        "per_family": dict(Counter(i["family"] for i in aprovadas if i["family"] not in ("", "n/a")).most_common(20)),
        "per_protein_class": dict(Counter(i["protein_class"] for i in aprovadas
                                            if i["protein_class"] not in ("", "n/a")).most_common(20)),
        "per_cdd_scope": dict(Counter(i["cdd_scope"] for i in aprovadas if i["cdd_scope"])),
    }
    return {"outdir": os.path.abspath(outdir), "items": itens, "summary": resumo,
            "classes": classes_padrao(),
            "custom_classes": estado.get("custom_classes", []),
            # Known families = those seen automatically plus the ones the researcher created
            "contigs": tamanhos,
            "families": sorted({i["family_auto"] for i in itens if i["family_auto"] not in ("", "n/a")}
                               | set(estado.get("custom_families", []))),
            "custom_families": estado.get("custom_families", []),
            "n_annotated": sum(1 for i in itens if i["class_manual"]),
            "n_families_annotated": sum(1 for i in itens if i["family_manual"]),
            # counted per locus, like the page's tile (a flag is also carried to the representative)
            "n_review": len({(i["species"], i["locus"] or i["orf"]) for i in itens if i["review"]}),
            "species_count": len({i["species"] for i in itens}),
            "cdsearch": tem_cdd,
            "pfam": tem_pfam,
            "sequence": {k: v for k, v in metricas.items() if k != "by_query"},
            # where each EVE falls on its viral protein (only the saved map: computing it
            # may read the DIAMOND database, so the page asks /api/protein-map for that)
            "protein_map": mapa_de_proteinas(outdir, itens, so_cache=True),
            # tblastn of the viral references over the genome (tblastn_scan.py), per species
            "tblastn_scan": _varredura_salva(outdir),
            "tblastn_removed": estado.get("tblastn_removed", []),
            # cross-species clusters with their tblastx pairs (multispecies_clustering.py)
            "cross_species": _comparacao_especies(outdir),
            # repeats and low complexity around the EVEs (repeats.py)
            "repeats": _repeticoes_salvas(outdir),
            "synteny_options": {"proteomes": proteomas_sintenia(),
                                "genomes": {sp: {"ok": bool(g["genome"]), "segments": g["segments"]}
                                            for sp, g in genomas_para_sintenia(outdir).items()}},
            "stamp": estado_dados(outdir)}


def _varredura_salva(outdir: str) -> list:
    try:
        from evetool.discovery import tblastn_scan
        return tblastn_scan.carregar(outdir)
    except Exception as erro:
        log(f"tblastn scan not loaded: {erro}")
        return []


def mapa_de_proteinas(outdir: str, itens: list = None, so_cache: bool = False, forcar: bool = False) -> dict:
    """Protein map (protein_map.py) of every candidate, with the families edited by hand."""
    from evetool.validation import protein_map
    if itens is None:
        itens = listar_candidatas(outdir)["items"]
    candidatos = [{"query": i["query"], "species": i["species"], "orf": i["orf"], "family": i["family"]}
                  for i in itens]
    try:
        return protein_map.mapa(outdir, candidatos=candidatos, so_cache=so_cache, forcar=forcar)
    except Exception as erro:
        return {"error": str(erro), "items": [], "references": {}, "no_alignment": []}


def _clusters_contexto(pasta: str) -> dict:
    """{contig: [{start, end, name}]} from the merged locus BED written by main.py."""
    regioes = {}
    dir_anotacao = os.path.join(pasta, "06_anotacao_genomica")
    caminho = next((p for p in (os.path.join(dir_anotacao, "eve_hits_diretos.merged.bed"),
                                os.path.join(dir_anotacao, "eve_clusters.bed"))
                    if os.path.exists(p) and os.path.getsize(p) > 0), None)
    if not caminho:
        return regioes
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if len(p) < 3 or linha.startswith(("#", "track")):
                continue
            try:
                a, b = int(p[1]), int(p[2])
            except ValueError:
                continue
            regioes.setdefault(p[0], []).append(
                {"start": a, "end": b, "name": p[3] if len(p) > 3 else f"{p[0]}:{a}-{b}"})
    return regioes


def _campos_pfam(r: dict, tem_pfam: bool) -> dict:
    """Pfam fields of one candidate; "Not searched" when the layer ran but not on this ORF."""
    if not r:
        return {"pfam_scope": "Not searched" if tem_pfam else "", "pfam_best": "",
                "pfam_best_desc": "", "pfam_domains": ""}
    return {"pfam_scope": r.get("pfam_scope", ""), "pfam_best": r.get("pfam_best", ""),
            "pfam_best_desc": r.get("pfam_best_description", ""),
            "pfam_domains": r.get("pfam_in_scope_domains", "")}


def _campos_cdd(r: dict, tem_cdd: bool) -> dict:
    """CD-Search fields of one row; empty scope when the layer has not run for it."""
    if not r:
        return {"cdd_scope": "Not searched" if tem_cdd else "", "cdd_best": "", "cdd_best_scope": "",
                "cdd_best_desc": "", "cdd_evalue": None, "cdd_bits": 0.0, "cdd_domains": ""}
    ev = r.get("cdd_best_evalue")
    return {"cdd_scope": r.get("cdd_scope", ""), "cdd_best": r.get("cdd_best", ""),
            "cdd_best_scope": r.get("cdd_best_scope", ""),
            "cdd_best_desc": r.get("cdd_best_description", ""),
            "cdd_evalue": float(ev) if ev not in (None, "") else None,
            "cdd_bits": _float(r.get("cdd_best_bitscore")),
            "cdd_domains": r.get("cdd_in_scope_domains", "")}


def _orfs_em_revisao(outdir: str) -> dict:
    """(species, orf) -> target whose genome type is "unknown" in RefSeq viral."""
    multi = os.path.isdir(os.path.join(outdir, "especies"))
    saida = {}
    for species, pasta in pastas_especies(outdir).items():
        caminho = next((c for c in (
            os.path.join(pasta, "02_alinhamento_diamond", "viral_screening_unknown_to_review.tsv"),
            os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_desconhecidos_revisar.tsv"))
            if os.path.exists(c)), None)
        if not caminho:
            continue
        with open(caminho, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                saida[(species if multi else "", r.get("orf", ""))] = {
                    "acc": r.get("acc", ""),
                    "product": r.get("product") or r.get("produto", ""),
                    "raw_genome": r.get("raw_genome") or r.get("genoma_bruto", "")}
    return saida


def metricas_sequencia(outdir: str) -> dict:
    """
    Per-EVE metrics of the overview charts, from the table dashboard.py writes
    (12_dashboard_web/dados/eves_dashboard.tsv). Empty when it has not been computed yet.
    """
    caminho = os.path.join(outdir, "12_dashboard_web", "dados", "eves_dashboard.tsv")
    ranks = ["realm", "kingdom", "phylum", "class", "order", "family", "genus"]
    janelas = [10, 25, 50]
    saida = {"available": False, "by_query": {}, "gc_genome": {}, "ranks": ranks, "windows": janelas,
             "no_hit": "No viral hit in the blastx (HMM only)"}
    if not os.path.exists(caminho):
        return saida
    num = lambda v: float(v) if v not in (None, "", "None", "nan") else None
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tax = {k: r.get(f"tax_{k}") for k in ranks if r.get(f"tax_{k}")}
            saida["by_query"][r["query"]] = {
                "len": num(r.get("length")), "gc": num(r.get("gc")), "dg": num(r.get("delta_gc_genome")),
                **{f"f{kb}": num(r.get(f"delta_gc_flank_{kb}kb")) for kb in janelas},
                "p": r.get("protein", ""), "pn": (r.get("product") or "")[:80], "d": r.get("domain", ""),
                "hm": [m for m in (r.get("hmm_models") or "").split("; ") if m],
                # no general hit at all -> null (a separate slice); a hit without lineage -> {}
                "t": tax if (tax or r.get("general_hit")) else None,
            }
            if r.get("gc_genome"):
                saida["gc_genome"][r.get("species", "")] = num(r["gc_genome"])
    saida["available"] = bool(saida["by_query"])
    return saida


def estado_dados(outdir: str) -> str:
    """
    Version stamp of every file the page reads: the page polls it and reloads when a
    run, another tab or a re-analysis changes the data on disk.
    """
    arquivos = [os.path.join(outdir, TABELA), caminho_selecao(outdir),
                os.path.join(outdir, "12_dashboard_web", "dados", "eves_dashboard.tsv"),
                os.path.join(outdir, SUBPASTA, "multiespecies", "cross_species_clusters.json"),
                os.path.join(outdir, SUBPASTA, "multiespecies_context", "cross_species_clusters.json"),
                os.path.join(outdir, SUBPASTA, "multiespecies_tblastn", "cross_species_clusters.json"),
                os.path.join(outdir, "11_sintenia_multiespecies", "cross_species_clusters.json"),
                os.path.join(outdir, SUBPASTA, "repeticoes", "repeats.json")]
    for _sp, pasta in pastas_especies(outdir).items():
        arquivos += [os.path.join(pasta, "04_candidatas_eves", "hmm_summary.tsv"),
                     os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta"),
                     os.path.join(pasta, "02_alinhamento_diamond", "viral_screening_unknown_to_review.tsv"),
                     os.path.join(pasta, "06_anotacao_genomica", "eve_hits_diretos.merged.bed"),
                     os.path.join(pasta, "04_candidatas_eves", "cdsearch_summary.tsv"),
                     os.path.join(pasta, "04_candidatas_eves", "pfam_summary.tsv")]
    try:
        from evetool.discovery import tblastn_scan
        arquivos += tblastn_scan.arquivos(outdir)
    except Exception:
        pass
    marcas = [f"{os.path.getmtime(a):.3f}" for a in arquivos if os.path.exists(a)]
    return f"{len(marcas)}-{max(marcas, default='0')}"


def gravar_snapshot(outdir: str) -> str:
    """
    Copy of the page data as a script (dados/snapshot.js), so the page opens with its data
    even straight from disk, without the server. The server rewrites it after every change.
    """
    import datetime
    dados = listar_candidatas(outdir)
    dados["saved_at"] = datetime.datetime.now().isoformat(timespec="seconds")
    destino_dir = os.path.join(outdir, "12_dashboard_web", "dados")
    os.makedirs(destino_dir, exist_ok=True)
    destino = os.path.join(destino_dir, "snapshot.js")
    corpo = json.dumps(dados, ensure_ascii=False).replace("</", "<\\/")
    with open(destino + ".tmp", "w", encoding="utf-8") as fh:
        fh.write("window.EVE_SNAPSHOT = " + corpo + ";\n")
    os.replace(destino + ".tmp", destino)
    try:
        exportar_html_autonomo(outdir, corpo)
    except Exception as erro:
        log(f"standalone page not written: {erro}")
    return destino


NOME_AUTONOMO = "curadoria_offline.html"


def exportar_html_autonomo(outdir: str, corpo_json: str = None) -> str:
    """
    One self-contained file (page + data + Plotly) in 12_dashboard_web: opens with a double
    click, without server.py and without the dados/ and assets/ folders, so it can also be
    copied or sent on its own. Read only: saving and the re-analysis still need the server.
    """
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
    with open(os.path.join(base, "curation.html"), encoding="utf-8") as fh:
        html = fh.read()
    if corpo_json is None:
        import datetime
        dados = listar_candidatas(outdir)
        dados["saved_at"] = datetime.datetime.now().isoformat(timespec="seconds")
        corpo_json = json.dumps(dados, ensure_ascii=False).replace("</", "<\\/")
    # A replacement function, so backslashes in the data are never read as regex escapes
    html = re.sub(r'<script src="dados/snapshot\.js"></script>',
                  lambda _m: "<script>window.EVE_SNAPSHOT = " + corpo_json + ";</script>", html, count=1)
    plotly = os.path.join(base, "assets", "plotly-2.24.1.min.js")
    if os.path.exists(plotly):
        with open(plotly, encoding="utf-8") as fh:
            js = fh.read().replace("</script", "<\\/script")
        html = re.sub(r'<script src="assets/plotly-[^"]+"></script>',
                      lambda _m: "<script>" + js + "</script>", html, count=1)
    destino = os.path.join(outdir, "12_dashboard_web", NOME_AUTONOMO)
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    with open(destino + ".tmp", "w", encoding="utf-8") as fh:
        fh.write(html)
    os.replace(destino + ".tmp", destino)
    return destino


def instalar_pagina(outdir: str) -> str:
    """
    Put the single page (curation + overview) in 12_dashboard_web as index.html and
    curation.html, with a fresh data snapshot, so it opens with its data from disk too.
    """
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
    origem = os.path.join(base, "curation.html")
    destino_dir = os.path.join(outdir, "12_dashboard_web")
    os.makedirs(destino_dir, exist_ok=True)
    for nome in ("index.html", "curation.html"):
        shutil.copy(origem, os.path.join(destino_dir, nome))
    velha = os.path.join(destino_dir, "revisao.html")   # the review now happens on the same page
    if os.path.exists(velha):
        os.remove(velha)
    plotly = os.path.join(base, "assets", "plotly-2.24.1.min.js")
    if os.path.exists(plotly):
        os.makedirs(os.path.join(destino_dir, "assets"), exist_ok=True)
        alvo = os.path.join(destino_dir, "assets", os.path.basename(plotly))
        if not os.path.exists(alvo):
            shutil.copy(plotly, alvo)
    gravar_snapshot(outdir)
    return os.path.join(destino_dir, "index.html")


def rodar_metricas(outdir: str, threads: int = 4) -> dict:
    """(Re)compute the sequence metrics of the overview charts with dashboard.py."""
    import sys
    from evetool.visualization.dashboard import gerar_dashboard, HMM_PADRAO
    genomas, nomes, hmm = {}, {}, None
    multi = os.path.isdir(os.path.join(outdir, "especies"))
    for species, pasta in pastas_especies(outdir).items():
        info = _info_pipeline(pasta) or _info_pipeline(outdir)
        if info.get("genome"):
            genomas[species if multi else ""] = info["genome"]
        hmm = hmm or info.get("hmm_db")
    gerar_dashboard(outdir, genomas=genomas, nomes=nomes, hmm_db=hmm or HMM_PADRAO, threads=threads)
    m = metricas_sequencia(outdir)
    return {"n_metrics": len(m["by_query"])}


def _fastas_por_especie(outdir: str, tipo: str) -> "OrderedDict[str, str]":
    nome = "eve_candidatas_unificadas_pep.fasta" if tipo == "pep" else "eve_candidatas_unificadas.fasta"
    mapa = OrderedDict()
    for species, pasta in pastas_especies(outdir).items():
        caminho = os.path.join(pasta, "04_candidatas_eves", nome)
        if os.path.exists(caminho):
            mapa[species] = caminho
    return mapa


def exportar_fasta(outdir: str, tipo: str = "nt", queries=None) -> dict:
    """
    Write a FASTA with the approved EVEs only (or with the given 'queries').
    Headers get the species prefix in multi-species mode.
    """
    if tipo not in ("nt", "pep"):
        raise ValueError("tipo must be 'nt' or 'pep'")

    dados = listar_candidatas(outdir)
    if dados.get("error"):
        raise RuntimeError(dados["error"])

    if queries:
        desejadas = set(queries)
    else:
        desejadas = {i["query"] for i in dados["items"] if not i["removed_flag"]}
    if not desejadas:
        raise RuntimeError("No EVE selected for export.")

    # query = "Species|ORF" in multi mode; map ORF -> query per species
    por_especie = {}
    for q in desejadas:
        species, _, orf = q.partition("|") if "|" in q else ("", "", q)
        por_especie.setdefault(species, {})[orf] = q

    destino_dir = os.path.join(outdir, SUBPASTA)
    os.makedirs(destino_dir, exist_ok=True)
    destino = os.path.join(destino_dir, f"eves_curadas_{tipo}.fasta")

    escritas, faltando = 0, []
    with open(destino, "w", encoding="utf-8") as saida:
        for species, caminho in _fastas_por_especie(outdir, tipo).items():
            alvo = por_especie.get(species) or por_especie.get("", {})
            if not alvo:
                continue
            gravar, encontrados = False, set()
            with open(caminho, encoding="utf-8", errors="replace") as fh:
                for linha in fh:
                    if linha.startswith(">"):
                        orf = linha[1:].split()[0]
                        gravar = orf in alvo
                        if gravar:
                            encontrados.add(orf)
                            escritas += 1
                            saida.write(f">{alvo[orf]}{linha[1 + len(orf):]}")
                            continue
                    if gravar:
                        saida.write(linha)
            faltando += [alvo[o] for o in set(alvo) - encontrados]

    log(f"{escritas} sequences ({tipo}) written to {destino}")
    return {"file": destino, "sequences": escritas, "not_found": faltando[:20],
            "total_not_found": len(faltando)}


def rodar_clusterizacao(outdir: str, threads: int = 4, queries=None) -> dict:
    """All-vs-all tblastx of the approved EVEs + hierarchical clustering (hierarchical_clustering.py)."""
    for exe in ("makeblastdb", "tblastx"):
        if shutil.which(exe) is None:
            raise RuntimeError(f"'{exe}' not found on the PATH.")

    fasta = exportar_fasta(outdir, "nt", queries=queries)["file"]
    n = sum(1 for l in open(fasta) if l.startswith(">"))
    if n < 3:
        raise RuntimeError(f"At least 3 sequences are needed to cluster (found {n}).")

    destino = os.path.join(outdir, SUBPASTA, "clusterizacao")
    os.makedirs(destino, exist_ok=True)
    db = os.path.join(destino, "eves_curadas_db")
    xml = os.path.join(destino, "tblastx_eves_curadas.xml")

    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    subprocess.run(["tblastx", "-query", fasta, "-db", db, "-out", xml,
                    "-evalue", "1e-3", "-num_threads", str(threads), "-outfmt", "5"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

    prefixo = os.path.join(destino, "eve_clusters_hierarch")
    import sys
    from evetool.comparative.hierarchical_clustering import run_clustering
    run_clustering(fmt5=xml, out=prefixo, no_r=True)

    imagens = sorted(glob.glob(prefixo + "*.png"))
    tabelas = sorted(glob.glob(prefixo + "*.tsv")) + sorted(glob.glob(prefixo + "*.csv"))
    clusters = _ler_clusters(tabelas)
    log(f"clustering finished: {len(imagens)} figure(s), {len(clusters)} sequences in clusters")
    return {"sequences": n, "xml": xml, "images": imagens, "tables": tabelas, "clusters": clusters}


def _ler_clusters(tabelas) -> dict:
    """Extract {sequence: cluster} from *_segmentos.tsv (or another table with those columns)."""
    # The distance matrix also carries 'cluster' in its name; the segments file comes first.
    ordenadas = sorted(tabelas, key=lambda t: (0 if t.endswith("_segmentos.tsv") else 1, t))
    for t in ordenadas:
        try:
            with open(t, encoding="utf-8", errors="replace") as fh:
                amostra = fh.readline()
                sep = "\t" if "\t" in amostra else ","
                fh.seek(0)
                leitor = csv.DictReader(fh, delimiter=sep)
                campos = [c.lower() for c in (leitor.fieldnames or [])]
                col_seq = next((c for c, lc in zip(leitor.fieldnames, campos)
                                if lc in ("segmento", "name", "seq", "sequencia", "label", "id", "query")), None)
                col_cl = next((c for c, lc in zip(leitor.fieldnames, campos)
                               if lc == "cluster" or lc in ("group", "group")), None)
                if col_seq and col_cl:
                    return {r[col_seq]: r[col_cl] for r in leitor if r.get(col_seq)}
        except Exception:
            continue
    return {}


def rodar_anotacao_hmm(outdir: str, threads: int = 4, hmm_db: str = None) -> dict:
    """Apply the calibrated HMM database to every candidate, as an extra layer."""
    import sys
    from evetool.annotation.hmm_domains import anotar, HMM_PADRAO
    info = _info_pipeline(outdir)
    banco = hmm_db or info.get("hmm_db") or HMM_PADRAO
    total = {"n_candidates": 0, "n_with_domain": 0, "n_domains": 0, "per_species": {}}
    for species, pasta in pastas_especies(outdir).items():
        if not os.path.exists(os.path.join(pasta, "04_candidatas_eves",
                                           "eve_candidatas_unificadas_pep.fasta")):
            continue
        r = anotar(outdir, hmm_db=banco, threads=threads, pasta_especie=pasta, forcar=True)
        total["per_species"][species or "single species"] = f"{r['n_com_dominio']}/{r['n_candidatas']}"
        # anotar() still returns the Portuguese keys; the page reads the English ones
        total["n_candidates"] += r["n_candidatas"]
        total["n_with_domain"] += r["n_com_dominio"]
        total["n_domains"] += r["n_domains"]
    return total


def rodar_cdsearch(outdir: str, evalue: float = 0.01, threads: int = 4, modo: str = "local") -> dict:
    """CD-Search over every candidate (local rpsblast+ by default); searched sequences are cached."""
    import sys
    from evetool.annotation.cdsearch import anotar
    return anotar(outdir, evalue=evalue, threads=threads, modo=modo)


def rodar_pfam(outdir: str, threads: int = 4) -> dict:
    """Pfam-A over every candidate (hmmscan --cut_ga); scanned sequences are cached."""
    import sys
    from evetool.annotation.pfam import anotar
    return anotar(outdir, threads=threads)


def rodar_varredura_tblastn(outdir: str, threads: int = 4, hits_por_orf: int = 1, evalue: float = 1e-5,
                            gap: int = 1000, forcar: bool = False) -> dict:
    """tblastn of the best viral hits over the genome, with CD-Search and Pfam on the new ORFs."""
    import sys
    from evetool.discovery.tblastn_scan import varrer
    r = varrer(outdir, hits_por_orf=hits_por_orf, evalue=evalue, gap=gap, threads=threads, forcar=forcar)
    r.pop("species", None)
    return r


def _genoma_da_especie(outdir: str, species: str = "") -> str:
    pastas = pastas_especies(outdir)
    pasta = pastas.get(species) or pastas.get("") or next(iter(pastas.values()), outdir)
    genoma = _info_pipeline(pasta).get("genome")
    if not genoma or not os.path.exists(genoma):
        raise RuntimeError(f"Genome of '{species or 'this run'}' not found (pipeline_info.json).")
    return genoma


def sequencia_regiao(outdir: str, species: str, contig: str, start: int, end: int, nome: str = "") -> dict:
    """Genome sequence of a region (1-based, inclusive), as used by 'copy' on the page."""
    import sys
    from evetool.discovery.ncbi_blastx import extrair_regiao
    start, end = int(start), int(end)
    if end < start or end - start > 5_000_000:
        raise RuntimeError("Invalid region (end before start, or longer than 5 Mb).")
    seq = extrair_regiao(_genoma_da_especie(outdir, species), contig, start, end)
    rotulo = f"{nome + ' ' if nome else ''}{species + '|' if species else ''}{contig}:{start}-{end}"
    return {"header": rotulo, "sequence": seq, "length": len(seq)}


def sequencias_regioes(outdir: str, regioes: list, flanco: int = 0) -> dict:
    """
    Multi-FASTA of several genome regions ([{sp, c, s, e, n, d}], 1-based inclusive; n names the
    record, d is free text for its header), each widened by `flanco` nt on both sides. One
    samtools faidx call per species. The + strand of the genome, as the single 'copy'.
    """
    import subprocess
    flanco = max(0, int(flanco or 0))
    por_sp, saida, erros, total = OrderedDict(), [], [], 0
    for k, r in enumerate(regioes):
        try:
            a, b = int(r["s"]), int(r["e"])
        except (KeyError, TypeError, ValueError):
            erros.append(f"region {k + 1}: invalid coordinates")
            continue
        if b < a:
            erros.append(f"{r.get('n') or k + 1}: end before start")
            continue
        a, b = max(1, a - flanco), b + flanco
        total += b - a + 1
        por_sp.setdefault(str(r.get("sp") or ""), []).append((k, str(r["c"]), a, b, r))
    if total > 50_000_000:
        raise RuntimeError(f"{total:,} nt requested; the limit is 50 Mb (fewer regions or a smaller flank).")
    registros = {}
    for sp, lista in por_sp.items():
        genoma = _genoma_da_especie(outdir, sp)
        p = subprocess.run(["samtools", "faidx", genoma, *[f"{c}:{a}-{b}" for _k, c, a, b, _r in lista]],
                           capture_output=True, text=True, timeout=600)
        if p.returncode != 0 and not p.stdout:
            raise RuntimeError(f"samtools faidx: {p.stderr.strip()[:300]}")
        seqs, atual = [], None
        for linha in p.stdout.splitlines():
            if linha.startswith(">"):
                atual = []
                seqs.append(atual)
            elif atual is not None:
                atual.append(linha.strip())
        if len(seqs) != len(lista):
            raise RuntimeError(f"samtools faidx returned {len(seqs)} of {len(lista)} regions: {p.stderr.strip()[:200]}")
        for (k, c, a, b, r), partes in zip(lista, seqs):
            seq = "".join(partes).upper()
            local = f"{sp + '|' if sp else ''}{c}:{a}-{b}"
            cab = " ".join(x for x in (str(r.get("n") or local), local if r.get("n") else "",
                                       f"len={len(seq)}", f"flank={flanco}" if flanco else "",
                                       str(r.get("d") or "")) if x)
            registros[k] = ">" + cab + "\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
    fasta = "".join(registros[k] for k in sorted(registros))
    return {"fasta": fasta, "n": len(registros), "nt": sum(len(x) for x in fasta.splitlines() if not x.startswith(">")),
            "errors": erros}


def _regioes_do_texto(texto: str, species_padrao: str = "") -> list:
    """'contig:start-end' (or 'Species|contig:start-end'), one per line or separated by spaces or ';'."""
    regioes = []
    for item in re.split(r"[\s;]+", texto or ""):
        if not item:
            continue
        m = re.match(r"^(?:([^|]+)\|)?(.+):([\d.]+)-([\d.]+)$", item.replace(",", ""))
        if not m:
            raise RuntimeError(f"Region not understood: '{item}' (use contig:start-end)")
        sp, contig, a, b = m.groups()
        regioes.append({"species": sp or species_padrao, "contig": contig,
                        "start": int(a.replace(".", "")), "end": int(b.replace(".", ""))})
    return regioes


def segmentos_para_blastx(outdir: str, queries=None, clusters=None, regioes: str = "") -> list:
    """
    Segments for the online BLASTX, always taken from the genome (+ strand, 1-based), so the
    hits map back to genome coordinates: the selected EVEs, the clusters sent by the page
    ([{name, species, contig, start, end}], start 0-based as in the page) and free regions.
    """
    import sys
    from evetool.discovery.ncbi_blastx import extrair_regiao
    segmentos, genomas = [], {}

    def genoma(sp):
        if sp not in genomas:
            genomas[sp] = _genoma_da_especie(outdir, sp)
        return genomas[sp]

    if queries:
        desejadas = set(queries)
        for species, pasta in pastas_especies(outdir).items():
            for contig, ini, fim, strand, orf in _todos_os_loci(pasta):
                q = f"{species}|{orf}" if species else orf
                if q in desejadas or orf in desejadas:
                    segmentos.append({"id": q, "origin": "EVE", "species": species, "contig": contig,
                                      "start": ini + 1, "end": fim, "strand": strand})
    for c in clusters or []:
        sp = c.get("species", "")
        segmentos.append({"id": f"{c.get('name') or 'cluster'}|{c['contig']}:{int(c['start']) + 1}-{int(c['end'])}",
                          "origin": "cluster", "species": sp, "contig": c["contig"],
                          "start": int(c["start"]) + 1, "end": int(c["end"])})
    padrao = next(iter(pastas_especies(outdir)), "")
    for r in _regioes_do_texto(regioes, padrao):
        segmentos.append({"id": f"{r['species'] + '|' if r['species'] else ''}{r['contig']}:{r['start']}-{r['end']}",
                          "origin": "region", **r})
    vistos, unicos = set(), []
    for sg in segmentos:
        if sg["id"] in vistos:
            continue
        vistos.add(sg["id"])
        sg["seq"] = extrair_regiao(genoma(sg["species"]), sg["contig"], sg["start"], sg["end"])
        unicos.append(sg)
    if not unicos:
        raise RuntimeError("Nothing to search: select EVEs, clusters or type regions.")
    return unicos


MODOS_UNIDADE = {"orfs": "multiespecies", "context": "multiespecies_context", "tblastn": "multiespecies_tblastn"}


def genomas_para_sintenia(outdir: str) -> dict:
    """
    {species: {"genome": path, "segments": bool}} used for the flanks. <results>/synteny_genomes.tsv
    ('species<TAB>genome') wins: a run made on pre-extracted segments points there to the whole
    genomes, and its segment names are mapped back onto them.
    """
    extra = {}
    cam = os.path.join(outdir, "synteny_genomes.tsv")
    if os.path.exists(cam):
        with open(cam, encoding="utf-8") as fh:
            for linha in fh:
                p = linha.split("#", 1)[0].rstrip("\n").split("\t")
                if len(p) >= 2 and p[0].strip() and p[1].strip():
                    extra[p[0].strip()] = p[1].strip()
    saida = {}
    for sp, pasta in pastas_especies(outdir).items():
        proprio = _info_pipeline(pasta).get("genome")
        g = extra.get(sp) or proprio
        saida[sp] = {"genome": g if g and os.path.exists(g) else None, "segments": sp in extra,
                     "own": proprio}
    return saida


def proteomas_sintenia() -> list:
    """Reference proteomes offered for the synteny labels: [{path, label}]."""
    import sys
    from evetool import config
    lista = []
    for cam in sorted(glob.glob(os.path.join(config.HMM_DIR, "proteomas_calibracao", "*.faa"))):
        rot = os.path.basename(cam)
        try:
            with open(cam, encoding="utf-8", errors="replace") as fh:
                m = re.search(r"\[([^\]]+)\]", fh.readline())
            if m:
                rot = f"{m.group(1)} ({os.path.basename(cam).split('_protein')[0]})"
        except Exception:
            pass
        lista.append({"path": cam, "label": rot})
    return lista


def rodar_multiespecies(outdir: str, threads: int = 4, queries=None, min_identity: float = 30.0,
                        min_coverage: float = 30.0, sintenia: bool = False,
                        proteoma_ref: str = None, renomear: str = None,
                        share_identity: float = 50.0, share_coverage: float = 50.0,
                        modo: str = "orfs", unidades: list = None, gap: int = None,
                        flanco: int = 20000, max_sitios: int = 3) -> dict:
    """
    Cross-species clustering (identity/coverage -> valleys) with a sharing heatmap.
    The units compared are chosen by the researcher: the approved ORFs ('orfs'), the
    genomic-context clusters ('context') or the tblastn segments ('tblastn'), the last two merged
    on the page at the distance they set and sent as regions [{sp, c, s, e, n, orfs}].
    Flank synteny of the shared clusters is optional.
    """
    import sys
    from argparse import Namespace
    from evetool.comparative.multispecies_clustering import executar
    from evetool.comparative.flank_synteny import origens_de_segmentos, calibrar_origens

    if modo not in MODOS_UNIDADE:
        raise RuntimeError(f"unknown unit type: {modo}")
    destino = os.path.join(outdir, SUBPASTA, MODOS_UNIDADE[modo])
    os.makedirs(destino, exist_ok=True)
    locais, membros = {}, {}
    if modo == "orfs":
        fasta = exportar_fasta(outdir, "nt", queries=queries)["file"]
        nomes = {l[1:].split()[0] for l in open(fasta) if l.startswith(">")}
        for sp, pasta in pastas_especies(outdir).items():
            for contig, ini, fim, _st, orf in _todos_os_loci(pasta):
                q = f"{sp}|{orf}" if sp else orf
                if q in nomes:
                    locais[q] = (sp, contig, ini, fim)
    else:
        unidades = [u for u in (unidades or []) if u.get("sp") is not None]
        if not unidades:
            raise RuntimeError("No units to compare: the page sends the clusters or segments it shows.")
        regs = []
        for k, u in enumerate(unidades):
            nome = re.sub(r"[^A-Za-z0-9_.-]", "_", str(u.get("n") or f"U{k + 1}"))
            q = f"{u['sp']}|{nome}" if u["sp"] else nome
            regs.append({"sp": u["sp"], "c": u["c"], "s": int(u["s"]), "e": int(u["e"]), "n": q})
            locais[q] = (u["sp"], u["c"], int(u["s"]) - 1, int(u["e"]))
            membros[q] = list(u.get("orfs") or [])
        res_fa = sequencias_regioes(outdir, regs)
        fasta = os.path.join(destino, f"units_{modo}.fasta")
        with open(fasta, "w", encoding="utf-8") as fh:
            fh.write(res_fa["fasta"])
    n_sp = len({l[1:].split("|")[0] for l in open(fasta) if l.startswith(">") and "|" in l})
    if n_sp < 2:
        raise RuntimeError(f"The cross-species comparison needs at least 2 species (found {n_sp}).")

    # A rename table belongs to one run (it fixes that run's labels): only the one saved in the
    # results folder is applied. One left next to the code would relabel every other run.
    if not renomear:
        padrao = os.path.join(outdir, "renomear_especies.tsv")
        renomear = padrao if os.path.exists(padrao) else None

    genomas, origens = {}, {}
    if sintenia:
        for sp, g in genomas_para_sintenia(outdir).items():
            genomas[sp] = g["genome"]
            if g["segments"] and g["own"]:
                origens[sp] = calibrar_origens(origens_de_segmentos(g["own"]), g["own"], g["genome"])
    cfg = Namespace(fasta=fasta, xml=None, outdir=destino, threads=threads,
                    min_identity=min_identity, min_coverage=min_coverage, filter_mode="max",
                    valley_ratio=0.5, max_bandwidth=0.05, min_pairs=15, max_depth=6,
                    w_cov=1.0, w_gap=1.0, max_overlap_nt=0, max_clusters_heatmap=40,
                    renomear=renomear, sintenia=sintenia, proteoma_ref=proteoma_ref,
                    flanco_sintenia=int(flanco or 20000), share_identity=share_identity,
                    share_coverage=share_coverage, locais=locais, genomas_sintenia=genomas,
                    origens=origens, unidades_modo=modo, unidades_gap=gap, unidades_membros=membros,
                    max_sitios_sintenia=int(max_sitios or 3))
    res = executar(cfg)
    resumo = res["resumo"]                 # multispecies_clustering.executar names its keys in Portuguese
    compartilhados = int((resumo["Total_Species"] >= 2).sum())
    imagens = [res["heatmap"]] if res.get("heatmap") else []
    sint = {}
    if sintenia and "synteny_status" in resumo:
        sint = resumo.loc[resumo["Total_Species"] >= 2, "synteny_status"].value_counts().to_dict()
    log(f"multi-species ({modo}): {len(resumo)} clusters, {compartilhados} shared")
    return {"clusters": len(resumo), "shared": compartilhados, "supported": res.get("n_suportados", 0),
            "species_count": res["n_especies"], "images": imagens, "directory": destino, "units": modo,
            "synteny": {str(k): int(v) for k, v in sint.items()},
            "tables": sorted(glob.glob(os.path.join(destino, "*.tsv"))
                              + glob.glob(os.path.join(destino, "*.csv")))}


def rodar_contexto(outdir: str, queries=None, flanco: int = 20000) -> dict:
    """Genomic context figures for the species of the approved EVEs."""
    import sys
    from evetool.visualization.plot_genomic_context import gerar_figuras_contexto_lote

    dados = listar_candidatas(outdir)
    if dados.get("error"):
        raise RuntimeError(dados["error"])
    alvo = set(queries) if queries else {i["query"] for i in dados["items"] if not i["removed_flag"]}
    especies_alvo = {q.split("|")[0] if "|" in q else "" for q in alvo}

    destino_base = os.path.join(outdir, SUBPASTA, "contexto")
    os.makedirs(destino_base, exist_ok=True)
    geradas, avisos = [], []

    for species, pasta in pastas_especies(outdir).items():
        if species not in especies_alvo:
            continue
        dir_anotacao = os.path.join(pasta, "06_anotacao_genomica")
        bed = next((p for p in (os.path.join(dir_anotacao, "eve_hits_diretos.merged.bed"),
                                os.path.join(dir_anotacao, "eve_clusters.bed"))
                    if os.path.exists(p) and os.path.getsize(p) > 0), None)
        bed_flanco = os.path.join(dir_anotacao, "eve_clusters_flanked.bed")
        if not bed:
            avisos.append(f"{species or 'single species'}: locus BED not found - run the full main.py.")
            continue
        if not (os.path.exists(bed_flanco) and os.path.getsize(bed_flanco) > 0):
            bed_flanco = bed
        depth = next(iter(glob.glob(os.path.join(pasta, "10_rnaseq_expressao", "*depth*.tsv"))), None)
        info = _info_pipeline(pasta) or _info_pipeline(outdir)
        gff = info.get("gff") if info.get("gff") and os.path.exists(info["gff"]) else \
            next(iter(glob.glob(os.path.join(pasta, "*.gff*")) + glob.glob(os.path.join(outdir, "*.gff*"))), None)
        destino = os.path.join(destino_base, species or "especie_unica")
        os.makedirs(destino, exist_ok=True)
        # Only the clusters holding the selected EVEs, not every cluster of the species
        orfs_sel = {q.split("|")[-1] for q in alvo if (q.split("|")[0] if "|" in q else "") == species}
        sel = [(c, a_, b_) for c, a_, b_, _st, o in _todos_os_loci(pasta) if o in orfs_sel]
        if sel:
            nomes_sel, linhas_sel = set(), []
            with open(bed, encoding="utf-8") as fh:
                for linha in fh:
                    c = linha.rstrip("\n").split("\t")
                    if len(c) >= 4 and any(c[0] == ct and int(c[1]) < b_ and int(c[2]) > a_ for ct, a_, b_ in sel):
                        nomes_sel.add(c[3]); linhas_sel.append(linha)
            bed = os.path.join(destino, "selected_clusters.bed")
            with open(bed, "w", encoding="utf-8") as fh:
                fh.writelines(linhas_sel)
            if bed_flanco.endswith("eve_clusters_flanked.bed"):
                bed_fl_sel = os.path.join(destino, "selected_clusters_flanked.bed")
                with open(bed_flanco, encoding="utf-8") as fh, open(bed_fl_sel, "w", encoding="utf-8") as out:
                    out.writelines(l for l in fh if len(l.split("\t")) > 3 and l.split("\t")[3] in nomes_sel)
                bed_flanco = bed_fl_sel
            for png in glob.glob(os.path.join(destino, "*.png")) + glob.glob(os.path.join(destino, "*.pdf")):
                os.remove(png)                 # figures of an earlier, different selection
        try:
            gerar_figuras_contexto_lote(bed_original=bed, bed_flanked=bed_flanco, gff_annotation=gff,
                                        depth_tsv=depth, output_dir=destino,
                                        species=species or None, flank=flanco)
            geradas += sorted(glob.glob(os.path.join(destino, "*.png")))
        except Exception as erro:
            avisos.append(f"{species or 'single species'}: {erro}")

    log(f"genomic context: {len(geradas)} figure(s), {len(avisos)} warning(s)")
    return {"images": geradas, "warnings": avisos}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--listar", action="store_true")
    ap.add_argument("--remover", nargs="*", default=None, help="Queries to mark as false positives")
    ap.add_argument("--exportar-fasta", choices=["nt", "pep"])
    ap.add_argument("--clusterizar", action="store_true")
    ap.add_argument("--multiespecies", action="store_true",
                    help="Cross-species clustering + heatmap over the approved EVEs")
    ap.add_argument("--sintenia", action="store_true", help="Include synteny (off by default)")
    ap.add_argument("--proteoma-ref", help="Reference proteome for the synteny step")
    ap.add_argument("--contexto", action="store_true")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--html-offline", action="store_true",
                    help=f"Refresh the page in 12_dashboard_web and write {NOME_AUTONOMO}, "
                         "a single file that opens with a double click (no server)")
    cfg = ap.parse_args()

    if cfg.remover is not None:
        salvar_selecao(cfg.outdir, cfg.remover)

    if cfg.listar or not any([cfg.remover is not None, cfg.exportar_fasta, cfg.clusterizar,
                              cfg.contexto, cfg.multiespecies, cfg.html_offline]):
        d = listar_candidatas(cfg.outdir)
        if d.get("error"):
            raise SystemExit(d["error"])
        r = d["summary"]
        print(f"\n{r['total']} candidates | approved: {r['approved']} | flagged as false positives: {r['rejected']}")
        for chave in ("per_group", "per_conservation", "per_species"):
            print(f"\n{chave}:")
            for k, v in sorted(r[chave].items(), key=lambda x: -x[1]):
                print(f"  {v:5}  {k}")

    if cfg.html_offline:
        instalar_pagina(cfg.outdir)
        print(f"Standalone page: {os.path.join(cfg.outdir, '12_dashboard_web', NOME_AUTONOMO)}")
    if cfg.exportar_fasta:
        res = exportar_fasta(cfg.outdir, cfg.exportar_fasta)
        print(f"FASTA: {res['arquivo']} ({res['sequencias']} sequences)")
    if cfg.clusterizar:
        res = rodar_clusterizacao(cfg.outdir, threads=cfg.threads)
        print(f"Clustering: {len(res['imagens'])} figure(s) in {os.path.dirname(res['xml'])}")
    if cfg.multiespecies:
        res = rodar_multiespecies(cfg.outdir, threads=cfg.threads, sintenia=cfg.sintenia,
                                  proteoma_ref=cfg.proteoma_ref)
        print(f"Multi-species: {res['clusters']} clusters, {res['compartilhados']} shared "
              f"across {res['especies']} species -> {res['diretorio']}")
    if cfg.contexto:
        res = rodar_contexto(cfg.outdir)
        print(f"Context: {len(res['imagens'])} figure(s)")
        for a in res["warnings"]:
            print(f"  [aviso] {a}")




# ---------------------------------------------------------------------------
# Inspecting one EVE: position, neighbourhood, GC composition and alignment
# ---------------------------------------------------------------------------
def _info_pipeline(pasta: str) -> dict:
    """Reads pipeline_info.json, accepting the older Portuguese keys as a fallback."""
    caminho = os.path.join(pasta, "pipeline_info.json")
    if os.path.exists(caminho):
        try:
            with open(caminho, encoding="utf-8") as fh:
                info = json.load(fh)
            for antigo, novo in (("genoma", "genome"), ("especie", "species"), ("flanco", "flank")):
                if antigo in info and novo not in info:
                    info[novo] = info[antigo]
            return info
        except Exception:
            pass
    return {}


def _gc(seq: str) -> float:
    s = seq.upper()
    validos = sum(s.count(b) for b in "ACGT")
    return round(100.0 * (s.count("G") + s.count("C")) / validos, 2) if validos else 0.0


def _perfil_do_genoma(genoma: str, cache_path: str) -> dict:
    """
    A single pass over the genome yields the GC and the length of every contig, plus the
    overall GC. Cached on disk, because reading a whole genome is expensive.
    """
    if os.path.exists(cache_path):
        try:
            cache = json.load(open(cache_path, encoding="utf-8"))
            if cache.get("genome") == genoma:
                return cache
        except Exception:
            pass
    perfil = {"genome": genoma, "contigs": {}, "genome_gc": 0.0, "tamanho_genoma": 0}
    if not genoma or not os.path.exists(genoma):
        return perfil
    atual, gc, at = None, 0, 0
    tot_gc = tot_at = 0

    def fechar():
        nonlocal tot_gc, tot_at
        if atual and (gc + at):
            perfil["contigs"][atual] = {"gc": round(100.0 * gc / (gc + at), 2), "tamanho": gc + at}
            tot_gc += gc
            tot_at += at

    with open(genoma, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                fechar()
                atual, gc, at = linha[1:].split()[0], 0, 0
            else:
                u = linha.strip().upper()
                gc += u.count("G") + u.count("C")
                at += u.count("A") + u.count("T")
    fechar()
    if tot_gc + tot_at:
        perfil["genome_gc"] = round(100.0 * tot_gc / (tot_gc + tot_at), 2)
        perfil["tamanho_genoma"] = tot_gc + tot_at
    try:
        json.dump(perfil, open(cache_path, "w", encoding="utf-8"))
    except Exception:
        pass
    return perfil


def _gc_das_janelas(genoma: str, contig: str, start: int, end: int, tam_contig: int,
                    janelas=(10000, 25000, 50000)) -> list:
    """
    Flank GC (excluding the EVE itself) for each requested window. A window only counts
    as valid when the contig offers at least half of the requested flank: on a short contig
    there is not enough surrounding sequence for the comparison to mean anything.
    """
    saida = []
    if not (genoma and os.path.exists(genoma) and shutil.which("samtools")):
        return saida
    for j in janelas:
        ini = max(0, start - j)
        termino = min(tam_contig, end + j) if tam_contig else end + j
        disponivel = (start - ini) + (termino - end)
        valida = disponivel >= j   # metade dos 2*j pedidos
        gc = 0.0
        if disponivel > 0:
            try:
                p = subprocess.run(["samtools", "faidx", genoma, f"{contig}:{ini + 1}-{termino}"],
                                   capture_output=True, text=True, timeout=180)
                regiao = "".join(l.strip() for l in p.stdout.splitlines() if not l.startswith(">"))
                if regiao:
                    rel = start - ini
                    gc = _gc(regiao[:rel] + regiao[rel + (end - start):])
            except Exception:
                pass
        saida.append({"window": j, "gc": gc, "nt_available": disponivel, "valid": bool(valida and gc)})
    return saida


def _perfil_gc_regiao(genoma: str, contig: str, start: int, end: int, tam_contig: int,
                      extensao: int = 50000, janela: int = 1000, passo: int = 500) -> dict:
    """
    GC in sliding windows along the region, with the EVE at the centre: it shows whether
    its composition stands out from what comes before and after, and how far away the
    surroundings return to the genome baseline.

    Output coordinates are relative to the CENTRE of the EVE (negative upstream, positive
    downstream), so the chart is symmetric around it.
    """
    vazio = {"points": [], "window": janela, "step": passo, "extent": extensao}
    if not (genoma and os.path.exists(genoma) and shutil.which("samtools")):
        return vazio
    ini = max(0, start - extensao)
    termino = min(tam_contig, end + extensao) if tam_contig else end + extensao
    try:
        p = subprocess.run(["samtools", "faidx", genoma, f"{contig}:{ini + 1}-{termino}"],
                           capture_output=True, text=True, timeout=240)
        seq = "".join(l.strip() for l in p.stdout.splitlines() if not l.startswith(">")).upper()
    except Exception:
        return vazio
    if not seq:
        return vazio

    centro = (start + end) / 2.0
    pontos = []
    for off in range(0, max(1, len(seq) - janela + 1), passo):
        jan_ini, jan_fim = ini + off, ini + off + janela
        # Windows touching the EVE are skipped: on a short EVE the window would be mostly
        # flank and would dilute its signal. The EVE GC enters as its own point, below.
        if jan_fim > start and jan_ini < end:
            continue
        trecho = seq[off:off + janela]
        validos = sum(trecho.count(b) for b in "ACGT")
        if validos < janela * 0.5:      # a window full of Ns carries no information
            continue
        pontos.append({"x": round((jan_ini + jan_fim) / 2.0 - centro, 1),
                       "gc": round(100.0 * (trecho.count("G") + trecho.count("C")) / validos, 2),
                       "zone": "upstream" if jan_fim <= start else "downstream"})

    trecho_eve = seq[max(0, start - ini):end - ini]
    validos_eve = sum(trecho_eve.count(b) for b in "ACGT")
    ponto_eve = {"x": 0.0, "gc": round(100.0 * (trecho_eve.count("G") + trecho_eve.count("C")) / validos_eve, 2),
                 "zone": "eve", "nt": validos_eve} if validos_eve else None
    return {"points": pontos, "eve_point": ponto_eve, "window": janela, "step": passo,
            "extent": extensao, "nt_upstream": start - ini, "nt_downstream": termino - end,
            "eve_half_width": round((end - start) / 2.0, 1)}


def _perfil_cpg(seq: str, ini_abs: int, centro: float, janela: int, passo: int) -> list:
    """
    Observed/expected CpG ratio per window: (nCpG * N) / (nC * nG).
    It sits near 1 when the dinucleotide appears at the expected frequency and drops well
    below that in methylated host DNA, which is useful for contrasting the EVE with its
    surroundings.
    """
    pontos = []
    for off in range(0, max(1, len(seq) - janela + 1), passo):
        t = seq[off:off + janela]
        nc, ng = t.count("C"), t.count("G")
        validos = sum(t.count(b) for b in "ACGT")
        if not (nc and ng and validos >= janela * 0.5):
            continue
        ncpg = t.count("CG")
        pontos.append({"x": round(ini_abs + off + janela / 2.0 - centro, 1),
                       "oe": round(ncpg * validos / (nc * ng), 3),
                       "cpg_pct": round(100.0 * ncpg / max(1, validos - 1), 2)})
    return pontos


def _genes_do_gff(gff: str, contig: str, ini: int, end: int, maximo: int = 60) -> list:
    """Host genes and exons overlapping the region, for the annotation track."""
    if not gff or not os.path.exists(gff):
        return []
    genes, exons = [], []
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith("#"):
                continue
            p = linha.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] != contig:
                continue
            try:
                a, b = int(p[3]) - 1, int(p[4])
            except ValueError:
                continue
            if b < ini or a > end:
                continue
            tipo = p[2]
            if tipo in ("gene", "pseudogene"):
                m = re.search(r"(?:Name|gene)=([^;\n]+)", p[8])
                genes.append({"start": a, "end": b, "strand": p[6],
                              "name": (m.group(1) if m else "")[:28], "type": tipo})
            elif tipo in ("exon", "CDS"):
                exons.append({"start": a, "end": b, "strand": p[6], "type": tipo})
    genes.sort(key=lambda g: g["start"])
    return [{"genes": genes[:maximo], "exons": exons[:400],
             "n_genes": len(genes), "n_exons": len(exons)}]


def _coordenadas(descricao: str):
    """Reuse the coordinate parser from main.py (orfipy and pre-extracted segments)."""
    import sys
    from evetool.main import _coordenadas_do_cabecalho
    return _coordenadas_do_cabecalho(descricao)


def _todos_os_loci(pasta: str) -> list:
    """[(contig, start, end, strand, orf)] for every candidate of the species."""
    caminho = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
    loci = []
    if not os.path.exists(caminho):
        return loci
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                c = _coordenadas(linha[1:].strip())
                if c:
                    loci.append((c[0], c[1], c[2], c[3], linha[1:].split()[0]))
    return loci


def _intervalos_mascarados(seq: str, tipo: str) -> list:
    """dustmasker (nt) / segmasker (aa) intervals of one sequence, 0-based inclusive [[a, b]]."""
    if not seq or not shutil.which("dustmasker" if tipo == "nt" else "segmasker"):
        return []
    exe = ["dustmasker", "-outfmt", "interval"] if tipo == "nt" else ["segmasker", "-outfmt", "interval"]
    p = subprocess.run(exe, input=f">s\n{seq}\n", capture_output=True, text=True, timeout=120)
    out = []
    for linha in p.stdout.splitlines():
        if " - " in linha:
            a, b = (int(x) for x in linha.split(" - "))
            out.append([a, b])
    return out


def _corridas_minusculas(seq: str, minimo: int = 20) -> list:
    """Runs of lower case (soft-masked repeats) as 0-based [[a, b)] on seq."""
    return [[m.start(), m.end()] for m in re.finditer(r"[a-z]{%d,}" % minimo, seq or "")]


def _repeticoes_salvas(outdir: str):
    try:
        import sys
        from evetool.annotation.repeats import carregar
        return carregar(outdir)
    except Exception:
        return None


def rodar_repeticoes(outdir: str, threads: int = 4, flanco: int = 20000, queries=None) -> dict:
    import sys
    from evetool.annotation.repeats import analisar
    r = analisar(outdir, flanco=flanco, threads=threads, queries=queries)
    # 'summary' is the page's own key in every server reply: the contexts go under another name
    return {"loci": r["loci"], "contexts": r["summary"], "directory": r["directory"]}


def _comparacao_especies(outdir: str):
    """The latest cross-species comparison: the page's run first, then the pipeline's."""
    cams = sorted((os.path.join(outdir, SUBPASTA, d, "cross_species_clusters.json") for d in MODOS_UNIDADE.values()),
                  key=lambda c: -os.path.getmtime(c) if os.path.exists(c) else 0)
    for cam in cams + [os.path.join(outdir, "11_sintenia_multiespecies", "cross_species_clusters.json")]:
        if os.path.exists(cam):
            try:
                with open(cam, encoding="utf-8") as fh:
                    d = json.load(fh)
                d["source"] = os.path.relpath(cam, outdir)
                d["mtime"] = os.path.getmtime(cam)
                return d
            except Exception:
                pass
    return None


def _anotacao_orfs(outdir: str, species: str, nomes: set) -> dict:
    """Class, best viral hit and curation state of the named ORFs of one species (hand edits win)."""
    caminho = os.path.join(outdir, TABELA)
    if not nomes or not os.path.exists(caminho):
        return {}
    estado = _ler_selecao(outdir)
    removidos = set(estado.get("removed", []))
    classes, familias = estado.get("proteins", {}), estado.get("families", {})
    saida = {}
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            nome = r.get("orf", "")
            if nome not in nomes or (species and r.get("species", "") not in ("", species)):
                continue
            q = r.get("query") or nome
            vazio = lambda v: "" if v in (None, "n/a", "NA") else v
            saida[nome] = {
                "protein_class": classes.get(q) or vazio(r.get("protein_class")),
                "family": familias.get(q) or vazio(r.get("viral_family")) or vazio(r.get("hmm_consensus_family")),
                "best_hit": vazio(r.get("viral_protein")), "product": vazio(r.get("viral_product")),
                "viral_bitscore": _float(r.get("viral_bitscore")), "group": r.get("group", ""),
                "hmm_model": vazio(r.get("hmm_best_model")), "query": q, "removed": q in removidos}
    return saida


def detalhar_eve(outdir: str, query: str, janela_vizinhos: int = 25000,
                 limite_contig_pequeno: int = 200000, extensao_perfil: int = 50000,
                 janela_gc: int = 2500, passo_gc: int = 500) -> dict:
    """
    Gather what is needed to judge an indistinguishable case: where the segment sits,
    whether other EVEs are nearby, how its GC compares with the flank and with the contig,
    and which viral alignment supports it.
    """
    species, _, orf = query.partition("|") if "|" in query else ("", "", query)
    pastas = pastas_especies(outdir)
    pasta = pastas.get(species) or pastas.get("") or outdir
    info = _info_pipeline(pasta)
    genoma = info.get("genome")
    flanco = int(info.get("flank") or 20000)

    # 1. ORF coordinates and sequence
    coord = seq_nt = None
    pep = []
    caminho_pep = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
    if os.path.exists(caminho_pep):
        with open(caminho_pep, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith(">"):
                    if coord:
                        break
                    if linha[1:].split()[0] == orf:
                        coord = _coordenadas(linha[1:].strip())
                elif coord:
                    pep.append(linha.strip())
    peptideo = "".join(pep).rstrip("*")
    caminho_nt = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas.fasta")
    if os.path.exists(caminho_nt):
        guardar, partes = False, []
        with open(caminho_nt, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith(">"):
                    if guardar:
                        break
                    guardar = linha[1:].split()[0] == orf
                elif guardar:
                    partes.append(linha.strip())
        seq_nt = "".join(partes) or None

    resultado = {"query": query, "species": species, "orf": orf, "genome": genoma}
    if not coord:
        resultado["error"] = "No coordinates found in the header of this ORF."
        return resultado
    contig, start, end, strand = coord
    resultado.update({"contig": contig, "start": start, "end": end, "strand": strand,
                      "length": end - start, "peptide_length": len(peptideo)})

    # 1b. Low complexity of the ORF: SEG on the peptide, DUST on its nucleotides (as aa positions)
    try:
        seg = _intervalos_mascarados(peptideo, "aa")
        dust = _intervalos_mascarados(seq_nt or "", "nt")
        n_aa = max(1, len(peptideo))
        resultado["low_complexity"] = {
            "seg": [[a + 1, b + 1] for a, b in seg],
            "dust": [[a // 3 + 1, min(n_aa, b // 3 + 1)] for a, b in dust],
            "seg_pct": round(100.0 * sum(b - a + 1 for a, b in seg) / n_aa, 1),
            "dust_pct": round(100.0 * sum(b - a + 1 for a, b in dust) / max(1, len(seq_nt or "")), 1),
        }
    except Exception as erro:
        resultado["low_complexity"] = {"error": str(erro)}
    # the saved repeat analysis of the locus this ORF belongs to
    rep_salvo = _repeticoes_salvas(outdir)
    if rep_salvo:
        resultado["repeats_locus"] = next((r for r in rep_salvo.get("loci", []) if query in r.get("orfs", [])), None)
        resultado["repeats_flank"] = rep_salvo.get("flank")

    # 2. GC composition: the EVE against its surroundings at several windows, and against the genome
    gc_eve = _gc(seq_nt) if seq_nt else 0.0
    cache = os.path.join(outdir, SUBPASTA, "cache_gc.json")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    perfil = _perfil_do_genoma(genoma, cache)
    dados_contig = perfil["contigs"].get(contig, {})
    tam_contig = dados_contig.get("tamanho", 0)
    gc_contig = dados_contig.get("gc", 0.0)
    gc_genoma = perfil.get("genome_gc", 0.0)

    # Short contig: the surroundings do not represent the genome, so only the global comparison holds
    contig_pequeno = bool(tam_contig) and tam_contig < limite_contig_pequeno
    janelas = [] if contig_pequeno else _gc_das_janelas(genoma, contig, start, end, tam_contig)

    comparacoes = [{"label": f"{j['window'] // 1000} kb", "ref_gc": j["gc"],
                    "delta": round(gc_eve - j["gc"], 2) if j["valid"] else None,
                    "valid": j["valid"], "nt_available": j["nt_available"]}
                   for j in janelas]
    if gc_contig and not contig_pequeno:
        comparacoes.append({"label": "contig", "ref_gc": gc_contig,
                            "delta": round(gc_eve - gc_contig, 2), "valid": True, "nt_available": tam_contig})
    if gc_genoma:
        comparacoes.append({"label": "genome", "ref_gc": gc_genoma,
                            "delta": round(gc_eve - gc_genoma, 2), "valid": True,
                            "nt_available": perfil["tamanho_genoma"]})

    resultado.update({
        "eve_gc": gc_eve, "contig_gc": gc_contig, "genome_gc": gc_genoma,
        "contig_length": tam_contig, "short_contig": contig_pequeno,
        "short_contig_limit": limite_contig_pequeno,
        "gc_comparisons": comparacoes,
        "gc_profile": _perfil_gc_regiao(genoma, contig, start, end, tam_contig,
                                       extensao=extensao_perfil, janela=janela_gc, passo=passo_gc),
        "query_sequence": seq_nt or "", "query_length": len(seq_nt) if seq_nt else 0,
    })

    # 3. Neighbouring EVEs on the same contig
    vizinhos = []
    for c, i, f, s, nome in _todos_os_loci(pasta):
        if c != contig or nome == orf:
            continue
        dist = 0 if (i < end and f > start) else (i - end if i >= end else start - f)
        if dist <= janela_vizinhos:
            vizinhos.append({"orf": nome, "start": i, "end": f, "strand": s, "distance": dist})
    vizinhos.sort(key=lambda v: v["distance"])
    resultado["neighbours"] = vizinhos[:12]
    resultado["n_neighbours"] = len(vizinhos)

    # 4. The cluster it belongs to, when the clustering has already run
    resultado["cluster"] = None
    for padrao in (os.path.join(pasta, "08_clusters_e_tensores", "*_segmentos.tsv"),
                   os.path.join(outdir, SUBPASTA, "multiespecies", "clusters_multispecies.tsv")):
        for caminho in glob.glob(padrao):
            with open(caminho, encoding="utf-8", errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    alvo = r.get("segmento") or r.get("eve") or ""
                    if orf in alvo or query == alvo:
                        resultado["cluster"] = r.get("cluster")
                        break
            if resultado["cluster"]:
                break

    # 4b. Cluster context: the ORFs that form it, CpG and the host annotation
    contexto = {"ok": False}
    try:
        flanco_ctx = min(max(2000, (end - start) * 6), 30000)
        ini_r, fim_r = max(0, start - flanco_ctx), min(tam_contig or (end + flanco_ctx), end + flanco_ctx)
        seq_reg = ""
        if genoma and os.path.exists(genoma) and shutil.which("samtools"):
            pr = subprocess.run(["samtools", "faidx", genoma, f"{contig}:{ini_r + 1}-{fim_r}"],
                                capture_output=True, text=True, timeout=180)
            seq_bruta = "".join(l.strip() for l in pr.stdout.splitlines() if not l.startswith(">"))
            seq_reg = seq_bruta.upper()
        else:
            seq_bruta = ""
        centro_r = (start + end) / 2.0
        jan_cpg = max(200, (fim_r - ini_r) // 120)
        # Candidate ORFs inside the window: these are the ones that make up the cluster
        orfs_regiao = [{"orf": nome, "start": i, "end": f, "strand": st,
                        "is_this": nome == orf}
                       for c2, i, f, st, nome in _todos_os_loci(pasta)
                       if c2 == contig and f >= ini_r and i <= fim_r]
        # Protein and best viral hit of each ORF, so close clusters can be told apart on the drawing
        anot_orfs = _anotacao_orfs(outdir, species, {o["orf"] for o in orfs_regiao})
        for o in orfs_regiao:
            o.update(anot_orfs.get(o["orf"], {}))
        orfs_regiao.sort(key=lambda o: o["start"])
        # Cluster extent according to the merged BED, when it exists
        cluster_reg = None
        for nome_bed in ("eve_hits_diretos.merged.bed", "eve_clusters.bed"):
            cam = os.path.join(pasta, "06_anotacao_genomica", nome_bed)
            if os.path.exists(cam):
                with open(cam, encoding="utf-8", errors="replace") as fh:
                    for linha in fh:
                        p = linha.split("\t")
                        if len(p) >= 3 and p[0] == contig:
                            a, b = int(p[1]), int(p[2])
                            if a <= start and b >= end:
                                cluster_reg = {"start": a, "end": b,
                                               "name": p[3] if len(p) > 3 else "cluster"}
                                break
                if cluster_reg:
                    break
        anot = _genes_do_gff(info.get("gff"), contig, ini_r, fim_r)
        contexto = {"ok": True, "start_region": ini_r, "end": fim_r, "centre": centro_r,
                    "cpg_window": jan_cpg,
                    "cpg": _perfil_cpg(seq_reg, ini_r, centro_r, jan_cpg, max(100, jan_cpg // 2)) if seq_reg else [],
                    # repeats along the window: soft-masked runs of the assembly and DUST low complexity
                    "softmask": [[ini_r + a, ini_r + b] for a, b in _corridas_minusculas(seq_bruta)],
                    "softmask_available": any(ch.islower() for ch in seq_bruta[:200000]) if seq_bruta else False,
                    "dust": [[ini_r + a, ini_r + b + 1] for a, b in _intervalos_mascarados(seq_reg, "nt")] if seq_reg else [],
                    "orfs": orfs_regiao, "cluster": cluster_reg,
                    "annotation": anot[0] if anot else None,
                    "has_gff": bool(info.get("gff") and os.path.exists(info.get("gff")))}
    except Exception as erro:
        contexto = {"ok": False, "error": str(erro)}
    resultado["cluster_context"] = contexto

    # 5. HMM domains: a layer independent of the alignment
    dominios = []
    cam_anot = os.path.join(pasta, "04_candidatas_eves", "hmm_annotation.tsv")
    if os.path.exists(cam_anot):
        with open(cam_anot, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf") == orf:
                    dominios.append({k: r.get(k, "") for k in
                                     ("model", "score", "ievalue", "model_coverage",
                                      "orf_from", "orf_to", "hmm_family", "hmm_genome", "source")})
    dominios.sort(key=lambda d: -float(d["score"] or 0))
    resultado["hmm_domains"] = dominios[:10]
    resultado["hmm_n_domains"] = len(dominios)
    resultado["hmm_summary"] = None
    cam_res = os.path.join(pasta, "04_candidatas_eves", "hmm_summary.tsv")
    if os.path.exists(cam_res):
        with open(cam_res, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf") == orf and int(r.get("n_domains") or 0):
                    resultado["hmm_summary"] = r
                    break

    # 5b. NCBI CD-Search domains (every scope, so the excluded ones can be checked too)
    resultado["cdd_domains"] = []
    cam_cdd = os.path.join(pasta, "04_candidatas_eves", "cdsearch_hits.tsv")
    if os.path.exists(cam_cdd):
        with open(cam_cdd, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf") == orf and r.get("hit_type") != "superfamily":
                    resultado["cdd_domains"].append({k: r.get(k, "") for k in (
                        "scope", "hit_type", "accession", "short_name", "from", "to", "evalue",
                        "bitscore", "incomplete", "description", "taxon")})
        resultado["cdd_domains"].sort(key=lambda d: -float(d["bitscore"] or 0))
    resultado["cdd_searched"] = os.path.exists(cam_cdd)

    # 5b'. Pfam-A domains (hmmscan with the Pfam gathering thresholds)
    resultado["pfam_domains"] = []
    cam_pfam = os.path.join(pasta, "04_candidatas_eves", "pfam_hits.tsv")
    if os.path.exists(cam_pfam):
        with open(cam_pfam, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("orf") == orf:
                    resultado["pfam_domains"].append({k: r.get(k, "") for k in (
                        "scope", "accession", "short_name", "from", "to", "model_coverage",
                        "ievalue", "bitscore", "description")})
        resultado["pfam_domains"].sort(key=lambda d: -float(d["bitscore"] or 0))
    resultado["pfam_searched"] = os.path.exists(cam_pfam)

    # 5c. Online BLASTX (NCBI) of this EVE or of any searched segment that covers it
    try:
        import sys
        from evetool.discovery.ncbi_blastx import hits_do_segmento
        resultado["blastx_ncbi"] = hits_do_segmento(outdir, contig, start + 1, end, query)[:40]
    except Exception as erro:
        resultado["blastx_ncbi"] = []
        log(f"[warning] online BLASTX hits not read: {erro}")

    # 6. The viral alignment that supports the candidate
    alinhamentos = []
    tsv = os.path.join(pasta, "02_alinhamento_diamond", "triagem_viral_diamond.tsv")
    if os.path.exists(tsv):
        with open(tsv, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                p = linha.rstrip("\n").split("\t")
                if len(p) > 8 and p[0] == orf:
                    campos = p[4].split("|")
                    alinhamentos.append({
                        "acc": p[2], "product": (campos[1].strip() if len(campos) > 1 else "")[:80],
                        "viral_genome": campos[2].strip() if len(campos) > 2 else "",
                        "family": campos[5].strip() if len(campos) > 5 else "",
                        "pident": float(p[5]), "evalue": p[6], "bitscore": float(p[7]),
                        "qcov": float(p[8]), **_posicao_no_orf(seq_nt, p[10] if len(p) > 10 else ""),
                    })
    alinhamentos.sort(key=lambda a: -a["bitscore"])
    resultado["alignments"] = alinhamentos[:8]
    resultado["n_alignments"] = len(alinhamentos)
    return resultado


def _posicao_no_orf(seq_nt: str, qseq: str) -> dict:
    """Amino-acid stretch of the ORF covered by an alignment, found by its nucleotide qseq."""
    trecho = (qseq or "").replace("-", "").upper()
    if not seq_nt or len(trecho) < 9:
        return {"orf_from": None, "orf_to": None}
    i = seq_nt.upper().find(trecho)
    if i < 0:
        return {"orf_from": None, "orf_to": None}
    return {"orf_from": i // 3 + 1, "orf_to": (i + len(trecho)) // 3}


# ---------------------------------------------------------------------------
# Chart export
# ---------------------------------------------------------------------------
def exportar_grafico(outdir: str, svg: str, nome: str = "grafico") -> dict:
    """
    Save the SVG sent by the page and, when a converter is available, produce the PDF too.
    The SVG arrives with the colours already resolved (the page does not use CSS variables
    on export), so the file looks the same outside the browser.
    """
    import re as _re
    seguro = _re.sub(r"[^A-Za-z0-9_.-]+", "_", nome).strip("_") or "grafico"
    destino_dir = os.path.join(outdir, SUBPASTA, "graficos")
    os.makedirs(destino_dir, exist_ok=True)
    caminho_svg = os.path.join(destino_dir, f"{seguro}.svg")
    with open(caminho_svg, "w", encoding="utf-8") as fh:
        fh.write(svg)

    caminho_pdf, erro = None, None
    conversor = shutil.which("rsvg-convert")
    if conversor:
        caminho_pdf = os.path.join(destino_dir, f"{seguro}.pdf")
        try:
            r = subprocess.run([conversor, "-f", "pdf", "-o", caminho_pdf, caminho_svg],
                               capture_output=True, text=True, timeout=180)
            if r.returncode != 0 or not os.path.exists(caminho_pdf) or os.path.getsize(caminho_pdf) == 0:
                erro = (r.stderr or "rsvg-convert falhou").strip()[:300]
                caminho_pdf = None
        except Exception as e:
            erro, caminho_pdf = str(e), None
    else:
        erro = "rsvg-convert not found on the PATH; only the SVG was written."

    log(f"chart exported: {caminho_svg}" + (f" and {caminho_pdf}" if caminho_pdf else ""))
    return {"svg": caminho_svg, "pdf": caminho_pdf, "error": erro}


if __name__ == "__main__":
    main()
