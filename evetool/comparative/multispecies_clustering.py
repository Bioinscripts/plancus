#!/usr/bin/env python3
"""
Two-stage clustering of EVEs across species.

Replaces the connected components of cross_species_comparator (which joined any pair at
e-value < 1e-5 and, through single linkage, merged distinct families) with the same
criterion used in the per-species analyses:

  Stage 1 - similarity split: only pairs with identity and coverage above the threshold
            (30% by default) become edges. The resulting components separate the groups
            that are distant from one another.
  Stage 2 - valleys: inside each component, recursive hierarchical splitting at the valley
            of the distance distribution (hierarchical_clustering.find_valley), with complete
            linkage, so every member of a cluster is similar to every other, not merely
            chained to it.

Outputs in the given directory:
  clusters_multispecies.tsv              one line per EVE, with its group and cluster
  cross_species_clusters_summary.tsv     one cluster per line, with evolutionary status
                                         and the synteny columns (empty when not requested)
  cross_species_presence_absence_matrix.csv
  cross_species_sharing_heatmap.png/.pdf the clusters whose sharing is supported (an
                                         inter-species pair above --share-identity and
                                         --share-coverage), species and clusters ordered
                                         by hierarchical clustering
  cross_species_clusters.json            members and tblastx pairs of every cluster, read by
                                         the "Cross-species" tab of the curation page
  cluster_multispecies_thresholds.tsv    threshold of each valley split

Usage:
  python -m evetool.comparative.multispecies_clustering -f eves_curadas_nt.fasta -o 11_sintenia_multiespecies
  python -m evetool.comparative.multispecies_clustering --xml tblastx.xml -o saida --min-identity 30 --min-coverage 30
"""
import os
import re
import sys
import csv
import shutil
import argparse
import subprocess
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from evetool.comparative.hierarchical_clustering import parse_pairs, build_matrices, split_recursive

# Single-hue sequential ramp (light -> dark). Never a rainbow for magnitude data.
CMAP_SEQ = LinearSegmentedColormap.from_list("eve_blue", ["#f4f7fb", "#bcd4ef", "#6fa3e0", "#2a78d6", "#14407a"])
COR_TEXTO = "#1d2430"
COR_GRID = "#d9dbd5"


def log(msg):
    print(f"[multi] {msg}", flush=True)


def especie_de(nome: str) -> str:
    """Headers shaped as 'Species|ORF'; without the prefix the species is empty."""
    return nome.split("|", 1)[0] if "|" in nome else ""


def carregar_renomeacao(caminho: str) -> dict:
    """Two-column TSV: current_name <TAB> correct_name."""
    mapa = {}
    if caminho and os.path.exists(caminho):
        with open(caminho, encoding="utf-8") as fh:
            for linha in fh:
                linha = linha.split("#", 1)[0].rstrip("\n")
                partes = [p.strip() for p in linha.split("\t")]
                if len(partes) >= 2 and partes[0] and partes[1]:
                    mapa[partes[0]] = partes[1]
        log(f"{len(mapa)} species will be renamed from {caminho}")
    return mapa


def rodar_tblastx(fasta: str, diretorio: str, threads: int = 4, evalue: float = 1e-3) -> str:
    for exe in ("makeblastdb", "tblastx"):
        if shutil.which(exe) is None:
            raise RuntimeError(f"'{exe}' not found on the PATH.")
    os.makedirs(diretorio, exist_ok=True)
    db = os.path.join(diretorio, "banco_multiespecies_eve")
    xml = os.path.join(diretorio, "all_vs_all_tblastx_multiespecies.xml")
    log("formatting the database and running all-vs-all tblastx...")
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    subprocess.run(["tblastx", "-query", fasta, "-db", db, "-out", xml, "-evalue", str(evalue),
                    "-num_threads", str(threads), "-outfmt", "5"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return xml


def separar_por_similaridade(names, df, min_identity: float, min_coverage: float,
                             modo: str = "max") -> list:
    """
    Stage 1. Connected components over the pairs above the thresholds only.
    The parse_pairs dataframe is already symmetric (one row per pair), so 'max' and 'mean'
    coincide; the parameter is kept for compatibility with the per-species flow.
    """
    import networkx as nx
    g = nx.Graph()
    g.add_nodes_from(names)
    n_arestas = 0
    for r in df.itertuples():
        ident = getattr(r, "identity_pct", np.nan)
        cov = getattr(r, "coverage_pct", np.nan)
        if np.isnan(ident) or np.isnan(cov):
            continue
        if ident >= min_identity and cov >= min_coverage:
            g.add_edge(r.seg_a, r.seg_b)
            n_arestas += 1
    comps = sorted(nx.connected_components(g), key=len, reverse=True)
    log(f"stage 1 (identity >= {min_identity}%, coverage >= {min_coverage}%): "
        f"{n_arestas} edges, {len(comps)} groups "
        f"({sum(1 for c in comps if len(c) == 1)} singletons)")
    return [sorted(c) for c in comps]


def dividir_por_vales(names, D, componentes, cfg) -> tuple:
    """Stage 2. Recursive splitting at the valley of the distance distribution, per component."""
    class _Args:
        max_depth = cfg.max_depth
        min_pairs = cfg.min_pairs
        valley_ratio = cfg.valley_ratio
        max_bandwidth = cfg.max_bandwidth

    idx_de = {n: i for i, n in enumerate(names)}
    cluster_de, splits, n_cluster = {}, [], 0
    for n_comp, membros in enumerate(componentes, start=1):
        idx = np.array([idx_de[m] for m in membros])
        if len(idx) < 3:
            n_cluster += 1
            for m in membros:
                cluster_de[m] = (f"G{n_comp}", f"C{n_cluster}")
            continue
        sub_splits = []
        partes = split_recursive(D, idx, "", 0, _Args, sub_splits)
        for s in sub_splits:
            s["similarity_group"] = f"G{n_comp}"
            splits.append({k: v for k, v in s.items() if k != "_info"})
        for rotulo, membros_idx in partes:
            n_cluster += 1
            for i in membros_idx:
                cluster_de[names[i]] = (f"G{n_comp}", f"C{n_cluster}")
    log(f"stage 2 (valleys): {n_cluster} final clusters from {len(componentes)} groups")
    return cluster_de, splits


def avaliar_sintenia(cluster_de: dict, cfg) -> dict:
    """
    Flank synteny of the shared clusters (flank_synteny.py): the nearest host genes on each
    side of every locus, labelled by homology to a reference proteome, compared between the
    loci of different species. Same neighbours on both sides = probably one ancestral
    insertion; different neighbours = independent integrations of related viruses.

    Off by default: --sintenia plus --proteoma-ref, and the coordinates of every unit
    (cfg.locais, given by the caller, which knows the genomes).
    """
    clusters = defaultdict(list)
    for eve, (_g, c) in cluster_de.items():
        clusters[c].append(eve)
    base = {"synteny_status": "not_requested", "synteny_orthologs": "", "synteny_n_shared": 0,
            "synteny_pairs": 0, "synteny_syntenic_pairs": 0}
    cfg._sintenia_detalhe = None
    if not cfg.sintenia:
        return {c: dict(base) for c in clusters}
    if not cfg.proteoma_ref:
        log("[warning] synteny requested without --proteoma-ref; the columns stay as 'no_reference_proteome'.")
        return {c: dict(base, synteny_status="no_reference_proteome") for c in clusters}
    locais = getattr(cfg, "locais", None)
    if not locais:
        log("[warning] synteny requested without the coordinates of the units; not evaluated.")
        return {c: dict(base, synteny_status="not_evaluated") for c in clusters}
    from evetool.comparative.flank_synteny import avaliar
    res = avaliar(dict(clusters), locais, getattr(cfg, "genomas_sintenia", {}) or {},
                  getattr(cfg, "origens", {}) or {}, cfg.proteoma_ref, int(cfg.flanco_sintenia or 20000),
                  os.path.join(cfg.outdir, "synteny"), threads=cfg.threads, especie_de=especie_de,
                  max_sitios=int(getattr(cfg, "max_sitios_sintenia", 3) or 3))
    cfg._sintenia_detalhe = res
    return res["clusters"]


def resumir(cluster_de: dict, mapa_nomes: dict, sintenia: dict = None) -> tuple:
    """Build the per-cluster summary and the presence/absence matrix per species."""
    membros_por_cluster = defaultdict(list)
    for eve, (group, cluster) in cluster_de.items():
        membros_por_cluster[(group, cluster)].append(eve)

    especies = sorted({mapa_nomes.get(especie_de(e), especie_de(e)) for e in cluster_de})
    resumo, matriz = [], []
    for (group, cluster), membros in sorted(membros_por_cluster.items(),
                                            key=lambda kv: -len(kv[1])):
        sp_membros = [mapa_nomes.get(especie_de(m), especie_de(m)) for m in membros]
        presentes = sorted(set(sp_membros))
        status = "Ancestral / Shared" if len(presentes) >= 2 else "Species-specific"
        linha_resumo = {
            "Cluster_ID": cluster, "Similarity_Group": group,
            "Total_Sequences": len(membros), "Total_Species": len(presentes),
            "Evolutionary_Status": status, "Especies": ";".join(presentes),
            "Example_Members": ";".join(sorted(membros)[:5]),
        }
        linha_resumo.update((sintenia or {}).get(cluster, {}))
        resumo.append(linha_resumo)
        linha = {"Cluster_ID": cluster, "Status": status, "Total_Copies": len(membros)}
        contagem = Counter(sp_membros)
        for sp in especies:
            linha[sp] = contagem.get(sp, 0)
        matriz.append(linha)
    return pd.DataFrame(resumo), pd.DataFrame(matriz), especies


def pares_entre_especies(cluster_de: dict, df: pd.DataFrame, mapa_nomes: dict) -> dict:
    """
    tblastx pairs inside each cluster: {cluster: [(a, b, identity, coverage, bitscore_pct)]}.
    A cluster is shared in a meaningful sense only when segments of different species align
    well to each other; those inter-species pairs are what tell it apart from a cluster of
    distant segments that the valleys kept together.
    """
    sp = lambda e: mapa_nomes.get(especie_de(e), especie_de(e))
    pares = defaultdict(list)
    for r in df.itertuples():
        ca, cb = cluster_de.get(r.seg_a), cluster_de.get(r.seg_b)
        if not ca or not cb or ca[1] != cb[1]:
            continue
        ident = 0.0 if np.isnan(r.identity_pct) else float(r.identity_pct)
        pares[ca[1]].append((r.seg_a, r.seg_b, ident, float(r.coverage_pct), float(r.bitscore_pct),
                             sp(r.seg_a) != sp(r.seg_b)))
    return pares


def estatisticas_compartilhamento(df_resumo: pd.DataFrame, pares: dict, ident_min: float, cov_min: float):
    """Inter-species identity/coverage per cluster and whether the sharing is supported."""
    linhas = []
    for c, st in zip(df_resumo["Cluster_ID"], df_resumo["Total_Species"]):
        inter = [p for p in pares.get(c, ()) if p[5]]
        fortes = [p for p in inter if p[2] >= ident_min and p[3] >= cov_min]
        linhas.append({
            "Inter_Species_Pairs": len(inter),
            "Inter_Identity_Max": round(max((p[2] for p in inter), default=0.0), 1),
            "Inter_Identity_Mean": round(float(np.mean([p[2] for p in inter])) if inter else 0.0, 1),
            "Inter_Coverage_Max": round(max((p[3] for p in inter), default=0.0), 1),
            "Inter_Coverage_Mean": round(float(np.mean([p[3] for p in inter])) if inter else 0.0, 1),
            "Inter_Supported_Pairs": len(fortes),
            "Sharing": ("supported" if fortes else "distant") if st >= 2 else "species-specific",
        })
    return pd.concat([df_resumo.reset_index(drop=True), pd.DataFrame(linhas)], axis=1)


def salvar_json(destino: str, cluster_de: dict, pares: dict, df_resumo: pd.DataFrame,
                especies: list, mapa_nomes: dict, cfg):
    """Everything the page needs to redraw the comparison with its own thresholds."""
    import json
    sp = lambda e: mapa_nomes.get(especie_de(e), especie_de(e))
    membros = defaultdict(list)
    for eve, (_g, c) in cluster_de.items():
        membros[c].append(eve)
    clusters = []
    for r in df_resumo.itertuples():
        ms = sorted(membros[r.Cluster_ID])
        ix = {m: i for i, m in enumerate(ms)}
        clusters.append({
            "id": r.Cluster_ID, "group": r.Similarity_Group,
            "members": [{"eve": m, "species": sp(m)} for m in ms],
            # [i, j, identity %, coverage %, bitscore % of the self-hit]
            "pairs": [[ix[a], ix[b], round(i_, 1), round(c_, 1), round(b_, 1)]
                      for a, b, i_, c_, b_, _x in pares.get(r.Cluster_ID, ())],
        })
    det = getattr(cfg, "_sintenia_detalhe", None)
    if det:
        for c in clusters:
            ms = [m["eve"] for m in c["members"]]
            ix = {m: i for i, m in enumerate(ms)}
            c["synteny"] = det["clusters"].get(c["id"])
            c["synteny_pairs"] = [[ix[a], ix[b], r["status"], r.get("left", []), r.get("right", []), r.get("orientation", "")]
                                  for a, b, r in det["pairs"].get(c["id"], ())]
            for m in c["members"]:
                f = det["flanks"].get(m["eve"])
                if f:
                    m["flank"] = {"L": f["L"], "R": f["R"], "Lp": f.get("Lp", []), "Rp": f.get("Rp", []),
                                  "where": f.get("where", ""), "status": f.get("status", "")}
    with open(destino, "w", encoding="utf-8") as fh:
        json.dump({"species": especies, "clusters": clusters,
                   "units": getattr(cfg, "unidades_modo", "orfs"), "units_gap": getattr(cfg, "unidades_gap", None),
                   "units_members": getattr(cfg, "unidades_membros", {}) or {},
                   "synteny": {"flank": getattr(cfg, "flanco_sintenia", None), "proteome": getattr(cfg, "proteoma_ref", None),
                               "names": det["names"] if det else {}, "excluded": det.get("excluded", {}),
                               "max_sites": getattr(cfg, "max_sitios_sintenia", 3)} if det else None,
                   "min_identity": cfg.min_identity, "min_coverage": cfg.min_coverage,
                   "share_identity": getattr(cfg, "share_identity", 50.0),
                   "share_coverage": getattr(cfg, "share_coverage", 50.0)}, fh, ensure_ascii=False)


def _ordem_hierarquica(M: np.ndarray, metrica: str = "jaccard") -> list:
    """Leaf order of an average-linkage tree over the rows of M (identity order when < 3 rows)."""
    if len(M) < 3:
        return list(range(len(M))), None
    from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
    from scipy.spatial.distance import pdist
    X = (M > 0).astype(bool) if metrica == "jaccard" else M
    d = np.nan_to_num(pdist(X, metric=metrica), nan=1.0)
    Z = linkage(d, method="average")
    try:
        Z = optimal_leaf_ordering(Z, d)
    except Exception:
        pass
    return list(leaves_list(Z)), Z


def heatmap(df_matriz: pd.DataFrame, especies: list, destino_png: str, max_clusters: int = 40,
            suportados: set = None):
    """
    Heatmap of the clusters present in >= 2 species only: those carry the evolutionary
    signal. Zero cells stay blank and unlabelled, because annotating every zero clutters

    """
    compartilhados = df_matriz[df_matriz["Status"].str.startswith("Ancestral")]
    n_distantes = 0
    if suportados is not None:
        n_distantes = int((~compartilhados["Cluster_ID"].isin(suportados)).sum())
        compartilhados = compartilhados[compartilhados["Cluster_ID"].isin(suportados)]
    if compartilhados.empty:
        log("no cluster with supported sharing between species; heatmap not generated.")
        return None

    m = compartilhados.set_index("Cluster_ID")[especies]
    m = m.assign(_n=(m > 0).sum(axis=1), _t=m.sum(axis=1)) \
         .sort_values(["_n", "_t"], ascending=False).drop(columns=["_n", "_t"])
    cortado = len(m) - max_clusters
    m = m.head(max_clusters)
    # Species with no copy among the displayed clusters drop off the axis
    m = m.loc[:, (m > 0).any(axis=0)]
    # Both axes ordered by average-linkage trees over the presence pattern (Jaccard), so the
    # species that share clusters sit side by side and so do the clusters with the same species
    ord_l, Z_l = _ordem_hierarquica(m.to_numpy(dtype=float))
    ord_c, Z_c = _ordem_hierarquica(m.to_numpy(dtype=float).T)
    m = m.iloc[ord_l, ord_c]

    n_l, n_c = m.shape
    from matplotlib import gridspec
    from scipy.cluster.hierarchy import dendrogram
    larg, alt = max(6.0, 1.05 * n_c + 3.6), max(3.8, 0.33 * n_l + 3.0)
    fig = plt.figure(figsize=(larg, alt))
    # the colour bar gets its own column, so the heatmap and both trees keep the same extent
    gs = gridspec.GridSpec(2, 4, width_ratios=[1.3, max(4.0, larg - 3.0), 0.9, 0.25],
                           height_ratios=[1.0, max(2.5, alt - 1.0)], wspace=0.02, hspace=0.01)
    ax = fig.add_subplot(gs[1, 1])
    for Z, pos, orient, n in ((Z_c, gs[0, 1], "top", n_c), (Z_l, gs[1, 0], "left", n_l)):
        axd = fig.add_subplot(pos)
        if Z is not None:
            dendrogram(Z, orientation=orient, ax=axd, no_labels=True, color_threshold=0,
                       above_threshold_color="#7b7a75", link_color_func=lambda _k: "#7b7a75")
        # scipy puts leaf k at 10k + 5: match the cell centres of the heatmap
        if orient == "top":
            axd.set_xlim(0, 10 * n)
        else:
            axd.set_ylim(10 * n, 0)
        axd.set_axis_off()
    dados = m.to_numpy(dtype=float)
    mascarado = np.ma.masked_where(dados == 0, dados)
    cmap = CMAP_SEQ.copy()
    cmap.set_bad("#ffffff")
    im = ax.imshow(mascarado, cmap=cmap, aspect="auto", vmin=1, vmax=max(2, dados.max()))

    ax.set_xticks(range(n_c))
    ax.set_xticklabels(m.columns, rotation=38, ha="right", fontsize=9.5, color=COR_TEXTO)
    ax.set_yticks(range(n_l))
    ax.set_yticklabels(m.index, fontsize=8.5, color=COR_TEXTO)
    ax.set_xticks(np.arange(-.5, n_c, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n_l, 1), minor=True)
    ax.grid(which="minor", color=COR_GRID, linewidth=0.6)
    ax.tick_params(which="both", length=0)
    for lado in ax.spines.values():
        lado.set_visible(False)

    limite = im.norm(dados.max()) if dados.max() else 1
    for i in range(n_l):
        for j in range(n_c):
            v = int(dados[i, j])
            if v:  # zero stays blank
                ax.text(j, i, v, ha="center", va="center", fontsize=8.5,
                        color="#ffffff" if im.norm(v) > 0.55 * limite else COR_TEXTO)

    ax.yaxis.tick_right()
    titulo = "EVE clusters shared between species"
    if cortado > 0:
        titulo += f" (the {max_clusters} largest of {len(compartilhados)})"
    if n_distantes:
        titulo += f"\n{n_distantes} shared cluster(s) of distant segments left out"
    fig.suptitle(titulo, fontsize=12.5, fontweight="bold", color=COR_TEXTO)
    ax.set_xlabel("Host species", fontsize=10.5, color=COR_TEXTO, labelpad=8)
    cax = fig.add_subplot(gs[1, 3])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("EVE copies", fontsize=9.5, color=COR_TEXTO)
    cb.outline.set_visible(False)
    fig.savefig(destino_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(destino_png.replace(".png", ".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"heatmap with {n_l} shared clusters x {n_c} species -> {destino_png}")
    return destino_png


def executar(cfg) -> dict:
    os.makedirs(cfg.outdir, exist_ok=True)
    xml = cfg.xml or rodar_tblastx(cfg.fasta, cfg.outdir, threads=cfg.threads)

    names, df = parse_pairs(xml, cfg.max_overlap_nt)
    if len(names) < 2:
        raise RuntimeError(f"Only {len(names)} segment(s) in the XML; nothing to compare.")
    mapa_nomes = carregar_renomeacao(cfg.renomear)

    componentes = separar_por_similaridade(names, df, cfg.min_identity, cfg.min_coverage, cfg.filter_mode)
    D, _Dgc, _BP = build_matrices(names, df, cfg.w_cov, cfg.w_gap)
    cluster_de, splits = dividir_por_vales(names, D, componentes, cfg)

    sintenia = avaliar_sintenia(cluster_de, cfg)
    df_resumo, df_matriz, especies = resumir(cluster_de, mapa_nomes, sintenia)
    ident_share = getattr(cfg, "share_identity", 50.0)
    cov_share = getattr(cfg, "share_coverage", 50.0)
    pares = pares_entre_especies(cluster_de, df, mapa_nomes)
    df_resumo = estatisticas_compartilhamento(df_resumo, pares, ident_share, cov_share)
    suportados = set(df_resumo.loc[df_resumo["Sharing"] == "supported", "Cluster_ID"])

    saida = lambda n: os.path.join(cfg.outdir, n)
    with open(saida("clusters_multispecies.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["eve", "species", "similarity_group", "cluster"])
        for eve, (group, cluster) in sorted(cluster_de.items()):
            sp = especie_de(eve)
            w.writerow([eve, mapa_nomes.get(sp, sp), group, cluster])
    df_resumo.to_csv(saida("cross_species_clusters_summary.tsv"), sep="\t", index=False)
    det = getattr(cfg, "_sintenia_detalhe", None)
    if det:
        with open(saida("synteny_pairs.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["cluster", "eve_a", "eve_b", "status", "orientation", "shared_left", "shared_right",
                        "flank_a", "flank_b"])
            for c, lista in det["pairs"].items():
                for a, b, r in lista:
                    nm = lambda gs: ";".join(det["names"].get(g, g) for g in gs)
                    w.writerow([c, a, b, r["status"], r.get("orientation", ""), nm(r.get("left", [])),
                                nm(r.get("right", [])), det["flanks"][a].get("where", ""), det["flanks"][b].get("where", "")])
    df_matriz.to_csv(saida("cross_species_presence_absence_matrix.csv"), index=False)
    if splits:
        pd.DataFrame(splits).to_csv(saida("cluster_multispecies_thresholds.tsv"), sep="\t", index=False)
    df.to_csv(saida("cluster_multispecies_pairs.tsv"), sep="\t", index=False)

    salvar_json(saida("cross_species_clusters.json"), cluster_de, pares, df_resumo, especies, mapa_nomes, cfg)
    png = heatmap(df_matriz, especies, saida("cross_species_sharing_heatmap.png"),
                  max_clusters=cfg.max_clusters_heatmap, suportados=suportados)

    comp = Counter(df_resumo["Evolutionary_Status"])
    log(f"{len(cluster_de)} EVEs | {len(df_resumo)} clusters | "
        f"{comp.get('Ancestral / Shared', 0)} shared, "
        f"{comp.get('Species-specific', 0)} species-specific; "
        f"{len(suportados)} with supported sharing (inter-species pair >= {ident_share}% identity "
        f"and >= {cov_share}% coverage)")
    return {"clusters": cluster_de, "resumo": df_resumo, "matriz": df_matriz,
            "heatmap": png, "xml": xml, "n_especies": len(especies), "n_suportados": len(suportados)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--fasta", help="Multi-species FASTA (headers 'Species|ORF')")
    ap.add_argument("--xml", help="Ready all-vs-all tblastx XML (skips the search)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--min-identity", type=float, default=30.0, help="Minimum identity for stage 1 (%%)")
    ap.add_argument("--min-coverage", type=float, default=30.0, help="Minimum coverage for stage 1 (%%)")
    ap.add_argument("--filter-mode", default="max", choices=["max", "mean", "both"])
    ap.add_argument("--valley-ratio", type=float, default=0.5)
    ap.add_argument("--max-bandwidth", type=float, default=0.05)
    ap.add_argument("--min-pairs", type=int, default=15)
    ap.add_argument("--max-depth", type=int, default=6)
    ap.add_argument("--w-cov", type=float, default=1.0)
    ap.add_argument("--w-gap", type=float, default=1.0)
    ap.add_argument("--max-overlap-nt", type=int, default=0)
    ap.add_argument("--max-clusters-heatmap", type=int, default=40)
    ap.add_argument("--share-identity", type=float, default=50.0,
                    help="A shared cluster enters the heatmap only with an inter-species tblastx pair "
                         "at or above this identity (%%, default 50)")
    ap.add_argument("--share-coverage", type=float, default=50.0,
                    help="... and at or above this coverage (%%, default 50)")
    ap.add_argument("--renomear", help="TSV 'current_name<TAB>correct_name' to fix species labels")
    ap.add_argument("--sintenia", action="store_true",
                    help="Evaluate flank synteny for the clusters (off by default). It can also be requested "
                         "per species, through the 'sintenia' column of the species TSV")
    ap.add_argument("--proteoma-ref",
                    help="Reference proteome (FASTA or .dmnd) used to label the flanking genes; required "
                         "when --sintenia is used")
    ap.add_argument("--flanco-sintenia", type=int, default=50000,
                    help="Flank size used by the synteny step (default: 50000 nt)")
    cfg = ap.parse_args()
    if not cfg.fasta and not cfg.xml:
        ap.error("provide --fasta or --xml.")
    executar(cfg)


if __name__ == "__main__":
    main()
