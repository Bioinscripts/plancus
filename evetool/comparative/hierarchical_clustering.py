#!/usr/bin/env python3
"""
Clusters de EVEs por bitscore (%), com limiares tirados dos "vales" entre picos
e figuras (heatmap em ggplot, dendrogramas em imagens separadas).

Entrada: BLAST tblastx todos-contra-todos em XML (outfmt 5), com o
auto-alinhamento de cada segmento presente (segmento x ele mesmo).

Métricas por par (HSPs não sobrepostos na query; direção com maior bitscore):
  bitscore_pct  soma do bitscore dos HSPs / MAIOR auto-score do par   (0-100)
  identity_pct  identidades / comprimento alinhado (aminoácidos)
  coverage_pct  nt alinhados na query / MAIOR comprimento do par
  gap_pct       % do trecho entre o 1o e o ultimo HSP que NAO esta alinhado
                (tblastx nao tem gaps no HSP; o gap e o espaco entre HSPs;
                 Hsp_gaps, se existir, tambem entra)

Clusterização
  distancia D = 1 - bitscore_pct/100, linkage completo: todos os pares de um
  cluster tem bitscore_pct >= limiar (conservados ENTRE SI).
  Limiar = vale da densidade (KDE) das distancias entre pares com alinhamento.
  Em cada cluster repete-se a busca de vale (recursivo) enquanto houver pares
  suficientes (--min-pairs) e um vale real (--valley-ratio). Sem vale no
  conjunto inicial, o corte e logo abaixo de 1 (so pares com alinhamento).

Gap e cobertura
  D_gc = sqrt((w_cov*(1-cov)^2 + w_gap*gap^2)/(w_cov+w_gap)); maior gap e menor
  cobertura = mais distante. Usada nos dendrogramas de clusters (distancia entre
  clusters) e nos dendrogramas internos de cada cluster, cada um em sua imagem.

Saidas (prefixo --out)
  _pares.tsv  _segmentos.tsv  _limiares.tsv  _clusters_resumo.tsv
  _clusters_distancias_gap_cobertura.tsv
  _picos_limiares.png  _dendrograma_bitscore.png
  _dendrograma_clusters_gap_cobertura.png  _cluster_<id>_dendrograma_gap_cobertura.png
  _heatmap_global.png  _heatmap_cluster_<id>.png   (R + ggplot2, eve_heatmap.R)
"""
import argparse
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# 1. XML -> métricas por par
# --------------------------------------------------------------------------
def wis(hsps, max_overlap_nt=0):
    """HSPs não sobrepostos (query) com máximo bitscore total (prog. dinâmica)."""
    if not hsps:
        return []
    s = sorted(hsps, key=lambda h: h["hi"])
    n = len(s)
    p = [-1] * n
    for i in range(n):
        for j in range(i - 1, -1, -1):
            if s[j]["hi"] - max_overlap_nt < s[i]["lo"]:
                p[i] = j
                break
    opt = [0.0] * (n + 1)
    for i in range(1, n + 1):
        incl = s[i - 1]["bitscore"] + (opt[p[i - 1] + 1] if p[i - 1] != -1 else 0.0)
        opt[i] = max(incl, opt[i - 1])
    out, cur = [], n
    while cur > 0:
        incl = s[cur - 1]["bitscore"] + (opt[p[cur - 1] + 1] if p[cur - 1] != -1 else 0.0)
        if incl >= opt[cur - 1]:
            out.append(s[cur - 1])
            cur = p[cur - 1] + 1 if p[cur - 1] != -1 else 0
        else:
            cur -= 1
    return out[::-1]


def _name(hid, hdef):
    hid = hid or ""
    if hid.startswith("lcl|"):
        return hid[4:]
    d = (hdef or "").split()
    return d[0] if d else hid


def parse_pairs(path, max_overlap_nt=0):
    print(f"[XML] lendo {path}...", file=sys.stderr)
    qlen, selfs, agg = {}, {}, {}
    for _, el in ET.iterparse(path):
        if el.tag != "Iteration":
            continue
        qd = el.findtext("Iteration_query-def", "") or el.findtext("Iteration_query-ID", "")
        q = (qd.split() or [""])[0]
        qlen[q] = int(el.findtext("Iteration_query-len", "0"))
        for hit in el.iter("Hit"):
            s = _name(hit.findtext("Hit_id", ""), hit.findtext("Hit_def", ""))
            hs = []
            for h in hit.iter("Hsp"):
                a, b = int(h.findtext("Hsp_query-from")), int(h.findtext("Hsp_query-to"))
                hs.append(dict(lo=min(a, b), hi=max(a, b),
                               bitscore=float(h.findtext("Hsp_bit-score", "0")),
                               ident=int(h.findtext("Hsp_identity", "0")),
                               aln=int(h.findtext("Hsp_align-len", "0")),
                               gaps=int(h.findtext("Hsp_gaps", "0"))))
            ch = wis(hs, max_overlap_nt)
            if not ch:
                continue
            bs = sum(h["bitscore"] for h in ch)
            if s == q:                                  # auto-alinhamento = máximo possível
                selfs[q] = max(selfs.get(q, 0.0), bs)
                continue
            cov = sum(h["hi"] - h["lo"] + 1 for h in ch)
            env = max(h["hi"] for h in ch) - min(h["lo"] for h in ch) + 1
            prev = agg.get((q, s))
            if prev is None or bs > prev["bs"]:
                agg[(q, s)] = dict(bs=bs, ident=sum(h["ident"] for h in ch),
                                   aln=sum(h["aln"] for h in ch), gaps=sum(h["gaps"] for h in ch),
                                   cov=cov, env=env)
        el.clear()
    names = sorted(set(qlen) | {s for _, s in agg})
    print(f"[XML] {len(names)} segmentos, {len(selfs)} com auto-alinhamento, "
          f"{len(agg)} pares dirigidos", file=sys.stderr)
    if len(selfs) < len(qlen):
        print(f"[aviso] {len(qlen) - len(selfs)} segmento(s) sem auto-alinhamento no XML: "
              f"pares sem nenhum auto-score disponivel sao descartados.", file=sys.stderr)

    best, dropped = {}, 0
    for (q, s), r in agg.items():
        den = max(selfs.get(q, 0.0), selfs.get(s, 0.0))
        if den <= 0:
            dropped += 1
            continue
        L = max(qlen.get(q, 0), qlen.get(s, 0)) or qlen[q]
        env = max(r["env"], 1)
        row = dict(bs=r["bs"],
                   bitscore_pct=min(100.0, 100.0 * r["bs"] / den),
                   identity_pct=min(100.0, 100.0 * r["ident"] / r["aln"]) if r["aln"] else np.nan,
                   coverage_pct=min(100.0, 100.0 * r["cov"] / L),
                   gap_pct=min(100.0, 100.0 * (max(env - r["cov"], 0) + 3 * r["gaps"]) / env))
        k = tuple(sorted((q, s)))
        if k not in best or row["bs"] > best[k]["bs"]:
            best[k] = row
    if dropped:
        print(f"[aviso] {dropped} par(es) dirigido(s) descartado(s) por falta de auto-score.", file=sys.stderr)
    df = pd.DataFrame([dict(seg_a=a, seg_b=b, **{k: round(v, 3) for k, v in r.items() if k != "bs"})
                       for (a, b), r in best.items()])
    return names, df


def build_matrices(names, df, w_cov, w_gap):
    n = len(names)
    ix = {x: i for i, x in enumerate(names)}
    BP = np.zeros((n, n))
    CV = np.zeros((n, n))
    GP = np.full((n, n), 100.0)
    for r in df.itertuples():
        i, j = ix[r.seg_a], ix[r.seg_b]
        BP[i, j] = BP[j, i] = r.bitscore_pct
        CV[i, j] = CV[j, i] = r.coverage_pct
        GP[i, j] = GP[j, i] = r.gap_pct
    D = 1.0 - BP / 100.0
    Dgc = np.sqrt((w_cov * (1 - CV / 100.0) ** 2 + w_gap * (GP / 100.0) ** 2) / (w_cov + w_gap))
    for M in (D, Dgc):
        np.fill_diagonal(M, 0.0)
    return D, Dgc, BP


# --------------------------------------------------------------------------
# 2. Vales entre picos
# --------------------------------------------------------------------------
def find_valley(vals, valley_ratio=0.5, min_pairs=15, prom=0.05, max_bw=0.05):
    """
    KDE das distancias; picos = maximos locais com proeminencia >= prom x pico maximo.
    Vale = minimo entre dois picos consecutivos; aceito se densidade_vale /
    min(densidade dos dois picos) <= valley_ratio e se cada lado tem pelo menos
    max(3, 5%) dos pares. Largura de banda = min(Scott, max_bw) (em unidades de
    distancia, 0.05 = 5 pontos de bitscore %): o Scott sozinho suaviza demais
    distribuicoes com varios grupos e esconde vales reais. Escolhe o mais fundo.
    Retorna dict (sempre com a curva) ou None se nao ha dados suficientes.
    """
    from scipy.signal import find_peaks
    from scipy.stats import gaussian_kde
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) < min_pairs or np.ptp(v) < 1e-4:
        return None
    if len(v) > 20000:
        v = np.random.default_rng(0).choice(v, 20000, replace=False)
    try:
        kde = gaussian_kde(v)
        sd = float(np.std(v, ddof=1))
        kde.set_bandwidth(min(kde.factor, max_bw / sd))
    except Exception:
        return None
    grid = np.linspace(-0.15, 1.15, 651)
    dens = kde(grid)
    pk, _ = find_peaks(dens, prominence=prom * dens.max())
    best = None
    for a, b in zip(pk[:-1], pk[1:]):
        m = a + int(np.argmin(dens[a:b + 1]))
        x = float(grid[m])
        if not 0.0 <= x <= 1.0:
            continue
        ratio = float(dens[m] / min(dens[a], dens[b]))
        side = max(3, int(0.05 * len(v)))
        if (v < x).sum() < side or (v >= x).sum() < side:
            continue
        if ratio <= valley_ratio and (best is None or ratio < best[1]):
            best = (x, ratio)
    return dict(values=v, grid=grid, dens=dens, peaks=grid[pk],
                thr=best[0] if best else None, ratio=best[1] if best else np.nan)


def split_recursive(D, idx, label, depth, a, splits):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    n = len(idx)
    npairs = n * (n - 1) // 2
    if n < 3 or depth >= a.max_depth or npairs < a.min_pairs:
        return [(label, idx)]
    sub = D[np.ix_(idx, idx)]
    vals = sub[np.triu_indices(n, 1)]
    aligned = vals[vals < 1 - 1e-9]
    info = find_valley(aligned, a.valley_ratio, a.min_pairs, max_bw=a.max_bandwidth)
    thr, why = (info["thr"], "vale") if info else (None, "")
    if thr is None and depth == 0:
        thr, why = 1.0 - 1e-6, "sem vale: so pares com alinhamento"
    if thr is None:
        return [(label, idx)]
    Z = linkage(squareform(sub, checks=False), method="complete")
    lab = fcluster(Z, t=thr, criterion="distance")
    groups = [np.where(lab == k)[0] for k in np.unique(lab)]
    if len(groups) == 1:
        return [(label, idx)]
    groups.sort(key=lambda g: -len(g))
    splits.append(dict(grupo=label or "raiz", n=n, limiar_dist=round(thr, 4),
                       limiar_bitscore_pct=round(100 * (1 - thr), 2), n_pares=len(aligned),
                       n_filhos=len(groups), criterio=why,
                       picos_bitscore_pct=";".join(f"{100 * (1 - p):.1f}" for p in
                                                   (info["peaks"] if info else [])
                                                   if 0 <= p <= 1),
                       _info=info))
    out = []
    for k, g in enumerate(groups, 1):
        out += split_recursive(D, idx[g], f"{label}.{k}" if label else str(k), depth + 1, a, splits)
    return out


# --------------------------------------------------------------------------
# 3. Figuras (matplotlib: picos e dendrogramas)
# --------------------------------------------------------------------------
def _colors(cluster_of):
    import matplotlib.pyplot as plt
    cm = plt.get_cmap("tab20")
    labs = [c for c in dict.fromkeys(cluster_of.values()) if c != "isolado"]
    dark_first = list(range(0, 20, 2)) + list(range(1, 20, 2))   # tons escuros antes dos claros
    cmap = {c: cm(dark_first[i % 20]) for i, c in enumerate(labs)}
    cmap["isolado"] = (0.55, 0.55, 0.55, 1.0)
    return cmap


def draw_dendro(Z, labels, cluster_of, path, title, xlabel, vlines=(), cmap=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    n = len(labels)
    fig, ax = plt.subplots(figsize=(9, max(4, min(60, 0.22 * n + 1.5))))
    dendrogram(Z, labels=labels, orientation="left", ax=ax, no_labels=n > 400,
               color_threshold=0, above_threshold_color="0.35", leaf_font_size=8)
    cmap = cmap or _colors(cluster_of)
    for t in ax.get_ymajorticklabels():
        c = cluster_of.get(t.get_text())
        if c is not None:
            t.set_color(cmap[c])
    for x, txt in vlines:
        ax.axvline(x, color="crimson", ls="--", lw=1)
        ax.text(x, ax.get_ylim()[1], f" {txt}", color="crimson", fontsize=8, va="bottom")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_valleys(splits, path, max_panels=12):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sp = [s for s in splits if s.get("_info")][:max_panels]
    if not sp:
        return
    cols = min(3, len(sp))
    rows = (len(sp) + cols - 1) // cols
    fig, axs = plt.subplots(rows, cols, figsize=(max(7.5, 4.6 * cols), 3.2 * rows + 0.4), squeeze=False)
    for ax, s in zip(axs.ravel(), sp):
        info = s["_info"]
        pct = 100 * (1 - info["grid"])
        ax.hist(100 * (1 - info["values"]), bins=40, density=True, color="0.8")
        ax.plot(pct, info["dens"] / 100.0, color="navy")
        for p in info["peaks"]:
            ax.axvline(100 * (1 - p), color="green", lw=0.8, ls=":")
        ax.axvline(s["limiar_bitscore_pct"], color="crimson", lw=1.5, ls="--")
        ax.set_xlim(0, 100)
        ax.set_title(f"grupo {s['grupo']} (n={s['n']}): limiar {s['limiar_bitscore_pct']}%", fontsize=9)
        ax.set_xlabel("bitscore % do par")
    for ax in axs.ravel()[len(sp):]:
        ax.axis("off")
    fig.suptitle("Pares com alinhamento: picos (verde) e vale usado como limiar (vermelho)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


# --------------------------------------------------------------------------
# 4. Pipeline
# --------------------------------------------------------------------------
def run_clustering(
    fmt5: str,
    out: str = "eve",
    max_overlap_nt: int = 0,
    valley_ratio: float = 0.5,
    max_bandwidth: float = 0.05,
    min_pairs: int = 15,
    max_depth: int = 6,
    w_cov: float = 1.0,
    w_gap: float = 1.0,
    linkage_gapcov: str = "average",
    max_cluster_images: int = 20,
    max_global_heatmap: int = 400,
    no_r: bool = False
):
    import sys
    import os
    import shutil
    import subprocess
    import numpy as np
    import pandas as pd
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform
    
    # Mocking args object to fit existing code easily
    class Args:
        pass
    a = Args()
    a.fmt5 = fmt5
    a.out = out
    a.max_overlap_nt = max_overlap_nt
    a.valley_ratio = valley_ratio
    a.max_bandwidth = max_bandwidth
    a.min_pairs = min_pairs
    a.max_depth = max_depth
    a.w_cov = w_cov
    a.w_gap = w_gap
    a.linkage_gapcov = linkage_gapcov
    a.max_cluster_images = max_cluster_images
    a.max_global_heatmap = max_global_heatmap
    a.no_r = no_r

    names, df = parse_pairs(a.fmt5, a.max_overlap_nt)
    if df.empty:
        print("[aviso] Nenhum par com alinhamento entre segmentos diferentes.", file=sys.stderr)
        return False
    names, df = parse_pairs(a.fmt5, a.max_overlap_nt)
    if df.empty:
        sys.exit("Nenhum par com alinhamento entre segmentos diferentes.")
    n = len(names)
    D, Dgc, BP = build_matrices(names, df, a.w_cov, a.w_gap)
    df.to_csv(f"{a.out}_pares.tsv", sep="\t", index=False)

    # ---- clusters por bitscore % com limiares de vale ----
    splits = []
    leaves = split_recursive(D, np.arange(n), "", 0, a, splits)
    leaves.sort(key=lambda x: -len(x[1]))
    cluster_of, rows_c, k = {}, [], 0
    for lab, idx in leaves:
        if len(idx) >= 2:
            k += 1
            lab = lab or "1"
            for i in idx:
                cluster_of[names[i]] = lab
            rows_c.append((lab, idx))
        else:
            cluster_of[names[idx[0]]] = "isolado"
    # ordem para o heatmap: clusters por tamanho, dentro de cada um pelo linkage de bitscore
    order = []
    for lab, idx in rows_c:
        if len(idx) >= 3:
            Zc = linkage(squareform(D[np.ix_(idx, idx)], checks=False), "complete")
            idx = idx[leaves_list(Zc)]
        order += list(idx)
    order += [i for i in range(n) if cluster_of[names[i]] == "isolado"]
    pos = {names[i]: r + 1 for r, i in enumerate(order)}
    size = pd.Series(cluster_of).value_counts()
    seg = pd.DataFrame(dict(segmento=names, cluster=[cluster_of[x] for x in names]))
    seg["n_cluster"] = seg["cluster"].map(lambda c: 1 if c == "isolado" else int(size[c]))
    seg["ordem"] = seg["segmento"].map(pos)
    seg.sort_values("ordem").to_csv(f"{a.out}_segmentos.tsv", sep="\t", index=False)

    lim = pd.DataFrame([{k_: v for k_, v in s.items() if k_ != "_info"} for s in splits])
    lim.to_csv(f"{a.out}_limiares.tsv", sep="\t", index=False)
    plot_valleys(splits, f"{a.out}_picos_limiares.png")

    # ---- resumo por cluster ----
    ix = {x: i for i, x in enumerate(names)}
    pk = df.assign(ca=df.seg_a.map(cluster_of), cb=df.seg_b.map(cluster_of))
    res = []
    for lab, idx in rows_c:
        sub = pk[(pk.ca == lab) & (pk.cb == lab)]
        m = len(idx)
        res.append(dict(cluster=lab, n=m, pares_possiveis=m * (m - 1) // 2, pares_com_alinhamento=len(sub),
                        bitscore_min=sub.bitscore_pct.min(), bitscore_medio=round(sub.bitscore_pct.mean(), 2),
                        identidade_media=round(sub.identity_pct.mean(), 2),
                        cobertura_media=round(sub.coverage_pct.mean(), 2),
                        gap_medio=round(sub.gap_pct.mean(), 2),
                        membros=";".join(names[i] for i in idx[:12]) + (";..." if m > 12 else "")))
    pd.DataFrame(res).to_csv(f"{a.out}_clusters_resumo.tsv", sep="\t", index=False)

    # ---- dendrogramas ----
    cmap = _colors(cluster_of)
    if n >= 3:
        Z = linkage(squareform(D, checks=False), "complete")
        vl = [(splits[0]["limiar_dist"], f"{splits[0]['limiar_bitscore_pct']}%")] if splits else []
        draw_dendro(Z, names, cluster_of, f"{a.out}_dendrograma_bitscore.png",
                    "Segmentos por bitscore % (linkage completo)", "1 - bitscore % do par", vl, cmap)

    labs = [lab for lab, _ in rows_c]
    if len(labs) >= 3:
        M = np.zeros((len(labs), len(labs)))
        for x, (la, ia) in enumerate(rows_c):
            for y, (lb, ib) in enumerate(rows_c):
                if x != y:
                    M[x, y] = Dgc[np.ix_(ia, ib)].mean()
        pd.DataFrame(M, index=labs, columns=labs).round(4).to_csv(
            f"{a.out}_clusters_distancias_gap_cobertura.tsv", sep="\t")
        Zk = linkage(squareform(M, checks=False), a.linkage_gapcov)
        lbl = [f"{lab} (n={len(i)})" for lab, i in rows_c]
        draw_dendro(Zk, lbl, {l: c for l, c in zip(lbl, labs)},
                    f"{a.out}_dendrograma_clusters_gap_cobertura.png",
                    "Distancia entre clusters por gap + cobertura (media dos pares)",
                    "distancia gap/cobertura", cmap=cmap)

    shown = 0
    for lab, idx in rows_c:
        if len(idx) < 3 or shown >= a.max_cluster_images:
            continue
        sub = Dgc[np.ix_(idx, idx)]
        Zc = linkage(squareform(sub, checks=False), a.linkage_gapcov)
        info = find_valley(sub[np.triu_indices(len(idx), 1)], a.valley_ratio, a.min_pairs, max_bw=a.max_bandwidth)
        vl = [(info["thr"], f"{info['thr']:.2f}")] if info and info["thr"] is not None else []
        nm = [names[i] for i in idx]
        draw_dendro(Zc, nm, {x: lab for x in nm}, f"{a.out}_cluster_{lab}_dendrograma_gap_cobertura.png",
                    f"Cluster {lab} (n={len(idx)}): gap + cobertura", "distancia gap/cobertura", vl, cmap)
        shown += 1

    # ---- heatmap (R / ggplot2) ----
    r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eve_heatmap.R")
    cmd = ["Rscript", r_script, a.out, str(a.max_cluster_images), str(a.max_global_heatmap)]
    if a.no_r:
        print("[heatmap] pulado (--no-r). Para gerar: " + " ".join(cmd), file=sys.stderr)
    elif shutil.which("Rscript") is None:
        print("[aviso] Rscript nao encontrado. Instale R + ggplot2 e rode: " + " ".join(cmd), file=sys.stderr)
    else:
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stderr.strip(), file=sys.stderr)
        if r.returncode:
            print("[aviso] o R falhou (veja acima); as tabelas e dendrogramas ja foram gerados.", file=sys.stderr)

    # ---- resumo no terminal ----
    print(f"\n[clusters] {n} segmentos | {len(rows_c)} clusters (>=2) | "
          f"{n - sum(len(i) for _, i in rows_c)} isolados | {len(splits)} separacoes por vale", file=sys.stderr)
    for s in splits[:10]:
        print(f"   grupo {s['grupo']:>6} n={s['n']:<5} limiar bitscore {s['limiar_bitscore_pct']:>6}% "
              f"-> {s['n_filhos']} grupos ({s['criterio']})", file=sys.stderr)
    if not splits:
        print("   nenhuma separacao: dados insuficientes ou sem vale entre picos "
              "(veja --min-pairs / --valley-ratio)", file=sys.stderr)


if __name__ == "__main__":
    main()
