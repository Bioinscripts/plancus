#!/usr/bin/env python3
"""
Host-genome figures of the curated EVEs: how many EVEs each genome holds, next to the size
and quality of the genome and to each filtering step of the pipeline.

Only species whose assembly is in NCBI are used by default (the own SPAdes assemblies are
much more fragmented, so their counts are not comparable); edit the "use" column of
13_curadoria/figuras/genomes.tsv to change that.

Figures (PDF, SVG, 300 dpi PNG, plus the table behind each):
  fig5_genome_features       per species: assembly size, predicted ORFs, ORFs kept by the
                             viral DIAMOND screen, approved EVE loci (single / paralogous),
                             in the layout of Fig. 1A of Conceição et al. 2026, Mobile DNA 17:8
  fig6_genome_correlations   approved EVE loci against assembly size, contig N50 and host
                             group, with Spearman's rho
  fig7_chromosome_map        position of every approved EVE locus on the chromosomes of the
                             chromosome-level assemblies

genomes.tsv is filled automatically: the assembly accession comes from the genome file name
or, for files named otherwise, from NCBI Datasets (one query per species, looked up by the
accession of the first sequence); assembly level, contig and scaffold N50 come from the
NCBI assembly report. Without internet the table keeps what it has and the figures use the
local .fai.

Usage:
  python -m evetool.curation.genome_figures -o results_eve
  eve-tool genome-figures -o results_eve --formats pdf,png
"""
import os
import re
import csv
import json
import argparse
import urllib.request
from collections import Counter, OrderedDict, defaultdict

from evetool.curation import figures as F

ARQ_GENOMAS = "genomes.tsv"
ARQ_CACHE = ".genome_stats_cache.json"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
ACC_MONTAGEM = re.compile(r"(GC[AF])_?(\d{9})\.?(\d+)")
ACC_SEQUENCIA = re.compile(r"^[A-Z]{1,6}_?\d{5,}\.\d+$")
COLS_GENOMAS = ["species", "use", "source", "accession", "assembly_level", "ncbi_length",
                "contig_n50", "scaffold_n50", "gc_percent", "organism"]


def log(msg):
    print(f"[genome-figures] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Assembly metadata
# ---------------------------------------------------------------------------
def _json_url(url, tempo=30):
    req = urllib.request.Request(url, headers={"User-Agent": "eve-tool genome-figures",
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=tempo) as r:
        return json.loads(r.read().decode("utf-8"))


def _primeira_sequencia(fai):
    try:
        with open(fai, encoding="utf-8", errors="replace") as fh:
            return fh.readline().split("\t", 1)[0].strip()
    except OSError:
        return ""


def _acesso_do_arquivo(genoma):
    m = ACC_MONTAGEM.search(os.path.basename(genoma))
    return f"{m.group(1)}_{m.group(2)}.{m.group(3)}" if m else ""


def _ler_genomas(caminho):
    linhas = OrderedDict()
    if os.path.exists(caminho):
        with open(caminho, encoding="utf-8") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                if r.get("species"):
                    linhas[r["species"].strip()] = {k: (r.get(k) or "").strip() for k in COLS_GENOMAS}
    return linhas


def metadados(outdir, especies, pasta_fig, online=True):
    """genomes.tsv: one row per species, filled from the file name, the sequence names and
    NCBI; values already in the file (also hand edits) are kept."""
    from evetool.curation import curation as cur
    caminho = os.path.join(pasta_fig, ARQ_GENOMAS)
    tab = _ler_genomas(caminho)
    pastas = cur.pastas_especies(outdir)
    mudou = not os.path.exists(caminho)
    for s in especies:
        r = tab.setdefault(s, {k: "" for k in COLS_GENOMAS})
        r["species"] = s
        genoma = cur._info_pipeline(pastas[s]).get("genome") or "" if s in pastas else ""
        primeira = _primeira_sequencia(genoma + ".fai")
        if not r["source"]:
            ncbi = bool(_acesso_do_arquivo(genoma)) or bool(ACC_SEQUENCIA.match(primeira))
            r["source"] = "NCBI" if ncbi else "own assembly"
            mudou = True
        if not r["use"]:
            r["use"] = "yes" if r["source"] == "NCBI" else "no"
            mudou = True
        if r["source"] != "NCBI" or not online:
            continue
        try:
            if not r["accession"]:
                acc = _acesso_do_arquivo(genoma)
                if not acc and primeira:
                    j = _json_url(f"{NCBI}/sequence_accession/{primeira}/sequence_assemblies")
                    acc = (j.get("accessions") or [""])[0]
                r["accession"] = acc
                mudou = True
            if r["accession"] and not r["assembly_level"]:
                j = _json_url(f"{NCBI}/accession/{r['accession']}/dataset_report")
                rep = (j.get("reports") or [{}])[0]
                st = rep.get("assembly_stats") or {}
                # an empty report = accession replaced by a newer version: leave it to the local files
                r.update({"assembly_level": (rep.get("assembly_info") or {}).get("assembly_level", "")
                          or "not in NCBI (replaced?)",
                          "ncbi_length": str(st.get("total_sequence_length", "")),
                          "contig_n50": str(st.get("contig_n50", "")),
                          "scaffold_n50": str(st.get("scaffold_n50", "")),
                          "gc_percent": str(st.get("gc_percent", "")),
                          "organism": (rep.get("organism") or {}).get("organism_name", "")})
                mudou = True
                log(f"{s}: {r['accession']} ({r['assembly_level']})")
        except Exception as erro:   # offline or NCBI down: keep what we have
            log(f"{s}: NCBI lookup failed ({erro}); using the local files only")
            online = False
    if mudou:
        with open(caminho, "w", encoding="utf-8", newline="") as fh:
            fh.write("# Filled automatically (file name, sequence names, NCBI Datasets). 'use' = yes puts the\n"
                     "# species in fig5-7; edit any value and run again, your edits are kept.\n")
            w = csv.DictWriter(fh, fieldnames=COLS_GENOMAS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for s in especies:
                w.writerow(tab[s])
    return tab


# ---------------------------------------------------------------------------
# Pipeline numbers per species
# ---------------------------------------------------------------------------
def _fai(fai):
    """[(name, length), ...] in file order."""
    out = []
    with open(fai, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            p = linha.split("\t", 2)
            if len(p) >= 2:
                try:
                    out.append((p[0], int(p[1])))
                except ValueError:
                    pass
    return out


def _n50(tamanhos):
    t = sorted(tamanhos, reverse=True)
    meio, soma = sum(t) / 2, 0
    for x in t:
        soma += x
        if soma >= meio:
            return x
    return 0


def _conta_fasta(caminho):
    n = 0
    with open(caminho, "rb") as fh:
        while True:
            bloco = fh.read(1 << 24)
            if not bloco:
                break
            n += bloco.count(b">")
    return n


CROMO = re.compile(r"chromosome:?\s+([\w.]+)")
N_GAP = re.compile(rb"[Nn]{10,}")


def _varre_genoma(genoma):
    """One pass over the FASTA: chromosome names (from the headers), contig N50 (scaffolds
    split at runs of >= 10 N, the same rule for every assembly) and GC %."""
    cromos, contigs, gc, acgt = {}, [], 0, 0

    def fecha(partes):
        nonlocal gc, acgt
        if not partes:
            return
        seq = b"".join(partes).upper()
        gc += seq.count(b"G") + seq.count(b"C")
        acgt += seq.count(b"A") + seq.count(b"T") + seq.count(b"G") + seq.count(b"C")
        contigs.extend(len(x) for x in N_GAP.split(seq) if x)

    partes = []
    with open(genoma, "rb") as fh:
        for linha in fh:
            if linha[:1] == b">":
                fecha(partes)
                partes = []
                t = linha[1:].decode("utf-8", "replace").strip()
                nome, _, desc = t.partition(" ")
                m = CROMO.search(desc)
                if m and not re.search(r"unlocalized|unplaced|organelle|mitochondri", desc):
                    cromos[nome] = m.group(1).rstrip(",")
            else:
                partes.append(linha.rstrip())
        fecha(partes)
    return {"chromosomes": cromos, "contig_n50": _n50(contigs), "contigs": len(contigs),
            "gc_percent": round(100 * gc / acgt, 2) if acgt else None}


def _arquivo_orfs(pasta):
    d = os.path.join(pasta, "01_orfs_cds")
    if not os.path.isdir(d):
        return ""
    for f in sorted(os.listdir(d)):
        if f.endswith("_orfs_pep.fasta") and "_nr_" not in f:
            return os.path.join(d, f)
    return ""


def _conta_mantidas(pasta):
    """ORFs kept by the viral DIAMOND screen (best hit to an RNA/ssDNA virus)."""
    arq = os.path.join(pasta, "02_alinhamento_diamond", "viral_screening_decision_per_orf.tsv")
    n = 0
    if os.path.exists(arq):
        with open(arq, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n += r.get("decision") == "kept"
    return n


def _segmentos(pasta):
    """Hierarchical tblastx clusters of the EVE regions of one genome:
    [(contig, start, end, cluster or None), ...]; cluster None = no paralogue."""
    arq = os.path.join(pasta, "08_clusters_e_tensores", "eve_clusters_hierarch_segmentos.tsv")
    out = []
    if os.path.exists(arq):
        with open(arq, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                reg = r.get("segmento", "").split("::", 1)[-1]
                contig, _, coords = reg.rpartition(":")
                a, _, b = coords.partition("-")
                try:
                    c = r.get("cluster", "")
                    out.append((contig, int(a), int(b), None if c in ("", "isolado") else c))
                except ValueError:
                    pass
    return out


def _coords(locus):
    contig, _, coords = locus.rpartition(":")
    a, _, b = coords.partition("-")
    try:
        return contig, int(a), int(b)
    except ValueError:
        return locus, 0, 0


def estatisticas(outdir, d, tab, pasta_fig):
    """One dict per species with the numbers of fig5/fig6; slow parts cached by file stamp."""
    from evetool.curation import curation as cur
    arq_cache = os.path.join(pasta_fig, ARQ_CACHE)
    try:
        with open(arq_cache, encoding="utf-8") as fh:
            cache = json.load(fh)
    except Exception:
        cache = {}
    pastas = cur.pastas_especies(outdir)
    loci_sp = defaultdict(list)
    for l in d["loci"]:
        loci_sp[l["species"]].append(l)
    est = OrderedDict()
    for s in d["species"]:
        pasta = pastas.get(s, "")
        genoma = cur._info_pipeline(pasta).get("genome") or ""
        fai = _fai(genoma + ".fai") if os.path.exists(genoma + ".fai") else []
        tam = [x for _, x in fai]
        orfs_arq = _arquivo_orfs(pasta)
        chave = {"orfs": orfs_arq, "genome": genoma}
        selo = {k: (os.path.getsize(v), int(os.path.getmtime(v))) if v and os.path.exists(v) else None
                for k, v in chave.items()}
        c = cache.get(s, {})
        if c.get("selo") != json.loads(json.dumps(selo)):
            log(f"{s}: counting ORFs and scanning the genome (once; cached)")
            c = {"selo": selo, "orfs": _conta_fasta(orfs_arq) if orfs_arq else 0}
            c.update(_varre_genoma(genoma) if os.path.exists(genoma) else {})
            cache[s] = c
        # paralogous = the EVE region of the locus shares a tblastx cluster with another region
        # of the same genome that also holds an approved locus
        segs = _segmentos(pasta)
        def seg_de(locus):
            ct, a, b = _coords(locus)
            for k, (sc, sa, sb, cl) in enumerate(segs):
                if sc == ct and a < sb and b > sa:
                    return k
            return None
        seg_locus = {l["locus"]: seg_de(l["locus"]) for l in loci_sp[s]}
        regioes_por_cluster = defaultdict(set)
        for k in seg_locus.values():
            if k is not None and segs[k][3] is not None:
                regioes_por_cluster[segs[k][3]].add(k)
        paralogos = sum(1 for k in seg_locus.values()
                        if k is not None and segs[k][3] is not None and len(regioes_por_cluster[segs[k][3]]) > 1)
        regioes = len({k for k in seg_locus.values() if k is not None}) + \
            sum(1 for k in seg_locus.values() if k is None)
        r = tab.get(s, {})
        n_loci = len(loci_sp[s])
        comp = sum(tam)
        est[s] = {
            "species": s, "host_group": d["host_group"].get(s, ""), "use": r.get("use") == "yes",
            "source": r.get("source", ""), "accession": r.get("accession", ""),
            "assembly_mb": comp / 1e6, "sequences": len(tam), "n50_local": _n50(tam),
            "assembly_level": (r.get("assembly_level") if not (r.get("assembly_level") or "n").startswith("n")
                               else ("Chromosome" if c.get("chromosomes") else
                                     "Scaffold" if (c.get("contigs") or 0) > len(tam) else "Contig")),
            "contig_n50": c.get("contig_n50"), "contigs": c.get("contigs"),
            "gc_percent": c.get("gc_percent"),
            "orfs": c.get("orfs", 0), "viral_orfs_kept": _conta_mantidas(pasta),
            "pipeline_candidates": _conta_fasta(os.path.join(pasta, "eve_candidatas_unificadas_pep.fasta"))
            if os.path.exists(os.path.join(pasta, "eve_candidatas_unificadas_pep.fasta")) else 0,
            "eve_loci": n_loci, "paralogous": paralogos, "single": n_loci - paralogos,
            "eve_regions": regioes,
            "loci_per_100mb": n_loci / (comp / 1e8) if comp else 0.0,
            "chromosomes": c.get("chromosomes", {}), "lengths": dict(fai),
        }
    with open(arq_cache, "w", encoding="utf-8") as fh:
        json.dump(cache, fh)
    return est


def tabela(est, pasta_fig):
    cols = ["species", "host_group", "use", "source", "accession", "assembly_level", "assembly_mb", "sequences",
            "contigs", "n50_local", "contig_n50", "gc_percent", "orfs", "viral_orfs_kept",
            "pipeline_candidates", "eve_loci", "single", "paralogous", "eve_regions", "loci_per_100mb"]
    F._gravar_tsv(os.path.join(pasta_fig, "fig5_genome_features.tsv"), cols,
                  [[round(e[k], 2) if isinstance(e[k], float) else ("" if e[k] is None else e[k]) for k in cols]
                   for e in est.values()])


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
MM = F.MM
COR_GENOMA, COR_PIPE = "#52514e", "#898781"
COR_UNICO, COR_PARALOGO = "#1c5cab", "#8dbbf1"
MARCADOR = {"Phlebotominae (New World)": ("D", True), "Phlebotominae (Old World)": ("o", True),
            "Psychodinae": ("^", False), "Other Diptera": ("o", False)}
NIVEL_CURTO = {"Chromosome": "chromosome", "Scaffold": "scaffold", "Contig": "contig",
               "Complete Genome": "complete"}


def _marcador(grupo):
    return MARCADOR.get(grupo, ("s", False))


def _abrev(s):
    """Lutzomyia_longipalpis -> L. lon ; Trichoceridae_sp -> Trichoceridae sp."""
    p = s.split("_")
    if len(p) < 2 or p[1] == "sp":
        return F._it(p[0]) + " sp."
    return F._it(f"{p[0][0]}. {p[1][:3]}")


def _usadas(est):
    return [e for e in est.values() if e["use"]]


def fig5(est, plt, pasta, formatos):
    """Per species: assembly size, predicted ORFs, ORFs kept by the viral screen, approved
    EVE loci (single / paralogous). Layout of Fig. 1A of Conceicao et al. 2026."""
    from matplotlib.patches import Rectangle, Patch
    usadas = _usadas(est)
    grupos = {e["species"]: e["host_group"] for e in usadas}
    linhas = F._linhas_com_grupos([e["species"] for e in usadas], grupos)
    ys, y = [], 0.0
    for tipo, _ in linhas:
        y -= 0.85 if tipo == "group" else 1.0
        ys.append(y)
    passo = 5.4 * MM
    alto = -ys[-1] * passo + 2 * MM
    W = 180 * MM
    altura = alto + 30 * MM
    fig = plt.figure(figsize=(W, altura))
    esq, larg, vao = 58, 25.5, 6.5
    paineis = [
        ("assembly_mb", "Assembly size (Mb)", "Genome", COR_GENOMA, False),
        ("orfs", "Predicted ORFs (millions)", "ORF prediction", COR_GENOMA, False),
        ("viral_orfs_kept", "ORFs with a viral hit", "Viral screen", COR_PIPE, True),
        ("eve_loci", "Approved EVE loci", "Manual curation", None, True),
    ]
    por = {e["species"]: e for e in usadas}
    eixos = []
    for k, (chave, xlabel, etapa, cor, rotulos) in enumerate(paineis):
        x0 = esq + k * (larg + vao)
        ax = fig.add_axes([x0 * MM / W, 12 * MM / altura, larg * MM / W, alto / altura])
        eixos.append(ax)
        escala = 1e6 if chave == "orfs" else 1
        vmax = max(e[chave] for e in usadas) / escala
        ax.set_xlim(0, vmax * (1.28 if rotulos else 1.05))
        ax.set_ylim(ys[-1] - 0.6, 0)
        for (tipo, nome), yy in zip(linhas, ys):
            if tipo != "species":
                continue
            e = por[nome]
            v = e[chave] / escala
            if chave == "eve_loci":
                ax.add_patch(Rectangle((0, yy - 0.34), e["single"], 0.68, facecolor=COR_UNICO,
                                       edgecolor="white", linewidth=0.7))
                ax.add_patch(Rectangle((e["single"], yy - 0.34), e["paralogous"], 0.68,
                                       facecolor=COR_PARALOGO, edgecolor="white", linewidth=0.7))
            else:
                ax.add_patch(Rectangle((0, yy - 0.34), v, 0.68, facecolor=cor, edgecolor="white", linewidth=0.7))
            if rotulos:
                ax.text(v + vmax * 0.03, yy, f"{e[chave]:,}", va="center", ha="left", fontsize=6.2, color=F.TINTA)
        ax.set_yticks([])
        for lado in ("top", "right", "left"):
            ax.spines[lado].set_visible(False)
        ax.xaxis.grid(True, color=F.GRADE, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_xlabel(xlabel, fontsize=6.5)
        ax.locator_params(axis="x", nbins=4)
        topo_ax = (12 * MM + alto) / altura
        fig.text((x0 + 0.5) * MM / W, topo_ax + 2.0 * MM / altura, f"{k + 1}  {etapa}", ha="left", va="bottom",
                 fontsize=6.8, color=F.TINTA, fontweight="bold")
        if k < 3:
            fig.text((x0 + larg + vao / 2) * MM / W, topo_ax + 2.6 * MM / altura, "→", ha="center",
                     va="bottom", fontsize=8, color=F.MUDO)
    ax0 = eixos[0]
    inv = fig.transFigure.inverted()
    for (tipo, nome), yy in zip(linhas, ys):
        fy = inv.transform(ax0.transData.transform((0, yy)))[1]
        if tipo == "group":
            fig.text(2 * MM / W, fy, nome, ha="left", va="center", fontsize=6.5, color=F.TINTA, fontweight="bold")
            continue
        e = por[nome]
        fig.text(4 * MM / W, fy, F.rotulo_especie(nome), ha="left", va="center", fontsize=6.5, color=F.TINTA)
        fig.text((esq - 2.5) * MM / W, fy, NIVEL_CURTO.get(e["assembly_level"], e["assembly_level"].lower()),
                 ha="right", va="center", fontsize=5.8, color=F.MUDO)
    topo = (12 * MM + alto) / altura
    fig.text((esq - 2.5) * MM / W, topo + 2.0 * MM / altura, "Assembly level", ha="right", va="bottom",
             fontsize=5.8, color=F.MUDO)
    alcas = [Patch(facecolor=COR_UNICO, edgecolor="none", label="single"),
             Patch(facecolor=COR_PARALOGO, edgecolor="none", label="paralogous")]
    fig.legend(handles=alcas, loc="lower left", ncol=2, frameon=False, handlelength=1.1, handleheight=0.9,
               columnspacing=1.0, borderaxespad=0, fontsize=6.2,
               bbox_to_anchor=((esq + 3 * (larg + vao)) * MM / W, topo + 6.5 * MM / altura))
    return F._salvar(fig, pasta, "fig5_genome_features", formatos)


def _rotula(ax, pontos, rotulos, fig, obstaculos=()):
    """Greedy label placement: each label goes to the first of 8 spots around its point
    that stays inside the axes and does not hit a point or an earlier label."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    cx = ax.get_window_extent(rend)
    disp = [ax.transData.transform(p) for p in pontos]
    r = 5 * fig.dpi / 72 / 2 + 1
    ocupado = [(x - r, y - r, x + r, y + r) for x, y in disp]
    for (xa, ya), (xb, yb) in obstaculos:      # segments in data coordinates (median ticks)
        (pa, qa), (pb, qb) = ax.transData.transform((xa, ya)), ax.transData.transform((xb, yb))
        ocupado.append((min(pa, pb), min(qa, qb) - 2, max(pa, pb), max(qa, qb) + 2))
    desloc = [(5, 0, "left", "center"), (-5, 0, "right", "center"), (0, 4.5, "center", "bottom"),
              (0, -4.5, "center", "top"), (4, 4, "left", "bottom"), (-4, 4, "right", "bottom"),
              (4, -4, "left", "top"), (-4, -4, "right", "top")]
    ordem = sorted(range(len(pontos)), key=lambda i: -disp[i][1])
    for i in ordem:
        melhor = None
        for dx, dy, ha, va in desloc:
            t = ax.annotate(rotulos[i], pontos[i], xycoords="data", xytext=(dx, dy), textcoords="offset points",
                            ha=ha, va=va, fontsize=5.6, color=F.TINTA_2, annotation_clip=False)
            bb = t.get_window_extent(rend)
            caixa = (bb.x0, bb.y0, bb.x1, bb.y1)
            dentro = bb.x0 >= cx.x0 and bb.x1 <= cx.x1 and bb.y0 >= cx.y0 and bb.y1 <= cx.y1
            sobrepoe = sum(1 for o in ocupado if not (caixa[2] < o[0] or caixa[0] > o[2]
                                                      or caixa[3] < o[1] or caixa[1] > o[3]))
            custo = sobrepoe * 10 + (0 if dentro else 5)
            if melhor is None or custo < melhor[0]:
                if melhor is not None:
                    melhor[1].remove()
                melhor = (custo, t, caixa)
            else:
                t.remove()
            if custo == 0:
                break
        ocupado.append(melhor[2])


def _spearman(x, y):
    from scipy import stats
    r = stats.spearmanr(x, y)
    return float(r.statistic if hasattr(r, "statistic") else r[0]), float(r.pvalue if hasattr(r, "pvalue") else r[1])


def _p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.3f}"


def correlacoes(est):
    """Spearman's rho between the genome measures and the EVE counts (species in use)."""
    usadas = _usadas(est)
    pares = [("assembly_mb", "eve_loci"), ("contig_n50", "eve_loci"), ("gc_percent", "eve_loci"),
             ("orfs", "eve_loci"), ("viral_orfs_kept", "eve_loci"), ("assembly_mb", "viral_orfs_kept"),
             ("assembly_mb", "loci_per_100mb"), ("contig_n50", "loci_per_100mb"), ("gc_percent", "loci_per_100mb")]
    out = []
    for a, b in pares:
        xs = [e[a] for e in usadas if e[a] is not None and e[b] is not None]
        ys = [e[b] for e in usadas if e[a] is not None and e[b] is not None]
        if len(xs) >= 4:
            rho, p = _spearman(xs, ys)
            out.append((a, b, len(xs), rho, p))
    return out


def fig6(est, plt, pasta, formatos):
    """Approved EVE loci against genome measures, one dot per species (NCBI assemblies)."""
    from matplotlib.lines import Line2D
    from scipy import stats
    usadas = _usadas(est)
    W, H = 180 * MM, 62 * MM
    fig = plt.figure(figsize=(W, H))
    larg, vao, esq = 34, 11.3, 10
    alt = 38
    paineis = [
        ("assembly_mb", "eve_loci", "Assembly size (Mb)", "Approved EVE loci", False),
        ("contig_n50", "eve_loci", "Contig N50 (Mb, log scale)", "Approved EVE loci", True),
        ("viral_orfs_kept", "eve_loci", "ORFs with a viral hit", "Approved EVE loci", False),
        ("grupo", "loci_per_100mb", "", "Approved EVE loci per 100 Mb", False),
    ]
    for k, (cx, cy, xl, yl, logx) in enumerate(paineis):
        ax = fig.add_axes([(esq + k * (larg + vao)) * MM / W, 12 * MM / H, larg * MM / W, alt * MM / H])
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        ax.grid(True, color=F.GRADE, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="both", length=2.5, labelsize=6)
        ax.set_ylabel(yl, fontsize=6.5)
        ax.text(-0.2, 1.06, "abcd"[k], transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        if cx == "grupo":
            fl = [e for e in usadas if "Phlebotominae" in e["host_group"]]
            ou = [e for e in usadas if "Phlebotominae" not in e["host_group"]]
            pts, rots, medianas = [], [], []
            for i, conj in enumerate((fl, ou)):
                vals = sorted(e[cy] for e in conj)
                for j, e in enumerate(sorted(conj, key=lambda e: e[cy])):
                    xj = i + (-0.16, 0.0, 0.16)[j % 3]
                    m, cheio = _marcador(e["host_group"])
                    ax.scatter([xj], [e[cy]], s=16, marker=m, facecolor=F.TINTA_2 if cheio else "white",
                               edgecolor=F.TINTA_2, linewidth=0.8, zorder=3)
                    pts.append((xj, e[cy]))
                    rots.append(_abrev(e["species"]))
                med = F._mediana(vals)
                ax.plot([i - 0.3, i + 0.3], [med, med], color=F.TINTA, linewidth=1.2, zorder=4, solid_capstyle="butt")
                medianas.append(((i - 0.3, med), (i + 0.3, med)))
            ax.set_xlim(-0.6, 1.9)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"Phlebotominae\n(n = {len(fl)})", f"Other\nDiptera\n(n = {len(ou)})"], fontsize=6)
            ax.tick_params(axis="x", length=0)
            ax.set_ylim(0, max(e[cy] for e in usadas) * 1.15)
            u = stats.mannwhitneyu([e[cy] for e in fl], [e[cy] for e in ou])
            ax.text(0.03, 0.97, f"Mann–Whitney {_p(float(u.pvalue))}", transform=ax.transAxes, ha="left", va="top",
                    fontsize=6, color=F.TINTA)
            _rotula(ax, pts, rots, fig, medianas)
            continue
        xs = [e[cx] / (1e6 if cx == "contig_n50" else 1) for e in usadas]
        ys = [e[cy] for e in usadas]
        for e, x, y in zip(usadas, xs, ys):
            m, cheio = _marcador(e["host_group"])
            ax.scatter([x], [y], s=16, marker=m, facecolor=F.TINTA_2 if cheio else "white", edgecolor=F.TINTA_2,
                       linewidth=0.8, zorder=3)
        if logx:
            from matplotlib.ticker import FixedLocator, NullLocator
            ax.set_xscale("log")
            ax.set_xlim(min(xs) / 2.5, max(xs) * 2.5)
            marcas = [m for m in (0.001, 0.01, 0.1, 1, 10, 100) if ax.get_xlim()[0] <= m <= ax.get_xlim()[1]]
            ax.xaxis.set_major_locator(FixedLocator(marcas))
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xticklabels([f"{m:g}" for m in marcas])
        else:
            ax.set_xlim(0, max(xs) * 1.12)
        ax.set_ylim(0, max(ys) * 1.15)
        ax.set_xlabel(xl, fontsize=6.5)
        rho, p = _spearman(xs, ys)
        ax.text(0.03, 0.97, f"ρ = {rho:.2f}, {_p(p)}", transform=ax.transAxes, ha="left", va="top",
                fontsize=6, color=F.TINTA)
        _rotula(ax, list(zip(xs, ys)), [_abrev(e["species"]) for e in usadas], fig)
    alcas = []
    for g in [g for g in MARCADOR if any(e["host_group"] == g for e in usadas)]:
        m, cheio = MARCADOR[g]
        alcas.append(Line2D([], [], linestyle="none", marker=m, markersize=4, markerfacecolor=F.TINTA_2 if cheio
                            else "white", markeredgecolor=F.TINTA_2, markeredgewidth=0.8, label=g))
    fig.legend(handles=alcas, loc="lower left", ncol=len(alcas), frameon=False, borderaxespad=0, fontsize=6.2,
               handletextpad=0.4, columnspacing=1.4, bbox_to_anchor=(esq * MM / W, (12 + alt + 8) * MM / H))
    fig.text(esq * MM / W, 0.5 * MM / H, f"One point per host species with an NCBI assembly (n = {len(usadas)}); "
             "ρ = Spearman's rank correlation.", ha="left", va="bottom", fontsize=5.8, color=F.MUDO)
    return F._salvar(fig, pasta, "fig6_genome_correlations", formatos)


def fig7(est, d, plt, pasta, formatos):
    """Approved EVE loci on the chromosomes of the chromosome-level assemblies."""
    from matplotlib.patches import FancyBboxPatch, Patch
    from matplotlib.lines import Line2D
    usadas = [e for e in _usadas(est) if e["chromosomes"]]
    if not usadas:
        return []
    loci_sp = defaultdict(list)
    for l in d["loci"]:
        loci_sp[l["species"]].append(l)

    def ordem_cromo(nome):
        try:
            return (0, float(nome))
        except ValueError:
            return (1, nome)

    linhas = []
    for e in usadas:
        cromos = sorted(e["chromosomes"].items(), key=lambda kv: ordem_cromo(kv[1]))
        linhas.append((e, [(seq, rot, e["lengths"].get(seq, 0)) for seq, rot in cromos]))
    vao_mb = 6.0
    maior = max(sum(c[2] for c in cs) / 1e6 + vao_mb * (len(cs) - 1) for _, cs in linhas)
    esq, larg = 40, 128
    passo = 14.0
    W = 180 * MM
    H = (len(linhas) * passo + 22) * MM
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([esq * MM / W, 12 * MM / H, larg * MM / W, len(linhas) * passo * MM / H])
    ax.set_xlim(0, maior)
    ax.set_ylim(-len(linhas) + 0.15, 0.78)
    ax.axis("off")
    mb_por_mm = maior / larg
    tipos_usados = set()
    for i, (e, cromos) in enumerate(linhas):
        y = -i
        fy = (12 * MM + (y - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0]) * len(linhas) * passo * MM) / H
        fig.text(2 * MM / W, fy, F.rotulo_especie(e["species"]), ha="left", va="center", fontsize=6.5, color=F.TINTA)
        x0, pos = 0.0, {}
        for seq, rot, comp in cromos:
            c = comp / 1e6
            ax.add_patch(FancyBboxPatch((x0, y - 0.09), c, 0.18, boxstyle="round,pad=0,rounding_size=0.09",
                                        mutation_aspect=1 / (mb_por_mm * 2.2), facecolor="#efeee9",
                                        edgecolor=F.BASE, linewidth=0.6))
            ax.text(x0 + c / 2, y - 0.13, rot, ha="center", va="top", fontsize=5.6, color=F.MUDO)
            pos[seq] = x0
            x0 += c + vao_mb
        fora = 0
        marcas = []
        for l in loci_sp[e["species"]]:
            ct, a, b = _coords(l["locus"])
            if ct not in pos:
                fora += 1
                continue
            marcas.append((pos[ct] + (a + b) / 2e6, l["genome_type"]))
        marcas.sort()
        grupos = []             # loci closer than 1.6 mm on the page form one stack
        for x, t in marcas:
            if grupos and x - grupos[-1][-1][0] < 1.6 * mb_por_mm:
                grupos[-1].append((x, t))
            else:
                grupos.append([(x, t)])
        for g in grupos:
            xm = sum(x for x, _ in g) / len(g)
            tipos = [t for t in F.TIPOS if any(tt == t for _, tt in g)]
            for nivel, t in enumerate(tipos):
                tipos_usados.add(t)
                ax.scatter([xm], [y + 0.17 + nivel * 0.1], s=9, marker="v", color=F.COR_TIPO.get(t, F.MUDO),
                           edgecolor="white", linewidth=0.4, zorder=3)
            if len(g) > 1:
                ax.text(xm, y + 0.17 + len(tipos) * 0.1 - 0.04, str(len(g)), ha="center", va="bottom",
                        fontsize=5.4, color=F.TINTA_2)
        tot = len(loci_sp[e["species"]])
        nota = f"{tot - fora} of {tot} loci on chromosomes" + (f"; {fora} on unplaced scaffolds" if fora else "")
        ax.text(x0 - vao_mb + 1.5 * mb_por_mm, y, nota, ha="left", va="center", fontsize=5.6, color=F.TINTA_2)
    # scale bar
    yb = -len(linhas) + 0.2
    ax.plot([0, 50], [yb, yb], color=F.TINTA_2, linewidth=0.8, solid_capstyle="butt")
    ax.text(25, yb - 0.08, "50 Mb", ha="center", va="top", fontsize=5.8, color=F.TINTA_2)
    alcas = [Line2D([], [], linestyle="none", marker="v", markersize=4, color=F.COR_TIPO[t], label=t)
             for t in F.TIPOS if t in tipos_usados]
    fig.legend(handles=alcas, loc="lower left", ncol=len(alcas), frameon=False, borderaxespad=0, fontsize=6.2,
               handletextpad=0.3, columnspacing=1.2, title=f"Viral genome type of the approved EVE loci (number above a stack = loci within "
                     f"{1.6 * mb_por_mm:.0f} Mb of each other)",
               title_fontsize=6.5, alignment="left",
               bbox_to_anchor=(esq * MM / W, (12 + len(linhas) * passo + 3) * MM / H))
    return F._salvar(fig, pasta, "fig7_chromosome_map", formatos)


# ---------------------------------------------------------------------------
def gerar(outdir, formatos=("pdf", "svg", "png"), online=True, d=None):
    d = d or F.carregar(outdir)
    pasta = os.path.join(outdir, F.SUBPASTA, F.PASTA_FIG)
    os.makedirs(pasta, exist_ok=True)
    tab = metadados(outdir, d["species"], pasta, online=online)
    est = estatisticas(outdir, d, tab, pasta)
    tabela(est, pasta)
    cor = correlacoes(est)
    F._gravar_tsv(os.path.join(pasta, "fig6_genome_correlations.tsv"),
                  ["x", "y", "n", "spearman_rho", "p_value"], [[a, b, n, round(r, 3), round(p, 4)] for a, b, n, r, p in cor])
    plt = F._estilo()
    saidas = [os.path.join(pasta, ARQ_GENOMAS), os.path.join(pasta, "fig5_genome_features.tsv"),
              os.path.join(pasta, "fig6_genome_correlations.tsv")]
    saidas += fig5(est, plt, pasta, formatos)
    plt.close("all")
    saidas += fig6(est, plt, pasta, formatos)
    plt.close("all")
    saidas += fig7(est, d, plt, pasta, formatos)
    plt.close("all")
    n = len(_usadas(est))
    log(f"{n} species in use ({', '.join(e['species'] for e in _usadas(est))})")
    for s in saidas:
        log(f"written: {s}")
    return {"folder": pasta, "files": saidas, "species": n}


def main():
    ap = argparse.ArgumentParser(description="Host-genome figures of the curated EVE loci.")
    ap.add_argument("-o", "--outdir", required=True, help="Results folder of the run")
    ap.add_argument("--formats", default="pdf,svg,png", help="Comma-separated: pdf,svg,png (default all)")
    ap.add_argument("--offline", action="store_true", help="Do not query NCBI Datasets")
    a = ap.parse_args()
    gerar(a.outdir, [f.strip() for f in a.formats.split(",") if f.strip()], online=not a.offline)


if __name__ == "__main__":
    main()
