#!/usr/bin/env python3
"""
Figures of the manually curated EVEs, for the paper and for sharing with the group.

The unit is the approved EVE locus: one locus representative that was not removed, the
same rows as the curation page's "one row per locus" view, so every number here matches
the page. The data come from the live curation state (curation_selection.json applied to
eves_classified.tsv), not from the exported TSV.

Figures (each as PDF, SVG and 300 dpi PNG, plus the table behind it):
  fig1_loci_per_species      approved loci per host species, split by viral genome type
  fig2_family_by_species     heatmap: viral family x host species (number of loci)
  fig3_family_by_protein     heatmap: viral family x protein class (number of loci)
  fig4_identity_by_family    amino-acid identity to the closest virus, per family
  fig5-7                     host genome (NCBI assemblies): see genome_figures.py
and figuras_curadoria.html, one self-contained interactive page with the same four charts.

The viral genome type comes from the NCBI lineage of the curated family (phylum/class),
then from the best viral hit when the family is not in the taxonomy cache. Loci left at
realm level (Riboviria, Viruses) or without a family stay "Unassigned".

Host species order and groups come from 13_curadoria/figuras/host_groups.tsv
(species <tab> group); a template is written on the first run, edit it and run again.

Usage:
  python -m evetool.curation.figures -o results_eve
  eve-tool figures -o results_eve --formats pdf,png
"""
import os
import csv
import json
import argparse
from collections import Counter, OrderedDict, defaultdict

SUBPASTA = "13_curadoria"
PASTA_FIG = "figuras"
ARQ_GRUPOS = "host_groups.tsv"

# Viral genome types, in the fixed colour order (colour follows the type, never its rank)
TIPOS = ["(−)ssRNA", "(+)ssRNA", "dsRNA", "ssDNA", "dsDNA", "RT", "Unassigned"]
COR_TIPO = {   # light-mode categorical slots 1-6 of the curation page; neutral grey = no type
    "(−)ssRNA": "#2a78d6", "(+)ssRNA": "#eb6834", "dsRNA": "#1baf7a", "ssDNA": "#eda100",
    "dsDNA": "#e87ba4", "RT": "#008300", "Unassigned": "#a8a7a0",
}
TINTA, TINTA_2, MUDO = "#0b0b0b", "#52514e", "#898781"
GRADE, BASE = "#e1e0d9", "#c3c2b7"
# Sequential blue ramp (steps 150-700) for counts; 0 stays blank
RAMPA = ["#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]
FAIXAS = [(1, 1), (2, 3), (4, 7), (8, 15), (16, None)]
CLASSES_NO_FIM = ("hypothetical", "other", "HMM model")


def log(msg):
    print(f"[figures] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Viral taxonomy -> genome type
# ---------------------------------------------------------------------------
class Taxonomia:
    """NCBI tree cached by the dashboard (pai/rank/nome), looked up by taxon name."""

    def __init__(self, caminho):
        self.pai, self.rank, self.nome, self.por_nome = {}, {}, {}, {}
        if caminho and os.path.exists(caminho):
            with open(caminho, encoding="utf-8") as fh:
                t = json.load(fh)
            self.pai, self.rank, self.nome = t.get("pai", {}), t.get("rank", {}), t.get("nome", {})
            for tid, n in self.nome.items():
                self.por_nome.setdefault(n, tid)

    def linhagem(self, nome):
        """[(rank, name), ...] from the root down to `nome`; [] when unknown."""
        tid, out, vistos = self.por_nome.get(nome), [], set()
        while tid and tid not in vistos:
            vistos.add(tid)
            out.append((self.rank.get(tid, ""), self.nome.get(tid, "")))
            if self.pai.get(tid) == tid:
                break
            tid = self.pai.get(tid)
        return out[::-1]


def _tipo_da_linhagem(nomes: set):
    if nomes & {"Negarnaviricota", "Ribozyviria"}:
        return "(−)ssRNA"
    if nomes & {"Duplornaviricota", "Duplopiviricetes"}:
        return "dsRNA"
    if nomes & {"Pisuviricota", "Kitrinoviricota", "Lenarviricota"}:
        return "(+)ssRNA"
    if nomes & {"Pararnavirae", "Artverviricota"}:
        return "RT"
    if nomes & {"Papovaviricetes", "Duplodnaviria", "Varidnaviria", "Adnaviria"}:
        return "dsDNA"
    if nomes & {"Monodnaviria", "Floreoviria", "Cressdnaviricota", "Quintoviricetes"}:
        return "ssDNA"
    return None


def _tipo_do_hit(texto: str):
    """Genome type field of the viral database title: ssRNA(-), dsRNA, ssDNA(+/-), ..."""
    t = (texto or "").replace(" ", "")
    if "RT" in t:
        return "RT"
    if t.startswith("ssRNA(-)"):
        return "(−)ssRNA"
    if t.startswith("ssRNA(+)") and "/" not in t:
        return "(+)ssRNA"
    if t.startswith("dsRNA"):
        return "dsRNA"
    if t.startswith("ssDNA"):
        return "ssDNA"
    if t.startswith("dsDNA"):
        return "dsDNA"
    return None   # 'RNA', 'unknown', ambisense


RANKS_ALTOS = {"superkingdom", "acellular root", "realm", "kingdom", "subkingdom", "no rank", "clade", ""}


def classificar_familia(familia: str, tax: Taxonomia, tipo_hit: str = ""):
    """
    (label, genome type) of a curated family name.
    Family names are kept; orders/classes become "unclassified X"; a name at realm level,
    "Viruses" or nothing stays Unassigned (the curator left it open, so do we).
    """
    f = (familia or "").strip()
    if f in ("", "n/a", "nan", "unclassified"):
        return "no assignment", "Unassigned"
    linha = tax.linhagem(f)
    ranks = {r: n for r, n in linha}
    rank_f = linha[-1][0] if linha else ""
    nomes = {n for _, n in linha}
    baixo = f.lower()
    if linha:
        tipo = _tipo_da_linhagem(nomes)
        if rank_f == "family":
            return f, tipo or _tipo_do_hit(tipo_hit) or "Unassigned"
        if "phylum" not in ranks and "class" not in ranks and "order" not in ranks:
            return f"unclassified {f}" if f != "Viruses" else "unclassified virus", "Unassigned"
        return f"unclassified {f}", tipo or "Unassigned"
    # Not in the taxonomy cache: the name itself, then the best viral hit
    if "negative-strand" in baixo:
        return "unclassified (−)ssRNA virus", "(−)ssRNA"
    if "positive-strand" in baixo:
        return "unclassified (+)ssRNA virus", "(+)ssRNA"
    if baixo.startswith("dsrna"):
        return "unclassified dsRNA virus", "dsRNA"
    if baixo.startswith("ssdna"):
        return "unclassified ssDNA virus", "ssDNA"
    tipo = _tipo_do_hit(tipo_hit) or "Unassigned"
    return (f if baixo.endswith("viridae") else f"unclassified {f}"), tipo


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def _ler_grupos(caminho):
    """species -> group, in file order; {} when there is no file."""
    grupos = OrderedDict()
    if os.path.exists(caminho):
        with open(caminho, encoding="utf-8") as fh:
            for r in csv.reader(fh, delimiter="\t"):
                if not r or r[0].startswith("#") or r[0] == "species":
                    continue
                grupos[r[0].strip()] = (r[1].strip() if len(r) > 1 else "")
    return grupos


def _tamanho_montagem(fai):
    """(total length, number of sequences) from a .fai; (0, 0) when missing."""
    total = n = 0
    if fai and os.path.exists(fai):
        with open(fai, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                p = linha.split("\t", 2)
                if len(p) >= 2:
                    try:
                        total += int(p[1])
                        n += 1
                    except ValueError:
                        pass
    return total, n


def carregar(outdir: str) -> dict:
    """Approved EVE loci with their curated family, genome type and protein class."""
    from evetool.curation import curation as cur
    dados = cur.listar_candidatas(outdir)
    if dados.get("error"):
        raise SystemExit(dados["error"])
    itens = [i for i in dados["items"] if not i.get("removed_flag") and i.get("locus_representative") is not False]

    tax = Taxonomia(os.path.join(outdir, "12_dashboard_web", "dados", ".cache_taxonomia_viral.json"))
    refs = {}
    mapa = os.path.join(outdir, SUBPASTA, "protein_map_cache.json")
    if os.path.exists(mapa):
        try:
            with open(mapa, encoding="utf-8") as fh:
                refs = (json.load(fh).get("result") or {}).get("references") or {}
        except Exception:
            refs = {}

    loci = []
    for i in itens:
        tipo_hit = (refs.get(i.get("protein") or "") or {}).get("genome_type", "")
        rotulo, tipo = classificar_familia(i.get("family", ""), tax, tipo_hit)
        pid = i.get("pident")
        try:
            pid = float(pid)
        except (TypeError, ValueError):
            pid = 0.0
        loci.append({
            "species": i.get("species", ""), "locus": i.get("locus", ""), "orf": i.get("orf", ""),
            "family": i.get("family", "") or "", "family_label": rotulo, "genome_type": tipo,
            "protein_class": i.get("protein_class", "") or "n/a",
            "pident": pid if pid > 0 else None, "qcov": i.get("qcov"),
            "group": i.get("group", ""), "product": i.get("product", "") or "",
            "protein": i.get("protein", "") or "",
            "family_manual": bool(i.get("family_manual")), "class_manual": bool(i.get("class_manual")),
        })

    # Host species: order and groups from host_groups.tsv, otherwise by number of loci
    pasta_fig = os.path.join(outdir, SUBPASTA, PASTA_FIG)
    os.makedirs(pasta_fig, exist_ok=True)
    arq_grupos = os.path.join(pasta_fig, ARQ_GRUPOS)
    grupos = _ler_grupos(arq_grupos)
    todas = list(cur.pastas_especies(outdir).keys())
    por_especie = Counter(l["species"] for l in loci)
    if not grupos:
        ordem = sorted(todas, key=lambda s: (-por_especie[s], s))
        with open(arq_grupos, "w", encoding="utf-8") as fh:
            fh.write("# species<TAB>group: the order of this file is the order of the figures. "
                     "Edit and run again.\nspecies\tgroup\n")
            for s in ordem:
                fh.write(f"{s}\t\n")
        log(f"host order/groups template written: {arq_grupos}")
        grupos = OrderedDict((s, "") for s in ordem)
    especies = [s for s in grupos if s in todas] + [s for s in todas if s not in grupos]

    montagem = {}
    for s, pasta in cur.pastas_especies(outdir).items():
        genoma = cur._info_pipeline(pasta).get("genome") or ""
        tam, n = _tamanho_montagem(genoma + ".fai")
        montagem[s] = {"length": tam, "contigs": n}

    # Families: grouped by genome type, proper families first (most loci first), then the
    # "unclassified X" rows of that type
    contagem_fam = Counter(l["family_label"] for l in loci)
    tipo_fam = {l["family_label"]: l["genome_type"] for l in loci}
    familias = sorted(contagem_fam, key=lambda f: (TIPOS.index(tipo_fam[f]),
                                                   f.startswith("unclassified") or f == "no assignment",
                                                   -contagem_fam[f], f.lower()))
    contagem_cls = Counter(l["protein_class"] for l in loci)
    classes = sorted(contagem_cls, key=lambda c: (c in CLASSES_NO_FIM, CLASSES_NO_FIM.index(c)
                                                  if c in CLASSES_NO_FIM else 0, -contagem_cls[c], c.lower()))
    return {"outdir": os.path.abspath(outdir), "loci": loci, "species": especies,
            "host_group": {s: grupos.get(s, "") for s in especies}, "assembly": montagem,
            "families": familias, "family_type": tipo_fam, "classes": classes,
            "types": [t for t in TIPOS if any(l["genome_type"] == t for l in loci)],
            "n_removed": sum(1 for i in dados["items"] if i.get("removed_flag")),
            "n_orfs": len([i for i in dados["items"] if not i.get("removed_flag")])}


# ---------------------------------------------------------------------------
# Tables behind the figures (the table view of each chart)
# ---------------------------------------------------------------------------
def _gravar_tsv(caminho, cabecalho, linhas):
    with open(caminho, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cabecalho)
        w.writerows(linhas)


def tabelas(d: dict, pasta: str):
    loci = d["loci"]
    _gravar_tsv(os.path.join(pasta, "curated_loci.tsv"),
                ["species", "host_group", "locus", "orf", "family_curated", "family_label", "genome_type",
                 "protein_class", "pident", "qcov", "evidence_group", "product", "protein",
                 "family_manual", "class_manual"],
                [[l["species"], d["host_group"].get(l["species"], ""), l["locus"], l["orf"], l["family"],
                  l["family_label"], l["genome_type"], l["protein_class"],
                  "" if l["pident"] is None else l["pident"], l["qcov"], l["group"], l["product"],
                  l["protein"], l["family_manual"], l["class_manual"]] for l in loci])
    por = Counter((l["species"], l["genome_type"]) for l in loci)
    _gravar_tsv(os.path.join(pasta, "fig1_loci_per_species.tsv"),
                ["species", "host_group", *d["types"], "total", "assembly_mb", "contigs"],
                [[s, d["host_group"].get(s, ""), *[por[(s, t)] for t in d["types"]],
                  sum(por[(s, t)] for t in d["types"]), round(d["assembly"][s]["length"] / 1e6, 1),
                  d["assembly"][s]["contigs"]] for s in d["species"]])
    m = Counter((l["family_label"], l["species"]) for l in loci)
    _gravar_tsv(os.path.join(pasta, "fig2_family_by_species.tsv"),
                ["family", "genome_type", *d["species"], "total"],
                [[f, d["family_type"][f], *[m[(f, s)] for s in d["species"]],
                  sum(m[(f, s)] for s in d["species"])] for f in d["families"]])
    c = Counter((l["family_label"], l["protein_class"]) for l in loci)
    _gravar_tsv(os.path.join(pasta, "fig3_family_by_protein.tsv"),
                ["family", "genome_type", *d["classes"], "total"],
                [[f, d["family_type"][f], *[c[(f, k)] for k in d["classes"]],
                  sum(c[(f, k)] for k in d["classes"])] for f in d["families"]])
    linhas = []
    for f in d["families"]:
        v = sorted(l["pident"] for l in loci if l["family_label"] == f and l["pident"] is not None)
        sem = sum(1 for l in loci if l["family_label"] == f and l["pident"] is None)
        if v:
            linhas.append([f, d["family_type"][f], len(v), sem, v[0], _mediana(v), v[-1]])
        else:
            linhas.append([f, d["family_type"][f], 0, sem, "", "", ""])
    _gravar_tsv(os.path.join(pasta, "fig4_identity_by_family.tsv"),
                ["family", "genome_type", "n_with_identity", "n_hmm_only", "min_pident",
                 "median_pident", "max_pident"], linhas)


def _mediana(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else round((v[n // 2 - 1] + v[n // 2]) / 2, 2)


# ---------------------------------------------------------------------------
# Static figures (matplotlib)
# ---------------------------------------------------------------------------
MM = 1 / 25.4


def _estilo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fonte = "Liberation Sans" if any("Liberation Sans" == f.name for f in font_manager.fontManager.ttflist) \
        else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": fonte, "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.edgecolor": BASE, "axes.labelcolor": TINTA_2, "xtick.color": MUDO, "ytick.color": MUDO,
        "text.color": TINTA, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 0, "svg.fonttype": "none", "pdf.fonttype": 42,
        "mathtext.fontset": "custom", "mathtext.rm": fonte, "mathtext.it": f"{fonte}:italic",
        "mathtext.bf": f"{fonte}:bold", "mathtext.sf": fonte, "mathtext.cal": fonte, "mathtext.tt": fonte,
        "figure.dpi": 150, "savefig.dpi": 300,
    })
    return plt


def _it(texto):
    """Italic for taxon names inside matplotlib text (mathtext, same sans face)."""
    t = texto.replace("\\", "").replace("$", "")
    return "$\\mathit{" + t.replace(" ", "\\ ").replace("(", "(").replace("_", "\\_") + "}$"


def rotulo_especie(s):
    nome = s.replace("_", " ").strip()
    if nome.endswith(" sp"):
        return _it(nome[:-3]) + " sp."
    return _it(nome)


ITALICO_SUFIXOS = ("viridae", "virinae", "virales", "viricetes", "viricota", "virae", "viria")


def rotulo_familia(f):
    if f.startswith("unclassified "):
        resto = f[len("unclassified "):]
        return "unclassified " + (_it(resto) if resto.lower().endswith(ITALICO_SUFIXOS) else resto)
    return _it(f) if f.lower().endswith(ITALICO_SUFIXOS) else f


def _salvar(fig, pasta, nome, formatos):
    saidas = []
    for fmt in formatos:
        caminho = os.path.join(pasta, f"{nome}.{fmt}")
        fig.savefig(caminho, bbox_inches="tight", pad_inches=0.04, facecolor="white")
        saidas.append(caminho)
    return saidas


def _faixa(n):
    for k, (a, b) in enumerate(FAIXAS):
        if n >= a and (b is None or n <= b):
            return k
    return None


def _tinta_na_celula(k):
    return "#ffffff" if k >= 2 else TINTA


def _linhas_com_grupos(especies, grupos):
    """[('group', name) | ('species', s), ...]: a header row before each run of a host group."""
    linhas, atual = [], None
    for s in especies:
        g = grupos.get(s, "")
        if g and g != atual:
            linhas.append(("group", g))
        atual = g
        linhas.append(("species", s))
    return linhas


def fig1(d, plt, pasta, formatos):
    """Horizontal stacked bars: loci per host species, by viral genome type."""
    from matplotlib.patches import Rectangle, Patch
    especies = d["species"]
    por = Counter((l["species"], l["genome_type"]) for l in d["loci"])
    total = {s: sum(por[(s, t)] for t in d["types"]) for s in especies}
    linhas = _linhas_com_grupos(especies, d["host_group"])
    # y of each row, top to bottom; a group header takes 0.8 of a species row
    ys, y = [], 0.0
    for tipo, _ in linhas:
        y -= 0.85 if tipo == "group" else 1.0
        ys.append(y)
    passo = 5.6 * MM
    alto_plot = -ys[-1] * passo + 2 * MM
    W = 180 * MM
    altura = alto_plot + 30 * MM
    fig = plt.figure(figsize=(W, altura))
    esq_mm, dir_mm = 50, 124
    ax = fig.add_axes([esq_mm * MM / W, 11 * MM / altura, (dir_mm - esq_mm) * MM / W, alto_plot / altura])
    xmax = max(total.values()) if total else 1
    ax.set_xlim(0, xmax * 1.1)        # limits first: the text columns below use the transform
    ax.set_ylim(ys[-1] - 0.6, 0)
    for (tipo, nome), y in zip(linhas, ys):
        if tipo == "group":
            ax.text(-0.02, y, nome, transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=6.5, color=TINTA, fontweight="bold")
            continue
        x = 0.0
        for t in d["types"]:
            v = por[(nome, t)]
            if v:
                ax.add_patch(Rectangle((x, y - 0.34), v, 0.68, facecolor=COR_TIPO[t],
                                       edgecolor="white", linewidth=0.7))
                x += v
        ax.text(x + xmax * 0.012, y, f"{total[nome]}", va="center", ha="left", fontsize=6.5, color=TINTA)
        ax.text(-0.02, y, rotulo_especie(nome), transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=6.5, color=TINTA)
        m = d["assembly"].get(nome, {})
        for x_mm, txt in ((142, f"{m.get('length', 0) / 1e6:,.0f}" if m.get("length") else "–"),
                          (170, f"{m.get('contigs', 0):,}" if m.get("contigs") else "–")):
            fy = fig.transFigure.inverted().transform(ax.transData.transform((0, y)))[1]
            fig.text(x_mm * MM / W, fy, txt, ha="right", va="center", fontsize=6.5, color=TINTA_2)
    ax.set_yticks([])
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.xaxis.grid(True, color=GRADE, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Approved EVE loci")
    topo = (11 * MM + alto_plot) / altura
    # assembly column headers: on the first group header row (empty on the right), else above
    if linhas and linhas[0][0] == "group":
        y_cab = fig.transFigure.inverted().transform(ax.transData.transform((0, ys[0])))[1]
        va = "center"
    else:
        y_cab, va = topo + 1.5 * MM / altura, "bottom"
    for x_mm, txt in ((142, "Assembly (Mb)"), (170, "Sequences")):
        fig.text(x_mm * MM / W, y_cab, txt, ha="right", va=va, fontsize=6.2, color=MUDO)
    alcas = [Patch(facecolor=COR_TIPO[t], edgecolor="none", label=t) for t in d["types"]]
    fig.legend(handles=alcas, loc="lower left", ncol=len(alcas), frameon=False,
               bbox_to_anchor=(esq_mm * MM / W, topo + 2 * MM / altura), handlelength=1.1,
               handleheight=0.9, columnspacing=1.3, borderaxespad=0, title="Viral genome type",
               title_fontsize=6.5, alignment="left")
    return _salvar(fig, pasta, "fig1_loci_per_species", formatos)


def plt_line(xs, ys, fig):
    from matplotlib.lines import Line2D
    return Line2D(xs, ys, transform=fig.transFigure, color=BASE, linewidth=0.8, solid_capstyle="butt")


def _runs_tipo(familias, tipo_fam):
    """[(genome type, first row, last row), ...] for consecutive families of the same type."""
    runs, atual, ini = [], None, 0
    for r, f in enumerate(list(familias) + [None]):
        t = tipo_fam.get(f) if f else None
        if t != atual:
            if atual:
                runs.append((atual, ini, r - 1))
            atual, ini = t, r
    return runs


def _margem_familias(familias, tipo_fam):
    """Left margin (mm): genome type band + the longest family name at 6.5 pt; a type with a
    single row has its label written across, so that row needs room for both."""
    margem = 12 + max((_larg_mm(f) for f in familias), default=12) + 3
    for t, a, b in _runs_tipo(familias, tipo_fam):
        if not _gira(t, a, b):
            margem = max(margem, 6.4 + _larg_mm(t) + 3 + max(_larg_mm(f) for f in familias[a:b + 1]) + 3)
    return margem


def _larg_mm(txt):
    """Approximate width of a 6.5 pt label, in mm."""
    return len(txt) * 1.18


def _gira(tipo, a, b, passo_mm=4.6):
    """Vertical genome type label only when it fits along its rows."""
    return (b - a + 1) * passo_mm >= _larg_mm(tipo) + 1


def _faixas_tipo(fig, ax, familias, tipo_fam, W, linha_sep=True):
    """Coloured band + vertical label per genome type, left of the family names, and a
    separator line between types."""
    from matplotlib.patches import Rectangle
    nl = len(familias)
    inv = fig.transFigure.inverted()
    for t, a, b in _runs_tipo(familias, tipo_fam):
        y_a = inv.transform(ax.transData.transform((0, nl - 1 - a + 0.42)))[1]
        y_b = inv.transform(ax.transData.transform((0, nl - 1 - b - 0.42)))[1]
        x = 4 * MM / W
        fig.add_artist(Rectangle((x, y_b), 1.4 * MM / W, y_a - y_b, transform=fig.transFigure,
                                 facecolor=COR_TIPO[t], edgecolor="none"))
        fig.text(x + 2.4 * MM / W, (y_a + y_b) / 2, t, ha="left", va="center",
                 rotation=90 if _gira(t, a, b) else 0, fontsize=6.5, color=TINTA_2)
        if linha_sep and b < nl - 1:
            ax.axhline(nl - 1 - b - 0.5, color=BASE, linewidth=0.8, zorder=3)


def _heatmap(d, plt, colunas, rotulos_col, contagem, grupos_col, nome, pasta, formatos, larg_col_mm=7.0):
    """Families (rows, grouped by genome type) x columns, count of loci per cell."""
    from matplotlib.patches import Rectangle, Patch
    familias = d["families"]
    nl, nc = len(familias), len(colunas)
    passo = 4.6 * MM
    esq_mm = _margem_familias(familias, d["family_type"])
    larg = esq_mm + nc * larg_col_mm + 12
    topo_mm = 40 if grupos_col else 32
    altura = nl * passo + (topo_mm + 18) * MM
    W = larg * MM
    fig = plt.figure(figsize=(W, altura))
    ax = fig.add_axes([esq_mm * MM / W, 16 * MM / altura, nc * larg_col_mm * MM / W, nl * passo / altura])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(-0.5, nl - 0.5)
    for r, f in enumerate(familias):
        y = nl - 1 - r
        for k, c in enumerate(colunas):
            v = contagem[(f, c)]
            if not v:
                continue
            b = _faixa(v)
            ax.add_patch(Rectangle((k - 0.5, y - 0.5), 1, 1, facecolor=RAMPA[b], edgecolor="white", linewidth=0.8))
            ax.text(k, y, str(v), ha="center", va="center", fontsize=5.8, color=_tinta_na_celula(b))
        tot = sum(contagem[(f, c)] for c in colunas)
        ax.text(nc - 0.5 + 0.35, y, str(tot), ha="left", va="center", fontsize=6, color=TINTA_2)
    for k in range(nc + 1):        # hairline grid, so empty cells still read as cells
        ax.axvline(k - 0.5, color=GRADE, linewidth=0.4, zorder=0)
    for r in range(nl + 1):
        ax.axhline(r - 0.5, color=GRADE, linewidth=0.4, zorder=0)
    ax.set_yticks(range(nl))
    ax.set_yticklabels([rotulo_familia(f) for f in familias][::-1], color=TINTA)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(rotulos_col, rotation=90, color=TINTA)
    ax.xaxis.tick_top()
    ax.tick_params(axis="both", length=0, pad=2)
    for lado in ax.spines.values():
        lado.set_visible(False)
    ax.text(nc - 0.5 + 0.35, nl - 0.5 + 0.15, "Total", ha="left", va="bottom", fontsize=6, color=MUDO)
    _faixas_tipo(fig, ax, familias, d["family_type"], W)
    runs = []
    if grupos_col:   # host groups: a separator between runs, a bracket + label above the names
        atual, ini = None, 0
        for k, c in enumerate(list(colunas) + [None]):
            g = grupos_col.get(c, "") if c else None
            if g != atual:
                if atual:
                    runs.append((atual, ini, k - 1))
                atual, ini = g, k
        for g, a, b in runs:
            if b < nc - 1:
                ax.axvline(b + 0.5, color=BASE, linewidth=0.8, zorder=3)
    rot = [f"{a}" if b == a else (f"{a}–{b}" if b else f"≥{a}") for a, b in FAIXAS]
    alcas = [Patch(facecolor=RAMPA[k], edgecolor="none", label=rot[k]) for k in range(len(FAIXAS))]
    fig.legend(handles=alcas, loc="upper left", ncol=len(alcas), frameon=False,
               bbox_to_anchor=(esq_mm * MM / W, 12 * MM / altura), handlelength=1.1, handleheight=0.9,
               columnspacing=1.2, borderaxespad=0, title="Approved EVE loci per cell", title_fontsize=6.5,
               alignment="left")
    if runs:
        fig.canvas.draw()
        inv = fig.transFigure.inverted()
        rend = fig.canvas.get_renderer()
        topo_lbl = max(t.get_window_extent(rend).y1 for t in ax.get_xticklabels())
        y_fig = inv.transform((0, topo_lbl))[1] + 1.5 * MM / altura
        for g, a, b in runs:
            xa = inv.transform(ax.transData.transform((a - 0.4, 0)))[0]
            xb = inv.transform(ax.transData.transform((b + 0.4, 0)))[0]
            fig.add_artist(plt_line([xa, xb], [y_fig, y_fig], fig))
            # "Phlebotominae (Old World)" -> two lines, so a narrow run still fits
            fig.text((xa + xb) / 2, y_fig + 1.0 * MM / altura, g.replace(" (", "\n("), ha="center",
                     va="bottom", fontsize=6.2, color=TINTA_2, linespacing=1.1)
    return _salvar(fig, pasta, nome, formatos)


def fig2(d, plt, pasta, formatos):
    m = Counter((l["family_label"], l["species"]) for l in d["loci"])
    return _heatmap(d, plt, d["species"], [rotulo_especie(s) for s in d["species"]], m,
                    d["host_group"], "fig2_family_by_species", pasta, formatos, larg_col_mm=6.2)


def fig3(d, plt, pasta, formatos):
    c = Counter((l["family_label"], l["protein_class"]) for l in d["loci"])
    return _heatmap(d, plt, d["classes"], d["classes"], c, None,
                    "fig3_family_by_protein", pasta, formatos, larg_col_mm=7.5)


def fig4(d, plt, pasta, formatos):
    """Strip plot: identity of each approved locus to its closest virus, per family.
    Families whose loci all lack a DIAMOND hit (HMM-only) have no identity and are left out."""
    import random
    com_id = {l["family_label"] for l in d["loci"] if l["pident"] is not None}
    familias = [f for f in d["families"] if f in com_id]
    fora = [f for f in d["families"] if f not in com_id]
    n_fora = sum(1 for l in d["loci"] if l["family_label"] in fora)
    nl = len(familias)
    passo = 4.6 * MM
    esq_mm = _margem_familias(familias, d["family_type"])
    W = (esq_mm + 118) * MM
    altura = nl * passo + 28 * MM
    fig = plt.figure(figsize=(W, altura))
    ax = fig.add_axes([esq_mm * MM / W, 14 * MM / altura, 88 * MM / W, nl * passo / altura])
    ax.set_xlim(15, 100)
    ax.set_ylim(-0.6, nl - 0.4)
    rnd = random.Random(7)
    for r, f in enumerate(familias):
        y = nl - 1 - r
        vals = [l["pident"] for l in d["loci"] if l["family_label"] == f and l["pident"] is not None]
        sem = sum(1 for l in d["loci"] if l["family_label"] == f and l["pident"] is None)
        ys = [y + rnd.uniform(-0.26, 0.26) for _ in vals]
        ax.scatter(vals, ys, s=7, color=COR_TIPO[d["family_type"][f]], edgecolors="white",
                   linewidths=0.4, zorder=3, alpha=0.9)
        if len(vals) >= 3:
            med = _mediana(vals)
            ax.plot([med, med], [y - 0.38, y + 0.38], color=TINTA, linewidth=1.2, zorder=4,
                    solid_capstyle="butt")
        nota = f"n = {len(vals)}" + (f"  ·  {sem} without hit" if sem else "")
        ax.text(1.02, y, nota, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.8,
                color=TINTA_2)
    ax.set_yticks(range(nl))
    ax.set_yticklabels([rotulo_familia(f) for f in familias][::-1], color=TINTA)
    ax.tick_params(axis="y", length=0, pad=3)
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.xaxis.grid(True, color=GRADE, linewidth=0.5)
    for r in range(nl):
        ax.axhline(r - 0.5, color=GRADE, linewidth=0.3, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel("Amino-acid identity to the closest viral protein (%)")
    _faixas_tipo(fig, ax, familias, d["family_type"], W)
    fig.text(esq_mm * MM / W, (14 * MM + nl * passo) / altura + 3 * MM / altura,
             "Each dot is one approved locus (best DIAMOND hit); black tick = median (n ≥ 3).",
             ha="left", va="bottom", fontsize=6.5, color=TINTA_2)
    if fora:
        fig.text(esq_mm * MM / W, 2 * MM / altura,
                 f"Not shown: {n_fora} loci in {len(fora)} families without a DIAMOND hit (HMM evidence only).",
                 ha="left", va="bottom", fontsize=6, color=MUDO)
    return _salvar(fig, pasta, "fig4_identity_by_family", formatos)


# ---------------------------------------------------------------------------
# Interactive page
# ---------------------------------------------------------------------------
def pagina(d: dict, pasta: str, destino: str = None, documento: bool = True) -> str:
    """
    figuras_curadoria.html: one self-contained page (data inline, no library, no server).
    documento=False writes only the content (title, style, body, script), for hosts that wrap
    the page in their own <html> skeleton.
    """
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
    with open(os.path.join(base, "figures.html"), encoding="utf-8") as fh:
        html = fh.read()
    if documento:
        cabeca, _, corpo_html = html.partition("</style>")
        html = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                + cabeca + "</style>\n</head>\n<body>" + corpo_html + "</body>\n</html>\n")
    import datetime
    dados = {
        "generated": datetime.datetime.now().isoformat(timespec="minutes"),
        "run": os.path.basename(d["outdir"].rstrip("/")),
        "species": d["species"], "host_group": d["host_group"], "assembly": d["assembly"],
        "families": d["families"], "family_type": d["family_type"], "classes": d["classes"],
        "types": d["types"], "n_removed": d["n_removed"], "n_orfs": d["n_orfs"],
        "loci": [{k: l[k] for k in ("species", "locus", "orf", "family_label", "genome_type",
                                    "protein_class", "pident", "product", "group")} for l in d["loci"]],
    }
    corpo = json.dumps(dados, ensure_ascii=False).replace("</", "<\\/")
    html = html.replace("/*__EVE_FIGURES__*/null", corpo, 1)
    destino = destino or os.path.join(pasta, "figuras_curadoria.html")
    with open(destino, "w", encoding="utf-8") as fh:
        fh.write(html)
    return destino


# ---------------------------------------------------------------------------
def gerar(outdir: str, formatos=("pdf", "svg", "png"), html=True, genomas=True, online=True) -> dict:
    d = carregar(outdir)
    pasta = os.path.join(outdir, SUBPASTA, PASTA_FIG)
    os.makedirs(pasta, exist_ok=True)
    log(f"{len(d['loci'])} approved loci in {len(d['species'])} species, {len(d['families'])} family rows")
    tabelas(d, pasta)
    plt = _estilo()
    saidas = []
    for f in (fig1, fig2, fig3, fig4):
        saidas += f(d, plt, pasta, formatos)
        plt.close("all")
    if html:
        saidas.append(pagina(d, pasta))
    for s in saidas:
        log(f"written: {s}")
    if genomas:     # fig5-7 (host genome); a failure there must not cost fig1-4
        try:
            from evetool.curation import genome_figures
            saidas += genome_figures.gerar(outdir, formatos, online=online, d=d)["files"]
        except Exception as erro:
            log(f"host-genome figures (fig5-7) skipped: {erro}")
    return {"folder": pasta, "files": saidas, "n_loci": len(d["loci"])}


def main():
    ap = argparse.ArgumentParser(description="Figures of the manually curated EVE loci.")
    ap.add_argument("-o", "--outdir", required=True, help="Results folder of the run")
    ap.add_argument("--formats", default="pdf,svg,png", help="Comma-separated: pdf,svg,png (default all)")
    ap.add_argument("--no-html", action="store_true", help="Skip the interactive page")
    ap.add_argument("--no-genomes", action="store_true", help="Skip fig5-7 (host-genome figures)")
    ap.add_argument("--offline", action="store_true", help="Do not query NCBI Datasets for fig5-7")
    a = ap.parse_args()
    gerar(a.outdir, [f.strip() for f in a.formats.split(",") if f.strip()], html=not a.no_html,
          genomas=not a.no_genomes, online=not a.offline)


if __name__ == "__main__":
    main()
