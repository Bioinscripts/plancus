import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Headless operation, no GUI
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D


def carregar_anotacoes_gff(gff_path: str, scaffold: str, start_pos: int, end_pos: int) -> list:
    """
    Extracts the genes and exons from the host genome GFF3/GTF that intersect the
    genomic-context window [start_pos, end_pos].
    Returns an empty list when gff_path is None or the file does not exist.
    """
    if not gff_path or not os.path.exists(gff_path):
        return []

    features = []
    with open(gff_path, "r") as f:
        for linha in f:
            if linha.startswith("#") or not linha.strip():
                continue
            partes = linha.strip().split("\t")
            if len(partes) < 9:
                continue

            seqid, source, ftype, fstart, fend, score, strand, frame, attrs = partes[:9]
            if seqid != scaffold:
                continue

            fstart, fend = int(fstart), int(fend)
            # Check the overlap with the window of interest
            if fend < start_pos or fstart > end_pos:
                continue

            # Extrai nome ou ID do atributo
            attr_dict = {}
            for par in attrs.split(";"):
                if "=" in par:
                    k, v = par.split("=", 1)
                    attr_dict[k.strip()] = v.strip()
                elif " " in par:
                    tokens = par.strip().split(" ", 1)
                    if len(tokens) == 2:
                        attr_dict[tokens[0].strip()] = tokens[1].strip().strip('"')

            nome = attr_dict.get("Name") or attr_dict.get("ID") or attr_dict.get("gene_id") or ftype
            features.append({
                "type": ftype,
                "start": fstart,
                "end": fend,
                "strand": strand,
                "name": nome,
                "parent": attr_dict.get("Parent") or attr_dict.get("gene_id") or ""
            })

    return features


def carregar_profundidade_locus(depth_tsv: str, cluster_id: str) -> pd.DataFrame:
    """
    Loads the positions and depths computed by samtools depth
    para o cluster especificado.
    """
    if not depth_tsv or not os.path.exists(depth_tsv):
        return None

    posicoes = []
    profundidades = []
    cid_alvo = cluster_id.split("(")[0]

    with open(depth_tsv, "r") as f:
        for linha in f:
            partes = linha.strip().split("\t")
            if len(partes) < 3:
                continue
            seq_id, pos, depth = partes[0], int(partes[1]), float(partes[2])
            if seq_id.split("(")[0] == cid_alvo:
                posicoes.append(pos)
                profundidades.append(depth)

    if not posicoes:
        return None

    return pd.DataFrame({"pos": posicoes, "depth": profundidades})


def plotar_locus_eve(
    scaffold: str,
    eve_start: int,
    eve_end: int,
    cluster_id: str,
    flank_start: int,
    flank_end: int,
    eve_strand: str = "+",
    gff_path: str = None,
    depth_df: pd.DataFrame = None,
    output_file: str = None,
    custom_title: str = None,
    species: str = None
):
    """
    Draws the high-resolution figure (PNG and PDF) of the EVE locus within its
    20,000 nt genomic neighbourhood, combining:
    1. An RNA-seq coverage track (when available).
    2. A genomic-architecture track with the host annotation (when available)
       e destaque proeminente da EVE integrada.
    """
    # 1. Load the host annotation
    anotacoes = carregar_anotacoes_gff(gff_path, scaffold, flank_start, flank_end)
    tem_rnaseq = depth_df is not None and not depth_df.empty

    # Build the figure layout
    if tem_rnaseq:
        fig, (ax_rna, ax_gene) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True,
            gridspec_kw={"height_ratios": [1.2, 1.0], "hspace": 0.12}
        )
    else:
        fig, ax_gene = plt.subplots(1, 1, figsize=(14, 4))
        ax_rna = None

    regiao_largura = flank_end - flank_start

    # ==========================
    # TRILHA 1: COBERTURA RNA-SEQ
    # ==========================
    if tem_rnaseq:
        # Map the relative positions of depth_df onto real genome coordinates
        # depth_df['pos'] vai de 1 a N
        pos_genomicas = flank_start + depth_df["pos"].values
        profundidades = depth_df["depth"].values

        ax_rna.plot(pos_genomicas, profundidades, color="#1f4e79", lw=1.2, label="RNA-seq depth")
        ax_rna.fill_between(pos_genomicas, profundidades, color="#4a90e2", alpha=0.35)

        # Highlight the EVE region on the RNA-seq track
        ax_rna.axvspan(eve_start, eve_end, color="#ff7675", alpha=0.25, label="EVE region")
        ax_rna.axvline(eve_start, color="#d63031", linestyle="--", lw=1.0)
        ax_rna.axvline(eve_end, color="#d63031", linestyle="--", lw=1.0)

        max_y = max(5.0, np.max(profundidades) * 1.15) if len(profundidades) > 0 else 5.0
        ax_rna.set_ylim(0, max_y)
        ax_rna.set_ylabel("Depth (reads)", fontsize=10, fontweight="bold")
        ax_rna.grid(True, linestyle=":", alpha=0.5, color="#cccccc")
        ax_rna.legend(loc="upper right", framealpha=0.85, fontsize=8)

        # Main title
        tit = custom_title or f"Locus {cluster_id} | Transcriptional activity & genomic context (20 kb)"
        if species:
            tit = f"[{species}] {tit}"
        ax_rna.set_title(tit, fontsize=12, fontweight="bold", pad=8)
    else:
        tit = custom_title or f"Genomic context of locus {cluster_id} (20 kb flanks)"
        if species:
            tit = f"[{species}] {tit}"
        ax_gene.set_title(tit, fontsize=12, fontweight="bold", pad=12)

    # ==================================
    # TRILHA 2: ARQUITETURA GENÔMICA
    # ==================================
    # Scaffold baseline (chromosomal backbone)
    y_backbone = 0.5
    ax_gene.plot([flank_start, flank_end], [y_backbone, y_backbone], color="#636e72", lw=3.0, zorder=1)

    # Marking of the 20 kb flanks
    ax_gene.axvline(eve_start, color="#d63031", linestyle=":", lw=1.2, alpha=0.8)
    ax_gene.axvline(eve_end, color="#d63031", linestyle=":", lw=1.2, alpha=0.8)

    # Host annotation (when present)
    genes = [f for f in anotacoes if f["type"] in ("gene", "mRNA")]
    exons = [f for f in anotacoes if f["type"] in ("exon", "CDS")]

    for g in genes:
        g_s = max(flank_start, g["start"])
        g_e = min(flank_end, g["end"])
        if g_e <= g_s:
            continue
        g_y = 0.75 if g["strand"] == "+" else 0.25
        cor_gene = "#0984e3" if g["strand"] == "+" else "#6c5ce7"

        # Linha do gene
        ax_gene.plot([g_s, g_e], [g_y, g_y], color=cor_gene, lw=1.5, zorder=2)
        # Gene label
        g_mid = (g_s + g_e) / 2.0
        offset_y = 0.12 if g["strand"] == "+" else -0.12
        ax_gene.text(g_mid, g_y + offset_y, g["name"], ha="center", va="center",
                     fontsize=7, color="#2d3436", style="italic")

    for ex in exons:
        e_s = max(flank_start, ex["start"])
        e_e = min(flank_end, ex["end"])
        if e_e <= e_s:
            continue
        e_y = 0.75 if ex["strand"] == "+" else 0.25
        h_rect = 0.14
        rect = patches.Rectangle(
            (e_s, e_y - h_rect / 2.0), e_e - e_s, h_rect,
            facecolor="#74b9ff", edgecolor="#0984e3", lw=0.8, zorder=3
        )
        ax_gene.add_patch(rect)

    # ==========================
    # ENDOGENOUS VIRAL ELEMENT (EVE)
    # ==========================
    eve_len = eve_end - eve_start
    h_eve = 0.28
    rect_eve = patches.Rectangle(
        (eve_start, y_backbone - h_eve / 2.0), eve_len, h_eve,
        facecolor="#d63031", edgecolor="#2d3436", lw=1.5, zorder=4
    )
    ax_gene.add_patch(rect_eve)

    # Arrow showing the strand orientation of the EVE
    sinal_fita = "→" if eve_strand == "+" else ("←" if eve_strand == "-" else "↔")
    texto_eve = f"EVE: {cluster_id} ({eve_len:,} bp) {sinal_fita}"
    ax_gene.text(
        (eve_start + eve_end) / 2.0, y_backbone + 0.22, texto_eve,
        ha="center", va="bottom", fontsize=9, fontweight="bold", color="#d63031",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#ffeaa7", edgecolor="#fdcb6e", lw=1.0)
    )

    # Text annotation of the flanks
    mid_up = (flank_start + eve_start) / 2.0
    mid_down = (eve_end + flank_end) / 2.0
    ax_gene.text(mid_up, y_backbone - 0.28, "5' flank (20 kb)", ha="center", va="top",
                 fontsize=8, color="#636e72")
    ax_gene.text(mid_down, y_backbone - 0.28, "3' flank (20 kb)", ha="center", va="top",
                 fontsize=8, color="#636e72")

    # Ajustes finais dos eixos
    ax_gene.set_xlim(flank_start, flank_end)
    ax_gene.set_ylim(-0.05, 1.15)
    ax_gene.set_yticks([])
    ax_gene.set_xlabel(f"Genomic coordinates on {scaffold} (bp)", fontsize=10, fontweight="bold")

    # X axis formatted with a thousands separator
    ax_gene.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
    ax_gene.tick_params(axis="x", labelsize=8)

    # Legenda customizada
    legend_elements = [
        patches.Patch(facecolor="#d63031", edgecolor="#2d3436", label="EVE insertion"),
        Line2D([0], [0], color="#636e72", lw=2.5, label="Host scaffold"),
    ]
    if anotacoes:
        legend_elements.extend([
            patches.Patch(facecolor="#74b9ff", edgecolor="#0984e3", label="Host exon (+/-)"),
            Line2D([0], [0], color="#0984e3", lw=1.5, label="Host gene"),
        ])
    if tem_rnaseq:
        legend_elements.append(
            patches.Patch(facecolor="#4a90e2", alpha=0.5, label="RNA-seq reads")
        )

    ax_gene.legend(handles=legend_elements, loc="lower right", framealpha=0.9, fontsize=8)

    try:
        plt.tight_layout()
    except Exception:
        pass

    # Export as PNG and PDF
    if output_file:
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        caminho_png = output_file if output_file.endswith(".png") else f"{output_file}.png"
        caminho_pdf = caminho_png[:-4] + ".pdf"
        plt.savefig(caminho_png, dpi=300, bbox_inches="tight")
        plt.savefig(caminho_pdf, bbox_inches="tight")
        plt.close()
        return caminho_png
    else:
        plt.close()
        return None


def gerar_figuras_contexto_lote(
    bed_original: str,
    bed_flanked: str,
    gff_annotation: str = None,
    depth_tsv: str = None,
    output_dir: str = "figuras_contexto",
    species: str = None,
    flank: int = 20000,
    max_figures: int = 100
) -> list:
    """
    Walks the EVE clusters found and draws a full figure for each locus,
    with 20 kb flanks, host annotation and RNA-seq coverage.
    """
    os.makedirs(output_dir, exist_ok=True)
    figuras_geradas = []

    if not os.path.exists(bed_original) or not os.path.exists(bed_flanked):
        print(f"[warning] BED files for plotting not found: {bed_original} or {bed_flanked}")
        return figuras_geradas

    # Original clusters
    df_orig = pd.read_csv(bed_original, sep="\t", header=None)
    # Flanked clusters
    df_flank = pd.read_csv(bed_flanked, sep="\t", header=None)

    mapa_flancos = {}
    for _, row in df_flank.iterrows():
        cid = str(row[3])
        mapa_flancos[cid] = (int(row[1]), int(row[2]))

    total = min(len(df_orig), max_figures)
    print(f"[Plot] Drawing genomic-context figures for {total} loci...")

    for i, (_, row) in enumerate(df_orig.iterrows()):
        if i >= max_figures:
            break

        scaff = str(row[0])
        eve_s = int(row[1])
        eve_e = int(row[2])
        cid = str(row[3])
        strand = str(row[5]) if len(row) > 5 else "+"

        if cid in mapa_flancos:
            flank_s, flank_e = mapa_flancos[cid]
        else:
            flank_s = max(0, eve_s - flank)
            flank_e = eve_e + flank

        # Load the depth data when depth_tsv exists
        df_depth = carregar_profundidade_locus(depth_tsv, cid) if depth_tsv else None

        caminho_fig = os.path.join(output_dir, f"{cid}_contexto_20kb.png")
        plotar_locus_eve(
            scaffold=scaff,
            eve_start=eve_s,
            eve_end=eve_e,
            cluster_id=cid,
            flank_start=flank_s,
            flank_end=flank_e,
            eve_strand=strand,
            gff_path=gff_annotation,
            depth_df=df_depth,
            output_file=caminho_fig,
            species=species
        )
        figuras_geradas.append(caminho_fig)

    print(f"[Plot] Done. {len(figuras_geradas)} figures saved to: {output_dir}")
    return figuras_geradas
