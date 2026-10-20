import os
import csv
import subprocess
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


from evetool.discovery.ncbi_genome_downloader import e_acesso_ncbi, baixar_genoma_ncbi
from evetool.annotation.rnaseq_coverage import is_bam


def papel_coluna(nome: str) -> str | None:
    """Role of a species-TSV header column, from keywords in its name (None = unknown)."""
    n = nome.strip().lower().replace("-", "_").replace(" ", "_")
    if "bam" in n:
        return "bam"
    if "gff" in n or "gtf" in n or "annot" in n:
        return "gff"
    if n.endswith(("r2", "_2")) or "read2" in n or "fastq_2" in n:
        return "r2"
    if n.endswith(("r1", "_1")) or "read1" in n or "fastq" in n or n in ("reads", "rnaseq"):
        return "r1"
    if "sint" in n or "synten" in n:
        return "sintenia"
    if n in ("sp", "species", "specie", "especie", "espécie", "especies", "organism", "organismo", "name", "nome") \
            or "species" in n or "especie" in n:
        return "species"
    if n in ("genome", "genoma", "accession", "acesso", "gca", "gcf", "fasta", "assembly", "code",
             "codigo", "código", "genome_fasta") or "genom" in n or "accession" in n or "assembly" in n:
        return "genome"
    return None


def carregar_lista_especies(tsv_path: str, pasta_genomas: str = None) -> list:
    """
    Read the TSV listing the species and their files.
    The column layout is detected automatically and supports:
      - a file with or without a header;
      - column 1 = GCA_/GCF_ (or a FASTA file) and column 2 = species name;
      - column 1 = species name and column 2 = GCA_/GCF_ (or a FASTA file);
      - optional extra columns: gff, rnaseq_r1, rnaseq_r2, sintenia (by position without a
        header; by name with one, in any order - see papel_coluna);
      - RNA-seq BAMs: a "bam" column, or any cell holding a BAM file (comma-separated for
        several samples), recognised by extension or content.
    """
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"Species TSV not found: {tsv_path}")

    pasta_base = pasta_genomas or os.path.dirname(tsv_path) or "."
    especies = []

    with open(tsv_path, "r", encoding="utf-8") as f:
        linhas = [l.strip() for l in f if l.strip() and not l.strip().startswith("#")]

    if not linhas:
        raise ValueError(f"The file '{tsv_path}' is empty.")

    # 1. Check whether the first line is a header
    primeira_linha = linhas[0]
    tokens_primeira = [t.strip() for t in primeira_linha.split("\t" if "\t" in primeira_linha else None)]
    tokens_primeira_lower = [t.lower() for t in tokens_primeira]

    # Each header column is given a role by keyword, so the columns can come in any order
    # and with free names ("code", "Specie", "Annotation gff", "BAM_file"...).
    papeis_cabecalho = [papel_coluna(t) for t in tokens_primeira_lower]
    tem_cabecalho = any(papeis_cabecalho)

    # A first column starting with GCA_ or GCF_ is definitely not a header
    if e_acesso_ncbi(tokens_primeira[0]):
        tem_cabecalho = False

    linhas_dados = linhas[1:] if tem_cabecalho else linhas

    # Column of each role when a header is present (first column wins)
    mapa_colunas = {}
    if tem_cabecalho:
        for idx, papel in enumerate(papeis_cabecalho):
            if papel and papel not in mapa_colunas:
                mapa_colunas[papel] = idx

    vazio = ("none", "na", "n/a", "", "-", "auto")

    for num_linha, linha in enumerate(linhas_dados, start=1 if not tem_cabecalho else 2):
        tokens = [t.strip() for t in linha.split("\t" if "\t" in linha else None)]
        if len(tokens) < 2:
            print(f"[warning] Line {num_linha} skipped: fewer than 2 columns: '{linha}'")
            continue

        def celula(papel, posicao):
            idx = mapa_colunas.get(papel, None if tem_cabecalho else posicao)
            if idx is None or idx >= len(tokens) or tokens[idx].lower() in vazio:
                return None
            return tokens[idx]

        tok0, tok1 = tokens[0], tokens[1]
        caminho_gff = celula("gff", 2)
        caminho_r1 = celula("r1", 3)
        caminho_r2 = celula("r2", 4)
        # Optional synteny switch ("sim"/"nao"), which the researcher can set line by line
        # in the TSV instead of for the whole analysis.
        tok_sint = (celula("sintenia", 5) or "").lower()
        sintenia = tok_sint in ("sim", "yes", "true", "1", "s", "y")
        # RNA-seq already aligned: the BAM column, plus any cell that is a BAM (by the .bam
        # extension or by its content, for files saved without extension); comma-separated
        # for several samples. A BAM never goes to the aligner as if it were FASTQ.
        def caminho_abs(c):
            return c if os.path.isabs(c) else os.path.join(pasta_base, c)
        celulas_bam = [celula("bam", None)] if "bam" in mapa_colunas else []
        celulas_bam += [t for i, t in enumerate(tokens[2:], 2) if i != mapa_colunas.get("bam")]
        caminhos_bam = []
        for c in celulas_bam:
            for b in (c or "").split(","):
                b = b.strip()
                if b and b.lower() not in vazio and (b.lower().endswith(".bam") or is_bam(caminho_abs(b))):
                    caminhos_bam.append(b)
        e_bam = lambda c: bool(c) and (c.lower().endswith(".bam") or is_bam(caminho_abs(c)))
        caminho_gff = None if e_bam(caminho_gff) else caminho_gff
        caminho_r1 = None if e_bam(caminho_r1) else caminho_r1
        caminho_r2 = None if e_bam(caminho_r2) else caminho_r2

        # 2. Heuristic to tell which column is the genome/accession and which is the species
        if e_acesso_ncbi(tok0) or tok0.endswith((".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz")):
            # Column 1 = accession/genome, column 2 = species
            gen_arq = tok0
            sp_nome = tok1
        elif e_acesso_ncbi(tok1) or tok1.endswith((".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz")):
            # Column 1 = species, column 2 = accession/genome
            sp_nome = tok0
            gen_arq = tok1
        elif tem_cabecalho:
            idx_sp = mapa_colunas.get("species", 0)
            idx_gen = mapa_colunas.get("genome", 1)
            sp_nome = tokens[idx_sp]
            gen_arq = tokens[idx_gen]
        else:
            # Default: column 1 = genome/accession, column 2 = species
            gen_arq = tok0
            sp_nome = tok1

        sp_nome = sp_nome.strip().replace(" ", "_")

        # 3. For an NCBI accession (GCA_ or GCF_), download it through NCBI Datasets
        if e_acesso_ncbi(gen_arq):
            res_ncbi = baixar_genoma_ncbi(gen_arq, pasta_destino=pasta_base, baixar_gff=True)
            caminho_gen = res_ncbi["genome"]
            if not caminho_gff and res_ncbi.get("gff"):
                caminho_gff = res_ncbi["gff"]
        else:
            caminho_gen = gen_arq if os.path.isabs(gen_arq) else os.path.join(pasta_base, gen_arq)
            if not os.path.exists(caminho_gen):
                candidato = os.path.join(pasta_base, os.path.basename(gen_arq))
                if os.path.exists(candidato):
                    caminho_gen = candidato
                else:
                    raise FileNotFoundError(f"Genome for '{sp_nome}' not found at: {caminho_gen}")

        if caminho_gff and not os.path.isabs(caminho_gff):
            caminho_gff = os.path.join(pasta_base, caminho_gff)

        if caminho_r1 and not os.path.isabs(caminho_r1):
            caminho_r1 = os.path.join(pasta_base, caminho_r1)

        if caminho_r2 and not os.path.isabs(caminho_r2):
            caminho_r2 = os.path.join(pasta_base, caminho_r2)

        caminhos_bam = list(dict.fromkeys(caminho_abs(b) for b in caminhos_bam))
        for b in caminhos_bam:
            if not os.path.exists(b):
                print(f"[warning] {sp_nome}: RNA-seq BAM not found: {b}")

        especies.append({
            "species": sp_nome,
            "genome": caminho_gen,
            "gff": caminho_gff,
            "rnaseq_r1": caminho_r1,
            "rnaseq_r2": caminho_r2,
            "rnaseq_bam": caminhos_bam,
            "sintenia": sintenia
        })

    print(f"[MultiSpecies] Loaded {len(especies)} species from: {tsv_path}")
    return especies


def consolidar_eves_multiespecies(mapa_resultados_especies: dict, diretorio_saida: str) -> str:
    """
    Merge the candidate EVEs found in each species into a single multi-species FASTA,
    adding the species prefix to the header (e.g. >Species|ID).
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    fasta_consolidado = os.path.join(diretorio_saida, "todas_especies_eves_unificadas.fasta")

    total_seqs = 0
    with open(fasta_consolidado, "w") as f_out:
        for sp_nome, arq_fasta in mapa_resultados_especies.items():
            if not arq_fasta or not os.path.exists(arq_fasta) or os.path.getsize(arq_fasta) == 0:
                continue
            for record in SeqIO.parse(arq_fasta, "fasta"):
                id_modificado = f"{sp_nome}|{record.id}"
                record.id = id_modificado
                record.description = ""
                SeqIO.write(record, f_out, "fasta")
                total_seqs += 1

    print(f"[MultiSpecies] {total_seqs} EVEs consolidated into: {fasta_consolidado}")
    return fasta_consolidado


def comparar_clusters_entre_especies(
    fasta_consolidado: str,
    diretorio_saida: str,
    threads: int = 4,
    evalue: float = 1e-5
) -> dict:
    """
    Run an all-vs-all tblastx search and a graph analysis (NetworkX) to group homologous
    EVEs shared between different species.
    It produces:
      1. a macrocluster summary table with ancestral vs species-specific status;
      2. a presence/absence matrix with the copy number per species;
      3. a comparative sharing heatmap.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    db_path = os.path.join(diretorio_saida, "banco_multiespecies_eve")
    tsv_out = os.path.join(diretorio_saida, "all_vs_all_tblastx_multiespecies.tsv")

    # 1. Build the BLAST database
    print("[MultiSpecies] Formatting the all-vs-all tblastx database...")
    cmd_db = ["makeblastdb", "-in", fasta_consolidado, "-dbtype", "nucl", "-out", db_path]
    subprocess.run(cmd_db, check=True, stdout=subprocess.DEVNULL)

    # 2. Run the all-vs-all tblastx
    print(f"[MultiSpecies] Running all-vs-all tblastx (threads={threads})...")
    cmd_blast = [
        "tblastx", "-query", fasta_consolidado, "-db", db_path,
        "-out", tsv_out, "-evalue", str(evalue), "-num_threads", str(threads),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore"
    ]
    subprocess.run(cmd_blast, check=True)

    # 3. Graph theory to define the multi-species homology clusters
    print("[MultiSpecies] Building the homology network with NetworkX...")
    grafo = nx.Graph()
    seqs = {seq.id: seq for seq in SeqIO.parse(fasta_consolidado, "fasta")}
    grafo.add_nodes_from(seqs.keys())

    with open(tsv_out, "r") as f:
        for linha in f:
            qseqid, sseqid, pident, length, ev, bs = linha.strip().split("\t")
            if qseqid != sseqid:
                grafo.add_edge(qseqid, sseqid, pident=float(pident), bitscore=float(bs))

    componentes = list(nx.connected_components(grafo))
    print(f"[MultiSpecies] Found {len(componentes)} homology clusters across species.")

    # 4. Compile the cluster metrics
    lista_resumo = []
    todas_especies = set()
    for seq_id in seqs.keys():
        sp = seq_id.split("|")[0]
        todas_especies.add(sp)

    lista_todas_especies = sorted(todas_especies)
    matriz_linhas = []

    for i, comp in enumerate(componentes, start=1):
        cid = f"MacroCluster_{i}"
        membros = list(comp)
        especies_presentes = set(m.split("|")[0] for m in membros)
        num_sp = len(especies_presentes)
        num_seqs = len(membros)

        status = "Ancestral / Shared" if num_sp >= 2 else "Species-specific"

        lista_resumo.append({
            "Cluster_ID": cid,
            "Total_Sequences": num_seqs,
            "Total_Species": num_sp,
            "Evolutionary_Status": status,
            "Especies": ";".join(sorted(especies_presentes)),
            "Example_Members": ";".join(membros[:5])
        })

        # Per-species counts for the presence/absence matrix
        linha_matriz = {"Cluster_ID": cid, "Status": status, "Total_Copies": num_seqs}
        for sp in lista_todas_especies:
            copias = sum(1 for m in membros if m.startswith(f"{sp}|"))
            linha_matriz[sp] = copias
        matriz_linhas.append(linha_matriz)

    df_resumo = pd.DataFrame(lista_resumo)
    tsv_resumo = os.path.join(diretorio_saida, "cross_species_clusters_summary.tsv")
    df_resumo.to_csv(tsv_resumo, sep="\t", index=False)

    df_matriz = pd.DataFrame(matriz_linhas)
    csv_matriz = os.path.join(diretorio_saida, "cross_species_presence_absence_matrix.csv")
    df_matriz.to_csv(csv_matriz, index=False)

    # 5. Sharing heatmap
    print("[MultiSpecies] Drawing the comparative heatmap across species...")
    heatmap_png = os.path.join(diretorio_saida, "cross_species_sharing_heatmap.png")
    heatmap_pdf = os.path.join(diretorio_saida, "cross_species_sharing_heatmap.pdf")

    # Pick the clusters to plot (shared and larger ones first)
    df_plot = df_matriz.set_index("Cluster_ID")[lista_todas_especies]
    # Sort by how many species carry the cluster
    df_plot["n_spp"] = (df_plot > 0).sum(axis=1)
    df_plot = df_plot.sort_values(by="n_spp", ascending=False).drop(columns=["n_spp"])

    # Cap at 50 clusters to stay legible when there are hundreds
    top_plot = df_plot.head(50)

    if not top_plot.empty:
        fig_h = max(6, len(top_plot) * 0.28)
        fig_w = max(8, len(lista_todas_especies) * 1.5)
        plt.figure(figsize=(fig_w, fig_h))

        cmap = sns.color_palette("YlOrRd", as_cmap=True)
        sns.heatmap(
            top_plot,
            annot=True,
            fmt="d",
            cmap=cmap,
            cbar_kws={"label": "EVE copy number"},
            linewidths=0.5,
            linecolor="#e0e0e0"
        )
        plt.title("EVE family sharing across host species", fontsize=14, fontweight="bold", pad=15)
        plt.xlabel("Host species", fontsize=11, fontweight="bold")
        plt.ylabel("EVE clusters", fontsize=11, fontweight="bold")
        plt.xticks(rotation=45, ha="right", fontsize=10)
        plt.tight_layout()
        plt.savefig(heatmap_png, dpi=300, bbox_inches="tight")
        plt.savefig(heatmap_pdf, bbox_inches="tight")
        plt.close()

    print(f"[MultiSpecies] Comparative analysis finished:")
    print(f"  - summary table: {tsv_resumo}")
    print(f"  - presence/absence matrix: {csv_matriz}")
    print(f"  - heatmap: {heatmap_png}")

    return {
        "resumo_tsv": tsv_resumo,
        "matriz_csv": csv_matriz,
        "heatmap_png": heatmap_png,
        "tblastx_tsv": tsv_out,
        "componentes": componentes
    }
