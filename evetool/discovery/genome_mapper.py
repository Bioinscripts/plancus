import os
import csv
import subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd


def preparar_banco_blast(genoma_fasta: str, diretorio_saida: str) -> str:
    """
    Format the host genome as a BLAST+ nucleotide database (makeblastdb). The same
    database is reused by the tblastn searches.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    nome_base = Path(genoma_fasta).stem
    db_path = os.path.join(diretorio_saida, f"{nome_base}_blastdb")

    comando_makeblastdb = [
        "makeblastdb",
        "-in", genoma_fasta,
        "-dbtype", "nucl",
        "-out", db_path
    ]

    try:
        subprocess.run(comando_makeblastdb, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(f"Could not format the BLAST database for {genoma_fasta}.\nDetails: {erro.stderr}")

    return db_path


def rodar_tblastn(
    fasta_query: str,
    db_path: str,
    diretorio_saida: str,
    nome_tabela_saida: str = "eve_tblastn_hits.tsv",
    threads: int = 4,
    evalue: float = 1e-5
) -> str:
    """
    Run tblastn of fasta_query (proteins) against the nucleotide genome. Exports the 11
    columns that feed the harmonic score: qseqid, sseqid, sstart, send, evalue, bitscore,
    pident, length, qlen, slen, sframe.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    tabela_saida = os.path.join(diretorio_saida, nome_tabela_saida)

    comando_tblastn = [
        "tblastn",
        "-query", fasta_query,
        "-db", db_path,
        "-out", tabela_saida,
        "-outfmt", "6 qseqid sseqid sstart send evalue bitscore pident length qlen slen sframe",
        "-evalue", str(evalue),
        "-num_threads", str(threads)
    ]

    try:
        subprocess.run(comando_tblastn, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(f"tblastn failed for {fasta_query}.\nDetails: {erro.stderr}")

    return tabela_saida


def mesclar_tabelas_tblastn(tabelas: list, arquivo_saida: str) -> str:
    """
    Concatenate several tblastn hit tables into one, skipping empty or missing files.
    """
    os.makedirs(os.path.dirname(arquivo_saida) or '.', exist_ok=True)
    linhas_escritas = 0

    with open(arquivo_saida, 'w') as f_out:
        for tabela in tabelas:
            if not tabela or not os.path.exists(tabela) or os.path.getsize(tabela) == 0:
                continue
            with open(tabela, 'r') as f_in:
                for linha in f_in:
                    if linha.strip():
                        f_out.write(linha)
                        linhas_escritas += 1

    if linhas_escritas == 0:
        raise RuntimeError("None of the tblastn tables held any hit; nothing to merge.")

    return arquivo_saida


def calcular_harmonic_score(pident: float, coverage: float) -> float:
    """Harmonic mean of percent identity and coverage."""
    if pident <= 0 or coverage <= 0:
        return 0.0
    return (2.0 * pident * coverage) / (pident + coverage)


def merge_intervals(frags: list, max_gap: int = 300) -> list:
    """
    Clustering level 1: merge adjacent or overlapping HSPs (gap <= max_gap, 300 nt by
    default) on the same scaffold into a single EVE.
    """
    frags = sorted(frags, key=lambda f: f["gmin"])
    groups = []
    if not frags:
        return groups

    current = [frags[0]]
    for fragment in frags[1:]:
        current_end = max(item["gmax"] for item in current)
        gap = fragment["gmin"] - current_end
        if gap <= max_gap:
            current.append(fragment)
        else:
            groups.append(current)
            current = [fragment]
    groups.append(current)
    return groups


def gerar_gff_clusters(
    tabela_tblastn: str,
    arquivo_gff_saida: str,
    merge_gap: int = 300,
    cluster_gap: int = 10000,
    min_pident: float = 25.0,
    min_length: int = 30,
    min_qcov: float = 20.0,
    min_score: float = 30.0,
    max_evalue: float = 1e-5
) -> str:
    """
    Filter HSPs by harmonic score and cluster in two levels:
      - level 1 (merge_gap=300 nt): merge adjacent HSPs into one EVE;
      - level 2 (cluster_gap=10000 nt): group nearby EVEs on the same scaffold into macro-clusters.
    Writes the result as GFF3 plus detailed TSV tables.
    """
    if not os.path.exists(tabela_tblastn) or os.path.getsize(tabela_tblastn) == 0:
        raise RuntimeError(f"The tblastn table '{tabela_tblastn}' is empty — no hit found in the genome.")

    os.makedirs(os.path.dirname(arquivo_gff_saida) or '.', exist_ok=True)
    diretorio_saida = os.path.dirname(arquivo_gff_saida) or '.'

    frags_by_contig = defaultdict(list)
    total_hsps = 0
    accepted_hsps = 0

    # 1. Read the HSPs and filter them by harmonic score
    with open(tabela_tblastn, 'r') as fh:
        for linha in fh:
            if not linha.strip() or linha.startswith("#"):
                continue
            total_hsps += 1
            partes = linha.strip().split("\t")
            if len(partes) < 8:
                continue

            qseqid = partes[0]
            sseqid = partes[1]
            sstart = int(partes[2])
            send = int(partes[3])
            evalue = float(partes[4])
            bitscore = float(partes[5])
            pident = float(partes[6])
            length = int(partes[7])

            # Supports 11-column tables (qlen, slen, sframe) or the 8-column fallback
            if len(partes) >= 11:
                qlen = int(partes[8])
                slen = int(partes[9])
                sframe = int(partes[10])
                strand = "+" if sframe > 0 else "-"
                qcov = (length / max(1, qlen)) * 100.0
            else:
                strand = "+" if sstart < send else "-"
                qcov = min(100.0, (length / 100.0) * 100.0)

            harmonic_score = calcular_harmonic_score(pident, qcov)

            # Strict quality filters
            if evalue > max_evalue:
                continue
            if pident < min_pident:
                continue
            if length < min_length:
                continue
            if qcov < min_qcov:
                continue
            if harmonic_score < min_score:
                continue

            accepted_hsps += 1
            gmin, gmax = min(sstart, send), max(sstart, send)

            frags_by_contig[sseqid].append({
                "gmin": gmin,
                "gmax": gmax,
                "strand": strand,
                "pident": pident,
                "length": length,
                "qcov": qcov,
                "harmonic_score": harmonic_score,
                "evalue": evalue,
                "bitscore": bitscore,
                "qseqid": qseqid
            })

    if not frags_by_contig:
        raise RuntimeError(
            f"No tblastn hit survived the strict filters "
            f"(total={total_hsps}, accepted=0, min_score={min_score}, min_pident={min_pident}%)."
        )

    print(f"[Clustering] HSPs evaluated: {total_hsps} | passed the harmonic score: {accepted_hsps}")

    # 2. Two-level processing (locus merge and macro-clusters)
    eve_rows = []
    cluster_rows = []
    gff_lines = ["##gff-version 3"]

    for contig, fragments in frags_by_contig.items():
        # Level 1: merge fragments of the same EVE with gap <= merge_gap
        grupos_hsp = merge_intervals(fragments, max_gap=merge_gap)
        eves_contig = []

        for grp in grupos_hsp:
            gmin = min(f["gmin"] for f in grp)
            gmax = max(f["gmax"] for f in grp)
            strands = [f["strand"] for f in grp]
            strand = max(set(strands), key=strands.count)
            queries = sorted(set(f["qseqid"] for f in grp))

            total_len = sum(f["length"] for f in grp)
            weighted_pident = sum(f["pident"] * f["length"] for f in grp) / max(1, total_len)
            weighted_harmonic = sum(f["harmonic_score"] * f["length"] for f in grp) / max(1, total_len)

            eves_contig.append({
                "contig": contig,
                "gmin": gmin,
                "gmax": gmax,
                "strand": strand,
                "n_fragments": len(grp),
                "queries": queries,
                "aligned_length": total_len,
                "pident_weighted": weighted_pident,
                "harmonic_score_weighted": weighted_harmonic,
                "best_evalue": min(f["evalue"] for f in grp),
                "best_bitscore": max(f["bitscore"] for f in grp)
            })

        eves_contig.sort(key=lambda e: e["gmin"])

        # Level 2: macro-clustering of EVEs with gap <= cluster_gap
        cluster_groups = []
        if eves_contig:
            current = [0]
            for i in range(1, len(eves_contig)):
                current_end = max(eves_contig[j]["gmax"] for j in current)
                if eves_contig[i]["gmin"] - current_end <= cluster_gap:
                    current.append(i)
                else:
                    cluster_groups.append(current)
                    current = [i]
            cluster_groups.append(current)

        for cluster_idx, indices in enumerate(cluster_groups, start=1):
            cluster_id = f"{contig}_cluster_{cluster_idx}"
            c_start = min(eves_contig[i]["gmin"] for i in indices)
            c_end = max(eves_contig[i]["gmax"] for i in indices)

            for i in indices:
                eves_contig[i]["cluster_id"] = cluster_id
                eves_contig[i]["cluster_size"] = len(indices)

            cluster_rows.append({
                "cluster_id": cluster_id,
                "contig": contig,
                "cluster_start": c_start,
                "cluster_end": c_end,
                "n_eves": len(indices)
            })

        for n, eve in enumerate(eves_contig, start=1):
            eve_id = f"{contig}_eve_{n}"
            eve_rows.append({
                "eve_id": eve_id,
                "cluster_id": eve["cluster_id"],
                "contig": contig,
                "start": eve["gmin"],
                "end": eve["gmax"],
                "strand": eve["strand"],
                "n_fragments": eve["n_fragments"],
                "viral_queries": ";".join(eve["queries"]),
                "aligned_length": eve["aligned_length"],
                "pident_weighted": round(eve["pident_weighted"], 2),
                "harmonic_score_weighted": round(eve["harmonic_score_weighted"], 2),
                "best_evalue": eve["best_evalue"],
                "best_bitscore": eve["best_bitscore"]
            })

            # Formatted GFF3 line
            attrs = (
                f"ID={eve['cluster_id']};Name={eve_id};eve_cluster={eve['cluster_id']};"
                f"harmonic_score={eve['harmonic_score_weighted']:.2f};pident={eve['pident_weighted']:.2f};"
                f"best_evalue={eve['best_evalue']:.2e};Members={','.join(eve['queries'][:5])}"
            )
            gff_lines.append("\t".join([
                contig, "EVE_Pipeline", "EVE_region", str(eve["gmin"]), str(eve["gmax"]),
                f"{eve['harmonic_score_weighted']:.2f}", eve["strand"], ".", attrs
            ]))

    # 3. Write the GFF3 and TSV files
    with open(arquivo_gff_saida, 'w') as f_gff:
        f_gff.write("\n".join(gff_lines) + "\n")

    caminho_eve_tsv = os.path.join(diretorio_saida, "eve_loci.tsv")
    pd.DataFrame(eve_rows).to_csv(caminho_eve_tsv, sep="\t", index=False)

    caminho_cluster_tsv = os.path.join(diretorio_saida, "eve_clusters.tsv")
    pd.DataFrame(cluster_rows).to_csv(caminho_cluster_tsv, sep="\t", index=False)

    print(f"[Clustering] Found {len(eve_rows)} EVE loci grouped into {len(cluster_rows)} clusters.")
    print(f"[Clustering] Files written: {arquivo_gff_saida}, {caminho_eve_tsv}")
    return arquivo_gff_saida


def converter_gff_para_bed(gff_path: str, bed_path: str) -> str:
    """
    Convert the EVE cluster GFF3 into the 0-based BED that bedtools slop and getfasta expect.
    """
    os.makedirs(os.path.dirname(bed_path) or '.', exist_ok=True)
    linhas_escritas = 0

    with open(gff_path, 'r') as f_gff, open(bed_path, 'w') as f_bed:
        for linha in f_gff:
            if linha.startswith("#") or not linha.strip():
                continue

            campos = linha.rstrip("\n").split("\t")
            if len(campos) < 9:
                continue

            scaffold, _fonte, _tipo, start, end, score, strand, _frame, atributos = campos[:9]

            id_cluster = "EVE_cluster"
            for par in atributos.split(";"):
                if par.startswith("ID="):
                    id_cluster = par[3:]
                    break

            inicio_bed = max(0, int(start) - 1)
            fim_bed = int(end)
            strand_bed = strand if strand in ("+", "-") else "."

            f_bed.write(f"{scaffold}\t{inicio_bed}\t{fim_bed}\t{id_cluster}\t{score}\t{strand_bed}\n")
            linhas_escritas += 1

    if linhas_escritas == 0:
        raise RuntimeError(f"No EVE cluster found in '{gff_path}'; nothing to convert to BED.")

    return bed_path
