import os
import re
import sys
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO


def gerar_tamanhos_cromossomos(genoma_fasta: str, arquivo_tamanhos: str) -> dict:
    """
    Write the chromosome sizes file required by bedtools slop. Uses samtools faidx when
    available, otherwise BioPython. Returns a {contig: length} dictionary.
    """
    tamanhos = {}
    fai_path = f"{genoma_fasta}.fai"
    
    # 1. Use or create the .fai index with samtools faidx
    if not os.path.exists(fai_path):
        try:
            subprocess.run(["samtools", "faidx", genoma_fasta], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except Exception:
            pass

    if os.path.exists(fai_path):
        with open(fai_path, "r") as fai:
            for line in fai:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    tamanhos[parts[0]] = int(parts[1])
    else:
        # BioPython fallback
        for rec in SeqIO.parse(genoma_fasta, "fasta"):
            tamanhos[rec.id.split()[0]] = len(rec.seq)

    with open(arquivo_tamanhos, "w") as out:
        for contig, size in tamanhos.items():
            out.write(f"{contig}\t{size}\n")

    return tamanhos


def agrupar_e_extrair_contexto(arquivo_bed: str, genoma_fasta: str, diretorio_saida: str, flanco: int = 20000) -> dict:
    """
    Extract the exact EVE cluster sequences and their flanking genomic context (e.g. 20 kb).
    
    Two modes are selected automatically:
      1. Whole-genome mode: uses samtools and bedtools (slop + getfasta) to slice the host
         chromosomes/scaffolds directly and take the real flanking sequence.
      2. Fragment-collection mode: when the input is a set of pre-extracted segments (headers
         like 'contig:start-end(...)'), the matching sequences are aggregated per cluster so
         the downstream steps can still run.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    
    tamanhos_path = os.path.join(diretorio_saida, "cromossomos.sizes")
    fasta_clusters_exact = os.path.join(diretorio_saida, "eve_clusters_exact.fasta")
    bed_flanqueado = os.path.join(diretorio_saida, "eve_clusters_flanked.bed")
    fasta_contexto = os.path.join(diretorio_saida, "eve_clusters_context.fasta")
    
    tamanhos_genoma = gerar_tamanhos_cromossomos(genoma_fasta, tamanhos_path)
    
    # Contigs present in the BED file
    contigs_bed = set()
    if os.path.exists(arquivo_bed):
        with open(arquivo_bed, "r") as f_bed:
            for line in f_bed:
                if line.strip() and not line.startswith("#"):
                    contigs_bed.add(line.strip().split("\t")[0])

    # Check whether the BED contigs match the FASTA contigs (whole-genome mode)
    contigs_em_comum = contigs_bed.intersection(tamanhos_genoma.keys())
    eh_genoma_completo = len(contigs_em_comum) > 0 and (len(contigs_em_comum) >= len(contigs_bed) * 0.5)

    if eh_genoma_completo:
        print(f"[Context] Whole genome detected ({len(contigs_em_comum)} matching contigs).")
        try:
            # 1. Extract the exact EVE sequences from the genome
            cmd_exact = [
                "bedtools", "getfasta",
                "-fi", genoma_fasta,
                "-bed", arquivo_bed,
                "-name",
                "-fo", fasta_clusters_exact
            ]
            subprocess.run(cmd_exact, check=True)

            # 2. Build the flanked BED with upstream and downstream margins
            cmd_slop = f"bedtools slop -i {arquivo_bed} -g {tamanhos_path} -b {flanco} > {bed_flanqueado}"
            subprocess.run(cmd_slop, shell=True, check=True)

            # 3. Extract the full genomic context sequences (EVE + flanks)
            cmd_context = [
                "bedtools", "getfasta",
                "-fi", genoma_fasta,
                "-bed", bed_flanqueado,
                "-name",
                "-fo", fasta_contexto
            ]
            subprocess.run(cmd_context, check=True)

            print(f"[Context] Exact sequences at: {fasta_clusters_exact}")
            print(f"[Context] Genomic context ({flanco} bp) at: {fasta_contexto}")

            return {
                "bed_final": arquivo_bed,
                "bed_flanqueado": bed_flanqueado,
                "fasta_contexto": fasta_contexto,
                "fasta_exato": fasta_clusters_exact,
                "tamanhos": tamanhos_path
            }
        except Exception as e:
            print(f"[warning] Direct bedtools extraction failed: {e}. Falling back...")

    # =========================================================================
    # Fragment-collection mode / fallback
    # =========================================================================
    print(f"[Context] Running in fragment mode / ORF aggregation...")

    seq_dict = {}
    for rec in SeqIO.parse(genoma_fasta, "fasta"):
        clean_id = rec.id.split()[0]
        seq_str = str(rec.seq)
        seq_dict[rec.id] = seq_str
        seq_dict[clean_id] = seq_str
        # Strand-free variants, for flexible matching
        if "(" in clean_id:
            seq_dict[clean_id.split("(")[0]] = seq_str

    gff_path = os.path.join(diretorio_saida, "eve_clusters.gff3")
    clusters_seqs = {}
    current_cluster = None

    if os.path.exists(gff_path):
        with open(gff_path, "r") as fgff:
            for line in fgff:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature = parts[2]
                attributes = parts[8]

                if feature == "region":
                    for attr in attributes.split(";"):
                        if attr.startswith("Name="):
                            current_cluster = attr.split("=")[1]
                            if current_cluster not in clusters_seqs:
                                clusters_seqs[current_cluster] = []
                elif feature == "CDS" and current_cluster:
                    orf_id = None
                    for attr in attributes.split(";"):
                        if attr.startswith("ID="):
                            orf_id = attr.split("=")[1]
                            break

                    if orf_id:
                        # Recover the base id by dropping the '_ORF.X' suffix orfipy adds
                        base_id = orf_id.rsplit("_ORF.", 1)[0]
                        if base_id in seq_dict:
                            clusters_seqs[current_cluster].append(seq_dict[base_id])
                        elif orf_id in seq_dict:
                            clusters_seqs[current_cluster].append(seq_dict[orf_id])
                        else:
                            # Fall back to prefix matching
                            achou = False
                            for k, v in seq_dict.items():
                                if orf_id.startswith(k) or k.startswith(base_id):
                                    clusters_seqs[current_cluster].append(v)
                                    achou = True
                                    break

    # No clusters from the GFF: fall back to the input BED
    if not clusters_seqs and os.path.exists(arquivo_bed):
        with open(arquivo_bed, "r") as f_bed:
            for line in f_bed:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                c_name = parts[3] if len(parts) >= 4 else "Cluster_1"
                contig_key = parts[0]
                if contig_key in seq_dict:
                    clusters_seqs.setdefault(c_name, []).append(seq_dict[contig_key])

    # Write the exact-cluster FASTA
    with open(fasta_clusters_exact, "w") as f_out:
        for c_name, seqs in clusters_seqs.items():
            merged_seq = "N".join(seqs) if seqs else "N"
            f_out.write(f">{c_name}\n{merged_seq}\n")

    # Make sure the flanked BED exists for the downstream tools
    if os.path.exists(arquivo_bed):
        shutil.copy(arquivo_bed, bed_flanqueado)
    else:
        with open(bed_flanqueado, "w") as bf:
            bf.write("")

    shutil.copy(fasta_clusters_exact, fasta_contexto)

    return {
        "bed_final": arquivo_bed,
        "bed_flanqueado": bed_flanqueado,
        "fasta_contexto": fasta_contexto,
        "fasta_exato": fasta_clusters_exact,
        "tamanhos": tamanhos_path
    }
