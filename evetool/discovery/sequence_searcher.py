import os
import subprocess
from glob import glob

def preparar_e_buscar(pasta_entrada: str, is_fasta: bool, orfs_input: str, diretorio_saida: str,
                      threads: int = 4, hits_por_query: int = 50) -> tuple:
    os.makedirs(diretorio_saida, exist_ok=True)
    
    if is_fasta:
        # 1. Find FASTAs in the folder and build the database on the fly
        fastas = glob(os.path.join(pasta_entrada, "*.fasta")) + glob(os.path.join(pasta_entrada, "*.faa"))
        if not fastas:
            raise FileNotFoundError("No FASTA file found in --fasta-dir.")
        
        # With several FASTAs the right move would be to concatenate them first;
        # here we assume the user supplied a single consolidated FASTA.
        db_path = os.path.join(diretorio_saida, "banco_customizado")
        comando_makedb = ["diamond", "makedb", "--in", fastas[0], "-d", db_path]
        subprocess.run(comando_makedb, check=True)
    else:
        # 2. Locate the existing .dmnd database. Accepts either a direct path to a
        # .dmnd file (the usual --dmnd-db case) or a folder holding one or more.
        if os.path.isfile(pasta_entrada):
            if not pasta_entrada.endswith(".dmnd"):
                raise ValueError(f"--dmnd-db must point to a .dmnd file: {pasta_entrada}")
            banco_dmnd = pasta_entrada
        elif os.path.isdir(pasta_entrada):
            bancos = glob(os.path.join(pasta_entrada, "*.dmnd"))
            if not bancos:
                raise FileNotFoundError(f"No .dmnd database found in {pasta_entrada}.")
            banco_dmnd = bancos[0]
        else:
            raise FileNotFoundError(f"DIAMOND database not found: {pasta_entrada}")

        # DIAMOND expects the path without the .dmnd extension in the search call
        db_path = banco_dmnd[:-len(".dmnd")]

    # 3. Run the screening
    saida_blast = os.path.join(diretorio_saida, "triagem_viral_diamond.tsv")
    fasta_nao_alinhados = os.path.join(diretorio_saida, "orfs_sem_hit_viral.fasta")
    
    comando_blast = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", orfs_input,
        "-o", saida_blast,
        # "qseq" at the end returns the aligned query sequence in the TSV
        "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "stitle", "pident", "evalue", "bitscore", "qcovhsp", "qtitle", "qseq",
        # "--un" exports the ORFs with no hit
        "--un", fasta_nao_alinhados,
        # 50 hits per ORF: the mix of genome types among the best hits is what separates a
        # genuine EVE from a host gene captured by a virus (see utils).
        "--max-target-seqs", str(hits_por_query),
        "--evalue", "1e-5",
        "--more-sensitive",
        "--threads", str(threads)
    ]
    
    subprocess.run(comando_blast, check=True)
    return saida_blast, fasta_nao_alinhados


def buscar_proteinas_referencia_diamond(
    fasta_query: str,
    banco_geral: str,
    diretorio_saida: str,
    threads: int = 4,
    hits_por_query: int = 100,
    top_n_total: int = None,
    evalue: float = 1e-3,
    export_xml: bool = False
) -> str:
    """
    Run DIAMOND blastx of the candidate EVEs against a general protein database (for
    example nr formatted as .dmnd), collecting up to 'hits_por_query' best hits per query.
    Captures coordinates, identity, positives, coverage and stitle, and pulls the subject
    protein sequences from the 'sseq' field. Can also export BLAST XML (outfmt 5).
    """
    os.makedirs(diretorio_saida, exist_ok=True)

    if os.path.isfile(banco_geral):
        if not banco_geral.endswith(".dmnd"):
            raise ValueError(f"--diamond-db-geral must point to a .dmnd file: {banco_geral}")
        db_path = banco_geral[:-len(".dmnd")]
    elif os.path.isdir(banco_geral):
        bancos = glob(os.path.join(banco_geral, "*.dmnd"))
        if not bancos:
            raise FileNotFoundError(f"No .dmnd database found in {banco_geral}.")
        db_path = bancos[0][:-len(".dmnd")]
    else:
        raise FileNotFoundError(f"General DIAMOND database not found: {banco_geral}")

    saida_tsv = os.path.join(diretorio_saida, "diamond_blastx_referencia_geral.tsv")

    # Extended outfmt 6 with 19 fields, including coordinates and sseq
    colunas_outfmt = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "positive", "ppos", "qcovhsp", "stitle", "sseq"
    ]

    comando_blast = [
        "diamond", "blastx",
        "-d", db_path,
        "-q", fasta_query,
        "-o", saida_tsv,
        "--outfmt", "6", *colunas_outfmt,
        "--max-target-seqs", str(hits_por_query),
        "--evalue", str(evalue),
        "--more-sensitive",
        "--threads", str(threads)
    ]

    try:
        subprocess.run(comando_blast, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(f"DIAMOND blastx against the general database {banco_geral} failed.\n"
                            f"Details: {erro.stderr}")

    if not os.path.exists(saida_tsv) or os.path.getsize(saida_tsv) == 0:
        raise RuntimeError(f"DIAMOND blastx found no hit against the general database {banco_geral}.")

    # Optional BLAST XML export (outfmt 5)
    if export_xml:
        saida_xml = os.path.join(diretorio_saida, "diamond_blastx_referencia_geral.xml")
        comando_xml = [
            "diamond", "blastx",
            "-d", db_path,
            "-q", fasta_query,
            "-o", saida_xml,
            "--outfmt", "5",
            "--max-target-seqs", str(hits_por_query),
            "--evalue", str(evalue),
            "--more-sensitive",
            "--threads", str(threads)
        ]
        try:
            subprocess.run(comando_xml, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(f"[run] Second BLASTX alignments also saved as XML: {saida_xml}")
        except Exception as erro_xml:
            print(f"[warning] Could not write the second BLASTX XML: {erro_xml}")

    # Keep only the best sequence per unique sseqid among all returned hits
    melhor_por_sseqid = {}
    with open(saida_tsv, 'r', encoding='utf-8', errors='replace') as f_in:
        for linha in f_in:
            campos = linha.rstrip("\n").split("\t")
            if len(campos) >= 19:
                # New extended format: sseqid=1, bitscore=11, sseq=18
                _qseqid, sseqid = campos[0], campos[1]
                bitscore_str = campos[11]
                sseq = campos[18]
            elif len(campos) >= 5:
                # Legacy 5-column format: sseqid=1, bitscore=2, sseq=4
                _qseqid, sseqid = campos[0], campos[1]
                bitscore_str = campos[2]
                sseq = campos[4]
            else:
                continue

            try:
                bitscore = float(bitscore_str)
            except ValueError:
                continue
            if sseqid not in melhor_por_sseqid or bitscore > melhor_por_sseqid[sseqid][0]:
                melhor_por_sseqid[sseqid] = (bitscore, sseq)

    # Sort by bitscore and apply the overall cap when requested
    hits_ordenados = sorted(melhor_por_sseqid.items(), key=lambda item: item[1][0], reverse=True)
    if top_n_total and top_n_total > 0:
        hits_ordenados = hits_ordenados[:top_n_total]

    if not hits_ordenados:
        raise RuntimeError("No valid hit extracted from the DIAMOND output against the general database.")

    fasta_saida = os.path.join(diretorio_saida, f"proteinas_referencia_diamond_top{len(hits_ordenados)}.fasta")
    with open(fasta_saida, 'w', encoding='utf-8') as f_out:
        for sseqid, (_bitscore, sseq) in hits_ordenados:
            sseq_limpo = sseq.replace("*", "")
            f_out.write(f">{sseqid}\n{sseq_limpo}\n")

    print(f"[run] {len(hits_ordenados)} viral reference proteins (local DIAMOND, up to {hits_por_query} hits/query) saved to: {fasta_saida}")
    return fasta_saida


def gerar_referencias_virais_filtradas(
    viral_tsv: str,
    geral_tsv: str,
    output_fasta: str,
    email_ncbi: str
) -> str:
    from Bio import Entrez, SeqIO
    import os
    if not email_ncbi:
        raise ValueError("a contact e-mail is required for NCBI Entrez "
                         "(--ncbi-email, EVE_NCBI_EMAIL or [ncbi] email in eve_config.ini)")
    Entrez.email = email_ncbi

    # 1. Load the viral hits (best hit per query)
    best_viral_hits = {}
    if os.path.exists(viral_tsv):
        with open(viral_tsv, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('	')
                qseqid = parts[0]
                # In triagem_viral_diamond.tsv, outfmt is "6 qseqid qlen sseqid slen stitle ..."
                # So sseqid is at index 2 if the second column is purely digits and qseqid matches ORF
                if len(parts) >= 11 and parts[1].isdigit():
                    sseqid = parts[2]
                else:
                    sseqid = parts[1]
                    
                if qseqid not in best_viral_hits:
                    best_viral_hits[qseqid] = sseqid

    # 2. Load the general hits and drop the disallowed queries
    queries_proibidas = set()
    if geral_tsv and os.path.exists(geral_tsv):
        with open(geral_tsv, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('	')
                qseqid = parts[0]
                sseqid = parts[1]
                
                stitle = ""
                if len(parts) >= 18:
                    # format from buscar_proteinas_referencia_diamond: outfmt 6 with 18+ cols
                    # title is near the end, index 16
                    stitle = parts[16]
                elif len(parts) >= 15:
                    stitle = parts[-1]

                title_lower = stitle.lower()
                
                from evetool.core.utils import is_hit_celular

                if is_hit_celular(sseqid):
                    queries_proibidas.add(qseqid)
                elif "zinc finger" in title_lower:
                    queries_proibidas.add(qseqid)
                elif "retrovir" in title_lower:
                    queries_proibidas.add(qseqid)
                elif "reverse transcriptase" in title_lower or "pol gag" in title_lower:
                    queries_proibidas.add(qseqid)

    # 3. Decide which viruses are valid
    virus_ids_validos = set()
    for qseqid, v_id in best_viral_hits.items():
        if qseqid not in queries_proibidas:
            # Clean up accession (e.g. ref|YP_123| -> YP_123)
            if '|' in v_id:
                parts = v_id.split('|')
                clean_vid = parts[-2] if parts[-1] == '' else parts[-1]
                # Sometimes it is like 'YP_010840456.1 |polymerase...' so split by space
                clean_vid = clean_vid.split()[0]
            else:
                clean_vid = v_id.split()[0]
            virus_ids_validos.add(clean_vid)

    if not virus_ids_validos:
        print("[warning] No valid virus left after the general filter.")
        return None
        
    print(f"[run] Downloading {len(virus_ids_validos)} validated viral proteins from the NCBI...")
    # 4. Download through EFetch
    with open(output_fasta, 'w') as f_out:
        batch_size = 50
        ids_list = list(virus_ids_validos)
        for i in range(0, len(ids_list), batch_size):
            batch = ids_list[i:i+batch_size]
            try:
                handle = Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text")
                records = list(SeqIO.parse(handle, "fasta"))
                SeqIO.write(records, f_out, "fasta")
            except Exception as e:
                print(f"[warning] Failed to download batch {i}: {e}")
                
    return output_fasta if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0 else None

