import os
import subprocess
from pathlib import Path

def remover_redundancia(pep_fasta: str, diretorio_saida: str, identidade: float = 0.9, cobertura: float = 0.9) -> str:
    """
    Run MMseqs2 to collapse highly similar proteins and return the path of the
    FASTA file holding the unique representatives.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    caminho_fasta = Path(pep_fasta)
    
    # Output prefix and the temporary folder MMseqs2 requires
    prefixo_saida = os.path.join(diretorio_saida, f"{caminho_fasta.stem}_nr")
    pasta_tmp = os.path.join(diretorio_saida, "tmp_mmseqs")
    
    comando_mmseqs = [
        "mmseqs", "easy-cluster",
        str(caminho_fasta),
        prefixo_saida,
        pasta_tmp,
        "--min-seq-id", str(identidade),
        "-c", str(cobertura),
        # cov-mode 0: coverage is required on BOTH sequences. In mode 1 (target only) a short
        # ORF is absorbed by a much larger one that merely contains it, and disappears
        # together with the representative when that one has no viral hit.
        "--cov-mode", "0"
    ]
    
    try:
        subprocess.run(
            comando_mmseqs,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # MMseqs2 appends '_rep_seq.fasta' to the given prefix
        fasta_representantes = f"{prefixo_saida}_rep_seq.fasta"
        return fasta_representantes
        
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(f"MMseqs2 failed for {caminho_fasta.name}.\nDetails: {erro.stderr}")

def caminho_tabela_clusters(fasta_representantes: str) -> str:
    """The 'representative<TAB>member' table MMseqs2 writes next to the representatives FASTA."""
    return fasta_representantes.replace("_rep_seq.fasta", "_cluster.tsv")


def resgatar_nucleotideos_representantes(fasta_rep_pep: str, fasta_nt_original: str, diretorio_saida: str) -> str:
    """
    Read the representative protein IDs chosen by MMseqs2 and extract the matching
    nucleotide ORFs from the original FASTA file.
    """
    caminho_pep = Path(fasta_rep_pep)
    arquivo_saida_nt = os.path.join(diretorio_saida, f"{caminho_pep.stem}_nt.fasta")
    
    # 1. Collect the cluster representative IDs
    ids_representantes = set()
    with open(fasta_rep_pep, 'r') as f_pep:
        for linha in f_pep:
            if linha.startswith(">"):
                # Take the exact ID, ignoring any description after the space
                id_seq = linha.strip().lstrip(">").split()[0]
                ids_representantes.add(id_seq)
                
    # 2. Filter the original nucleotide file
    with open(fasta_nt_original, 'r') as f_nt_in, open(arquivo_saida_nt, 'w') as f_nt_out:
        gravar_sequencia = False
        for linha in f_nt_in:
            if linha.startswith(">"):
                id_seq = linha.strip().lstrip(">").split()[0]
                # Check whether this ORF survived the redundancy filter
                if id_seq in ids_representantes:
                    gravar_sequencia = True
                    f_nt_out.write(linha)
                else:
                    gravar_sequencia = False
            elif gravar_sequencia:
                f_nt_out.write(linha)
                
    return arquivo_saida_nt
