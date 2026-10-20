import os
import subprocess
from pathlib import Path

def extrair_orfs(genoma_fasta: str, diretorio_saida: str) -> dict:
    """
    Run orfipy to predict ORFs and write the results to disk.
    Returns a dictionary with the absolute paths of the generated files.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    caminho_genoma = Path(genoma_fasta)
    nome_base = caminho_genoma.stem
    
    # Absolute output paths
    saida_nt = os.path.join(diretorio_saida, f"{nome_base}_orfs_nt.fasta")
    saida_pep = os.path.join(diretorio_saida, f"{nome_base}_orfs_pep.fasta")
    
    # Building the command as a list lets subprocess escape special characters for us
    comando_orfipy = [
        "orfipy", str(caminho_genoma),
        "--min", "100",
        "--max", "10000",
        "--strand", "b",
        "--ignore-case",
        "--include-stop",
        "--partial-3",
        "--partial-5",
        # NOTE: orfipy already joins --outdir with the names given here, so we pass
        # only the file name, not the full path (otherwise it would try to write
        # to "<outdir>/<outdir>/file.fasta").
        "--dna", os.path.basename(saida_nt),
        "--pep", os.path.basename(saida_pep),
        "--outdir", diretorio_saida
    ]
    
    try:
        # Run quietly, but capture the error if orfipy fails
        subprocess.run(
            comando_orfipy, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True
        )
        return {"nt_fasta": saida_nt, "pep_fasta": saida_pep}
        
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(f"orfipy failed for {nome_base}.\nDetails: {erro.stderr}")
