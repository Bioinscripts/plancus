import os
import shutil
import zipfile
import subprocess
from glob import glob
from pathlib import Path


def verificar_datasets_cli():
    """Checks whether the NCBI 'datasets' binary is installed on the PATH."""
    if shutil.which("datasets") is None:
        raise RuntimeError(
            "The NCBI 'datasets' utility was not found on the PATH.\n"
            "Instale-o via conda/mamba com: 'conda install -c conda-forge ncbi-datasets-cli' "
            "ou via Pixi adicionando 'ncbi-datasets-cli' ao pixi.toml."
        )


def e_acesso_ncbi(texto: str) -> bool:
    """Verifica se uma string possui formato de accession do NCBI GenBank ou RefSeq (GCA_ ou GCF_)."""
    if not texto or not isinstance(texto, str):
        return False
    texto_limpo = texto.strip().upper()
    return texto_limpo.startswith("GCA_") or texto_limpo.startswith("GCF_")


def baixar_genoma_ncbi(
    accession: str,
    pasta_destino: str = "genomas_ncbi",
    baixar_gff: bool = True,
    sobrescrever: bool = False
) -> dict:
    """
    Faz o download automatizado de um genoma pelo NCBI Datasets CLI (GCA_ ou GCF_).
    Extracts the genomic FASTA and, when available, the GFF3 annotation.
    Returns a dictionary with the local paths of the files:
      {"accession": accession, "genome": caminho_fna, "gff": caminho_gff}
    """
    verificar_datasets_cli()
    os.makedirs(pasta_destino, exist_ok=True)
    accession = accession.strip()

    # 1. Check whether the files were already downloaded
    fna_existentes = glob(os.path.join(pasta_destino, f"*{accession}*.fna")) + \
                     glob(os.path.join(pasta_destino, f"*{accession}*.fasta"))
    gff_existentes = glob(os.path.join(pasta_destino, f"*{accession}*.gff3")) + \
                     glob(os.path.join(pasta_destino, f"*{accession}*.gff"))

    if fna_existentes and not sobrescrever:
        caminho_fna = fna_existentes[0]
        caminho_gff = gff_existentes[0] if gff_existentes else None
        print(f"[NCBI Datasets] Genome {accession} already present locally: {caminho_fna}")
        return {
            "accession": accession,
            "genome": caminho_fna,
            "gff": caminho_gff
        }

    # 2. Executa o download via NCBI Datasets CLI
    print(f"\n[NCBI Datasets] Baixando genoma para o accession: {accession}...")
    arquivo_zip = os.path.join(pasta_destino, f"{accession}_dataset.zip")
    pasta_tmp = os.path.join(pasta_destino, f"tmp_{accession}")

    include_tipos = "genome,gff3" if baixar_gff else "genome"
    cmd_download = [
        "datasets", "download", "genome", "accession", accession,
        "--include", include_tipos,
        "--filename", arquivo_zip,
        "--no-progressbar"
    ]

    try:
        subprocess.run(cmd_download, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as erro:
        raise RuntimeError(
            f"[NCBI Datasets] Falha ao baixar o accession '{accession}'.\n"
            f"Erro detalhado: {erro.stderr.strip()}"
        )

    if not os.path.exists(arquivo_zip) or os.path.getsize(arquivo_zip) == 0:
        raise RuntimeError(f"[NCBI Datasets] Download failed: no ZIP file was produced for {accession}.")

    # 3. Descompacta o arquivo ZIP
    print(f"[NCBI Datasets] Descompactando pacote de dados de {accession}...")
    os.makedirs(pasta_tmp, exist_ok=True)
    with zipfile.ZipFile(arquivo_zip, "r") as z:
        z.extractall(pasta_tmp)

    # 4. Localiza o FASTA (.fna) e o GFF (.gff) descompactados
    fnas_encontrados = glob(os.path.join(pasta_tmp, "**", "*.fna"), recursive=True)
    if not fnas_encontrados:
        raise RuntimeError(f"[NCBI Datasets] Nenhum arquivo .fna encontrado no pacote de {accession}.")

    fna_orig = fnas_encontrados[0]
    nome_fna_final = os.path.basename(fna_orig)
    caminho_fna_final = os.path.join(pasta_destino, nome_fna_final)
    shutil.copy2(fna_orig, caminho_fna_final)

    # Localiza o GFF3, se existir
    gffs_encontrados = glob(os.path.join(pasta_tmp, "**", "genomic.gff"), recursive=True) + \
                      glob(os.path.join(pasta_tmp, "**", "*.gff3"), recursive=True)
    caminho_gff_final = None
    if gffs_encontrados:
        gff_orig = gffs_encontrados[0]
        caminho_gff_final = os.path.join(pasta_destino, f"{accession}_genomic.gff3")
        shutil.copy2(gff_orig, caminho_gff_final)
        print(f"[NCBI Datasets] GFF3 annotation found and saved to: {caminho_gff_final}")
    else:
        print(f"[NCBI Datasets] No GFF3 annotation available for {accession} (only the genome was saved).")

    # 5. Clean up the temporary files
    try:
        os.remove(arquivo_zip)
        shutil.rmtree(pasta_tmp)
    except Exception:
        pass

    print(f"[NCBI Datasets] Genoma FASTA salvo com sucesso em: {caminho_fna_final}")
    return {
        "accession": accession,
        "genome": caminho_fna_final,
        "gff": caminho_gff_final
    }


def resolver_genoma_ou_accession(identificador: str, pasta_genomas: str = "genomas") -> dict:
    """
    Decides whether the identifier is an accession (GCA_/GCF_) or a local file.
    An accession is downloaded automatically.
    Se for arquivo local existente, retorna os caminhos locais.
    """
    if e_acesso_ncbi(identificador):
        return baixar_genoma_ncbi(identificador, pasta_destino=pasta_genomas)

    # Caminho de arquivo local
    caminho_arquivo = identificador if os.path.isabs(identificador) else os.path.join(pasta_genomas, identificador)
    if not os.path.exists(caminho_arquivo):
        # Try the current directory
        if os.path.exists(identificador):
            caminho_arquivo = identificador
        else:
            raise FileNotFoundError(f"Local genome file not found: {identificador}")

    return {
        "accession": None,
        "genome": caminho_arquivo,
        "gff": None
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Automated download of genomes and GFF3 annotations through NCBI Datasets.")
    parser.add_argument("-a", "--accession", required=True, help="Accession do NCBI (ex: GCA_058968885.1 ou GCF_000466445.1)")
    parser.add_argument("-o", "--outdir", default="genomas_ncbi", help="Pasta de destino para o FASTA e GFF3")
    parser.add_argument("--no-gff", action="store_true", help="Do not download the GFF3 annotation")
    parser.add_argument("--force", action="store_true", help="Force a new download even when the file already exists")

    args = parser.parse_args()
    resultado = baixar_genoma_ncbi(
        accession=args.accession,
        pasta_destino=args.outdir,
        baixar_gff=not args.no_gff,
        sobrescrever=args.force
    )
    print("\n[Resultado do Download]")
    print(f"  Genoma FASTA : {resultado['genome']}")
    print(f"  GFF3 annotation: {resultado['gff']}")
