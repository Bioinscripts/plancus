import os
import shutil
import subprocess
from glob import glob
from pathlib import Path


def resgatar_peptideos_nao_alinhados(fasta_un_nt: str, fasta_pep_original: str, diretorio_saida: str) -> str:
    """
    Read the IDs of the nucleotide ORFs that failed DIAMOND and extract their
    amino-acid versions for the HMMER search.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    caminho_un = Path(fasta_un_nt)
    saida_pep_un = os.path.join(diretorio_saida, f"{caminho_un.stem}_pep.fasta")
    
    # 1. Collect the IDs of the unaligned sequences
    ids_nao_alinhados = set()
    if os.path.exists(fasta_un_nt) and os.path.getsize(fasta_un_nt) > 0:
        with open(fasta_un_nt, 'r') as f_in:
            for linha in f_in:
                if linha.startswith(">"):
                    id_seq = linha.strip().lstrip(">").split()[0]
                    ids_nao_alinhados.add(id_seq)
                
    # 2. Scan the original peptide file and keep only the matching records
    with open(fasta_pep_original, 'r') as f_pep, open(saida_pep_un, 'w') as f_out:
        gravar = False
        for linha in f_pep:
            if linha.startswith(">"):
                id_seq = linha.strip().lstrip(">").split()[0]
                if id_seq in ids_nao_alinhados:
                    gravar = True
                    f_out.write(linha)
                else:
                    gravar = False
            elif gravar:
                f_out.write(linha)
                
    return saida_pep_un


def _escanear_com_pyhmmer(fasta_pep_un: str, diretorio_saida: str, hmm_dir: str = None,
                          hmm_db: str = None, threads: int = 4, lote_modelos: int = 200,
                          usar_tc: bool = None, evalue: float = 1e-3) -> str:
    """
    Run the deep HMM domain search with pyhmmer (direct HMMER3/Easel C bindings):
    no manual hmmpress needed, and the multithreaded search releases the GIL.
    """
    import pyhmmer

    os.makedirs(diretorio_saida, exist_ok=True)
    saida_tabela_hmm = os.path.join(diretorio_saida, "triagem_profunda_hmm.tsv")

    # 1. Validate the peptide file
    if not os.path.exists(fasta_pep_un) or os.path.getsize(fasta_pep_un) == 0:
        print("[HMM Scanner/pyhmmer] No unaligned peptides to screen. Writing an empty table.")
        with open(saida_tabela_hmm, "w") as f_out:
            f_out.write("# pyhmmer: no unaligned peptide found to scan.\n")
        return saida_tabela_hmm

    alphabet = pyhmmer.easel.Alphabet.amino()
    try:
        with pyhmmer.easel.SequenceFile(fasta_pep_un, digital=True, alphabet=alphabet) as seq_file:
            targets = seq_file.read_block()
    except (EOFError, ValueError) as err:
        print(f"[HMM Scanner/pyhmmer] Warning: empty or invalid peptide file ({err}). Writing an empty table.")
        with open(saida_tabela_hmm, "w") as f_out:
            f_out.write("# pyhmmer: empty peptide file or no valid sequences.\n")
        return saida_tabela_hmm

    if len(targets) == 0:
        print("[HMM Scanner/pyhmmer] Zero sequences read from the peptide FASTA.")
        with open(saida_tabela_hmm, "w") as f_out:
            f_out.write("# pyhmmer: zero sequences to scan.\n")
        return saida_tabela_hmm

    # 2. Load the HMM models
    fontes_hmm = []
    if hmm_dir:
        arquivos_hmm = sorted(glob(os.path.join(hmm_dir, "*.hmm")))
        if not arquivos_hmm:
            raise FileNotFoundError(f"No .hmm file found in '{hmm_dir}'.")

        # Concatenate into a single file in the output directory, for the record and reproducibility
        banco_hmm_final = os.path.join(diretorio_saida, "banco_viral_customizado.hmm")
        with open(banco_hmm_final, 'w') as out_db:
            for hmm_file in arquivos_hmm:
                with open(hmm_file, 'r') as in_hmm:
                    out_db.write(in_hmm.read())

        fontes_hmm = arquivos_hmm

    elif hmm_db:
        # If it points to a .h3m file, drop the suffix to open the base .hmm
        caminho_base = hmm_db[:-4] if hmm_db.endswith(".h3m") else hmm_db
        if not os.path.exists(caminho_base) and os.path.exists(hmm_db):
            caminho_base = hmm_db

        fontes_hmm = [caminho_base]
    else:
        raise ValueError("Either --hmm-dir or --hmm-db must be provided.")

    def _modelos_em_lotes(caminhos, tamanho):
        """Iterate the models in batches, without holding the whole database in memory."""
        lote = []
        for caminho in caminhos:
            with pyhmmer.plan7.HMMFile(caminho) as hf:
                source = hf.optimized_profiles() if hf.is_pressed() else hf
                for model in source:
                    lote.append(model)
                    if len(lote) >= tamanho:
                        yield lote
                        lote = []
        if lote:
            yield lote

    n_modelos = 0
    for caminho in fontes_hmm:
        with pyhmmer.plan7.HMMFile(caminho) as hf:
            n_modelos += sum(1 for _ in (hf.optimized_profiles() if hf.is_pressed() else hf))

    if not n_modelos:
        raise RuntimeError("No valid HMM profile was loaded for the search.")

    print(f"[HMM Scanner/pyhmmer] Running native C hmmsearch "
          f"(HMM models: {n_modelos}, peptides: {len(targets)}, CPUs: {threads}, "
          f"batches of {lote_modelos} models)...")

    # 3. HMMER search through pyhmmer (hmmsearch in C)
    cpus_count = max(1, int(threads))
    first = True
    total_hits = 0

    criterion = None
    processados = 0
    with open(saida_tabela_hmm, "wb") as f_out:
        for lote in _modelos_em_lotes(fontes_hmm, lote_modelos):
            if criterion is None:
                # A calibrated database (TC stored in each model) uses the per-model cutoff;
                # with usar_tc=False we force an e-value so weak hits are recorded too.
                tem_tc = all(getattr(h.cutoffs, "trusted", None) is not None for h in lote)
                aplicar_tc = tem_tc if usar_tc is None else bool(usar_tc and tem_tc)
                criterion = ({"bit_cutoffs": "trusted"} if aplicar_tc
                            else {"E": evalue, "domE": evalue})
                print("[HMM Scanner/pyhmmer] " + ("Banco calibrado: limiar TC por model."
                      if aplicar_tc else f"Permissive mode: e-value <= {evalue} (TC recorded separately)."))
            for top_hits in pyhmmer.hmmer.hmmsearch(lote, targets, cpus=cpus_count, **criterion):
                if len(top_hits) > 0:
                    total_hits += len(top_hits)
                top_hits.write(f_out, format="domains", header=first)
                first = False
            processados += len(lote)
            if n_modelos > lote_modelos:
                print(f"[HMM Scanner/pyhmmer] {processados}/{n_modelos} models scanned...", flush=True)
            del lote

    if not os.path.exists(saida_tabela_hmm) or os.path.getsize(saida_tabela_hmm) == 0:
        with open(saida_tabela_hmm, "w") as f_out:
            f_out.write("# pyhmmer: no hit found at E-value <= 1e-3.\n")

    print(f"[HMM Scanner/pyhmmer] Scan finished. Profiles with hits: {total_hits} | domtblout table: {saida_tabela_hmm}")
    return saida_tabela_hmm


def _escanear_com_subprocess(fasta_pep_un: str, diretorio_saida: str, hmm_dir: str = None, hmm_db: str = None, threads: int = 4) -> str:
    """
    Legacy fallback: prepare the HMM database with hmmpress (when needed) and run the
    system hmmsearch binary through subprocess.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    banco_hmm_final = None
    
    if hmm_dir:
        arquivos_hmm = glob(os.path.join(hmm_dir, "*.hmm"))
        if not arquivos_hmm:
            raise FileNotFoundError(f"No .hmm file found in {hmm_dir}.")
        
        banco_hmm_final = os.path.join(diretorio_saida, "banco_viral_customizado.hmm")
        with open(banco_hmm_final, 'w') as out_db:
            for hmm_file in arquivos_hmm:
                with open(hmm_file, 'r') as in_hmm:
                    out_db.write(in_hmm.read())
        
        if shutil.which("hmmpress"):
            subprocess.run(["hmmpress", banco_hmm_final], check=True, stdout=subprocess.DEVNULL)
        
    elif hmm_db:
        banco_hmm_final = hmm_db
        if not os.path.exists(f"{banco_hmm_final}.h3m") and shutil.which("hmmpress"):
            subprocess.run(["hmmpress", banco_hmm_final], check=True, stdout=subprocess.DEVNULL)
    else:
        raise ValueError("Either --hmm-dir or --hmm-db must be provided.")

    saida_tabela_hmm = os.path.join(diretorio_saida, "triagem_profunda_hmm.tsv")
    
    comando_hmmsearch = [
        "hmmsearch",
        "--domtblout", saida_tabela_hmm,
        "--cpu", str(threads),
        "-E", "1e-3",
        "--domE", "1e-3",
        banco_hmm_final,
        fasta_pep_un
    ]
    
    subprocess.run(comando_hmmsearch, check=True, stdout=subprocess.DEVNULL)
    return saida_tabela_hmm


def escanear_hmm(fasta_pep_un: str, diretorio_saida: str, hmm_dir: str = None, hmm_db: str = None,
                 threads: int = 4, lote_modelos: int = 200, usar_tc: bool = None,
                 evalue: float = 1e-3) -> str:
    """
    Single entry point for the deep HMM scan. Prefers pyhmmer (fast, direct C bindings)
    and falls back automatically to the hmmsearch binary when pyhmmer is missing or fails.
    """
    try:
        import pyhmmer
        return _escanear_com_pyhmmer(
            fasta_pep_un=fasta_pep_un,
            diretorio_saida=diretorio_saida,
            hmm_dir=hmm_dir,
            hmm_db=hmm_db,
            threads=threads,
            lote_modelos=lote_modelos,
            usar_tc=usar_tc,
            evalue=evalue
        )
    except ImportError:
        print("[HMM Scanner] pyhmmer not found in this environment. Falling back to the system hmmsearch CLI...")
        return _escanear_com_subprocess(
            fasta_pep_un=fasta_pep_un,
            diretorio_saida=diretorio_saida,
            hmm_dir=hmm_dir,
            hmm_db=hmm_db,
            threads=threads
        )
    except Exception as err_pyhmmer:
        print(f"[HMM Scanner] pyhmmer raised an unexpected error ({err_pyhmmer}). Falling back to the hmmsearch CLI...")
        return _escanear_com_subprocess(
            fasta_pep_un=fasta_pep_un,
            diretorio_saida=diretorio_saida,
            hmm_dir=hmm_dir,
            hmm_db=hmm_db,
            threads=threads
        )
