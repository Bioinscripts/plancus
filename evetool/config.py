"""
Where the pipeline finds its databases and writes its results, kept out of the code.

Each path is resolved in this order:
  1. environment variable (EVE_DATA_DIR, EVE_HMM_DIR, EVE_DIAMOND_DIR, EVE_RESULTS_DIR, EVE_PFAM_HMM)
  2. [paths] section of the config file: $EVE_CONFIG, otherwise eve_config.ini next to this file
  3. default under data_dir (~/eve_data unless set):
       hmm_dir     = <data_dir>/Data_base_hmm_eves
       diamond_dir = <data_dir>/DIAMOND
       results_dir = <data_dir>/Resultados
       pfam_hmm    = <data_dir>/Pfam/Pfam-A.hmm

Copy eve_config.example.ini to eve_config.ini and edit it, or run
`./eve-tool config` to see the paths in use.
"""
import os
import sys
import configparser

# Repository root (eve_config.ini lives there, next to README.md)
PASTA_CODIGO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARQUIVO_CONFIG = os.environ.get("EVE_CONFIG") or os.path.join(PASTA_CODIGO, "eve_config.ini")

_ini = configparser.ConfigParser()
_ini.read(ARQUIVO_CONFIG, encoding="utf-8")


def _valor(chave: str, padrao: str) -> str:
    """Environment variable > config file > default, with ~ and $VARS expanded."""
    v = os.environ.get(f"EVE_{chave.upper()}") or _ini.get("paths", chave, fallback="") or padrao
    return os.path.abspath(os.path.expandvars(os.path.expanduser(v)))


DATA_DIR = _valor("data_dir", "~/eve_data")
HMM_DIR = _valor("hmm_dir", os.path.join(DATA_DIR, "Data_base_hmm_eves"))
DIAMOND_DIR = _valor("diamond_dir", os.path.join(DATA_DIR, "DIAMOND"))
RESULTS_DIR = _valor("results_dir", os.path.join(DATA_DIR, "Resultados"))
PFAM_HMM = _valor("pfam_hmm", os.path.join(DATA_DIR, "Pfam", "Pfam-A.hmm"))


# Standard database files inside those folders (written by prepare.py)
VIRAL_DMND = os.path.join(DIAMOND_DIR, "mydb.dmnd")
GERAL_DMND = os.path.join(DIAMOND_DIR, "geral", "uniprot_sprot.dmnd")
HMM_DB = os.path.join(HMM_DIR, "COMBINED_DATABASE", "combined_viral_calibrado.hmm")
CDD_RPSDB = os.path.join(HMM_DIR, "CDD", "rpsdb", "Cdd")
TAXDUMP = os.path.join(HMM_DIR, "RVDB", "taxdump.tar.gz")
PROTEOMAS_DIR = os.path.join(HMM_DIR, "proteomas_calibracao")


def hmm_path(*partes: str) -> str:
    """Path inside the HMM/annotation database folder (Data_base_hmm_eves)."""
    return os.path.join(HMM_DIR, *partes)


def diamond_path(*partes: str) -> str:
    """Path inside the DIAMOND database folder."""
    return os.path.join(DIAMOND_DIR, *partes)


if __name__ == "__main__":
    caminhos = {"DATA_DIR": DATA_DIR, "HMM_DIR": HMM_DIR, "DIAMOND_DIR": DIAMOND_DIR,
                "RESULTS_DIR": RESULTS_DIR, "PFAM_HMM": PFAM_HMM, "VIRAL_DMND": VIRAL_DMND,
                "GERAL_DMND": GERAL_DMND, "HMM_DB": HMM_DB, "CDD_RPSDB": CDD_RPSDB}
    if "--shell" in sys.argv:
        # For the shell scripts: eval "$(./eve-tool config --shell)"
        import shlex
        for nome, valor in caminhos.items():
            print(f"EVE_{nome}={shlex.quote(valor)}")
    else:
        lido = ARQUIVO_CONFIG if os.path.exists(ARQUIVO_CONFIG) else f"{ARQUIVO_CONFIG} (not found, using defaults)"
        print(f"Config: {lido}")
        for nome, valor in caminhos.items():
            existe = os.path.exists(valor) or os.path.exists(valor + ".pal") or os.path.exists(valor + ".rps")
            print(f"  {nome:12s} {valor}{'' if existe else '   [missing]'}")
        print("Full check of tools and databases: ./eve-tool check")
