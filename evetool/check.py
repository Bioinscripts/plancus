"""
Checks that this machine can run the pipeline: command-line tools, Python libraries and databases.

  ./eve-tool check

Each missing item comes with how to get it. Exit code 1 when something required is missing.
"""
import importlib
import os
import shutil
import sys

from evetool import config

OK, FALTA, OPC = "ok", "MISSING", "optional, missing"

# (command, required, what it is used for)
FERRAMENTAS = [
    ("orfipy", True, "ORF prediction"),
    ("mmseqs", True, "redundancy removal"),
    ("diamond", True, "viral and cellular screening"),
    ("makeblastdb", True, "BLAST+ databases"),
    ("tblastn", True, "genome anchoring and tblastn scan"),
    ("tblastx", True, "cross-species clustering"),
    ("bedtools", True, "loci and flanks"),
    ("samtools", True, "sequence extraction"),
    ("dustmasker", False, "low complexity (nucleotide)"),
    ("segmasker", False, "low complexity (protein)"),
    ("rpsblast", False, "CD-Search, local mode"),
    ("hmmscan", False, "Pfam annotation"),
    ("datasets", False, "NCBI genome download from GCA_/GCF_ accessions"),
    ("hisat2", False, "RNA-seq coverage (short reads)"),
    ("minimap2", False, "RNA-seq coverage (long reads)"),
    ("zenity", False, "folder dialog of open_curation.sh"),
]
BIBLIOTECAS = [("Bio", True), ("pandas", True), ("numpy", True), ("scipy", True), ("matplotlib", True),
               ("seaborn", True), ("networkx", True), ("requests", True), ("pyhmmer", True), ("pysam", False)]
# (label, path, required, part of prepare.py)
BANCOS = [
    ("viral DIAMOND", config.VIRAL_DMND, True, "bundle"),
    ("viral HMMs (calibrated)", config.HMM_DB, True, "bundle"),
    ("HMM annotations", config.hmm_path("COMBINED_DATABASE", "combined_viral_annotations.tsv"), True, "bundle"),
    ("RVDB keywords", config.hmm_path("RVDB", "rvdb_base_annotations.tsv"), False, "bundle"),
    ("Swiss-Prot DIAMOND", config.GERAL_DMND, True, "swissprot"),
    ("NCBI taxonomy", config.TAXDUMP, True, "taxdump"),
    ("synteny proteome(s)", config.PROTEOMAS_DIR, False, "proteome"),
    ("Pfam-A", config.PFAM_HMM, False, "pfam"),
    ("CDD (rpsblast)", config.CDD_RPSDB, False, "cdd"),
]


def existe(caminho: str) -> bool:
    if os.path.isdir(caminho):
        return any(n.endswith(".faa") for n in os.listdir(caminho))
    return any(os.path.exists(caminho + s) for s in ("", ".pal", ".rps"))


def linha(estado: str, nome: str, extra: str):
    marca = {"ok": "  ok ", "MISSING": " !!  ", "optional, missing": "  -  "}[estado]
    print(f"{marca} {nome:26s} {extra}")


def main():
    erros = 0
    print(f"Python {sys.version.split()[0]}  ({sys.executable})")
    print(f"Config: {config.ARQUIVO_CONFIG if os.path.exists(config.ARQUIVO_CONFIG) else 'defaults (no eve_config.ini)'}")
    print(f"Data folder: {config.DATA_DIR}\n\nTools")
    for cmd, obrig, uso in FERRAMENTAS:
        achou = shutil.which(cmd) or (cmd == "rpsblast" and shutil.which("rpsblast+"))
        estado = OK if achou else (FALTA if obrig else OPC)
        erros += estado == FALTA
        linha(estado, cmd, uso)
    print("\nPython libraries")
    for mod, obrig in BIBLIOTECAS:
        try:
            importlib.import_module(mod)
            estado = OK
        except Exception:
            estado = FALTA if obrig else OPC
        erros += estado == FALTA
        linha(estado, mod, "")
    print("\nDatabases")
    faltam = set()
    for rot, cam, obrig, parte in BANCOS:
        estado = OK if existe(cam) else (FALTA if obrig else OPC)
        erros += estado == FALTA
        if estado != OK:
            faltam.add(parte)
        linha(estado, rot, cam)

    print()
    if erros:
        print(f"{erros} required item(s) missing.")
        print("  tools/libraries:  conda env create -f environment.yml && conda activate eve_env")
    if faltam:
        opc = [p for p in ("pfam", "cdd") if p in faltam]
        base = sorted(faltam - {"pfam", "cdd"})
        if base:
            print(f"  databases:        ./eve-tool databases --only {' '.join(base)}")
        if opc:
            print(f"  optional:         ./eve-tool databases --only {' '.join(opc)}")
    if not erros:
        print("Ready to run. Example: ./run_pipeline.sh example -i genome.fna")
    sys.exit(1 if erros else 0)


if __name__ == "__main__":
    main()
