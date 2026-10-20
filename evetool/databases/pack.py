"""
Packs the curated viral databases of this machine into one file to publish (e.g. on Zenodo),
so that `./eve-tool databases --bundle-url <address>` installs them anywhere.

  ./eve-tool pack-databases -o eve_tool_db_v1.tar
  ./eve-tool pack-databases -o new.tar --hmm-dir <build>/hmm --diamond-dir <build>/diamond   # a rebuilt database

Contents (paths as in config.py; large text files gzipped inside the archive):
  diamond/sequences.fasta                                  viral proteins (DIAMOND is rebuilt on install)
  hmm/COMBINED_DATABASE/combined_viral_calibrado.hmm       calibrated RVDB + eggNOG models
  hmm/COMBINED_DATABASE/combined_viral_calibrado_tc.json   their trusted cut-offs
  hmm/COMBINED_DATABASE/combined_viral_annotations.tsv     taxonomy and genome type of each model
  hmm/COMBINED_DATABASE/build_manifest.json                source versions and recipe (when present)
  hmm/RVDB/rvdb_base_annotations.tsv                       RVDB functional keywords
  hmm/RVDB/rvdb_family_genome_mapping.tsv                  viral family -> genome type
  hmm/CDD/rpsdb/Cdd_domains.tsv                            source taxon of every CDD domain

Public databases (Swiss-Prot, taxonomy, Pfam, CDD, proteomes) are not packed: prepare.py
downloads them from their sources. Prints the MD5 to paste into BUNDLE_MD5 of prepare.py.
"""
import argparse
import gzip
import hashlib
import io
import os
import shutil
import sys
import tarfile
import tempfile
import time

from evetool import config

ARQUIVOS = [
    ("diamond", "sequences.fasta", True),
    ("hmm", "COMBINED_DATABASE/combined_viral_calibrado.hmm", True),
    ("hmm", "COMBINED_DATABASE/combined_viral_calibrado_tc.json", False),
    ("hmm", "COMBINED_DATABASE/combined_viral_annotations.tsv", False),
    ("hmm", "RVDB/rvdb_base_annotations.tsv", False),
    ("hmm", "RVDB/rvdb_family_genome_mapping.tsv", False),
    ("hmm", "CDD/rpsdb/Cdd_domains.tsv", False),
]
OPCIONAIS = [("hmm", "COMBINED_DATABASE/build_manifest.json", False)]


def empacotar(saida: str, hmm_dir: str = None, diamond_dir: str = None) -> str:
    raiz = {"hmm": hmm_dir or config.HMM_DIR, "diamond": diamond_dir or config.DIAMOND_DIR}
    faltam = [os.path.join(raiz[r], c) for r, c, _ in ARQUIVOS if not os.path.exists(os.path.join(raiz[r], c))]
    if faltam:
        sys.exit("Missing:\n  " + "\n  ".join(faltam))
    itens = ARQUIVOS + [x for x in OPCIONAIS if os.path.exists(os.path.join(raiz[x[0]], x[1]))]

    with tarfile.open(saida, "w") as tar, tempfile.TemporaryDirectory() as tmp:
        for r, rel, comprimir in itens:
            origem = os.path.join(raiz[r], rel)
            nome = f"{r}/{rel}"
            if comprimir:
                gz = os.path.join(tmp, os.path.basename(rel) + ".gz")
                print(f"compressing {origem} ...", flush=True)
                with open(origem, "rb") as src, gzip.open(gz, "wb", compresslevel=6) as dst:
                    shutil.copyfileobj(src, dst, length=16 << 20)
                tar.add(gz, arcname=nome + ".gz")
                os.remove(gz)
            else:
                tar.add(origem, arcname=nome)
        leia = (b"EVE tool viral databases. Install with: ./eve-tool databases --bundle <this file>\n")
        info = tarfile.TarInfo("LEIAME.txt")
        info.size = len(leia)
        info.mtime = int(time.time())
        tar.addfile(info, io.BytesIO(leia))

    h = hashlib.md5()
    with open(saida, "rb") as fh:
        for bloco in iter(lambda: fh.read(8 << 20), b""):
            h.update(bloco)
    print(f"\n{saida}: {os.path.getsize(saida) / 1e6:.0f} MB  md5 {h.hexdigest()}")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="eve_tool_db_v1.tar")
    ap.add_argument("--hmm-dir", help="HMM folder to pack (default: the installed one)")
    ap.add_argument("--diamond-dir", help="DIAMOND folder to pack (default: the installed one)")
    a = ap.parse_args()
    empacotar(a.output, a.hmm_dir, a.diamond_dir)
    print("After uploading, set BUNDLE_URL and BUNDLE_MD5 in prepare.py.")


if __name__ == "__main__":
    main()
