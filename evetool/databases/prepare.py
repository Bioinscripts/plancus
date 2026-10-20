"""
Downloads and installs the databases the pipeline uses, in the folders set by config.py
(eve_config.ini or EVE_* variables; by default under ~/eve_data).

  ./eve-tool databases                 # essentials: viral bundle, Swiss-Prot, taxonomy, synteny proteome
  ./eve-tool databases --pfam --cdd    # also Pfam-A (~300 MB) and the local CDD (~4.5 GB)
  ./eve-tool databases --all
  ./eve-tool databases --bundle eve_tool_db_v1.tar   # viral bundle from a local file

Parts:
  bundle     the curated viral databases of this tool (viral protein FASTA for DIAMOND, calibrated
             RVDB + eggNOG HMMs and their annotations, CDD domain taxa); built by pack.py
  swissprot  UniProtKB/Swiss-Prot, the cellular reference for the classification (DIAMOND)
  taxdump    NCBI taxonomy, for the viral family of each hit
  proteome   host proteome that labels the flanking genes in the synteny step (default P. papatasi)
  pfam       Pfam-A.hmm, pressed (domain annotation, transposon domains of the repeat analysis)
  cdd        NCBI CDD for local rpsblast+ (CD-Search annotation)

Every part is skipped when already installed; --force reinstalls it. Check with
`./eve-tool check`.

The bundle is a fixed release (cite its version). To rebuild the viral databases from the current
RVDB, eggNOG and NCBI Virus releases with the same criteria, see `./eve-tool build-databases -h`.
"""
import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import zipfile

from evetool import config

# Published viral bundle (pack.py). Set it when the Zenodo record exists, or pass --bundle-url.
BUNDLE_URL = os.environ.get("EVE_BUNDLE_URL", "")
BUNDLE_MD5 = ""

URLS = {
    "swissprot": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    "taxdump": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    "pfam": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    "cdd": "https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz",
}
PROTEOMA_PADRAO = "GCF_024763615.1"     # Phlebotomus papatasi, RefSeq

# Bundle members: prefix "hmm/" goes to hmm_dir, "diamond/" to diamond_dir
DESTINOS = {"hmm": lambda: config.HMM_DIR, "diamond": lambda: config.DIAMOND_DIR}


def log(msg):
    print(f"[databases] {msg}", flush=True)


def baixar(url: str, destino: str):
    """Download with progress to destino (via a .part file, so a broken download is never used)."""
    os.makedirs(os.path.dirname(destino), exist_ok=True)
    parte = destino + ".part"
    log(f"downloading {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "eve-tool/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(parte, "wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        feito, ultimo = 0, -1
        while True:
            bloco = r.read(4 << 20)
            if not bloco:
                break
            fh.write(bloco)
            feito += len(bloco)
            pct = int(100 * feito / total) if total else -1
            if pct != ultimo and (pct % 5 == 0 or pct < 0):
                ultimo = pct
                txt = f"{pct}%" if pct >= 0 else f"{feito >> 20} MB"
                print(f"\r  {txt} of {total >> 20} MB" if total else f"\r  {txt}", end="", flush=True)
    print()
    os.replace(parte, destino)


def md5(caminho: str) -> str:
    h = hashlib.md5()
    with open(caminho, "rb") as fh:
        for bloco in iter(lambda: fh.read(8 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def gunzip(origem: str, destino: str):
    with gzip.open(origem, "rb") as src, open(destino + ".part", "wb") as dst:
        shutil.copyfileobj(src, dst, length=16 << 20)
    os.replace(destino + ".part", destino)


def rodar(cmd: list):
    log(" ".join(cmd))
    subprocess.run(cmd, check=True)


def precisa(nome: str, alvo: str, forcar: bool) -> bool:
    if os.path.exists(alvo) and not forcar:
        log(f"{nome}: already installed ({alvo})")
        return False
    return True


def ferramenta(nome: str):
    if not shutil.which(nome):
        sys.exit(f"'{nome}' is not on the PATH. Activate the environment first: conda activate eve_env")


# ---------------------------------------------------------------------------
# Parts
# ---------------------------------------------------------------------------
def instalar_bundle(arquivo: str = None, url: str = None, forcar: bool = False):
    alvo = config.HMM_DB
    if not precisa("viral bundle", alvo, forcar) and os.path.exists(config.VIRAL_DMND):
        return
    url = url or BUNDLE_URL
    if not arquivo:
        if not url:
            sys.exit("The viral bundle has no download address yet. Pass --bundle <file.tar> or "
                     "--bundle-url <address> (or set EVE_BUNDLE_URL).")
        arquivo = os.path.join(config.DATA_DIR, "downloads", os.path.basename(url.split("?")[0]) or "bundle.tar.gz")
        if not os.path.exists(arquivo):
            baixar(url, arquivo)
    if BUNDLE_MD5 and md5(arquivo) != BUNDLE_MD5:
        sys.exit(f"Checksum mismatch for {arquivo}; delete it and download again.")
    log(f"extracting {arquivo}")
    with tarfile.open(arquivo) as tar:
        for m in tar.getmembers():
            if not m.isfile():
                continue
            raiz, _, resto = m.name.partition("/")
            if raiz not in DESTINOS or not resto or ".." in resto.split("/"):
                continue
            destino = os.path.join(DESTINOS[raiz](), resto)
            os.makedirs(os.path.dirname(destino), exist_ok=True)
            with tar.extractfile(m) as src, open(destino + ".part", "wb") as dst:
                shutil.copyfileobj(src, dst, length=16 << 20)
            os.replace(destino + ".part", destino)
            if destino.endswith(".gz") and not destino.endswith(".tar.gz"):
                gunzip(destino, destino[:-3])
                os.remove(destino)
    # The viral DIAMOND database is rebuilt here, so it always matches the local DIAMOND version
    fasta = os.path.join(config.DIAMOND_DIR, "sequences.fasta")
    if os.path.exists(fasta) and (forcar or not os.path.exists(config.VIRAL_DMND)):
        ferramenta("diamond")
        rodar(["diamond", "makedb", "--in", fasta, "-d", config.VIRAL_DMND[:-5], "--quiet"])
    log("viral bundle installed")


def instalar_swissprot(forcar: bool = False):
    alvo = config.GERAL_DMND
    if not precisa("Swiss-Prot", alvo, forcar):
        return
    ferramenta("diamond")
    pasta = os.path.dirname(alvo)
    gz = os.path.join(pasta, "uniprot_sprot.fasta.gz")
    baixar(URLS["swissprot"], gz)
    fasta = gz[:-3]
    gunzip(gz, fasta)
    os.remove(gz)
    rodar(["diamond", "makedb", "--in", fasta, "-d", alvo[:-5], "--quiet"])


def instalar_taxdump(forcar: bool = False):
    if precisa("NCBI taxonomy", config.TAXDUMP, forcar):
        baixar(URLS["taxdump"], config.TAXDUMP)


def instalar_proteoma(acesso: str = PROTEOMA_PADRAO, forcar: bool = False, pasta: str = None):
    pasta = pasta or config.PROTEOMAS_DIR
    alvo = os.path.join(pasta, f"{acesso}_protein.faa")
    if not precisa(f"proteome {acesso}", alvo, forcar):
        return
    os.makedirs(pasta, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=pasta) as tmp:
        zipado = os.path.join(tmp, "proteome.zip")
        if shutil.which("datasets"):
            rodar(["datasets", "download", "genome", "accession", acesso, "--include", "protein",
                   "--filename", zipado, "--no-progressbar"])
        else:
            baixar(f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acesso}/download"
                   f"?include_annotation_type=PROT_FASTA", zipado)
        with zipfile.ZipFile(zipado) as z:
            nomes = [n for n in z.namelist() if n.endswith("protein.faa")]
            if not nomes:
                sys.exit(f"{acesso} has no annotated proteome at NCBI (choose a RefSeq GCF_ accession).")
            with z.open(nomes[0]) as src, open(alvo + ".part", "wb") as dst:
                shutil.copyfileobj(src, dst)
    os.replace(alvo + ".part", alvo)
    log(f"proteome {acesso} installed")


def instalar_pfam(forcar: bool = False):
    alvo = config.PFAM_HMM
    if not precisa("Pfam-A", alvo + ".h3m", forcar):
        return
    gz = alvo + ".gz"
    baixar(URLS["pfam"], gz)
    gunzip(gz, alvo)
    os.remove(gz)
    if shutil.which("hmmpress"):
        rodar(["hmmpress", "-f", alvo])
    else:
        import pyhmmer
        log("hmmpress (pyhmmer)")
        with pyhmmer.plan7.HMMFile(alvo) as fh:
            pyhmmer.hmmer.hmmpress(fh, alvo)


def instalar_cdd(forcar: bool = False):
    pasta = os.path.dirname(config.CDD_RPSDB)
    if not precisa("CDD", config.CDD_RPSDB + ".pal", forcar) and \
            not precisa("CDD", config.CDD_RPSDB + ".rps", forcar):
        return
    tgz = os.path.join(config.DATA_DIR, "downloads", "Cdd_LE.tar.gz")
    if not os.path.exists(tgz):
        baixar(URLS["cdd"], tgz)
    log(f"extracting the CDD into {pasta}")
    os.makedirs(pasta, exist_ok=True)
    with tarfile.open(tgz) as tar:
        tar.extractall(pasta, filter="data") if sys.version_info >= (3, 12) else tar.extractall(pasta)
    os.remove(tgz)
    if not os.path.exists(os.path.join(pasta, "Cdd_domains.tsv")):
        log("note: Cdd_domains.tsv (domain taxa, from the viral bundle) is missing; "
            "the CD-Search annotation will use the domain descriptions only")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", help="Viral bundle already downloaded (.tar from pack.py)")
    ap.add_argument("--bundle-url", help="Download address of the viral bundle")
    ap.add_argument("--proteome", default=PROTEOMA_PADRAO,
                    help="RefSeq accession of the synteny reference proteome (default %(default)s)")
    ap.add_argument("--pfam", action="store_true", help="Also install Pfam-A")
    ap.add_argument("--cdd", action="store_true", help="Also install the local CDD (~4.5 GB)")
    ap.add_argument("--all", action="store_true", help="Everything")
    ap.add_argument("--only", nargs="+", choices=["bundle", "swissprot", "taxdump", "proteome", "pfam", "cdd"],
                    help="Install only these parts")
    ap.add_argument("--force", action="store_true", help="Reinstall what already exists")
    a = ap.parse_args()

    partes = a.only or (["bundle", "swissprot", "taxdump", "proteome"]
                        + (["pfam"] if a.pfam or a.all else []) + (["cdd"] if a.cdd or a.all else []))
    log(f"data folder: {config.DATA_DIR}  (change it in eve_config.ini or with EVE_DATA_DIR)")
    for p in partes:
        if p == "bundle":
            instalar_bundle(a.bundle, a.bundle_url, a.force)
        elif p == "swissprot":
            instalar_swissprot(a.force)
        elif p == "taxdump":
            instalar_taxdump(a.force)
        elif p == "proteome":
            instalar_proteoma(a.proteome, a.force)
        elif p == "pfam":
            instalar_pfam(a.force)
        elif p == "cdd":
            instalar_cdd(a.force)
    log("done. Check everything with: ./eve-tool check")


if __name__ == "__main__":
    main()
