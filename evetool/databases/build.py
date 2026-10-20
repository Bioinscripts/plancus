"""
Rebuilds the viral HMM databases of the EVE tool from a newer RVDB release, with the same
selection criteria as the published bundle (recipe_viral_db.json), so a newer database keeps the
same biases. Optional: the published bundle stays the reference for reproducibility.

The RefSeq viral proteins and Swiss-Prot are NOT downloaded: the build takes the installed files
(or the ones given with --refseq-fasta / --swissprot), so only the HMM sources change.

  ./eve-tool build-databases                       # everything, into <data_dir>/viral_db_builds/<date>
  ./eve-tool build-databases --rvdb-version 32.0   # a given RVDB release instead of the latest
  ./eve-tool build-databases --until merge         # stop before the calibration (quick look)
  ./eve-tool build-databases --out <folder> -t 8   # resume a build: finished steps are skipped
  ./eve-tool build-databases --refseq-fasta my_ncbi_virus.fasta --swissprot uniprot_sprot.fasta.gz

Steps (each one is skipped when its output exists, so an interrupted build resumes):
  refseq     the RefSeq viral proteins (a fixed file, checked for the bundle defline:
             accession |title|molecule type|species|genus|family|length)
  rvdb       RVDB models of viral families, without dsDNA, dsDNA-RT and ssRNA-RT; models with no
             family are kept when their seed sequences are RNA viruses (see the recipe)
  eggnog     eggNOG 5.0 viral models of the RNA (and ssDNA) groups of the recipe
  merge      RVDB + eggNOG into combined_viral_curated.hmm and its annotation table
  inputs     calibration negatives: Swiss-Prot (a fixed file), NCBI taxonomy and the dipteran proteomes
  diamond    DIAMOND database of the RefSeq viral proteins
  calibrate  trusted cut-off of every model (calibrate_hmm.py; the longest step, ~1 day)
  report     build_manifest.json (versions, counts, checksums) and build_report.md, which
             compares the new database with the installed one
  pack       eve_tool_db_<date>.tar, installed with ./eve-tool databases --bundle <file> --force
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
import lzma
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import urllib.parse
import urllib.request

from evetool import config
from evetool.databases.prepare import URLS, baixar, ferramenta, md5

RECEITA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "recipe_viral_db.json")
PASSOS = ["refseq", "rvdb", "eggnog", "merge", "inputs", "diamond", "calibrate", "report", "pack"]
SUFIXOS_VIRAIS = ("virus", "virales", "viridae", "viricetes", "viria", "viricota")


def log(msg):
    print(f"[build {datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


def feito(caminho: str) -> bool:
    if os.path.exists(caminho) and os.path.getsize(caminho) > 0:
        log(f"  already done: {caminho}")
        return True
    return False


class Pastas:
    def __init__(self, raiz):
        self.raiz = raiz
        self.downloads = os.path.join(raiz, "downloads")
        self.inputs = os.path.join(raiz, "inputs")
        self.diamond = os.path.join(raiz, "diamond")
        self.hmm = os.path.join(raiz, "hmm")
        self.rvdb = os.path.join(self.hmm, "RVDB")
        self.eggnog = os.path.join(self.hmm, "EGGNOG")
        self.comb = os.path.join(self.hmm, "COMBINED_DATABASE")
        for p in (self.downloads, self.inputs, self.diamond, self.rvdb, self.eggnog, self.comb):
            os.makedirs(p, exist_ok=True)
        self.fasta = os.path.join(self.diamond, "sequences.fasta")
        self.dmnd = os.path.join(self.diamond, "mydb.dmnd")
        self.taxdump = os.path.join(self.inputs, "taxdump.tar.gz")
        self.swissprot = os.path.join(self.inputs, "uniprot_sprot.fasta")
        self.proteomas = os.path.join(self.inputs, "proteomes")
        self.estado = os.path.join(raiz, "build_state.json")

    def ler_estado(self) -> dict:
        try:
            with open(self.estado) as fh:
                return json.load(fh)
        except Exception:
            return {}

    def gravar_estado(self, **kw):
        e = self.ler_estado()
        e.update(kw)
        with open(self.estado, "w") as fh:
            json.dump(e, fh, indent=1)


# ---------------------------------------------------------------------------
# NCBI taxonomy (taxdump), without ete3
# ---------------------------------------------------------------------------
class Taxonomia:
    def __init__(self, taxdump: str):
        self.pai, self.rank, self.nome, self.por_nome, self.fundidos = {}, {}, {}, {}, {}
        with tarfile.open(taxdump) as tar:
            for l in tar.extractfile("nodes.dmp"):
                p = l.decode().split("\t|\t")
                self.pai[p[0]], self.rank[p[0]] = p[1], p[2]
            for l in tar.extractfile("names.dmp"):
                p = l.decode().split("\t|\t")
                if p[3].startswith("scientific"):
                    self.nome[p[0]] = p[1]
                    self.por_nome[p[1].lower()] = p[0]
                else:
                    self.por_nome.setdefault(p[1].lower(), p[0])
            for l in tar.extractfile("merged.dmp"):
                p = l.decode().rstrip("\t|\n").split("\t|\t")
                self.fundidos[p[0]] = p[1]

    def resolver(self, taxid) -> str:
        t = str(int(taxid)) if str(taxid).strip().lstrip("-").isdigit() else ""
        return self.fundidos.get(t, t)

    def linhagem(self, taxid) -> list:
        t, saida = self.resolver(taxid), []
        while t in self.pai and t != "1" and len(saida) < 100:
            saida.append(t)
            t = self.pai[t]
        return saida

    def familia(self, taxid) -> tuple:
        """(name, family) as the first curation did with ete3: the family-rank ancestor, else the name."""
        t = self.resolver(taxid)
        if not t or t == "0":
            return "Unknown", "Unknown"
        if t not in self.pai:
            return "Not Found", "Not Found"
        nome = self.nome.get(t, "Not Found")
        fam = [self.nome[x] for x in self.linhagem(t) if self.rank.get(x) == "family" and x in self.nome]
        return nome, (fam[0] if fam else nome)


# ---------------------------------------------------------------------------
# Fixed inputs: the RefSeq viral proteins and Swiss-Prot are never downloaded by the build.
# They come from the installed database or from files the researcher gives, so a rebuild changes
# only the HMM sources (RVDB) and keeps the same viral and cellular references.
# ---------------------------------------------------------------------------
def _copiar_fixo(origem: str, destino: str, rotulo: str) -> dict:
    origem = os.path.abspath(os.path.expanduser(origem or ""))
    if not os.path.exists(origem):
        sys.exit(f"{rotulo}: {origem} not found. Pass the file with the option of the build "
                 f"(see ./eve-tool build-databases -h).")
    log(f"  {rotulo}: copying {origem}")
    abrir = gzip.open if origem.endswith(".gz") else open
    with abrir(origem, "rb") as src, open(destino + ".part", "wb") as dst:
        shutil.copyfileobj(src, dst, length=16 << 20)
    os.replace(destino + ".part", destino)
    with open(destino) as fh:
        n = sum(1 for l in fh if l.startswith(">"))
    return {"file": origem, "file_date": f"{datetime.date.fromtimestamp(os.path.getmtime(origem))}",
            "sequences": n, "md5": md5(destino)}


def passo_refseq(p: Pastas, receita: dict, origem: str):
    if feito(p.fasta):
        return
    info = _copiar_fixo(origem, p.fasta, "RefSeq viral proteins")
    # The screening reads the genome type and taxonomy from the defline: check it is the bundle's
    with open(p.fasta) as fh:
        cabecalhos = [l for l in (next(fh, "") for _ in range(4000)) if l.startswith(">")][:200]
    ruins = [c for c in cabecalhos if len(c.split("|")) < 7]
    if not cabecalhos or len(ruins) > len(cabecalhos) // 10:
        os.remove(p.fasta)
        campos = " | ".join(receita["refseq_viral"]["defline"])
        sys.exit(f"{origem} does not have the defline of the bundle ({campos}); "
                 f"e.g. {(ruins or [''])[0].strip()[:120]}\n"
                 f"How to download it from NCBI Virus: docs/refseq_viral_proteins.md")
    p.gravar_estado(refseq=info)
    log(f"  {info['sequences']:,} proteins -> {p.fasta}")


def genomas_por_familia_ncbi(fasta: str) -> dict:
    """Most common NCBI Virus molecule type of each family (fallback for the curated table)."""
    c = collections.defaultdict(collections.Counter)
    with open(fasta) as fh:
        for l in fh:
            if l.startswith(">"):
                partes = l.rstrip("\n").split("|")
                if len(partes) >= 7 and partes[5] and partes[2]:
                    c[partes[5]][partes[2]] += 1
    return {f: k.most_common(1)[0][0] for f, k in c.items()}


# ---------------------------------------------------------------------------
# rvdb
# ---------------------------------------------------------------------------
def versao_rvdb(url: str, pedida: str) -> str:
    if pedida and pedida != "latest":
        return pedida
    with urllib.request.urlopen(url.replace("/files/", "/"), timeout=60) as resp:
        html = resp.read().decode(errors="replace")
    versoes = sorted({v for v in re.findall(r"U-RVDBv([\d.]+)-prot\.hmm\.xz", html)},
                     key=lambda v: [int(x) for x in v.split(".")])
    if not versoes:
        sys.exit(f"No RVDB release found at {url}")
    return versoes[-1]


def _xz(origem: str, destino: str):
    with lzma.open(origem) as src, open(destino + ".part", "wb") as dst:
        shutil.copyfileobj(src, dst, length=16 << 20)
    os.replace(destino + ".part", destino)


def passo_rvdb(p: Pastas, receita: dict, versao: str):
    r = receita["rvdb"]
    saida_hmm = os.path.join(p.rvdb, f"U-RVDBv{versao}-prot-curated.hmm")
    if feito(saida_hmm):
        return
    arquivos = [f.format(v=versao) for f in r["files"]]
    for f in arquivos:
        destino = os.path.join(p.downloads, f)
        if not os.path.exists(destino):
            baixar(r["url"] + f, destino)
    hmm_bruto = os.path.join(p.rvdb, f"U-RVDBv{versao}-prot.hmm")
    sqlite_ = os.path.join(p.rvdb, f"U-RVDBv{versao}-prot-hmm.sqlite")
    annot = os.path.join(p.rvdb, "annot")
    if not os.path.exists(hmm_bruto):
        log("  decompressing the RVDB HMMs")
        _xz(os.path.join(p.downloads, arquivos[0]), hmm_bruto)
    if not os.path.exists(sqlite_):
        _xz(os.path.join(p.downloads, arquivos[2]), sqlite_)
    if not os.path.isdir(annot):
        log("  extracting the RVDB model annotations")
        os.makedirs(annot + ".part", exist_ok=True)
        with tarfile.open(os.path.join(p.downloads, arquivos[1]), "r:xz") as tar:
            for m in tar.getmembers():
                if m.isfile() and "/annot/" in "/" + m.name and m.name.endswith(".txt") and ".." not in m.name:
                    with tar.extractfile(m) as src, open(os.path.join(annot + ".part", os.path.basename(m.name)), "wb") as dst:
                        shutil.copyfileobj(src, dst)
        if not os.listdir(annot + ".part"):
            sys.exit(f"No model annotation (annot/*.txt) found in {arquivos[1]}")
        os.replace(annot + ".part", annot)

    # 1. base annotations: model, LCA taxid and the keywords of its seed names
    con = sqlite3.connect(sqlite_)
    tabela = r["keywords_table"]
    kw = collections.defaultdict(list)
    for fam, s in con.execute(f"SELECT t.famId, k.str FROM {tabela} t JOIN keyword k ON k.id = t.kwId ORDER BY t.famId, t.rowid"):
        kw[fam].append(s)
    base = [(f"FAM{fid:06d}", lca, ",".join(kw.get(fid, []))) for fid, lca in con.execute("SELECT id, LCAtaxid FROM family ORDER BY id")]
    con.close()
    with open(os.path.join(p.rvdb, "rvdb_base_annotations.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["HMM_ID", "LCAtaxid", "Protein_Keywords"])
        w.writerows(base)

    # 2. taxonomy and genome type
    tax = Taxonomia(p.taxdump)
    mapa_curado = {}
    with open(os.path.join(p.rvdb, "rvdb_family_genome_mapping.tsv")) as fh:
        for linha in csv.DictReader(fh, delimiter="\t"):
            mapa_curado[linha["Family"]] = linha["Genome_Type"]
    mapa_ncbi = genomas_por_familia_ncbi(p.fasta)

    def genoma(familia):
        return mapa_curado.get(familia) or mapa_ncbi.get(familia) or "Unknown"

    resgate = r["rescue_without_family"]
    rna, retro = set(resgate["rna_taxa"]), set(resgate["retro_taxa"])

    def sementes(hmm_id):
        """Taxonomic class and family genome type of the seed organisms of a model."""
        caminho = os.path.join(annot, f"{int(hmm_id[3:])}.txt")
        c, genomas = collections.Counter(), collections.Counter()
        if not os.path.exists(caminho):
            return c, genomas
        with open(caminho, errors="replace") as fh:
            bloco = re.split(r"^SEQUENCES:\s*$", fh.read(), maxsplit=1, flags=re.MULTILINE)
        for titulo in (bloco[1].splitlines() if len(bloco) > 1 else []):
            m = re.search(r"\[([^\]]+)\]\s*$", titulo)
            t = tax.por_nome.get(m.group(1).lower()) if m else None
            if not t:
                c["unknown"] += 1
                continue
            lin = set(tax.linhagem(t))
            classe = "retro" if lin & retro else "rna" if lin & rna else "virus" if "10239" in lin else "cellular"
            c[classe] += 1
            if classe == "rna":
                fam = [tax.nome[x] for x in tax.linhagem(t) if tax.rank.get(x) == "family" and x in tax.nome]
                if fam and genoma(fam[0]) != "Unknown":
                    genomas[genoma(fam[0])] += 1
        return c, genomas

    linhas, motivos = [], collections.Counter()
    for hmm_id, lca, palavras in base:
        nome, familia = tax.familia(lca)
        g = genoma(familia)
        if str(familia).lower().endswith(SUFIXOS_VIRAIS):
            if g in r["drop_genome_types"]:
                motivos["drop: genome type"] += 1
                continue
            motivos["keep: viral family"] += 1
        else:
            c, gs = sementes(hmm_id)
            classificadas = c["rna"] + c["retro"] + c["virus"] + c["cellular"]
            if (c["rna"] and c["retro"] <= resgate["max_retro_seeds"]
                    and c["rna"] / classificadas >= resgate["min_rna_seed_fraction"]):
                g = gs.most_common(1)[0][0] if gs else resgate["default_genome_type"]
                motivos["keep: RNA-virus seeds (no family)"] += 1
            else:
                motivos["drop: no viral family"] += 1
                continue
        linhas.append({"HMM_ID": hmm_id, "LCAtaxid": lca, "Tax_Name": nome, "Family": familia,
                       "Genome_Type": g, "Protein_Keywords": palavras})
    log("  RVDB selection: " + ", ".join(f"{v} {k}" for k, v in motivos.items()))
    tsv = os.path.join(p.rvdb, "RVDB_Full_Annotations.tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t")
        w.writeheader()
        w.writerows(linhas)
    manter = {l["HMM_ID"] for l in linhas}
    with open(os.path.join(p.rvdb, "modelos_curados_ids.txt"), "w") as fh:
        fh.writelines(f"{l['HMM_ID']}\n" for l in linhas)

    # 3. the selected models
    import pyhmmer
    n = 0
    with pyhmmer.plan7.HMMFile(hmm_bruto) as hf, open(saida_hmm + ".part", "wb") as out:
        for hmm in hf:
            nome = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
            if nome in manter:
                hmm.write(out)
                n += 1
    os.replace(saida_hmm + ".part", saida_hmm)
    p.gravar_estado(rvdb={"version": versao, "models_total": len(base), "models_kept": n, "selection": dict(motivos)})
    log(f"  {n} RVDB models -> {saida_hmm}")


# ---------------------------------------------------------------------------
# eggnog
# ---------------------------------------------------------------------------
def passo_eggnog(p: Pastas, receita: dict):
    e = receita["eggnog"]
    saida = os.path.join(p.eggnog, "eggnog_virus_db.hmm")
    if feito(saida):
        return
    anotacoes = []
    with open(saida + ".part", "w") as out:
        for taxid in sorted(e["groups"]):
            nome_tax, rank, genoma = e["groups"][taxid]
            tgz = os.path.join(p.downloads, f"{taxid}_hmms.tar.gz")
            if not os.path.exists(tgz):
                baixar(e["url"].format(taxid=taxid), tgz)
            with tarfile.open(tgz) as tar:
                for m in sorted((m for m in tar.getmembers() if m.isfile() and ".hmm" in m.name), key=lambda m: m.name):
                    texto = tar.extractfile(m).read()
                    if m.name.endswith(".gz"):
                        texto = gzip.decompress(texto)
                    modelo, original, leng, nseq = None, None, None, None
                    for linha in texto.decode(errors="ignore").splitlines(keepends=True):
                        if linha.startswith("NAME"):
                            original = linha.split()[1]
                            modelo = original.split(".")[0]
                            linha = f"NAME  {modelo}\n"
                        elif linha.startswith("LENG"):
                            leng = int(linha.split()[1])
                        elif linha.startswith("NSEQ"):
                            nseq = int(linha.split()[1])
                        out.write(linha)
                    anotacoes.append({"HMM_ID": modelo, "LCAtaxid": taxid, "Tax_Name": nome_tax, "Tax_Rank": rank,
                                      "Family": nome_tax, "Genome_Type": genoma, "Length": leng, "NSEQ": nseq,
                                      "Original_NAME": original, "Original_File": f"{taxid}_hmms/{m.name}"})
    with open(os.path.join(p.eggnog, "eggnog_virus_annotations.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(anotacoes[0]), delimiter="\t")
        w.writeheader()
        w.writerows(anotacoes)
    os.replace(saida + ".part", saida)
    p.gravar_estado(eggnog={"version": e["version"], "models": len(anotacoes)})
    log(f"  {len(anotacoes)} eggNOG models -> {saida}")


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------
COLUNAS_COMBINADAS = ["HMM_ID", "Source_DB", "LCAtaxid", "Tax_Name", "Family", "Genome_Type", "Protein_Keywords",
                      "Length", "Original_File"]


def passo_merge(p: Pastas, versao: str):
    saida = os.path.join(p.comb, "combined_viral_curated.hmm")
    if feito(saida):
        return
    fontes = [("RVDB", os.path.join(p.rvdb, f"U-RVDBv{versao}-prot-curated.hmm"), os.path.join(p.rvdb, "RVDB_Full_Annotations.tsv")),
              ("eggNOG", os.path.join(p.eggnog, "eggnog_virus_db.hmm"), os.path.join(p.eggnog, "eggnog_virus_annotations.tsv"))]
    linhas = []
    with open(saida + ".part", "wb") as out:
        for fonte, hmm, tsv in fontes:
            with open(hmm, "rb") as src:
                shutil.copyfileobj(src, out, length=64 << 20)
            with open(tsv) as fh:
                for l in csv.DictReader(fh, delimiter="\t"):
                    linhas.append({c: (fonte if c == "Source_DB" else l.get(c, "")) for c in COLUNAS_COMBINADAS})
    with open(os.path.join(p.comb, "combined_viral_annotations.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUNAS_COMBINADAS, delimiter="\t")
        w.writeheader()
        w.writerows(linhas)
    with open(os.path.join(p.comb, "modelos_unificados_ids.txt"), "w") as fh:
        fh.writelines(f"{l['HMM_ID']}\n" for l in linhas)
    os.replace(saida + ".part", saida)
    log(f"  {len(linhas)} models -> {saida}")


# ---------------------------------------------------------------------------
# inputs and diamond
# ---------------------------------------------------------------------------
def passo_taxdump(p: Pastas):
    if not feito(p.taxdump):
        baixar(URLS["taxdump"], p.taxdump)


def passo_inputs(p: Pastas, receita: dict, swissprot: str, versao_sp: str = ""):
    passo_taxdump(p)
    if not feito(p.swissprot):
        info = _copiar_fixo(swissprot, p.swissprot, "Swiss-Prot")
        info["release"] = versao_sp or f"file {os.path.basename(info['file'])} of {info['file_date']}"
        p.gravar_estado(swissprot=info)
    from evetool.databases.prepare import instalar_proteoma
    for acc in receita["calibration"]["negatives"]["proteomes"]:
        instalar_proteoma(acc, pasta=p.proteomas)


def passo_diamond(p: Pastas):
    if feito(p.dmnd):
        return
    ferramenta("diamond")
    subprocess.run(["diamond", "makedb", "--in", p.fasta, "-d", p.dmnd[:-5], "--quiet"], check=True)


# ---------------------------------------------------------------------------
# calibrate
# ---------------------------------------------------------------------------
def passo_calibrate(p: Pastas, receita: dict, threads: int, teste: int = 0):
    final = os.path.join(p.comb, "combined_viral_calibrado.hmm")
    tabela = os.path.join(p.comb, "calibracao", "4_model_calibration.tsv")
    if feito(final) and os.path.exists(tabela):
        return
    par = dict(receita["calibration"]["params"])
    par["threads"] = threads or par["threads"]
    cmd = [sys.executable, "-m", "evetool.databases.calibrate_hmm",
           "--hmm", os.path.join(p.comb, "combined_viral_curated.hmm"),
           "--anotacoes", os.path.join(p.comb, "combined_viral_annotations.tsv"),
           "--rvdb-annot-dir", os.path.join(p.rvdb, "annot"), "--taxdump", p.taxdump,
           "--proteomas-dir", p.proteomas, "--swissprot", p.swissprot,
           "--refseq-viral-fasta", p.fasta, "--refseq-viral-dmnd", p.dmnd,
           "--saida", os.path.join(p.comb, "calibracao"),
           "--threads", str(par["threads"]), "--lote", str(par["lote"]), "--t-report", str(par["t_report"]),
           "--margem", str(par["margem"]), "--piso", str(par["piso"]), "--min-pos", str(par["min_pos"]),
           "--id-negativos", str(par["id_negativos"]), "--chunk", str(par["chunk"])]
    ambiente = dict(os.environ, PYTHONPATH=config.PASTA_CODIGO + os.pathsep + os.environ.get("PYTHONPATH", ""))
    if teste:                             # a few models only, to check the inputs; nothing is written as final
        subprocess.run(cmd[:-2] + ["--max-modelos", str(teste)], env=ambiente, check=True)
        log(f"  test calibration of {teste} models: {os.path.join(p.comb, 'calibracao', '4_model_calibration.tsv.parcial')}")
        return
    # One process per chunk of models: the memory of a single long process keeps growing
    for bloco in range(1, 2001):
        log(f"  calibration chunk {bloco}")
        rc = subprocess.run(cmd, env=ambiente).returncode
        if os.path.exists(tabela) and os.path.exists(final):
            log(f"  calibration finished -> {final}")
            return
        if rc != 0:
            log(f"  chunk ended with code {rc}; trying again")
    sys.exit("The calibration did not finish; run the build again to resume it.")


# ---------------------------------------------------------------------------
# report and pack
# ---------------------------------------------------------------------------
def _ids_hmm(caminho: str) -> dict:
    """Model -> trusted cut-off (sequence) of a calibrated .hmm."""
    import pyhmmer
    saida = {}
    if not os.path.exists(caminho):
        return saida
    with pyhmmer.plan7.HMMFile(caminho) as hf:
        for hmm in hf:
            nome = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
            tc = hmm.cutoffs.trusted
            saida[nome] = round(tc[0], 1) if tc else None
    return saida


def _fasta_resumo(caminho: str) -> tuple:
    acc, genomas = set(), collections.Counter()
    if not os.path.exists(caminho):
        return acc, genomas
    with open(caminho) as fh:
        for l in fh:
            if l.startswith(">"):
                partes = l[1:].split("|")
                acc.add(partes[0].split()[0])
                genomas[partes[2] if len(partes) > 2 else "?"] += 1
    return acc, genomas


def _anotacoes(caminho: str) -> dict:
    if not os.path.exists(caminho):
        return {}
    with open(caminho) as fh:
        return {l["HMM_ID"]: l for l in csv.DictReader(fh, delimiter="\t")}


def passo_report(p: Pastas, receita: dict):
    novo_hmm = os.path.join(p.comb, "combined_viral_calibrado.hmm")
    instalado_hmm = config.HMM_DB
    estado = p.ler_estado()
    log("  comparing with the installed database")
    tc_novo, tc_velho = _ids_hmm(novo_hmm), _ids_hmm(instalado_hmm)
    an_novo = _anotacoes(os.path.join(p.comb, "combined_viral_annotations.tsv"))
    an_velho = _anotacoes(os.path.join(config.HMM_DIR, "COMBINED_DATABASE", "combined_viral_annotations.tsv"))
    acc_novo, gen_novo = _fasta_resumo(p.fasta)
    acc_velho, gen_velho = _fasta_resumo(config.diamond_path("sequences.fasta"))

    entraram = sorted(set(tc_novo) - set(tc_velho))
    sairam = sorted(set(tc_velho) - set(tc_novo))
    comuns = set(tc_novo) & set(tc_velho)
    mudou_tc = [m for m in comuns if tc_novo[m] is not None and tc_velho[m] is not None and abs(tc_novo[m] - tc_velho[m]) >= 1]
    with open(os.path.join(p.raiz, "build_model_changes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["HMM_ID", "change", "source", "family", "genome_type", "TC_installed", "TC_new"])
        for m in entraram + sairam + sorted(mudou_tc):
            a = an_novo.get(m) or an_velho.get(m) or {}
            mud = "added" if m in entraram else "removed" if m in sairam else "TC changed"
            w.writerow([m, mud, a.get("Source_DB", ""), a.get("Family", ""), a.get("Genome_Type", ""),
                        tc_velho.get(m, ""), tc_novo.get(m, "")])

    def por(campo, ids, an):
        return collections.Counter((an.get(m) or {}).get(campo, "?") for m in ids)

    linhas = [f"# Viral database build {os.path.basename(p.raiz)}", "",
              f"Built {datetime.datetime.now():%Y-%m-%d %H:%M} with recipe v{receita['recipe_version']} "
              f"(evetool/databases/recipe_viral_db.json). Compared with the installed database in {config.HMM_DIR}.", "",
              "## Sources", "",
              f"- RVDB {estado.get('rvdb', {}).get('version', '?')}: {estado.get('rvdb', {}).get('models_kept', '?')} of "
              f"{estado.get('rvdb', {}).get('models_total', '?')} models selected",
              f"- eggNOG {receita['eggnog']['version']}: {estado.get('eggnog', {}).get('models', '?')} models",
              f"- RefSeq viral proteins (fixed file {estado.get('refseq', {}).get('file', '?')}, "
              f"{estado.get('refseq', {}).get('file_date', '?')}): {len(acc_novo):,}",
              f"- Swiss-Prot (fixed file): {estado.get('swissprot', {}).get('release', '(not recorded)')}", "",
              "## HMM models after calibration", "",
              "| | installed | new |", "|---|---|---|",
              f"| models | {len(tc_velho)} | {len(tc_novo)} |",
              f"| added / removed | | +{len(entraram)} / −{len(sairam)} |",
              f"| trusted cut-off changed by ≥ 1 bit | | {len(mudou_tc)} of {len(comuns)} shared |", "",
              "| genome type | installed | new |", "|---|---|---|"]
    g_v, g_n = por("Genome_Type", tc_velho, an_velho), por("Genome_Type", tc_novo, an_novo)
    linhas += [f"| {g} | {g_v.get(g, 0)} | {g_n.get(g, 0)} |" for g in sorted(set(g_v) | set(g_n))]
    linhas += ["", "## RefSeq viral proteins (DIAMOND)", "",
               "| | installed | new |", "|---|---|---|",
               f"| proteins | {len(acc_velho):,} | {len(acc_novo):,} |",
               f"| added / removed | | +{len(acc_novo - acc_velho):,} / −{len(acc_velho - acc_novo):,} |"]
    for g in sorted(set(gen_velho) | set(gen_novo), key=lambda g: -(gen_novo.get(g, 0))):
        if gen_velho.get(g, 0) + gen_novo.get(g, 0) >= 50:
            linhas.append(f"| {g} | {gen_velho.get(g, 0):,} | {gen_novo.get(g, 0):,} |")
    linhas += ["", "Every added, removed or re-calibrated model: build_model_changes.tsv. "
               "A different database finds a different set of EVEs: record the database version with the results."]
    with open(os.path.join(p.raiz, "build_report.md"), "w") as fh:
        fh.write("\n".join(linhas) + "\n")

    manifesto = {"build": os.path.basename(p.raiz), "date": f"{datetime.datetime.now():%Y-%m-%d %H:%M}",
                 "recipe": receita, "sources": estado,
                 "counts": {"hmm_models": len(tc_novo), "refseq_proteins": len(acc_novo)},
                 "md5": {os.path.relpath(f, p.raiz): md5(f) for f in (novo_hmm, p.fasta) if os.path.exists(f)}}
    with open(os.path.join(p.comb, "build_manifest.json"), "w") as fh:
        json.dump(manifesto, fh, indent=1)
    if os.path.exists(novo_hmm):
        from evetool.annotation.hmm_domains import carregar_tc
        carregar_tc(novo_hmm)                 # writes combined_viral_calibrado_tc.json, packed with the bundle
    log(f"  report -> {os.path.join(p.raiz, 'build_report.md')}")


def passo_pack(p: Pastas):
    from evetool.databases.pack import empacotar
    cdd = os.path.join(p.hmm, "CDD", "rpsdb", "Cdd_domains.tsv")
    if not os.path.exists(cdd) and os.path.exists(config.hmm_path("CDD", "rpsdb", "Cdd_domains.tsv")):
        os.makedirs(os.path.dirname(cdd), exist_ok=True)
        shutil.copy2(config.hmm_path("CDD", "rpsdb", "Cdd_domains.tsv"), cdd)
    saida = os.path.join(p.raiz, f"eve_tool_db_{os.path.basename(p.raiz)}.tar")
    empacotar(saida, p.hmm, p.diamond)
    log(f"  install it with: ./eve-tool databases --bundle {saida} --force")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="Build folder (default: <data_dir>/viral_db_builds/<date>); an existing one is resumed")
    ap.add_argument("--rvdb-version", default="", help="RVDB release, e.g. 32.0 (default: the recipe's, i.e. the latest)")
    ap.add_argument("--steps", nargs="+", choices=PASSOS, help="Run only these steps")
    ap.add_argument("--until", choices=PASSOS, help="Stop after this step")
    ap.add_argument("-t", "--threads", type=int, default=0, help="Threads for the calibration (default: the recipe's)")
    ap.add_argument("--calibration-test", type=int, default=0, metavar="N",
                    help="Calibrate only N models and stop (checks the inputs in minutes instead of a day)")
    ap.add_argument("--refseq-fasta", default=config.diamond_path("sequences.fasta"),
                    help="RefSeq viral proteins with the bundle defline (not downloaded: see "
                         "docs/refseq_viral_proteins.md; default: the installed %(default)s)")
    ap.add_argument("--swissprot", default=config.diamond_path("geral", "uniprot_sprot.fasta"),
                    help="Swiss-Prot FASTA (.fasta or .fasta.gz) for the calibration negatives (not downloaded; "
                         "default: the installed %(default)s)")
    ap.add_argument("--swissprot-release", default="", help="Release label recorded in the manifest, e.g. 2026_02")
    a = ap.parse_args()

    with open(RECEITA) as fh:
        receita = json.load(fh)
    raiz = os.path.abspath(a.out or os.path.join(config.DATA_DIR, "viral_db_builds", f"{datetime.date.today()}"))
    p = Pastas(raiz)
    versao = p.ler_estado().get("rvdb", {}).get("version") or versao_rvdb(receita["rvdb"]["url"], a.rvdb_version or receita["rvdb"]["version"])
    passos = a.steps or (PASSOS[:PASSOS.index(a.until) + 1] if a.until else PASSOS)
    log(f"build folder {raiz} · RVDB {versao} · eggNOG {receita['eggnog']['version']} · steps: {' '.join(passos)}")

    # The curated family -> genome table is part of the recipe: it starts from the installed one
    mapa = os.path.join(p.rvdb, "rvdb_family_genome_mapping.tsv")
    if not os.path.exists(mapa):
        origem = config.hmm_path("RVDB", "rvdb_family_genome_mapping.tsv")
        if not os.path.exists(origem):
            sys.exit(f"Missing {origem}: install the viral bundle first (./eve-tool databases).")
        shutil.copy2(origem, mapa)

    for passo in passos:
        log(f"step {passo}")
        if passo == "refseq":
            passo_refseq(p, receita, a.refseq_fasta)
        elif passo == "rvdb":
            passo_taxdump(p)
            passo_rvdb(p, receita, versao)
        elif passo == "eggnog":
            passo_eggnog(p, receita)
        elif passo == "merge":
            passo_merge(p, versao)
        elif passo == "inputs":
            passo_inputs(p, receita, a.swissprot, a.swissprot_release)
        elif passo == "diamond":
            passo_diamond(p)
        elif passo == "calibrate":
            passo_calibrate(p, receita, a.threads, a.calibration_test)
            if a.calibration_test:
                break
        elif passo == "report":
            passo_report(p, receita)
        elif passo == "pack":
            passo_pack(p)
    log("done")


if __name__ == "__main__":
    main()
