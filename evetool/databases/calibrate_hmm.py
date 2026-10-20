#!/usr/bin/env python3
"""
Calibration of the viral HMM models (COMBINED_DATABASE: curated RVDB + viral eggNOG).

Stages (each one reuses the previous output when it already exists):
  1. selection : keep the models of non-retroviral RNA viruses and of ssDNA viruses.
                 Drop dsDNA viruses, Anelloviridae and retroelements (Pararnavirae/Ortervirales,
                 or a retroelement-taxon seed) and models with >= 50% cellular seeds.
  2. negatives : Swiss-Prot without viral proteins + dipteran proteomes. Dipteran proteins
                 whose best RefSeq viral hit is a target virus (possibly domesticated EVEs) are
                 set aside and judged separately, so they do not inflate the cutoff.
  3. positives : RefSeq viral proteins of non-retroviral RNA and ssDNA viruses (no Anelloviridae).
  4. calibrate : for each model, the highest score on a cellular protein sets the TC cutoff
                 (stored in the .hmm). Models that detect no viral protein above their own TC
                 are removed.

Usage:
  python -m evetool.databases.calibrate_hmm --threads 7
  python -m evetool.databases.calibrate_hmm --threads 7 --max-modelos 50   # quick test
"""
import os
import re
import sys
import csv
import gzip
import time
import tarfile
import argparse
import subprocess
from collections import defaultdict

import pandas as pd
import pyhmmer
from evetool import config

BASE_HMM = config.HMM_DIR
PADROES = {
    "hmm": f"{BASE_HMM}/COMBINED_DATABASE/combined_viral_curated.hmm",
    "anotacoes": f"{BASE_HMM}/COMBINED_DATABASE/combined_viral_annotations.tsv",
    "rvdb_annot_dir": f"{BASE_HMM}/RVDB/U-RVDBv32.0-prot-hmm-txt/annot",
    "taxdump": f"{BASE_HMM}/RVDB/taxdump.tar.gz",
    "proteomas_dir": f"{BASE_HMM}/proteomas_calibracao",
    "swissprot": config.diamond_path("geral", "uniprot_sprot.fasta"),
    "refseq_viral_fasta": config.diamond_path("sequences.fasta"),
    "refseq_viral_dmnd": config.diamond_path("mydb.dmnd"),
    "saida": f"{BASE_HMM}/COMBINED_DATABASE/calibracao",
}

REINOS_DSDNA = ("Duplodnaviria", "Varidnaviria", "Adnaviria")
FAMILIAS_EXCLUIDAS = ("Anelloviridae",)
TAXONS_RETRO = ("Pararnavirae", "Ortervirales", "Metaviridae", "Pseudoviridae", "Retroviridae",
                "Belpaoviridae", "Caulimoviridae", "Hepadnaviridae")
# Organismo da semente pertence a um taxon de retroelemento. "reverse transcriptase-like" sozinho
# NOT included: many RNA-virus RdRps (partiti, noda, tombus) are annotated that way in GenBank.
REGEX_ORGANISMO_RETRO = re.compile(
    r"errantivir|metavir|pseudovir|semotivir|belpaovir|hemivir|chromovir|retrovir|"
    r"\bTED virus|\bgypsy|\bcopia\b|\bTy[13]\b",
    re.IGNORECASE,
)
REGEX_PRODUTO_RETRO = re.compile(r"integrase|retrotranspos|transposon|gag-pol polyprotein from", re.IGNORECASE)
REGEX_ORGANISMO_VIRAL = re.compile(r"vir(us|al|idae|ales|icota|icetes|ia)|phage|viroid|\bsp\.\]", re.IGNORECASE)
FRACAO_SEMENTES_CELULARES = 0.5
GENOMAS_ALVO = {"ssRNA(+)", "ssRNA(-)", "ssRNA(+/-)", "ssRNA", "dsRNA", "RNA",
                "ssDNA", "ssDNA(+)", "ssDNA(-)", "ssDNA(+/-)"}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def nome(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def eh_genoma_alvo(genoma) -> bool:
    """Non-retroviral RNA or ssDNA (every part of a mixed label must be on target)."""
    if not isinstance(genoma, str):
        return False
    partes = [p.strip() for p in genoma.split(";")]
    return all(p in GENOMAS_ALVO for p in partes)


# ---------------------------------------------------------------------------
# 1. Model selection
# ---------------------------------------------------------------------------
def ler_annot_rvdb(annot_dir: str, hmm_id: str) -> dict:
    """Read the LCA, keywords and seed titles of an RVDB model (FAM000469 -> 469.txt)."""
    caminho = os.path.join(annot_dir, f"{int(hmm_id[3:])}.txt")
    info = {"lca": "", "n_sementes": 0, "n_retro": 0, "n_celulares": 0}
    if not os.path.exists(caminho):
        return info
    with open(caminho, errors="replace") as fh:
        conteudo = fh.read()
    m = re.search(r"^LCA\t(.*)$", conteudo, re.MULTILINE)
    info["lca"] = m.group(1) if m else ""
    # Final "SEQUENCES:" block (not to be confused with "KEYWORDS FROM SEQUENCES:")
    bloco = re.split(r"^SEQUENCES:\s*$", conteudo, maxsplit=1, flags=re.MULTILINE)
    sementes = [l for l in bloco[1].splitlines() if l.strip()] if len(bloco) > 1 else []
    for titulo in sementes:
        m_org = re.search(r"\[([^\]]+)\]\s*$", titulo)
        organismo = m_org.group(1) if m_org else ""
        if REGEX_ORGANISMO_RETRO.search(organismo) or REGEX_PRODUTO_RETRO.search(titulo):
            info["n_retro"] += 1
        elif organismo and not REGEX_ORGANISMO_VIRAL.search(f"[{organismo}]"):
            info["n_celulares"] += 1
    info["n_sementes"] = len(sementes)
    return info


def selecionar_modelos(cfg) -> pd.DataFrame:
    destino = os.path.join(cfg.saida, "1_model_selection.tsv")
    if os.path.exists(destino):
        log(f"Selection already exists: {destino}")
        return pd.read_csv(destino, sep="\t")
    # A selection written before the columns were renamed is still usable: the column
    # names and the keep/drop verbs are mapped on the way in, so a calibration started
    # with the old layout can be resumed instead of redone.
    antigo = os.path.join(cfg.saida, "1_selecao_modelos.tsv")
    if os.path.exists(antigo):
        log(f"Reusing the selection written with the previous column names: {antigo}")
        sel = pd.read_csv(antigo, sep="\t").rename(columns={
            "genoma": "genome", "n_sementes": "n_seeds", "sementes_retro": "retro_seeds",
            "sementes_celulares": "cellular_seeds", "selecao": "selection"})
        sel["selection"] = (sel["selection"].astype(str)
                            .str.replace(r"^manter", "keep", regex=True)
                            .str.replace(r"^remover", "drop", regex=True))
        sel.to_csv(destino, sep="\t", index=False)
        return sel

    df = pd.read_csv(cfg.anotacoes, sep="\t")
    linhas = []
    for _, r in df.iterrows():
        hmm_id, source = r["HMM_ID"], r["Source_DB"]
        genoma, family = r["Genome_Type"], str(r["Family"])
        lca, reason = "", None
        info = {"n_sementes": 0, "n_retro": 0, "n_celulares": 0}

        if source == "RVDB":
            info = ler_annot_rvdb(cfg.rvdb_annot_dir, hmm_id)
            lca = info["lca"]
            if any(f in lca or f == family for f in FAMILIAS_EXCLUIDAS):
                reason = "drop:excluded_family"
            elif any(t in lca for t in REINOS_DSDNA) or "dsDNA" in str(genoma):
                reason = "drop:dsDNA_virus"
            elif any(t in lca for t in TAXONS_RETRO):
                reason = "drop:retroelement_LCA"
            elif info["n_retro"] > 0:
                reason = "drop:retroelement_seeds"
            elif info["n_sementes"] and info["n_celulares"] / info["n_sementes"] >= FRACAO_SEMENTES_CELULARES:
                reason = "drop:cellular_seeds"
            elif "Orthornavirae" in lca or ("Monodnaviria" in lca and eh_genoma_alvo(genoma)):
                reason = "keep"
            elif eh_genoma_alvo(genoma):
                reason = "keep:broad_LCA"
            else:
                reason = "drop:genome_out_of_scope"
        else:  # eggNOG: no keywords, so decide by the genome type of the taxonomic group
            if eh_genoma_alvo(genoma):
                amplo = family in ("ssRNA viruses", "dsRNA viruses") or "unclassified" in family
                reason = "keep:broad_LCA" if amplo else "keep"
            else:
                reason = "drop:genome_out_of_scope"

        linhas.append({"HMM_ID": hmm_id, "source": source, "family": family, "genome": genoma,
                       "n_seeds": info["n_sementes"], "retro_seeds": info["n_retro"],
                       "cellular_seeds": info["n_celulares"], "lca": lca, "selection": reason})

    sel = pd.DataFrame(linhas)
    sel.to_csv(destino, sep="\t", index=False)
    log("Model selection:\n" + sel["selection"].value_counts().to_string())
    return sel


# ---------------------------------------------------------------------------
# 2. Conjunto negativo
# ---------------------------------------------------------------------------
def taxids_virais(taxdump: str) -> set:
    filhos = defaultdict(list)
    with tarfile.open(taxdump) as tar:
        fh = tar.extractfile("nodes.dmp")
        for linha in fh:
            p = linha.decode().split("\t|\t")
            filhos[p[1]].append(p[0])
    virais, pilha = set(), ["10239"]
    while pilha:
        t = pilha.pop()
        virais.add(t)
        pilha.extend(filhos.get(t, []))
    return virais


def escrever_fasta_filtrado(entrada, saida_fh, manter, prefixo=""):
    """Copy records from 'entrada' to 'saida_fh' when manter(header) is True."""
    n, gravar = 0, False
    with open(entrada) as fh:
        for linha in fh:
            if linha.startswith(">"):
                gravar = manter(linha)
                if gravar:
                    n += 1
                    saida_fh.write(">" + prefixo + linha[1:])
            elif gravar:
                saida_fh.write(linha)
    return n


def montar_negativos(cfg) -> tuple:
    pct = int(cfg.id_negativos * 100)
    neg_final = os.path.join(cfg.saida, f"2_negativos_nr{pct}.faa")
    eve_like = os.path.join(cfg.saida, "2_dipteros_hit_virus_alvo.faa")
    if os.path.exists(neg_final) and os.path.exists(eve_like):
        log(f"Negatives already exist: {neg_final}")
        return neg_final, eve_like

    # 2a. Dipteran proteins whose best RefSeq viral hit is an RNA virus
    dipteros = os.path.join(cfg.saida, "2_dipteros_todos.faa")
    proteomas = sorted(f for f in os.listdir(cfg.proteomas_dir) if f.endswith("_protein.faa"))
    with open(dipteros, "w") as out:
        for f in proteomas:
            acc = f.replace("_protein.faa", "")
            escrever_fasta_filtrado(os.path.join(cfg.proteomas_dir, f), out, lambda h: True, prefixo=f"{acc}|")

    tsv = os.path.join(cfg.saida, "2_dipteros_vs_refseq_viral.tsv")
    if not os.path.exists(tsv):
        log("DIAMOND blastp: dipteran proteomes vs RefSeq viral...")
        subprocess.run(["diamond", "blastp", "-q", dipteros, "-d", cfg.refseq_viral_dmnd, "-o", tsv,
                        "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "stitle",
                        "--max-target-seqs", "1", "--evalue", "1e-10", "--threads", str(cfg.threads)],
                       check=True, stderr=subprocess.DEVNULL)
    ids_eve_like = set()
    with open(tsv) as fh:
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            campos = p[4].split("|")
            genoma = campos[2] if len(campos) > 2 else ""
            family = campos[5] if len(campos) > 5 else ""
            if eh_genoma_alvo(genoma) and family not in FAMILIAS_EXCLUIDAS:
                ids_eve_like.add(p[0])
    log(f"Dipteran proteins whose best hit is a target virus (set aside): {len(ids_eve_like)}")

    with open(eve_like, "w") as fe:
        escrever_fasta_filtrado(dipteros, fe, lambda h: h[1:].split()[0] in ids_eve_like)

    # 2b. Swiss-Prot without viruses + dipterans minus the EVE-like proteins
    log("Filtering Swiss-Prot (dropping viral proteins by taxonomy)...")
    virais = taxids_virais(cfg.taxdump)
    neg_bruto = os.path.join(cfg.saida, "2_negativos_bruto.faa")
    with open(neg_bruto, "w") as out:
        def sp_nao_viral(h):
            m = re.search(r"OX=(\d+)", h)
            return bool(m) and m.group(1) not in virais
        n_sp = escrever_fasta_filtrado(cfg.swissprot, out, sp_nao_viral)
        n_dip = escrever_fasta_filtrado(dipteros, out, lambda h: h[1:].split()[0] not in ids_eve_like)
    log(f"Negatives: non-viral Swiss-Prot={n_sp}, dipterans={n_dip}")

    # 2c. Redundancy reduction (it barely changes the maximum score)
    log(f"MMseqs2 linclust at {pct}% on the negatives...")
    tmp = os.path.join(cfg.saida, "tmp_mmseqs")
    prefixo = os.path.join(cfg.saida, "2_negativos")
    subprocess.run(["mmseqs", "easy-linclust", neg_bruto, prefixo, tmp, "--min-seq-id", str(cfg.id_negativos),
                    "-c", "0.9", "--threads", str(cfg.threads)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.replace(f"{prefixo}_rep_seq.fasta", neg_final)
    for sufixo in ("_all_seqs.fasta", "_cluster.tsv"):
        if os.path.exists(prefixo + sufixo):
            os.remove(prefixo + sufixo)
    n = sum(1 for l in open(neg_final) if l.startswith(">"))
    log(f"Negatives after reduction: {n}")
    return neg_final, eve_like


# ---------------------------------------------------------------------------
# 3. Conjunto positivo
# ---------------------------------------------------------------------------
def montar_positivos(cfg) -> str:
    destino = os.path.join(cfg.saida, "3_positivos_refseq_RNA_ssDNA.faa")
    if os.path.exists(destino):
        log(f"Positives already exist: {destino}")
        return destino

    def alvo(h):
        campos = h.split("|")
        family = campos[5].strip() if len(campos) > 5 else ""
        return len(campos) > 2 and eh_genoma_alvo(campos[2].strip()) and family not in FAMILIAS_EXCLUIDAS

    with open(destino, "w") as out:
        n = escrever_fasta_filtrado(cfg.refseq_viral_fasta, out, alvo)
    log(f"Positives (RefSeq viral RNA without RT + ssDNA): {n}")
    return destino


# ---------------------------------------------------------------------------
# 4. Calibration
# ---------------------------------------------------------------------------
def carregar_bloco(fasta):
    alfabeto = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(fasta, digital=True, alphabet=alfabeto) as sf:
        return sf.read_block()


def maximos_por_modelo(hmms, bloco, cpus, t_report, guardar_top=0):
    """Return {model: (max_seq, max_dom, n_hits, top_names)} for hits >= t_report bits."""
    res = {}
    for top in pyhmmer.hmmer.hmmsearch(hmms, bloco, cpus=cpus, T=t_report, domT=t_report):
        q = nome(top.query.name)
        hits = sorted(top, key=lambda h: h.score, reverse=True)
        if not hits:
            res[q] = (0.0, 0.0, 0, "")
            continue
        max_dom = max((d.score for h in hits for d in h.domains), default=0.0)
        tops = ";".join(f"{nome(h.name)}:{h.score:.1f}" for h in hits[:guardar_top])
        res[q] = (hits[0].score, max_dom, len(hits), tops)
    return res


def calibrar(cfg, sel, negativos, eve_like, positivos):
    tabela = os.path.join(cfg.saida, "4_model_calibration.tsv")
    hmm_saida = os.path.join(os.path.dirname(cfg.saida), "combined_viral_calibrado.hmm")
    parcial = tabela + ".parcial"

    manter = set(sel.loc[sel["selection"].str.startswith("keep"), "HMM_ID"])
    info_sel = sel.set_index("HMM_ID")
    feitos = set()
    if os.path.exists(parcial):
        feitos = set(pd.read_csv(parcial, sep="\t")["HMM_ID"])
        log(f"Resuming: {len(feitos)} models already calibrated")

    log("Loading the sequence sets...")
    b_neg, b_eve, b_pos = carregar_bloco(negativos), carregar_bloco(eve_like), carregar_bloco(positivos)
    log(f"Negatives={len(b_neg)} | dipterans with a viral hit={len(b_eve)} | positives={len(b_pos)}")

    colunas = ["HMM_ID", "source", "family", "selection", "max_neg_seq", "max_neg_dom", "n_neg",
               "top_negatives", "max_dipterans_viral_hit", "TC_seq", "TC_dom", "max_pos", "n_pos_above_TC",
               "decision"]
    novo = not os.path.exists(parcial)
    fh_tab = open(parcial, "a", newline="")
    escritor = csv.DictWriter(fh_tab, fieldnames=colunas, delimiter="\t")
    if novo:
        escritor.writeheader()

    total_alvo = len(manter) if not cfg.max_modelos else min(cfg.max_modelos, len(manter))
    lote, processados, start = [], len(feitos), time.time()

    nonlocal_novos = [0]

    def processar(lote):
        nonlocal processados
        neg = maximos_por_modelo(lote, b_neg, cfg.threads, cfg.t_report, guardar_top=3)
        eve = maximos_por_modelo(lote, b_eve, cfg.threads, cfg.t_report)
        pos_res = {}
        for top in pyhmmer.hmmer.hmmsearch(lote, b_pos, cpus=cfg.threads, T=cfg.t_report, domT=cfg.t_report):
            pos_res[nome(top.query.name)] = [h.score for h in top]
        for hmm in lote:
            q = nome(hmm.name)
            mx_seq, mx_dom, n_neg, tops = neg.get(q, (0.0, 0.0, 0, ""))
            tc_seq = round(max(mx_seq + cfg.margem, cfg.piso), 1)
            tc_dom = round(max(mx_dom + cfg.margem, cfg.piso), 1)
            escores_pos = pos_res.get(q, [])
            n_pos = sum(1 for s in escores_pos if s >= tc_seq)
            decision = "keep" if n_pos >= cfg.min_pos else "drop:no_discriminating_power"
            escritor.writerow({
                "HMM_ID": q, "source": info_sel.at[q, "source"], "family": info_sel.at[q, "family"],
                "selection": info_sel.at[q, "selection"], "max_neg_seq": round(mx_seq, 1),
                "max_neg_dom": round(mx_dom, 1), "n_neg": n_neg, "top_negatives": tops,
                "max_dipterans_viral_hit": round(eve.get(q, (0.0,))[0], 1), "TC_seq": tc_seq, "TC_dom": tc_dom,
                "max_pos": round(max(escores_pos, default=0.0), 1), "n_pos_above_TC": n_pos,
                "decision": decision,
            })
        fh_tab.flush()
        nonlocal_novos[0] += len(lote)
        processados += len(lote)
        dt = time.time() - start
        feitos_sessao = processados - len(feitos)
        eta = dt / max(1, feitos_sessao) * (total_alvo - processados)
        log(f"Calibrados {processados}/{total_alvo} | decorrido {dt/60:.1f} min | restante ~{eta/60:.0f} min")

    log(f"Calibrando {total_alvo} modelos (lotes de {cfg.lote})...")
    with pyhmmer.plan7.HMMFile(cfg.hmm) as hf:
        for hmm in hf:
            q = nome(hmm.name)
            if q not in manter or q in feitos:
                continue
            if cfg.max_modelos and processados + len(lote) >= cfg.max_modelos:
                break
            if cfg.chunk and nonlocal_novos[0] + len(lote) >= cfg.chunk:
                break
            lote.append(hmm)
            if len(lote) >= cfg.lote:
                processar(lote)
                lote = []
    if lote:
        processar(lote)
    fh_tab.close()

    if cfg.max_modelos:
        log(f"Test mode finished. Partial table: {parcial}")
        return
    if processados < len(manter):
        log(f"Chunk finished ({processados}/{len(manter)}). Run again to continue.")
        return

    # Grava o .hmm calibrado com TC em cada modelo mantido
    res = pd.read_csv(parcial, sep="\t").drop_duplicates("HMM_ID", keep="last")
    res.to_csv(tabela, sep="\t", index=False)
    cortes = res[res["decision"] == "keep"].set_index("HMM_ID")[["TC_seq", "TC_dom"]].to_dict("index")
    log(f"Gravando {len(cortes)} modelos calibrados em {hmm_saida} ...")
    with pyhmmer.plan7.HMMFile(cfg.hmm) as hf, open(hmm_saida, "wb") as out:
        for hmm in hf:
            q = nome(hmm.name)
            if q in cortes:
                c = cortes[q]
                hmm.cutoffs.trusted = (float(c["TC_seq"]), float(c["TC_dom"]))
                hmm.write(out)
    log("Resumo final:\n" + res["decision"].value_counts().to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for chave, valor in PADROES.items():
        ap.add_argument(f"--{chave.replace('_', '-')}", dest=chave, default=valor)
    ap.add_argument("--threads", type=int, default=5)
    ap.add_argument("--lote", type=int, default=25, help="Modelos por lote de hmmsearch")
    ap.add_argument("--t-report", type=float, default=10.0, help="Minimum score (bits) recorded in the searches")
    ap.add_argument("--margem", type=float, default=1.0, help="Bits added to the highest negative score to set the TC")
    ap.add_argument("--piso", type=float, default=25.0, help="Minimum TC (bits) for models with no negative hit")
    ap.add_argument("--min-pos", type=int, default=1, help="Minimum number of viral proteins above the TC for the model to be kept")
    ap.add_argument("--max-modelos", type=int, default=0, help="Calibrate only N models (test run)")
    ap.add_argument("--chunk", type=int, default=0,
                    help="Calibrate at most N pending models and exit (works around the memory growth "
                         "of the process; use a shell loop to resume until done)")
    ap.add_argument("--id-negativos", type=float, default=0.7,
                    help="Clustering identity for the negatives; lower is faster (default: 0.7)")
    cfg = ap.parse_args()
    os.makedirs(cfg.saida, exist_ok=True)

    sel = selecionar_modelos(cfg)
    negativos, eve_like = montar_negativos(cfg)
    positivos = montar_positivos(cfg)
    calibrar(cfg, sel, negativos, eve_like, positivos)


if __name__ == "__main__":
    main()
