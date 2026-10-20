#!/usr/bin/env python3
"""
Flank synteny of EVE loci across species.

Two loci of different species in one cross-species cluster may be the same ancestral
insertion (orthologous) or two independent integrations of related viruses. Only the first
keeps the same host neighbours. For every locus the flanks are taken from its genome, the
host genes in them are labelled by homology (DIAMOND blastx against a reference proteome,
isoforms collapsed with MMseqs2), and each inter-species pair compares the nearest labelled
genes of its left and right flanks, in both orientations (assemblies may be reversed).

Pair status:
  syntenic            shared host genes on both sides
  one_side            shared genes on one side only (e.g. a short contig or a rearrangement)
  not_syntenic        both loci have labelled genes on both sides, none shared
  insufficient_flank  a locus has no labelled gene on one side: nothing can be said

Pre-extracted segments (contigs named 'contig_start-end', the original location kept in the
FASTA header as 'contig:start-end(strand)') are mapped back onto the whole genome given in
<results>/synteny_genomes.tsv ('species<TAB>genome.fasta'), so their real flanks are used.

Usage (normally called by multispecies_clustering.py / the curation page):
  python -m evetool.comparative.flank_synteny --loci loci.tsv --clusters clusters.tsv --proteome ref.faa -o out
"""
import os
import re
import csv
import json
import shutil
import hashlib
import argparse
import subprocess
from collections import defaultdict


def log(msg):
    print(f"[synteny] {msg}", flush=True)


ORDEM_STATUS = {"syntenic": 4, "one_side": 3, "not_syntenic": 2, "insufficient_flank": 1, "no_genome": 0}


# ---------------------------------------------------------------------------
# Reference proteome: DIAMOND database and isoform collapse, cached by content
# ---------------------------------------------------------------------------
def _cache_dir() -> str:
    d = os.path.join(os.path.expanduser("~"), ".cache", "eve_tool", "synteny")
    os.makedirs(d, exist_ok=True)
    return d


def preparar_proteoma(proteoma: str, threads: int = 4) -> dict:
    """{'dmnd': path, 'grupo': {protein: representative}, 'nome': {protein: description}}."""
    if not proteoma or not os.path.exists(proteoma):
        raise RuntimeError(f"Reference proteome not found: {proteoma}")
    st = os.stat(proteoma)
    chave = hashlib.md5(f"{os.path.abspath(proteoma)}|{st.st_size}|{st.st_mtime}".encode()).hexdigest()[:12]
    base = os.path.join(_cache_dir(), chave)
    dmnd, mapa = base + ".dmnd", base + "_groups.json"
    nomes = {}
    with open(proteoma, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                p = linha[1:].strip().split(None, 1)
                nomes[p[0]] = p[1] if len(p) > 1 else p[0]
    if not os.path.exists(dmnd):
        log(f"formatting the reference proteome ({len(nomes):,} proteins)...")
        subprocess.run(["diamond", "makedb", "--in", proteoma, "-d", base, "--threads", str(threads), "--quiet"],
                       check=True)
    if os.path.exists(mapa):
        with open(mapa, encoding="utf-8") as fh:
            grupo = json.load(fh)
    else:
        grupo = {p: p for p in nomes}
        if shutil.which("mmseqs"):
            log("collapsing isoforms of the reference proteome (MMseqs2, 90% identity, 80% coverage)...")
            tmp = base + "_mmseqs"
            r = subprocess.run(["mmseqs", "easy-cluster", proteoma, tmp, tmp + "_tmp", "--min-seq-id", "0.9",
                                "-c", "0.8", "--cov-mode", "1", "--threads", str(threads), "-v", "1"],
                               capture_output=True, text=True)
            tsv = tmp + "_cluster.tsv"
            if r.returncode == 0 and os.path.exists(tsv):
                with open(tsv) as fh:
                    for linha in fh:
                        rep, membro = linha.rstrip("\n").split("\t")[:2]
                        grupo[membro] = rep
            for sufixo in ("_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta"):
                if os.path.exists(tmp + sufixo):
                    os.remove(tmp + sufixo)
            shutil.rmtree(tmp + "_tmp", ignore_errors=True)
        with open(mapa, "w", encoding="utf-8") as fh:
            json.dump(grupo, fh)
    return {"dmnd": dmnd, "grupo": grupo, "nome": nomes}


def descrever(nome: str) -> str:
    """'uncharacterized protein LOC1 [Phlebotomus papatasi]' -> 'uncharacterized protein LOC1'."""
    return re.sub(r"\s*\[[^\]]+\]\s*$", "", nome or "").strip()


# ---------------------------------------------------------------------------
# Genomes and coordinates
# ---------------------------------------------------------------------------
def _fai(genoma: str) -> dict:
    fai = genoma + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", genoma], check=True, capture_output=True)
    tam = {}
    with open(fai) as fh:
        for linha in fh:
            p = linha.split("\t")
            tam[p[0]] = int(p[1])
    return tam


def origens_de_segmentos(fasta: str) -> dict:
    """{segment contig: (original contig, start1, end1, strand)} from 'contig:start-end(strand)' headers."""
    origem = {}
    if not fasta or not os.path.exists(fasta):
        return origem
    rx = re.compile(r"(\S+?):(\d+)-(\d+)(?:\(([^)]*)\))?$")
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if not linha.startswith(">"):
                continue
            partes = linha[1:].split()
            nome = partes[0]
            for campo in partes[1:][::-1]:
                m = rx.match(campo)
                if m:
                    fita = "-" if (m.group(4) or "+").startswith("-") else "+"
                    origem[nome] = (m.group(1), int(m.group(2)), int(m.group(3)), fita)
                    break
            else:
                m = re.match(r"^(.+?)_(\d+)-(\d+)$", nome)
                if m:
                    origem[nome] = (m.group(1), int(m.group(2)), int(m.group(3)), "+")
    return origem


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _faidx(fasta: str, regioes: list) -> list:
    seqs = []
    for i in range(0, len(regioes), 500):
        p = subprocess.run(["samtools", "faidx", fasta, *regioes[i:i + 500]], capture_output=True, text=True)
        atual = None
        for linha in p.stdout.splitlines():
            if linha.startswith(">"):
                atual = []
                seqs.append(atual)
            elif atual is not None:
                atual.append(linha.strip())
    return ["".join(x).upper() for x in seqs]


def calibrar_origens(origens: dict, fasta_segmentos: str, genoma: str) -> dict:
    """
    Fix each segment's offset and orientation by comparing its sequence with the genome: the
    header coordinates may be 0- or 1-based and its strand tag is not always the orientation of
    the stored sequence. Returns {segment: (contig, start0, end, strand)}; unmatched segments are left out.
    """
    if not origens or not genoma or not os.path.exists(genoma):
        return {}
    tam = _fai(genoma)
    _fai(fasta_segmentos)
    nomes = [n for n, (c, a, b, _f) in origens.items() if c in tam]
    segs = _faidx(fasta_segmentos, nomes)
    regs = [f"{origens[n][0]}:{max(1, origens[n][1])}-{min(tam[origens[n][0]], origens[n][2] + 1)}" for n in nomes]
    gens = _faidx(genoma, regs)
    saida = {}
    for n, sq, g in zip(nomes, segs, gens):
        c, a, b, _f = origens[n]
        base = max(1, a)                     # 1-based position of g[0]
        for ini1 in (a + 1, a):              # 0-based header first (BED-like), then 1-based
            k = ini1 - base
            alvo = g[k:k + len(sq)] if k >= 0 else ""
            if len(alvo) != len(sq) or not sq:
                continue
            if alvo == sq:
                saida[n] = (c, ini1 - 1, ini1 - 1 + len(sq), "+")
                break
            if alvo.translate(_RC)[::-1] == sq:
                saida[n] = (c, ini1 - 1, ini1 - 1 + len(sq), "-")
                break
    if len(saida) < len(origens):
        log(f"{len(origens) - len(saida)} of {len(origens)} segment(s) could not be placed on {os.path.basename(genoma)}")
    return saida


def mapear_para_genoma(contig, ini0, fim, origem: dict):
    """0-based [ini0, fim) on a segment -> the same span on the original contig (and the strand flip)."""
    o = origem.get(contig)
    if not o:
        return contig, ini0, fim, False
    c, a0, b0, fita = o                      # calibrated: 0-based start, end, orientation
    if fita == "-":
        return c, b0 - fim, b0 - ini0, True
    return c, a0 + ini0, a0 + fim, False


# ---------------------------------------------------------------------------
# Flank labelling
# ---------------------------------------------------------------------------
def rotular_flancos(loci: dict, genomas: dict, origens: dict, ref: dict, flanco: int, workdir: str,
                    threads: int = 4, max_genes: int = 6) -> dict:
    """
    loci: {name: (species, contig, ini0, fim)} on the analysed sequences.
    Returns {name: {"L": [labels nearest first], "R": [...], "where": "contig:a-b", "status": ...}}.
    Left/right follow the genome's + strand.
    """
    os.makedirs(workdir, exist_ok=True)
    cache_p = os.path.join(workdir, "flank_labels_cache.json")
    try:
        with open(cache_p, encoding="utf-8") as fh:
            cache = json.load(fh)
    except Exception:
        cache = {}
    resultado, pedidos, tamanhos = {}, defaultdict(list), {}
    for nome, (sp, contig, ini0, fim) in loci.items():
        genoma = genomas.get(sp)
        if not genoma or not os.path.exists(genoma):
            resultado[nome] = {"L": [], "R": [], "where": "", "status": "no_genome"}
            continue
        if genoma not in tamanhos:
            tamanhos[genoma] = _fai(genoma)
        c, a, b, _inv = mapear_para_genoma(contig, ini0, fim, origens.get(sp, {}))
        if c not in tamanhos[genoma]:
            resultado[nome] = {"L": [], "R": [], "where": "", "status": "no_genome"}
            continue
        L = tamanhos[genoma][c]
        chave = f"v2|{os.path.abspath(genoma)}|{c}|{a}|{b}|{flanco}|{ref['dmnd']}"
        resultado[nome] = {"where": f"{c}:{a + 1}-{b}", "contig_len": L, "chave": chave}
        if chave in cache:
            resultado[nome].update(cache[chave])
            continue
        esq = (max(0, a - flanco), a)
        dir_ = (b, min(L, b + flanco))
        pedidos[genoma].append((nome, c, esq, dir_))

    if pedidos:
        fasta = os.path.join(workdir, "flanks.fasta")
        n = 0
        with open(fasta, "w") as fo:
            for genoma, lista in pedidos.items():
                regs, rotulos = [], []
                for nome, c, esq, dir_ in lista:
                    for lado, (x, y) in (("L", esq), ("R", dir_)):
                        if y - x >= 300:
                            regs.append(f"{c}:{x + 1}-{y}")
                            rotulos.append((nome, lado, x, y))
                for i in range(0, len(regs), 500):
                    lote, lrot = regs[i:i + 500], rotulos[i:i + 500]
                    p = subprocess.run(["samtools", "faidx", genoma, *lote], capture_output=True, text=True)
                    seqs, atual = [], None
                    for linha in p.stdout.splitlines():
                        if linha.startswith(">"):
                            atual = []
                            seqs.append(atual)
                        elif atual is not None:
                            atual.append(linha.strip())
                    for (nome, lado, x, y), partes in zip(lrot, seqs):
                        fo.write(f">{n}|{lado}|{x}|{y}\n{''.join(partes)}\n")
                        resultado[nome].setdefault("_ids", {})[str(n)] = lado
                        n += 1
        hits = defaultdict(list)
        if n:
            log(f"labelling {n} flanks of {sum(len(v) for v in pedidos.values())} loci "
                f"(±{flanco:,} nt) against the reference proteome...")
            tsv = os.path.join(workdir, "flanks_vs_proteome.tsv")
            subprocess.run(["diamond", "blastx", "-q", fasta, "-d", ref["dmnd"], "-o", tsv, "--evalue", "1e-10",
                            "-k", "1", "-F", "15", "--range-culling", "--sensitive", "--threads", str(threads), "--quiet",
                            "--outfmt", "6", "qseqid", "sseqid", "qstart", "qend", "bitscore"], check=True)
            with open(tsv) as fh:
                for linha in fh:
                    q, s, qs, qe, bits = linha.rstrip("\n").split("\t")
                    hits[q].append((min(int(qs), int(qe)), max(int(qs), int(qe)), s, float(bits)))
        idx_nome = {i: (nome, lado) for nome, r in resultado.items() for i, lado in r.get("_ids", {}).items()}
        brutos = defaultdict(lambda: {"L": [], "R": []})
        for q, lista in hits.items():
            i, lado, x, y = q.split("|")
            nome, _ = idx_nome[i]
            comp = int(y) - int(x)
            for qs, qe, s, bits in lista:
                # distance to the EVE: from the flank end next to it (near and far edge of the hit)
                perto, longe = (comp - qe, comp - qs) if lado == "L" else (qs, qe)
                brutos[nome][lado].append((perto, longe, ref["grupo"].get(s, s)))
        for nome, r in resultado.items():
            if "_ids" not in r:
                continue
            for lado in ("L", "R"):
                vistos, ordem, pos = set(), [], []
                for perto, longe, g in sorted(brutos[nome][lado]):
                    if g not in vistos:
                        vistos.add(g)
                        ordem.append(g)
                        pos.append([g, perto, longe])
                r[lado] = ordem[:max_genes]
                r[lado + "p"] = pos[:max_genes]          # [label, near edge, far edge] in nt from the EVE
            r.pop("_ids", None)
            cache[r["chave"]] = {"L": r["L"], "R": r["R"], "Lp": r["Lp"], "Rp": r["Rp"]}
        with open(cache_p, "w", encoding="utf-8") as fh:
            json.dump(cache, fh)
    for r in resultado.values():
        for k in ("L", "R", "Lp", "Rp"):
            r.setdefault(k, [])
        r.pop("chave", None)
        if "status" not in r:
            r["status"] = "ok" if (r["L"] or r["R"]) else "no_gene"
    return resultado


# ---------------------------------------------------------------------------
# Pairwise comparison
# ---------------------------------------------------------------------------
def comparar_par(a: dict, b: dict) -> dict:
    if a.get("status") == "no_genome" or b.get("status") == "no_genome":
        return {"status": "no_genome", "shared": [], "orientation": ""}
    melhor = None
    for orient, (bl, br) in (("same", (b["L"], b["R"])), ("reversed", (b["R"], b["L"]))):
        esq = [g for g in a["L"] if g in bl]
        dir_ = [g for g in a["R"] if g in br]
        nota = (bool(esq) + bool(dir_), len(esq) + len(dir_))
        if melhor is None or nota > melhor[0]:
            melhor = (nota, orient, esq, dir_)
    (lados, _n), orient, esq, dir_ = melhor
    if lados == 2:
        st = "syntenic"
    elif lados == 1:
        st = "one_side"
    elif a["L"] and a["R"] and b["L"] and b["R"]:
        st = "not_syntenic"
    else:
        st = "insufficient_flank"
    return {"status": st, "left": esq, "right": dir_, "shared": esq + dir_,
            "orientation": orient if lados else ""}


# Neighbours that say nothing about orthology: mobile elements, and genes found next to many
# distinct sites of one species (repeats, EVE-rich arrays copied around the genome)
RX_MOVEL = re.compile(r"transpos|reverse transcriptase|\bRNA-directed DNA polymerase|integrase|\bgag\b|\bpol\b|"
                      r"retrovirus-related|jockey|piggyBac|mariner|tigger|pogo|helitron|endonuclease-reverse", re.I)


def vizinhos_repetitivos(flancos: dict, especie_de, nomes: dict, flanco: int, max_sitios: int = 3) -> dict:
    """{label: reason} for labels to ignore: mobile-element descriptions, or more than max_sitios
    distinct sites (loci further apart than two flanks) in one species."""
    sitios = defaultdict(lambda: defaultdict(list))
    for nome, f in flancos.items():
        w = f.get("where") or ""
        if ":" not in w:
            continue
        c, ab = w.rsplit(":", 1)
        a, b = (int(x) for x in ab.split("-"))
        for g in set(f["L"] + f["R"]):
            sitios[g][especie_de(nome)].append((c, (a + b) // 2))
    fora = {}
    for g, por_sp in sitios.items():
        if RX_MOVEL.search(nomes.get(g, "")):
            fora[g] = "mobile element"
            continue
        for sp, lista in por_sp.items():
            n, ult = 0, None
            for c, x in sorted(lista):
                if not ult or c != ult[0] or x - ult[1] > 2 * flanco:
                    n += 1
                ult = (c, x)
            if n > max_sitios:
                fora[g] = f"{n} sites in {sp}"
                break
    return fora


def avaliar(clusters: dict, loci: dict, genomas: dict, origens: dict, proteoma: str, flanco: int,
            workdir: str, threads: int = 4, especie_de=lambda n: n.split("|", 1)[0], max_sitios: int = 3) -> dict:
    """
    clusters: {cluster: [unit names]}; only clusters with two or more species are evaluated.
    Returns {"clusters": {c: summary}, "pairs": {c: [(a, b, result)]}, "flanks": {name: labels},
             "names": {label: description}}.
    """
    alvo = {c: ms for c, ms in clusters.items() if len({especie_de(m) for m in ms}) >= 2}
    vazio = {"synteny_status": "not_shared", "synteny_orthologs": "", "synteny_n_shared": 0,
             "synteny_pairs": 0, "synteny_syntenic_pairs": 0}
    if not alvo:
        return {"clusters": {c: dict(vazio) for c in clusters}, "pairs": {}, "flanks": {}, "names": {}}
    ref = preparar_proteoma(proteoma, threads)
    usados = {m for ms in alvo.values() for m in ms}
    flancos = rotular_flancos({m: loci[m] for m in usados if m in loci}, genomas, origens, ref, flanco,
                              workdir, threads)
    for m in usados:
        flancos.setdefault(m, {"L": [], "R": [], "where": "", "status": "no_genome"})
    nomes_ref = {g: descrever(ref["nome"].get(g, g)) for r in flancos.values() for g in r["L"] + r["R"]}
    excluidos = vizinhos_repetitivos(flancos, especie_de, nomes_ref, flanco, max_sitios)
    if excluidos:
        log(f"{len(excluidos)} neighbour gene(s) ignored as repetitive or mobile (more than {max_sitios} sites in a species)")
    uteis = {m: dict(f, L=[g for g in f["L"] if g not in excluidos], R=[g for g in f["R"] if g not in excluidos])
             for m, f in flancos.items()}
    resumo, pares = {}, {}
    for c, ms in clusters.items():
        if c not in alvo:
            resumo[c] = dict(vazio)
            continue
        lista = []
        for i, a in enumerate(ms):
            for b in ms[i + 1:]:
                if especie_de(a) != especie_de(b):
                    lista.append((a, b, comparar_par(uteis[a], uteis[b])))
        pares[c] = lista
        melhor = max((r["status"] for _a, _b, r in lista), key=lambda s: ORDEM_STATUS[s], default="insufficient_flank")
        genes = defaultdict(int)
        for _a, _b, r in lista:
            for g in set(r["shared"]):
                genes[g] += 1
        n_sint = sum(1 for _a, _b, r in lista if r["status"] == "syntenic")
        resumo[c] = {"synteny_status": melhor,
                     "synteny_orthologs": ";".join(descrever(ref["nome"].get(g, g)) or g
                                                   for g, _n in sorted(genes.items(), key=lambda x: -x[1])[:6]),
                     "synteny_n_shared": len(genes), "synteny_pairs": len(lista),
                     "synteny_syntenic_pairs": n_sint}
    nomes = nomes_ref
    log(f"{len(alvo)} shared cluster(s): " + ", ".join(
        f"{sum(1 for v in resumo.values() if v['synteny_status'] == s)} {s}" for s in ORDEM_STATUS))
    return {"clusters": resumo, "pairs": pares, "flanks": flancos, "names": nomes, "excluded": excluidos}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True, help="TSV: name, species, contig, start0, end, genome")
    ap.add_argument("--clusters", required=True, help="TSV: name, cluster")
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--flank", type=int, default=20000)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-t", "--threads", type=int, default=4)
    a = ap.parse_args()
    loci, genomas, clusters = {}, {}, defaultdict(list)
    with open(a.loci) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            loci[r[0]] = (r[1], r[2], int(r[3]), int(r[4]))
            genomas[r[1]] = r[5]
    with open(a.clusters) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            clusters[r[1]].append(r[0])
    res = avaliar(clusters, loci, genomas, {}, a.proteome, a.flank, a.outdir, a.threads)
    with open(os.path.join(a.outdir, "synteny.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=list)


if __name__ == "__main__":
    main()
