#!/usr/bin/env python3
"""
Repeats and low-complexity around the approved EVEs.

For every EVE locus (overlapping ORFs of one insertion count once):

  1. Low complexity of the EVE itself: dustmasker on its nucleotide sequence and segmasker on
     the peptide of its representative ORF. A candidate that is mostly low complexity is
     probably a compositional artefact, not a viral fragment.
  2. Repeat content from the genome's soft-masking (NCBI assemblies are soft-masked: repeats in
     lower case): masked share of the EVE and of each flank, against the genome-wide share.
  3. Transposable elements next to the EVE: ORFs of the flanks (six frames, >= 100 aa) searched
     with the transposon domains of Pfam-A (pyhmmer, gathering cut-offs) and grouped into
     LTR retrotransposon, retrotransposon RT, non-LTR (endonuclease), DNA transposon and
     Helitron; the nearest hit of each class and its distance to the EVE.

Context of each locus: 'low complexity' (artefact flag), 'inside TE' (a TE domain overlaps the
EVE, or the EVE is mostly soft-masked), 'TE nearby' (a TE domain within the flank), 'repeat-rich
flank' (flanks masked above twice the genome share) or 'clean'.

Runs on pre-extracted segments use the whole genomes of <results>/synteny_genomes.tsv for the
soft-masking and the flanks, like the synteny step.

Usage:
  python -m evetool.annotation.repeats -o <results> [--flank 20000] [-t 8]
Outputs in <results>/13_curadoria/repeticoes/: repeats_per_locus.tsv, repeats_te_hits.tsv, repeats.json
"""
import os
import re
import csv
import json
import glob
import argparse
import subprocess
from collections import defaultdict, Counter

SUBPASTA = "13_curadoria"

# Pfam transposon domains -> class. Viral RdRps are not here: only reverse transcriptases,
# integrases, gag and transposases, so a viral EVE never counts as its own neighbour.
CLASSES_TE = [
    ("LTR retrotransposon", r"^(rve|rve_2|rve_3|Integrase_.*|gag-asp_proteas|gag_pre-integrs|Peptidase_A17|Gypsy|"
                            r"Retrotrans_gag|RT_RNaseH|RT_RNaseH_2|RVT_2|DUF4817)$"),
    ("retrotransposon RT",  r"^(RVT_1|RVT_3|RVT_connect|RVT_thumb|RVT_N)$"),
    ("non-LTR (endonuclease)", r"^(Exo_endo_phos_2)$"),
    ("Helitron",            r"^(Helitron_like_N)$"),
    ("DNA transposon",      r"^(DDE_.*|Dimer_Tnp_hAT|hAT-like_RNase-H|MULE|Transposase_.*|Tnp_.*|HTH_Tnp_.*|"
                            r"TnpB_IS66|Transp_Tc5_C)$"),
]


def log(msg):
    print(f"[repeats] {msg}", flush=True)


def classe_te(modelo: str) -> str:
    for nome, rx in CLASSES_TE:
        if re.match(rx, modelo):
            return nome
    return ""


# ---------------------------------------------------------------------------
# Pfam subset of transposon domains (cached)
# ---------------------------------------------------------------------------
def hmms_te(pfam: str):
    import pyhmmer
    cache = os.path.join(os.path.expanduser("~"), ".cache", "eve_tool", "pfam_te.hmm")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    if not os.path.exists(cache) or os.path.getmtime(cache) < os.path.getmtime(pfam):
        log("extracting the transposon domains of Pfam-A (once)...")
        n = 0
        with pyhmmer.plan7.HMMFile(pfam) as fh, open(cache + ".tmp", "wb") as fo:
            for hmm in fh:
                nome = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                if classe_te(nome):
                    hmm.write(fo)
                    n += 1
        os.replace(cache + ".tmp", cache)
        log(f"{n} transposon models kept")
    with pyhmmer.plan7.HMMFile(cache) as fh:
        return list(fh)


def _pfam_padrao() -> str:
    try:
        from evetool import config
        cand = [getattr(config, "PFAM_HMM", None)]
    except Exception:
        cand = []
    cand += [os.environ.get("EVE_PFAM_HMM")]
    for c in cand:
        if c and os.path.exists(c):
            return c
    raise RuntimeError("Pfam-A.hmm not found: run `./eve-tool databases --pfam`, "
                       "set pfam_hmm in eve_config.ini or pass --pfam.")


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------
def _faidx(fasta: str, regioes: list, upper: bool = False) -> list:
    seqs = []
    for i in range(0, len(regioes), 400):
        p = subprocess.run(["samtools", "faidx", fasta, *regioes[i:i + 400]], capture_output=True, text=True)
        atual = None
        for linha in p.stdout.splitlines():
            if linha.startswith(">"):
                atual = []
                seqs.append(atual)
            elif atual is not None:
                atual.append(linha.strip())
    seqs = ["".join(x) for x in seqs]
    return [s.upper() for s in seqs] if upper else seqs


def fracao_genoma_mascarada(genoma: str) -> float:
    """Lower-case share of a genome (soft-masked repeats), cached next to the user's cache."""
    st = os.stat(genoma)
    cache = os.path.join(os.path.expanduser("~"), ".cache", "eve_tool", "softmask.json")
    try:
        with open(cache) as fh:
            d = json.load(fh)
    except Exception:
        d = {}
    chave = f"{os.path.abspath(genoma)}|{st.st_size}|{st.st_mtime}"
    if chave in d:
        return d[chave]
    baixa = total = 0
    minus = bytes(range(ord("a"), ord("z") + 1))
    with open(genoma, "rb") as fh:
        for linha in fh:
            if linha.startswith(b">"):
                continue
            linha = linha.rstrip()
            total += len(linha)
            baixa += len(linha) - len(linha.translate(None, minus))
    d[chave] = round(100.0 * baixa / total, 2) if total else 0.0
    with open(cache, "w") as fh:
        json.dump(d, fh)
    return d[chave]


def _mascarado(intervalos_txt: str) -> dict:
    """dustmasker/segmasker '-outfmt interval' -> {id: masked positions}."""
    out, atual = defaultdict(int), None
    for linha in intervalos_txt.splitlines():
        if linha.startswith(">"):
            atual = linha[1:].split()[0]
            out[atual] += 0
        elif " - " in linha and atual is not None:
            a, b = (int(x) for x in linha.split(" - "))
            out[atual] += b - a + 1
    return out


def baixa_complexidade(seqs: dict, tipo: str, workdir: str) -> dict:
    """{id: % masked} with dustmasker (nt) or segmasker (aa)."""
    if not seqs:
        return {}
    fa = os.path.join(workdir, f"lowcomplexity_{tipo}.fasta")
    with open(fa, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")
    exe = ["dustmasker", "-in", fa, "-outfmt", "interval"] if tipo == "nt" else \
          ["segmasker", "-in", fa, "-outfmt", "interval"]
    p = subprocess.run(exe, capture_output=True, text=True)
    if p.returncode != 0:
        log(f"[warning] {exe[0]} failed: {p.stderr.strip()[:200]}")
        return {}
    m = _mascarado(p.stdout)
    return {k: round(100.0 * m.get(k, 0) / max(1, len(s.replace("*", ""))), 1) for k, s in seqs.items()}


_CODONS = None


def _orfs_seis_quadros(seq: str, min_aa: int = 100):
    """[(aa, start0, end, strand)] of stop-to-stop ORFs on both strands of seq (upper case)."""
    from Bio.Seq import Seq
    L, saida = len(seq), []
    for fita, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for f in range(3):
            trecho = s[f:]
            trecho = trecho[:len(trecho) - len(trecho) % 3]
            aa = str(Seq(trecho).translate())
            pos = 0
            for pedaco in aa.split("*"):
                if len(pedaco) >= min_aa and pedaco.count("X") < len(pedaco) * 0.2:
                    a_nt, b_nt = f + 3 * pos, f + 3 * (pos + len(pedaco))
                    if fita == "+":
                        saida.append((pedaco, a_nt, b_nt, "+"))
                    else:
                        saida.append((pedaco, L - b_nt, L - a_nt, "-"))
                pos += len(pedaco) + 1
    return saida


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def analisar(outdir: str, flanco: int = 20000, threads: int = 4, pfam: str = None, queries=None) -> dict:
    import sys
    from evetool.curation import curation as cur
    from evetool.comparative.flank_synteny import origens_de_segmentos, calibrar_origens, mapear_para_genoma, _fai
    import pyhmmer

    dados = cur.listar_candidatas(outdir)
    if dados.get("error"):
        raise RuntimeError(dados["error"])
    alvo = set(queries) if queries else {i["query"] for i in dados["items"] if not i["removed_flag"]}
    workdir = os.path.join(outdir, SUBPASTA, "repeticoes")
    os.makedirs(workdir, exist_ok=True)

    # 1. loci: approved ORFs of one insertion together
    loci = {}
    for i in dados["items"]:
        if i["query"] not in alvo or not i.get("contig") or not i.get("end"):
            continue
        k = f"{i['species']}|{i.get('locus') or i['query']}"
        l = loci.setdefault(k, {"id": k, "species": i["species"], "contig": i["contig"], "start": i["start"],
                                "end": i["end"], "queries": [], "rep": i})
        l["start"], l["end"] = min(l["start"], i["start"]), max(l["end"], i["end"])
        l["queries"].append(i["query"])
        if i.get("locus_representative") and not l["rep"].get("locus_representative"):
            l["rep"] = i
    if not loci:
        raise RuntimeError("No approved EVE with coordinates.")
    log(f"{len(loci)} EVE loci from {len(alvo)} approved ORFs")

    # 2. genomes: the analysed one for the EVE itself, the whole one for masking and flanks
    gen_sint = cur.genomas_para_sintenia(outdir)
    origens, base_mask = {}, {}
    for sp, g in gen_sint.items():
        if g["segments"] and g["own"] and g["genome"]:
            origens[sp] = calibrar_origens(origens_de_segmentos(g["own"]), g["own"], g["genome"])
    for sp in {l["species"] for l in loci.values()}:
        g = (gen_sint.get(sp) or {}).get("genome")
        if g:
            base_mask[sp] = fracao_genoma_mascarada(g)
            if base_mask[sp] < 1.0:          # not soft-masked: masking says nothing here
                log(f"[warning] {sp}: genome not soft-masked ({base_mask[sp]}% lower case); masking columns left empty")
                base_mask[sp] = None

    # 3. EVE sequences (nt from the analysed genome; peptide of the representative ORF)
    por_sp = defaultdict(list)
    for l in loci.values():
        por_sp[l["species"]].append(l)
    nt = {}
    for sp, lista in por_sp.items():
        genoma = cur._genoma_da_especie(outdir, sp)
        seqs = _faidx(genoma, [f"{l['contig']}:{l['start'] + 1}-{l['end']}" for l in lista], upper=True)
        for l, s in zip(lista, seqs):
            nt[l["id"]] = s
    pep = {}
    for sp, pasta in cur.pastas_especies(outdir).items():
        cam = os.path.join(pasta, "04_candidatas_eves", "eve_candidatas_unificadas_pep.fasta")
        if not os.path.exists(cam):
            continue
        reps = {l["rep"]["orf"]: l["id"] for l in loci.values() if l["species"] == sp}
        atual = None
        with open(cam, encoding="utf-8", errors="replace") as fh:
            for linha in fh:
                if linha.startswith(">"):
                    atual = reps.get(linha[1:].split()[0])
                    if atual:
                        pep[atual] = ""
                elif atual:
                    pep[atual] += linha.strip()
    dust = baixa_complexidade(nt, "nt", workdir)
    seg = baixa_complexidade({k: v.rstrip("*") for k, v in pep.items() if v}, "aa", workdir)

    # 4. soft-masking and flanks on the whole genome
    regs_por_gen = defaultdict(list)
    for l in loci.values():
        g = (gen_sint.get(l["species"]) or {}).get("genome")
        l["genome"] = g
        if not g:
            continue
        c, a, b, _inv = mapear_para_genoma(l["contig"], l["start"], l["end"], origens.get(l["species"], {}))
        tam = _fai(g).get(c)
        if tam is None:
            l["genome"] = None
            continue
        l["g_contig"], l["g_start"], l["g_end"], l["contig_len"] = c, a, b, tam
        l["esq"] = (max(0, a - flanco), a)
        l["dir"] = (b, min(tam, b + flanco))
        regs_por_gen[g].append(l)
    flancos_fa = {}
    for g, lista in regs_por_gen.items():
        regs = []
        for l in lista:
            regs.append(f"{l['g_contig']}:{l['g_start'] + 1}-{l['g_end']}")
            for x, y in (l["esq"], l["dir"]):
                regs.append(f"{l['g_contig']}:{x + 1}-{max(x + 1, y)}")
        seqs = _faidx(g, regs)
        for k, l in enumerate(lista):
            eve_s, esq_s, dir_s = seqs[3 * k:3 * k + 3]
            com_mask = base_mask.get(l["species"]) is not None
            frac = lambda s: round(100.0 * sum(ch.islower() for ch in s) / len(s), 1) if s and com_mask else None
            l["flank_nt"] = (l["esq"][1] - l["esq"][0]) + (l["dir"][1] - l["dir"][0])
            l["mask_eve"] = frac(eve_s)
            l["mask_left"] = frac(esq_s) if l["esq"][1] > l["esq"][0] else None
            l["mask_right"] = frac(dir_s) if l["dir"][1] > l["dir"][0] else None
            flancos_fa[l["id"]] = {"L": esq_s.upper(), "R": dir_s.upper(), "E": eve_s.upper()}

    # 5. TE domains in the ORFs of the flanks and of the EVE span
    hmms = hmms_te(pfam or _pfam_padrao())
    alfa = pyhmmer.easel.Alphabet.amino()
    seqs, meta = [], {}
    for lid, fl in flancos_fa.items():
        l = loci[lid]
        for lado, s in fl.items():
            for n, (aa, a, b, fita) in enumerate(_orfs_seis_quadros(s)):
                nome = f"{len(meta)}"
                meta[nome] = (lid, lado, a, b, fita, len(s))
                seqs.append(pyhmmer.easel.TextSequence(name=nome.encode(), sequence=aa).digitize(alfa))
    log(f"{len(seqs)} ORFs (>= 100 aa) in the EVEs and flanks of {len(flancos_fa)} loci; "
        f"searching {len(hmms)} transposon domains...")
    te_hits = defaultdict(list)
    if seqs:
        for top in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=threads, bit_cutoffs="gathering"):
            modelo = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
            for hit in top:
                if not hit.included:
                    continue
                nome = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                lid, lado, a, b, fita, comp = meta[nome]
                for d in hit.domains.included:
                    ia, ib = d.env_from - 1, d.env_to          # aa on the ORF
                    if fita == "+":
                        na, nb = a + 3 * ia, a + 3 * ib
                    else:
                        na, nb = b - 3 * ib, b - 3 * ia
                    # distance to the EVE (0 = inside it)
                    dist = 0 if lado == "E" else (comp - nb if lado == "L" else na)
                    te_hits[lid].append({"side": lado, "model": modelo, "class": classe_te(modelo),
                                         "distance": max(0, dist), "score": round(d.score, 1),
                                         "evalue": float(f"{d.i_evalue:.2e}")})

    # 6. per-locus verdict
    linhas = []
    for lid, l in loci.items():
        hits = sorted(te_hits.get(lid, []), key=lambda h: h["distance"])
        dentro = [h for h in hits if h["side"] == "E"]
        perto = [h for h in hits if h["side"] != "E"]
        base = base_mask.get(l["species"])
        lc_nt, lc_aa = dust.get(lid), seg.get(lid)
        mflank = [x for x in (l.get("mask_left"), l.get("mask_right")) if x is not None]
        if (lc_nt or 0) >= 50 or (lc_aa or 0) >= 50:
            ctx = "low complexity"
        elif dentro or (l.get("mask_eve") or 0) >= 50:
            ctx = "inside TE"
        elif perto:
            ctx = "TE nearby"
        elif base and mflank and sum(mflank) / len(mflank) >= max(2 * base, base + 20):
            ctx = "repeat-rich flank"
        elif not l.get("genome"):
            ctx = "no genome"
        elif (l.get("flank_nt") or 0) < min(2000, flanco):
            ctx = "short flanks"             # nothing around to judge (short contig)
        else:
            ctx = "clean"
        classes = Counter(h["class"] for h in perto)
        prox = {}
        for h in perto:
            prox.setdefault(h["class"], h)
        linhas.append({
            "locus": lid, "species": l["species"], "contig": l["contig"], "start": l["start"] + 1, "end": l["end"],
            "genome_location": f"{l['g_contig']}:{l['g_start'] + 1}-{l['g_end']}" if l.get("g_contig") else "",
            "orfs": l["queries"], "representative": l["rep"]["query"],
            "family": l["rep"].get("family", ""), "protein_class": l["rep"].get("protein_class", ""),
            "lowcomplexity_nt_pct": lc_nt, "lowcomplexity_aa_pct": lc_aa,
            "softmask_eve_pct": l.get("mask_eve"), "softmask_left_pct": l.get("mask_left"),
            "softmask_right_pct": l.get("mask_right"), "softmask_genome_pct": base,
            "te_inside": sorted({f"{h['model']} ({h['class']})" for h in dentro}),
            "te_classes_nearby": dict(classes),
            "te_nearest": [{"class": c, "model": h["model"], "distance": h["distance"], "side": h["side"]}
                           for c, h in sorted(prox.items(), key=lambda kv: kv[1]["distance"])],
            "te_nearest_distance": perto[0]["distance"] if perto else None,
            "flank_nt": l.get("flank_nt"),
            "context": ctx,
        })
    linhas.sort(key=lambda r: (r["species"], r["contig"], r["start"]))

    # 7. outputs
    cols = ["locus", "species", "contig", "start", "end", "genome_location", "family", "protein_class",
            "lowcomplexity_nt_pct", "lowcomplexity_aa_pct", "softmask_eve_pct", "softmask_left_pct",
            "softmask_right_pct", "softmask_genome_pct", "flank_nt", "context", "te_inside", "te_classes_nearby",
            "te_nearest_distance", "te_nearest", "orfs"]
    with open(os.path.join(workdir, "repeats_per_locus.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in linhas:
            w.writerow([";".join(r[c]) if c in ("te_inside", "orfs") else
                        ";".join(f"{k}:{v}" for k, v in r[c].items()) if c == "te_classes_nearby" else
                        ";".join(f"{x['class']}:{x['model']}@{x['distance']}{x['side']}" for x in r[c]) if c == "te_nearest" else
                        ("" if r[c] is None else r[c]) for c in cols])
    with open(os.path.join(workdir, "repeats_te_hits.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["locus", "side", "model", "class", "distance_to_eve", "score", "evalue"])
        for lid, hs in te_hits.items():
            for h in sorted(hs, key=lambda h: h["distance"]):
                w.writerow([lid, {"L": "left", "R": "right", "E": "inside"}[h["side"]], h["model"], h["class"],
                            h["distance"], h["score"], h["evalue"]])
    resumo = Counter(r["context"] for r in linhas)
    res = {"flank": flanco, "loci": linhas, "summary": dict(resumo),
           "genome_softmask": base_mask, "pfam_models": len(hmms)}
    with open(os.path.join(workdir, "repeats.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False)
    log(", ".join(f"{v} {k}" for k, v in resumo.most_common()))
    return {"loci": len(linhas), "summary": dict(resumo), "directory": workdir}


def carregar(outdir: str):
    cam = os.path.join(outdir, SUBPASTA, "repeticoes", "repeats.json")
    if os.path.exists(cam):
        try:
            with open(cam, encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            pass
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--flank", type=int, default=20000)
    ap.add_argument("--pfam", default=None)
    ap.add_argument("-t", "--threads", type=int, default=4)
    a = ap.parse_args()
    print(json.dumps(analisar(a.outdir, a.flank, a.threads, a.pfam), indent=1))


if __name__ == "__main__":
    main()
