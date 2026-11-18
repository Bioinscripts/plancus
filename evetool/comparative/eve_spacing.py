#!/usr/bin/env python3
"""
How far apart are the EVEs in the genome, and is there a pattern?

For every approved EVE locus (the overlapping ORFs of one insertion together), the distance to the
nearest other approved locus of the same genome, species by species: the gap between the two loci,
edge to edge (0 when they touch or overlap), on the original chromosome or contig. A run made on
extracted segments is placed back on the whole genome first, so two EVEs of two segments of one
chromosome are neighbours. A locus alone on its contig has no neighbour: it is counted apart, as
it is in the random genomes.

Random placement (the null model): the same number of loci, with the same lengths, dropped at random
on the assembly - a contig drawn in proportion to its length, a uniform position on it - `reps`
times, and measured the same way. Per species:
  * median gap observed and at random, and their ratio (< 1: closer together than chance);
  * share of the loci with another EVE within 10 kb, observed and at random;
  * p_clustered = share of the random genomes whose median gap is <= the observed one (small: the
    EVEs sit closer together than chance), p_spread = share >= (small: more evenly spread);
  * Kolmogorov-Smirnov D between the observed gaps and the pooled random gaps, with its asymptotic p;
  * densities of log10(gap + 1) on one grid (Gaussian kernel, Silverman bandwidth), for the
    ridgeline of the page.
The null ignores what the genome is made of (repeats, gaps, heterochromatin): a departure from it
says the EVEs are not dropped uniformly, not why. A species whose sequence is not a whole genome
(a segment run without the genome in synteny_genomes.tsv: under 5 Mb, or more than 1% of it EVEs)
is reported as unavailable instead of being compared with a meaningless random placement.

With --merge N (the page's "merge EVEs within N bp"), loci closer than N bp are first joined into one
insertion - the ORFs of one insertion broken by stop codons - so the distances are between insertions,
and the random genomes get the joined insertions.

    python -m evetool.comparative.eve_spacing -o <results> [--reps 200] [--merge 0] [--tsv]
"""
import argparse
import json
import math
import os
from collections import OrderedDict

import numpy as np

GRADE = np.round(np.arange(0.0, 9.0001, 0.02), 3)      # log10(gap + 1), 1 bp to 1 Gb
PERTO = 10_000
MIN_GENOMA = 5_000_000          # below this the "genome" of a species is its extracted segments, not an assembly
_cache = OrderedDict()


def fmt_bp(v: float) -> str:
    return f"{v / 1e6:.1f} Mb" if v >= 1e6 else f"{v / 1e3:.1f} kb" if v >= 1e3 else f"{int(v)} bp"


def _tamanhos(genoma: str) -> dict:
    fai = genoma + ".fai"
    if os.path.exists(fai):
        with open(fai, encoding="utf-8") as fh:
            return {c[0]: int(c[1]) for c in (l.split("\t") for l in fh) if len(c) > 1}
    import pysam
    with pysam.FastaFile(genoma) as fa:
        return dict(zip(fa.references, fa.lengths))


def _vizinhos(contigs, ini, fim):
    """Nearest-neighbour gap of every interval among the intervals of its own contig (nan: alone).
    Sorted by start, the closest interval after one is the next start, the closest before it is
    the largest end among the earlier ones; overlapping or touching intervals give 0."""
    contigs, ini, fim = np.asarray(contigs), np.asarray(ini, dtype=float), np.asarray(fim, dtype=float)
    gap, viz = np.full(len(ini), np.nan), np.full(len(ini), -1)
    ordem = np.lexsort((ini, contigs))
    c, a, b = contigs[ordem], ini[ordem], fim[ordem]
    corte = np.flatnonzero(c[1:] != c[:-1]) + 1
    for bloco in np.split(np.arange(len(ordem)), corte):
        if len(bloco) < 2:
            continue
        aa, bb = a[bloco], b[bloco]
        fim_ant = np.maximum.accumulate(bb)                     # largest end up to each one
        quem_ant = np.zeros(len(bloco), dtype=int)
        for k in range(1, len(bloco)):
            quem_ant[k] = k if bb[k] >= bb[quem_ant[k - 1]] else quem_ant[k - 1]
        antes = np.full(len(bloco), np.inf)
        antes[1:] = np.maximum(0, aa[1:] - fim_ant[:-1])
        depois = np.full(len(bloco), np.inf)
        depois[:-1] = np.maximum(0, aa[1:] - bb[:-1])
        g = np.minimum(antes, depois)
        gap[ordem[bloco]] = g
        v = np.where(antes <= depois, np.r_[-1, quem_ant[:-1]], np.r_[np.arange(1, len(bloco)), -1])
        viz[ordem[bloco]] = [ordem[bloco][x] if x >= 0 else -1 for x in v]
    return gap, viz


def _kde(valores, bw):
    if len(valores) == 0:
        return np.zeros(len(GRADE))
    z = (GRADE[:, None] - np.asarray(valores)[None, :]) / bw
    return np.exp(-0.5 * z * z).sum(axis=1) / (len(valores) * bw * math.sqrt(2 * math.pi))


def _silverman(v):
    v = np.asarray(v)
    if len(v) < 2:
        return 0.25
    dp = v.std(ddof=1) or 0.5
    q75, q25 = np.percentile(v, [75, 25])
    iqr = (q75 - q25) / 1.34 or dp
    # capped: gaps of 0 bp and of megabases in one species would otherwise flatten its ridge
    return max(0.08, min(0.35, 0.9 * min(dp, iqr) * len(v) ** -0.2))


def _juntar(pos, distancia):
    """Loci on one contig within `distancia` bp (edge to edge) become one insertion: [(contig, start, end, [indices])]."""
    ordem = sorted(range(len(pos)), key=lambda k: (pos[k][0], pos[k][1]))
    saida = []
    for k in ordem:
        c, a, b = pos[k]
        if saida and saida[-1][0] == c and a - saida[-1][2] <= distancia:
            saida[-1][2] = max(saida[-1][2], b)
            saida[-1][3].append(k)
        else:
            saida.append([c, a, b, [k]])
    return [tuple(u) for u in saida]


def _ks(x, y):
    """Two-sample Kolmogorov-Smirnov D and its asymptotic p."""
    x, y = np.sort(x), np.sort(y)
    if not len(x) or not len(y):
        return None, None
    tudo = np.concatenate([x, y])
    d = float(np.max(np.abs(np.searchsorted(x, tudo, side="right") / len(x) - np.searchsorted(y, tudo, side="right") / len(y))))
    en = math.sqrt(len(x) * len(y) / (len(x) + len(y)))
    lam = (en + 0.12 + 0.11 / en) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    return round(d, 4), round(min(1.0, max(0.0, p)), 6)


def _nulo(contigs_tam: dict, comprimentos, reps: int, rng):
    nomes = np.array(list(contigs_tam))
    tam = np.array([contigs_tam[c] for c in nomes], dtype=float)
    prob = tam / tam.sum()
    comp = np.asarray(comprimentos, dtype=float)
    gaps, medianas, com_viz, perto = [], [], [], []
    for _ in range(reps):
        k = rng.choice(len(nomes), size=len(comp), p=prob)
        ini = np.floor(rng.random(len(comp)) * np.maximum(1, tam[k] - comp + 1))
        g, _v = _vizinhos(nomes[k], ini, ini + comp)
        ok = g[~np.isnan(g)]
        gaps.append(ok)
        com_viz.append(len(ok) / len(comp))
        perto.append(float((ok <= PERTO).sum()) / len(comp))
        medianas.append(float(np.median(ok)) if len(ok) else np.nan)
    return np.concatenate(gaps) if gaps else np.array([]), np.array(medianas), float(np.mean(com_viz)), float(np.mean(perto))


def distancias(outdir: str, reps: int = 200, semente: int = 1, juntar: int = 0) -> dict:
    """Nearest-EVE gaps per species, observed and under random placement (see the module docstring).
    `juntar` > 0 first merges the loci within that many bp into one insertion (e.g. the ORFs of one
    insertion broken by stop codons), and the random genomes get the merged insertions."""
    from evetool.curation.curation import estado_dados, loci_aprovados, genomas_para_sintenia
    from evetool.annotation.rnaseq_tracks import _origem
    from evetool.comparative.flank_synteny import mapear_para_genoma
    juntar = max(0, int(juntar or 0))
    chave = (outdir, estado_dados(outdir), int(reps), int(semente), juntar)
    if chave in _cache:
        return _cache[chave]
    rng = np.random.default_rng(semente)
    genomas = genomas_para_sintenia(outdir)
    por_sp = OrderedDict()
    for l in loci_aprovados(outdir):
        por_sp.setdefault(l["species"], []).append(l)
    especies = []
    for sp, loci in por_sp.items():
        item = {"species": sp, "n_loci": len(loci), "loci": []}
        g = genomas.get(sp) or {}
        try:
            origens = _origem(outdir, sp) if g.get("segments") else {}
            pos = []
            for l in loci:
                c, a, b, _inv = mapear_para_genoma(l["contig"], l["start"], l["end"], origens)
                pos.append((c, a, b))
            unidades = _juntar(pos, juntar) if juntar else [(c, a, b, [k]) for k, (c, a, b) in enumerate(pos)]
            nome = lambda u: loci[u[3][0]]["orf"] + (f" (+{len(u[3]) - 1})" if len(u[3]) > 1 else "")
            gap, viz = _vizinhos([u[0] for u in unidades], [u[1] for u in unidades], [u[2] for u in unidades])
            for k, u in enumerate(unidades):
                l = loci[u[3][0]]
                item["loci"].append({"orf": nome(u), "query": l["query"], "locus": l["locus"], "family": l["family"],
                                     "members": [loci[x]["orf"] for x in u[3]],
                                     "genome_contig": u[0], "genome_start": int(u[1]), "genome_end": int(u[2]),
                                     "gap": None if np.isnan(gap[k]) else int(gap[k]),
                                     "nearest": nome(unidades[viz[k]]) if viz[k] >= 0 else None})
            item["n_units"] = len(unidades)
            pos = [(u[0], u[1], u[2]) for u in unidades]
            obs = gap[~np.isnan(gap)]
            item.update(with_neighbour=int(len(obs)), contigs_with_eves=len({p[0] for p in pos}))
            tam = _tamanhos(g["genome"]) if g.get("genome") else {}
            if not tam:
                raise RuntimeError("genome not found: no contig lengths for the random placement")
            comp = np.sort(np.fromiter(tam.values(), dtype=float))[::-1]
            total = float(comp.sum())
            item.update(genome_size=int(total), genome_contigs=len(tam),
                        genome_n50=int(comp[np.searchsorted(np.cumsum(comp), total / 2)]))
            # a species of a segment run without its whole genome: what is there is the segments themselves
            soma = sum(p[2] - p[1] for p in pos)
            if total < MIN_GENOMA or soma / total > 0.01:
                item["unavailable"] = (f"the sequence of this species is not a whole genome ({fmt_bp(total)} in {len(tam):,} "
                                       f"sequence(s), {100 * soma / total:.1f}% of it EVEs): add the genome to synteny_genomes.tsv")
                especies.append(item)
                continue
            nul, medianas, com_viz_nul, perto_nul = _nulo(tam, [p[2] - p[1] for p in pos], reps, rng)
            lo = np.log10(obs + 1)
            bw = _silverman(lo)
            med = float(np.median(obs)) if len(obs) else None
            validas = medianas[~np.isnan(medianas)]
            d, p = _ks(obs, nul) if len(obs) >= 3 else (None, None)
            item.update(
                reps=reps, bandwidth=round(bw, 3),
                median_gap=med, mean_log10=round(float(lo.mean()), 3) if len(lo) else None,
                median_gap_random=float(np.median(nul)) if len(nul) else None,
                with_neighbour_random=round(com_viz_nul, 4),
                near_share=round(float((obs <= PERTO).sum()) / len(pos), 4), near_share_random=round(perto_nul, 4),
                p_clustered=round(float((validas <= med).mean()), 4) if med is not None and len(validas) else None,
                p_spread=round(float((validas >= med).mean()), 4) if med is not None and len(validas) else None,
                ks_d=d, ks_p=p,
                density=[round(float(x), 5) for x in _kde(lo, bw)],
                density_random=[round(float(x), 5) for x in _kde(np.log10(nul + 1), bw)] if len(nul) else [])
            if item["median_gap_random"]:
                item["median_ratio"] = round(med / item["median_gap_random"], 4) if med is not None else None
        except Exception as erro:
            item["error"] = str(erro)
        especies.append(item)
    resposta = {"grid": [float(x) for x in GRADE], "near": PERTO, "reps": reps, "species": especies,
                "stamp": chave[1], "merge": juntar}
    _cache[chave] = resposta
    while len(_cache) > 4:
        _cache.popitem(last=False)
    return resposta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--merge", type=int, default=0, help="merge the loci within this many bp into one insertion first")
    ap.add_argument("--tsv", action="store_true", help="write eve_spacing_loci.tsv and eve_spacing_species.tsv in 13_curadoria")
    a = ap.parse_args()
    r = distancias(a.outdir, a.reps, a.seed, a.merge)
    campos = ["species", "n_loci", "n_units", "with_neighbour", "median_gap", "median_gap_random", "median_ratio", "near_share",
              "near_share_random", "p_clustered", "p_spread", "ks_d", "ks_p", "genome_size", "genome_contigs", "genome_n50",
              "unavailable", "error"]
    for s in r["species"]:
        print("\t".join(str(s.get(k, "")) for k in campos))
    if a.tsv:
        pasta = os.path.join(a.outdir, "13_curadoria")
        os.makedirs(pasta, exist_ok=True)
        with open(os.path.join(pasta, "eve_spacing_species.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\t".join(campos) + "\n")
            for s in r["species"]:
                fh.write("\t".join("" if s.get(k) is None else str(s.get(k)) for k in campos) + "\n")
        with open(os.path.join(pasta, "eve_spacing_loci.tsv"), "w", encoding="utf-8") as fh:
            fh.write("species\torf\tfamily\tgenome_contig\tgenome_start\tgenome_end\tnearest_eve\tgap_bp\n")
            for s in r["species"]:
                for l in s.get("loci", []):
                    fh.write("\t".join(str(x) for x in (s["species"], l["orf"], l["family"], l["genome_contig"], l["genome_start"] + 1,
                                                        l["genome_end"], l["nearest"] or "", "" if l["gap"] is None else l["gap"])) + "\n")
        print(json.dumps({"written": pasta}))


if __name__ == "__main__":
    main()
