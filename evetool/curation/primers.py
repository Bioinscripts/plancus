#!/usr/bin/env python3
"""
Primer suggestions for an EVE (Primer3, through primer3-py).

Three kinds of pair are designed on the genome sequence of the EVE and its flanks:

  internal        both primers inside the inserted sequence: amplifies the EVE itself
                  (detection, RT-PCR / qPCR of its transcript);
  junction_left   forward primer in the host flank before the EVE, reverse primer inside it;
  junction_right  forward primer inside the EVE, reverse primer in the host flank after it.

A junction product exists only if the viral sequence is integrated at that point of the host
genome, so it tells an EVE from an exogenous virus or a contaminant. "Left" and "right" are on
the + strand of the contig, whatever the strand of the EVE.

Soft-masked (lower-case) bases of the genome are kept: Primer3 refuses a primer whose 3' end
sits on one, which keeps the external primers out of annotated repeats when the assembly is
soft-masked. Coordinates in the answer are 1-based, inclusive, on the contig.

  python -m evetool.curation.primers -o <results> --query "Species|ORF_id"
  python -m evetool.curation.primers -o <results> --all          # every approved locus -> TSV
"""
import argparse
import csv
import json
import os

SUBPASTA = "13_curadoria"
TIPOS = ("internal", "junction_left", "junction_right")
ROTULOS = {"internal": "internal (inside the EVE)",
           "junction_left": "junction, left flank → EVE",
           "junction_right": "junction, EVE → right flank"}

PADRAO = {
    "flank": 400,                 # host sequence offered on each side for the external primer
    "n": 3,                       # pairs returned per kind
    "internal_size": (90, 300),   # product size range of the internal pair
    "junction_size": (150, 700),  # product size range of the junction pairs
    "tm": (57.0, 60.0, 63.0),     # min, optimum, max
    "size": (18, 20, 25),
    "gc": (35.0, 50.0, 65.0),
    "max_tm_diff": 3.0,
    "avoid_masked": True,
}


def _opcoes(opcoes: dict = None) -> dict:
    o = dict(PADRAO)
    for k, v in (opcoes or {}).items():
        if k in o and v not in (None, ""):
            o[k] = v
    for k in ("internal_size", "junction_size"):
        a, b = o[k]
        o[k] = (max(40, int(a)), max(int(a) + 10, int(b)))
    o["flank"] = max(60, min(int(o["flank"]), 20000))
    o["n"] = max(1, min(int(o["n"]), 10))
    return o


def _globais(o: dict, faixa: tuple) -> dict:
    tmin, topt, tmax = o["tm"]
    smin, sopt, smax = o["size"]
    gmin, gopt, gmax = o["gc"]
    return {
        "PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_NUM_RETURN": o["n"],
        "PRIMER_MIN_SIZE": smin, "PRIMER_OPT_SIZE": sopt, "PRIMER_MAX_SIZE": smax,
        "PRIMER_MIN_TM": tmin, "PRIMER_OPT_TM": topt, "PRIMER_MAX_TM": tmax,
        "PRIMER_PAIR_MAX_DIFF_TM": o["max_tm_diff"],
        "PRIMER_MIN_GC": gmin, "PRIMER_OPT_GC_PERCENT": gopt, "PRIMER_MAX_GC": gmax,
        "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_NS_ACCEPTED": 0, "PRIMER_GC_CLAMP": 1,
        "PRIMER_MAX_END_GC": 3,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_MAX_SELF_ANY_TH": 45.0, "PRIMER_MAX_SELF_END_TH": 35.0, "PRIMER_MAX_HAIRPIN_TH": 24.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0, "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
        "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_SALT_DIVALENT": 1.5, "PRIMER_DNTP_CONC": 0.6,
        "PRIMER_DNA_CONC": 50.0, "PRIMER_TM_FORMULA": 1, "PRIMER_SALT_CORRECTIONS": 1,
        "PRIMER_LOWERCASE_MASKING": 1 if o["avoid_masked"] else 0,
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": 3, "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": 3,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(faixa)],
    }


def _oligo(res: dict, lado: str, k: int, t0: int, seq: str) -> dict:
    pos, comp = res[f"PRIMER_{lado}_{k}"]
    ini = pos if lado == "LEFT" else pos - comp + 1          # 0-based on the template
    return {
        "sequence": res[f"PRIMER_{lado}_{k}_SEQUENCE"],
        "start": t0 + ini + 1, "end": t0 + ini + comp,       # 1-based inclusive on the contig
        "strand": "+" if lado == "LEFT" else "-",
        "length": comp,
        "tm": round(res[f"PRIMER_{lado}_{k}_TM"], 1),
        "gc": round(res[f"PRIMER_{lado}_{k}_GC_PERCENT"], 1),
        "self_any": round(res.get(f"PRIMER_{lado}_{k}_SELF_ANY_TH", 0.0), 1),
        "self_end": round(res.get(f"PRIMER_{lado}_{k}_SELF_END_TH", 0.0), 1),
        "hairpin": round(res.get(f"PRIMER_{lado}_{k}_HAIRPIN_TH", 0.0), 1),
        "masked": any(c.islower() for c in seq[ini:ini + comp]),
    }


def _por_que_nao(res: dict) -> str:
    partes = [f"{lado.lower()}: {res[ch]}" for lado, ch in (("LEFT", "PRIMER_LEFT_EXPLAIN"), ("RIGHT", "PRIMER_RIGHT_EXPLAIN"),
                                                         ("PAIR", "PRIMER_PAIR_EXPLAIN")) if res.get(ch)]
    return "; ".join(partes)


def desenhar(seq: str, eve_ini: int, eve_fim: int, t0: int = 0, opcoes: dict = None, tipos=TIPOS) -> dict:
    """
    Primer pairs for the EVE at seq[eve_ini:eve_fim] (0-based, half-open on the template, whose
    first base is position t0 + 1 of the contig). Returns {"pairs": [...], "notes": [...]}.
    """
    import primer3
    o = _opcoes(opcoes)
    n, comp_eve = len(seq), eve_fim - eve_ini
    pares, notas = [], []
    esq, dir_ = eve_ini, n - eve_fim                       # flank available on each side
    smin = o["size"][0]

    for tipo in tipos:
        if tipo == "internal":
            if comp_eve < 2 * smin + 10:
                notas.append(f"internal: the EVE is too short ({comp_eve} bp)")
                continue
            a, b = o["internal_size"]
            faixa = (min(a, comp_eve), min(b, comp_eve))
            args = {"SEQUENCE_INCLUDED_REGION": [eve_ini, comp_eve]}
        else:
            flanco = esq if tipo == "junction_left" else dir_
            if flanco < smin + 5:
                notas.append(f"{tipo}: no host flank on this side (contig end at {flanco} bp)")
                continue
            if comp_eve < smin:
                notas.append(f"{tipo}: the EVE is too short ({comp_eve} bp)")
                continue
            faixa = o["junction_size"]
            regioes = [0, esq, eve_ini, comp_eve] if tipo == "junction_left" else [eve_ini, comp_eve, eve_fim, dir_]
            args = {"SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [regioes]}
        try:
            res = primer3.bindings.design_primers({"SEQUENCE_ID": tipo, "SEQUENCE_TEMPLATE": seq, **args},
                                                  _globais(o, faixa))
        except Exception as erro:
            notas.append(f"{tipo}: Primer3 refused the input ({erro})")
            continue
        achados = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
        if not achados:
            notas.append(f"{tipo}: no pair within the constraints ({_por_que_nao(res)})")
        for k in range(achados):
            f, r = _oligo(res, "LEFT", k, t0, seq), _oligo(res, "RIGHT", k, t0, seq)
            pares.append({
                "type": tipo, "label": ROTULOS[tipo], "rank": k + 1,
                "forward": f, "reverse": r,
                "product_size": int(res[f"PRIMER_PAIR_{k}_PRODUCT_SIZE"]),
                "product_start": f["start"], "product_end": r["end"],
                "penalty": round(res[f"PRIMER_PAIR_{k}_PENALTY"], 2),
                "pair_compl_any": round(res.get(f"PRIMER_PAIR_{k}_COMPL_ANY_TH", 0.0), 1),
                "pair_compl_end": round(res.get(f"PRIMER_PAIR_{k}_COMPL_END_TH", 0.0), 1),
                # share of the product that is EVE: a junction product carries host and virus
                "eve_bases": max(0, min(r["end"], t0 + eve_fim) - max(f["start"] - 1, t0 + eve_ini)),
            })
    return {"pairs": pares, "notes": notas, "options": {k: (list(v) if isinstance(v, tuple) else v) for k, v in o.items()}}


def _molde(genoma: str, contig: str, ini: int, fim: int, flanco: int):
    """Genome sequence of [ini - flanco, fim + flanco) clipped to the contig, case kept."""
    import pysam
    with pysam.FastaFile(genoma) as fa:
        if contig not in fa.references:
            raise RuntimeError(f"contig {contig} is not in {os.path.basename(genoma)}")
        tam = fa.get_reference_length(contig)
        t0, t1 = max(0, ini - flanco), min(tam, fim + flanco)
        return fa.fetch(contig, t0, t1), t0, tam


def para_regiao(outdir: str, species: str, contig: str, ini: int, fim: int, nome: str = "",
                opcoes: dict = None, tipos=TIPOS) -> dict:
    """Primers for contig:[ini, fim) (0-based, half-open, as the coordinates of the page)."""
    from evetool.curation.curation import _genoma_da_especie, genomas_para_sintenia
    ini, fim = int(ini), int(fim)
    if fim <= ini:
        raise RuntimeError("empty region")
    if fim - ini > 200_000:
        raise RuntimeError(f"region too long for primer design ({fim - ini:,} bp)")
    o = _opcoes(opcoes)
    genoma, colocado = _genoma_da_especie(outdir, species), None
    g = genomas_para_sintenia(outdir).get(species) or {}
    if g.get("segments") and g.get("genome"):
        # run made on extracted segments: the flanks are on the whole genome, so design there
        from evetool.annotation.rnaseq_tracks import _origem
        from evetool.comparative.flank_synteny import mapear_para_genoma
        c2, a2, b2, invertido = mapear_para_genoma(contig, ini, fim, _origem(outdir, species))
        if c2 != contig:
            colocado = {"contig": contig, "start": ini + 1, "end": fim, "flipped": invertido}
            genoma, contig, ini, fim = g["genome"], c2, a2, b2
    seq, t0, tam = _molde(genoma, contig, ini, fim, o["flank"])
    r = desenhar(seq, ini - t0, fim - t0, t0, o, tipos)
    mascarado = sum(1 for c in seq if c.islower())
    if o["avoid_masked"] and mascarado > 0.9 * len(seq):
        # a fully lower-case assembly is not soft-masking: design again ignoring the case
        r = desenhar(seq.upper(), ini - t0, fim - t0, t0, {**o, "avoid_masked": False}, tipos)
        r["notes"].append("the genome is all lower-case here, so soft-masking was ignored")
    if colocado:
        r["placed_from"] = colocado
        r["notes"].append(f"designed on the whole genome: segment {colocado['contig']}:{colocado['start']}-{colocado['end']} "
                          f"is {contig}:{ini + 1}-{fim}{' (reversed)' if colocado['flipped'] else ''}")
    r.update({"name": nome, "species": species, "contig": contig, "start": ini + 1, "end": fim,
              "length": fim - ini, "contig_length": tam,
              "template_start": t0 + 1, "template_end": t0 + len(seq),
              "soft_masked_pct": round(100 * mascarado / max(1, len(seq)), 1)})
    return r


COLUNAS = ["species", "eve", "contig", "eve_start", "eve_end", "eve_length", "family", "protein_class",
           "pair_type", "rank", "forward", "forward_start", "forward_end", "forward_tm", "forward_gc",
           "reverse", "reverse_start", "reverse_end", "reverse_tm", "reverse_gc",
           "product_size", "product_start", "product_end", "eve_bases_in_product", "penalty", "notes"]


def para_aprovadas(outdir: str, queries=None, opcoes: dict = None, tipos=TIPOS) -> dict:
    """
    Primers for every approved locus (its representative ORF; `queries` restricts the set).
    Writes 13_curadoria/primers/eve_primers.tsv, one row per pair, and a row with the reason
    for each kind of pair that could not be designed.
    """
    from evetool.curation.curation import listar_candidatas
    itens = [i for i in listar_candidatas(outdir)["items"] if not i["removed_flag"] and i["contig"] and i["end"]]
    if queries is not None:
        quais = set(queries)
        itens = [i for i in itens if i["query"] in quais]
    # one design per locus: the span of all its approved reading frames
    loci = {}
    for i in itens:
        k = (i["species"], i["locus"] or i["query"])
        l = loci.setdefault(k, {"rep": i, "start": i["start"], "end": i["end"]})
        l["start"], l["end"] = min(l["start"], i["start"]), max(l["end"], i["end"])
        if i["locus_representative"] and not l["rep"]["locus_representative"]:
            l["rep"] = i
    pasta = os.path.join(outdir, SUBPASTA, "primers")
    os.makedirs(pasta, exist_ok=True)
    destino = os.path.join(pasta, "eve_primers.tsv")
    n_pares = com = 0
    with open(destino, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLUNAS)
        for l in sorted(loci.values(), key=lambda x: (x["rep"]["species"], x["rep"]["contig"], x["start"])):
            i = l["rep"]
            base = [i["species"], i["orf"], i["contig"], l["start"] + 1, l["end"], l["end"] - l["start"],
                    i["family"], i["protein_class"]]
            try:
                r = para_regiao(outdir, i["species"], i["contig"], l["start"], l["end"], i["orf"], opcoes, tipos)
            except Exception as erro:
                w.writerow(base + ["", ""] + [""] * 15 + [str(erro)])
                continue
            com += bool(r["pairs"])
            # the coordinates the primers are in (the whole genome when the run was made on segments)
            base[2:6] = [r["contig"], r["start"], r["end"], r["length"]]
            for p in r["pairs"]:
                n_pares += 1
                f, rv = p["forward"], p["reverse"]
                w.writerow(base + [p["type"], p["rank"], f["sequence"], f["start"], f["end"], f["tm"], f["gc"],
                                   rv["sequence"], rv["start"], rv["end"], rv["tm"], rv["gc"],
                                   p["product_size"], p["product_start"], p["product_end"], p["eve_bases"],
                                   p["penalty"], ""])
            for nota in r["notes"]:
                if nota.split(":")[0] in TIPOS:          # a kind of pair that could not be designed
                    w.writerow(base + [nota.split(":")[0], ""] + [""] * 15 + [nota])
    return {"file": destino, "loci": len(loci), "loci_with_primers": com, "pairs": n_pares}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--query", help='"Species|ORF" (or the ORF alone in a single-genome run)')
    ap.add_argument("--all", action="store_true", help="every approved locus, to a TSV")
    ap.add_argument("--flank", type=int, default=PADRAO["flank"])
    a = ap.parse_args()
    if a.all:
        print(json.dumps(para_aprovadas(a.outdir, opcoes={"flank": a.flank}), indent=1))
        return
    from evetool.curation.curation import listar_candidatas
    it = next((i for i in listar_candidatas(a.outdir)["items"] if i["query"] == a.query or i["orf"] == a.query), None)
    if not it:
        raise SystemExit(f"{a.query}: not among the candidates")
    r = para_regiao(a.outdir, it["species"], it["contig"], it["start"], it["end"], it["orf"], {"flank": a.flank})
    for p in r["pairs"]:
        print(f'{p["type"]:15s} #{p["rank"]} {p["forward"]["sequence"]} ({p["forward"]["tm"]}°C)  '
              f'{p["reverse"]["sequence"]} ({p["reverse"]["tm"]}°C)  {p["product_size"]} bp')
    for nota in r["notes"]:
        print("note:", nota)


if __name__ == "__main__":
    main()
