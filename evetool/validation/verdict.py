#!/usr/bin/env python3
"""
Verdict of each candidate: is it an EVE of an RNA (non-retroviral) or ssDNA virus, or
something else (host gene, transposon, virus outside the scope)?

Two things are reported, on purpose kept apart:

  detection   HOW the candidate was found: which of the three independent lines of viral
              evidence it has (alignment to a viral protein, viral HMM, viral conserved domain)
  verdict     WHAT it most likely is, from explicit rules over that evidence plus the
              evidence against (cellular hit, transposon domains, virus outside the scope)

Verdicts, strongest first:

  EVE (high confidence)       alignment + (HMM or in-scope domain), nothing against
  EVE (probable)              one alignment, or HMM + in-scope domain; or high with a caveat
  Ambiguous                   a single weak line (HMM only), or viral ~ cellular
  Possible ancient insertion  the cellular hit wins, but it is an uncharacterised protein:
                              an old EVE annotated as a host gene looks exactly like this
  No viral evidence           neither alignment, HMM nor in-scope domain
  Host gene                   a characterised cellular protein wins and there is no in-scope
                              viral domain. Labelled only: host genes are removed by hand, from
                              the list of the Host genes tab
  Transposon / Out of scope   explicitly not an RNA/ssDNA-virus insertion: the page moves
                              these to the Removed box by itself (the researcher can restore)

Mapped reads, when there are any, add to this (see _com_reads): DNA reads of the assembly that
confirm the integration and small RNA production support an EVE; assembly reads that do not
cover it, or pile up far above the flanks, count against.

Poly-A, poly-anything and short tandem repeats inside the segment, and the low-complexity and
soft-masking figures of the repeats step, also count against it (see _com_repeticoes).

Only the last line is ever removed automatically. A candidate that aligns only to host
proteins is NOT removed for that: uncharacterised host proteins may be the EVE itself.

The rules work per ORF; the automatic removal is decided per locus (all overlapping reading
frames), and a locus is kept when any of its frames has in-scope viral support.
"""

import csv
import os
import re

from evetool.core.utils import classe_genoma_viral

AUTO_PREFIXO = "auto:"
VEREDITOS = ["EVE (high confidence)", "EVE (probable)", "Ambiguous", "Possible ancient insertion",
             "No viral evidence", "Host gene", "Transposon", "Out of scope"]
# "Host gene" is a verdict but is NOT removed by itself: a host gene is not always a non-EVE, so
# the researcher removes host genes by hand from the list of the Host genes tab.
NAO_VIRAIS = ("Host gene", "Transposon", "Out of scope")
REMOVIDOS_SOZINHOS = ("Transposon", "Out of scope")
ESCOPOS_ALVO = ("RNA virus", "ssDNA virus")

# Bacteriophage families (ssDNA, ssRNA and dsRNA ones included: in a eukaryotic assembly a
# phage is bacterial contamination, not an insertion) and reverse-transcribing families.
FAGOS = {"inoviridae", "microviridae", "plectroviridae", "paulinoviridae", "finnlakeviridae",
         "fiersviridae", "leviviridae", "steitzviridae", "atkinsviridae", "duinviridae",
         "blumeviridae", "solspiviridae", "cystoviridae", "tectiviridae", "corticoviridae",
         "plasmaviridae", "myoviridae", "siphoviridae", "podoviridae", "autographiviridae",
         "herelleviridae", "ackermannviridae", "demerecviridae", "drexlerviridae",
         "straboviridae", "kyanoviridae", "peduoviridae", "casjensviridae", "zobellviridae"}
RX_FAGO = re.compile(r"\b(bacterio)?phage\b|caudovir|leviviricetes", re.I)
RETRO = {"retroviridae", "metaviridae", "pseudoviridae", "belpaoviridae", "ortervirales",
         "caulimoviridae", "hepadnaviridae", "nudnaviridae"}
RX_NAO_CARACTERIZADA = re.compile(
    r"uncharacteri[sz]ed|hypothetical|unknown function|unnamed protein|predicted protein|"
    r"putative protein\b|\bDUF\d+|\bUPF\d+|\bORF\d*\b|\bLOC\d+", re.I)

_GENOMA_DA_FAMILIA = None


def _genoma_da_familia(family: str) -> str:
    """Genome type of a viral family, from the table shipped with the viral HMMs."""
    global _GENOMA_DA_FAMILIA
    if _GENOMA_DA_FAMILIA is None:
        _GENOMA_DA_FAMILIA = {}
        try:
            from evetool import config
            caminho = config.hmm_path("RVDB", "rvdb_family_genome_mapping.tsv")
            if os.path.exists(caminho):
                with open(caminho, encoding="utf-8", errors="replace") as fh:
                    for r in csv.reader(fh, delimiter="\t"):
                        if len(r) >= 2 and r[0] != "Family":
                            _GENOMA_DA_FAMILIA[r[0].strip().lower()] = r[1].strip()
        except Exception:
            pass
    return _GENOMA_DA_FAMILIA.get((family or "").strip().lower(), "")


def escopo_viral(family: str, product: str = "", genoma_hmm: str = "") -> str:
    """
    Scope of the viral evidence: "RNA", "ssDNA" (in scope), "phage", "RT", "dsDNA" (out of
    scope) or "" when it cannot be told. The family of the best viral hit decides; the genome
    type of the best HMM model is the fallback.
    """
    f = (family or "").strip().lower()
    if f in FAGOS or RX_FAGO.search(product or ""):
        return "phage"
    if f in RETRO:
        return "RT"
    classe = classe_genoma_viral(_genoma_da_familia(f))
    if classe == "unknown":
        if f in ("riboviria", "orthornavirae") or f.endswith(("viricota", "viricetes")) and "rna" in f:
            classe = "RNA"
        else:
            classe = classe_genoma_viral((genoma_hmm or "").split(";")[0])
    return "" if classe == "unknown" else classe


def avaliar(i: dict, te_dentro=()) -> dict:
    """
    Detection class, verdict and reasons of one candidate (an item of listar_candidatas).
    te_dentro: transposon domains found inside its locus by the repeats step.
    """
    aln = i.get("viral_bitscore", 0) > 0
    hmm = i.get("hmm_n_domains", 0) > 0
    escopos_dom = {i.get("cdd_scope", ""), i.get("pfam_scope", "")}
    dom = bool(escopos_dom & set(ESCOPOS_ALVO))
    linhas = [nome for nome, tem in (("Alignment", aln), ("HMM", hmm), ("Domain", dom)) if tem]
    deteccao = " + ".join(linhas) if linhas else "No evidence"

    grupo = i.get("group", "")
    escopo = escopo_viral(i.get("family_auto") or i.get("family", ""), i.get("product", ""),
                          i.get("hmm_model_genome", "")) if (aln or hmm) else ""
    fora = escopo in ("phage", "RT", "dsDNA")
    retro_dom = "Retro (excluded)" in escopos_dom
    nao_viral_dom = "Non-viral" in escopos_dom and not dom
    outro_virus_dom = "Other virus" in escopos_dom and not dom
    bits_c = i.get("cellular_bitscore", 0)
    titulo_c = i.get("cellular_title", "") or ""
    celular_caracterizada = bits_c > 0 and titulo_c not in ("", "n/a") and not RX_NAO_CARACTERIZADA.search(titulo_c)
    celular_domina = grupo.startswith("4_")

    # an alignment where the virus clearly wins, to a virus inside the scope
    apoio_aln = aln and grupo[:1] in ("1", "2") and not fora

    motivos, contra = [], []
    if aln:
        motivos.append(f"viral alignment ({i.get('viral_bitscore', 0):.0f} bits, {i.get('family') or 'no family'})")
    if hmm:
        motivos.append(f"{i.get('hmm_n_domains')} viral HMM domain(s)")
    if dom:
        motivos.append("RNA/ssDNA-virus conserved domain")
    if bits_c > 0:
        contra.append(f"cellular hit {bits_c:.0f} bits" + ("" if celular_caracterizada else " (uncharacterised)"))
    if nao_viral_dom:
        contra.append("only non-viral conserved domains")
    if outro_virus_dom:
        contra.append("domain of a virus outside the scope")
    if i.get("review"):
        contra.append("viral target of unknown genome type")
    if retro_dom:
        contra.append("reverse-transcriptase / integrase / gag domain")
    if te_dentro:
        contra.append("transposon domain inside the locus: " + ", ".join(sorted(te_dentro)[:3]))

    # explicitly not an RNA/ssDNA-virus insertion (removed automatically)
    if fora and not dom:
        veredito = "Out of scope"
        contra.insert(0, {"phage": "best viral evidence is a bacteriophage",
                          "RT": "best viral evidence is a reverse-transcribing virus or retroelement",
                          "dsDNA": "best viral evidence is a dsDNA virus"}[escopo])
    elif (retro_dom or te_dentro) and not dom and not apoio_aln:
        veredito = "Transposon"
    elif not dom and celular_caracterizada and (celular_domina or (not aln and not hmm)
                                                 or (nao_viral_dom and not aln)):
        veredito = "Host gene"
        contra.insert(0, "characterised cellular protein wins and no in-scope viral domain")
    # the researcher decides
    elif celular_domina:
        veredito = "Possible ancient insertion"
    elif not linhas:
        veredito = "No viral evidence"
    elif aln and (hmm or dom) and grupo[:1] in ("1", "2"):
        veredito = "EVE (high confidence)" if not contra else "EVE (probable)"
    elif (aln and grupo[:1] in ("1", "2")) or (hmm and dom):
        veredito = "EVE (probable)" if not (nao_viral_dom or outro_virus_dom) else "Ambiguous"
    else:
        veredito = "Ambiguous"

    return {"detection": deteccao, "n_evidence": len(linhas), "viral_scope": escopo,
            "verdict": veredito, "verdict_for": motivos, "verdict_against": contra,
            "in_scope_support": bool(dom or apoio_aln),
            "auto_reject": veredito if veredito in REMOVIDOS_SOZINHOS else ""}


def _com_reads(x: dict, ev: dict) -> None:
    """
    Add what the mapped reads say (validation/read_evidence.py) to the verdict of one frame.
    Reads never make a non-viral sequence viral and never remove anything by themselves:
      - integration confirmed or small RNA production lifts a clean "probable" to "high";
      - both together lift an "Ambiguous" (a single weak line, e.g. HMM only) to "probable";
      - assembly reads that do not cover the locus, or pile up far above the flanks, hold a
        verdict back one step; a junction no read crosses is a caveat on "high".
    """
    dna, srna = ev.get("dna"), ev.get("srna")
    apoio = []
    limpo = not x["verdict_against"]
    if dna:
        if dna["state"] == "integrated":
            x["verdict_for"].append(dna["text"])
            apoio.append("DNA")
        elif dna["state"] != "partial":
            x["verdict_against"].append(dna["text"])
    if srna:
        x["verdict_for"].append(srna["text"])
        apoio.append("sRNA")
    x["read_support"] = apoio
    v = x["verdict"]
    if v in NAO_VIRAIS or v == "No viral evidence":
        return
    if dna and dna["state"] in ("not covered", "excess"):
        x["verdict"] = {"EVE (high confidence)": "EVE (probable)", "EVE (probable)": "Ambiguous"}.get(v, v)
    elif dna and dna["state"] == "junction open":
        x["verdict"] = "EVE (probable)" if v == "EVE (high confidence)" else v
    elif v == "EVE (probable)" and apoio and limpo:
        x["verdict"] = "EVE (high confidence)"
    elif v == "Ambiguous" and len(apoio) == 2:
        x["verdict"] = "EVE (probable)"


def _com_hospedeiro(x: dict, h: dict) -> None:
    """
    Best hit of the frame in the host proteomes (validation/host_genes.py). It is listed against
    the frame only when it could explain it away: a characterised host protein scoring at least
    as much as the best viral hit, and not the annotation of this very ORF in its own genome
    (same species, >= 95 % identity). Such a hit holds "high confidence" back to "probable";
    removing the gene stays a decision of the researcher.
    """
    if not h:
        return
    especie = (x.get("species") or "").replace("_", " ").lower()
    proprio = bool(especie) and especie == (h.get("organism") or "").lower() and (h.get("pident") or 0) >= 95
    x["host_hit"] = {"description": h.get("description", ""), "organism": h.get("organism", ""),
                     "bitscore": h["bitscore"], "pident": h.get("pident"), "self": proprio,
                     "uncharacterised": bool(h.get("uncharacterised"))}
    if proprio or h.get("uncharacterised") or h["bitscore"] < (x.get("viral_bitscore") or 0):
        return
    x["verdict_against"].append(f"host proteome hit as strong as the viral one: {h.get('description') or 'protein'} "
                                f"({h.get('organism') or 'host'}, {h['bitscore']:.0f} bits)")
    if x["verdict"] == "EVE (high confidence)":
        x["verdict"] = "EVE (probable)"


def _com_repeticoes(x: dict, simples: dict, rep: dict) -> None:
    """
    Homopolymers and simple repeats in the segment (validation/low_complexity.py), plus the
    low-complexity and soft-masking figures of the repeats step when it ran, count against it:
    they are always listed; from 5 % of the segment (or one homopolymer of 20 nt) they hold
    "high confidence" back to "probable"; a segment that is mostly simple
    sequence (>= 30 %) is held back one more step, to "Ambiguous".
    """
    linhas = []
    forte = False
    pesa = False          # enough to change the verdict, not only to be listed
    if simples and simples.get("flags"):
        linhas.append(f"simple sequence in the segment ({simples['simple_frac']:.0%}): " + "; ".join(simples["flags"]))
        forte = bool(simples.get("strong"))
        # one short run in a long ORF is listed but changes nothing: it weighs from 5 % of the
        # segment, or from a single homopolymer of 20 nt
        pesa = simples["simple_frac"] >= 0.05 or simples.get("longest_homopolymer", 0) >= 20
    if rep:
        if (rep.get("lowcomplexity_nt_pct") or 0) >= 20:
            linhas.append(f"{rep['lowcomplexity_nt_pct']:.0f}% of the segment is low complexity (DUST)")
            forte = forte or rep["lowcomplexity_nt_pct"] >= 50
            pesa = True
        if (rep.get("softmask_eve_pct") or 0) >= 50:
            linhas.append(f"{rep['softmask_eve_pct']:.0f}% of the segment is a soft-masked repeat of the genome")
            pesa = True
    x["simple_seq"] = {"flags": (simples or {}).get("flags", []), "frac": (simples or {}).get("simple_frac", 0.0),
                       "longest_homopolymer": (simples or {}).get("longest_homopolymer", 0), "strong": forte}
    x["simple_sequence"] = "; ".join(x["simple_seq"]["flags"])
    if not linhas:
        return
    x["verdict_against"].extend(linhas)
    if not pesa:
        return
    if x["verdict"] == "EVE (high confidence)":
        x["verdict"] = "EVE (probable)"
    if forte and x["verdict"] == "EVE (probable)":
        x["verdict"] = "Ambiguous"


def aplicar(itens: list, te_por_locus: dict = None, reads: dict = None, hospedeiro: dict = None,
            simples: dict = None, repeticoes: dict = None) -> list:
    """
    Fill the verdict fields of every item, then settle the automatic removal per locus:
    a frame that is explicitly non-viral does not take the locus down when another frame of
    the same insertion has in-scope viral support. Returns the queries to remove.
    """
    te_por_locus = te_por_locus or {}
    por_locus = {}
    for i in itens:
        chave = (i.get("species", ""), i.get("locus") or i["query"])
        i.update(avaliar(i, te_por_locus.get(chave, ())))
        _com_reads(i, (reads or {}).get(chave, {}))
        _com_hospedeiro(i, (hospedeiro or {}).get(i["query"]))
        _com_repeticoes(i, (simples or {}).get(i["query"]), (repeticoes or {}).get(chave))
        por_locus.setdefault(chave, []).append(i)
    remover = []
    for grupo in por_locus.values():
        if any(x["in_scope_support"] for x in grupo):
            for x in grupo:
                x["auto_reject"] = ""
            continue
        motivo = next((x["auto_reject"] for x in sorted(grupo, key=lambda x: not x.get("locus_representative", True))
                       if x["auto_reject"]), "")
        # only when the representative frame (the row the researcher sees) is itself explicit
        rep = next((x for x in grupo if x.get("locus_representative", True)), grupo[0])
        if not rep["auto_reject"]:
            motivo = ""
        for x in grupo:
            x["auto_reject"] = motivo
            if motivo:
                remover.append(x["query"])
    return remover
