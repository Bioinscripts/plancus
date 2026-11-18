#!/usr/bin/env python3
"""
Homopolymers and simple repeats inside a candidate segment, counted against it in the verdict.

A viral ORF hit that is mostly poly-A, poly-anything or a short tandem repeat is usually a
low-complexity artefact: such stretches translate into runs of one or two amino acids that
align to anything of similar composition. A long poly-A / poly-T also marks retroposed
sequence. Measured straight on the nucleotide sequence of the candidate, so it needs neither
the repeats step nor any external tool:

  homopolymer     a run of one base of at least HOMO_MIN nt (reported as poly-A, poly-T, ...)
  tandem repeat   a unit of 2-6 nt repeated at least 4 times in a row, at least TANDEM_MIN nt
  simple_frac     fraction of the segment inside either kind

`flags` is what the page shows; `strong` is set when a large share of the segment is simple
(simple_frac >= FORTE) and holds the verdict back further.
"""

import os
import re

HOMO_MIN = 10
TANDEM_MIN = 16
FORTE = 0.30
RX_HOMO = re.compile(r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % ((HOMO_MIN,) * 4))
RX_TANDEM = re.compile(r"([ACGT]{2,6}?)\1{3,}")

_cache = {}


def analisar(seq: str) -> dict:
    s = (seq or "").upper()
    n = len(s)
    if not n:
        return {"flags": [], "simple_frac": 0.0, "longest_homopolymer": 0, "strong": False}
    coberto = bytearray(n)
    homos, tandems = [], []
    for m in RX_HOMO.finditer(s):
        homos.append((m.group(0)[0], m.start(), m.end() - m.start()))
        coberto[m.start():m.end()] = b"\x01" * (m.end() - m.start())
    for m in RX_TANDEM.finditer(s):
        unidade = m.group(1)
        if len(set(unidade)) == 1 or m.end() - m.start() < TANDEM_MIN:
            continue                                     # a homopolymer, counted above, or too short
        tandems.append((unidade, m.start(), m.end() - m.start(), (m.end() - m.start()) // len(unidade)))
        coberto[m.start():m.end()] = b"\x01" * (m.end() - m.start())
    frac = sum(coberto) / n
    flags = []
    por_base = {}
    for base, _ini, tam in homos:
        por_base[base] = max(por_base.get(base, 0), tam)
    for base, tam in sorted(por_base.items(), key=lambda x: -x[1]):
        k = sum(1 for b, _i, _t in homos if b == base)
        flags.append(f"poly-{base} of {tam} nt" + (f" (+{k - 1} more)" if k > 1 else ""))
    for unidade, _ini, tam, copias in sorted(tandems, key=lambda x: -x[2])[:2]:
        flags.append(f"tandem repeat ({unidade})x{copias}, {tam} nt")
    if len(tandems) > 2:
        flags.append(f"{len(tandems) - 2} more tandem repeat(s)")
    return {"flags": flags, "simple_frac": round(frac, 3),
            "longest_homopolymer": max((t for _b, _i, t in homos), default=0),
            "n_homopolymers": len(homos), "n_tandem": len(tandems),
            "strong": frac >= FORTE}


def do_fasta(caminho: str) -> dict:
    """{sequence id: analisar(sequence)} of a nucleotide FASTA, cached on its modification time."""
    try:
        chave = (caminho, os.path.getmtime(caminho))
    except OSError:
        return {}
    if _cache.get(caminho, (None,))[0] == chave:
        return _cache[caminho][1]
    saida, nome, partes = {}, None, []
    with open(caminho, encoding="utf-8", errors="replace") as fh:
        for linha in fh:
            if linha.startswith(">"):
                if nome is not None:
                    saida[nome] = analisar("".join(partes))
                nome, partes = linha[1:].split()[0], []
            else:
                partes.append(linha.strip())
    if nome is not None:
        saida[nome] = analisar("".join(partes))
    _cache[caminho] = (chave, saida)
    return saida
