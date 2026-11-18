#!/usr/bin/env python3
"""
viral_positional_analyzer.py
============================
Módulo de Validação Posicional e Evolutiva de Elementos Virais Endógenos (EVEs):

1. Triagem e Alinhamento Posicional:
   - Parsing de alinhamentos em formato tabular estendido (outfmt 6) e BLAST XML (outfmt 5).
   - tBLASTn recíproco de proteínas virais (aminoácidos) contra os segmentos nucleotídicos
     de EVEs candidatas, anotando posições que não se sobreponham via Weighted Interval Scheduling.

2. Positional Viral Consensus Score:
   - Mapeamento base por base de múltiplos vírus alinhando na mesma região da EVE.
   - Pontuação de consenso amplificada quando múltiplos vírus corroboram a mesma coordenada.

3. 2º BLASTX (Banco Geral - Top 100 Best Hits):
   - Captura precisa das posições de alinhamento na query (qstart, qend), identidade (% pident),
     positivos (% ppos / positive) e cobertura de segmento alinhado (% qcovhsp).
   - Classificação taxonômica automática dos hits:
       * VIRAL: linhagens virais de referência.
       * NÃO-VIRAL - PERSISTÊNCIA EM OUTRAS ESPÉCIES: hits em outros genomas eucariotos/artrópodes
         anotados como proteínas hipotéticas, não caracterizadas ou EVEs ortólogas (persistência evolutiva).
       * NÃO-VIRAL - SOBREPOSIÇÃO HOSPEDEIRO: genes celulares estruturais/metabólicos do hospedeiro
         (actina, tubulina, ribossomal, quinase - potenciais falsos positivos).

4. Comparação Posicional Comparativa (Viral vs. Não-Viral):
   - Comparação intervalo por intervalo entre coordenadas virais e não-virais.
   - Avaliação de melhor posição/score, persistência em outras espécies e integridade do elemento.

5. Exportação de Relatórios e Formatos:
   - Suporte completo a TSV tradicional estendido e BLAST XML.
   - Geração de tabelas de sumário, intervalos detalhados, perfil de cobertura e JSON estruturado.
"""

import os
import sys
import re
import csv
import bisect
import shutil
import argparse
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Any, Set
from collections import defaultdict
from Bio import SeqIO
from evetool import config


@dataclass
class AlignmentHsp:
    """Representa um High-Scoring Segment Pair (HSP) normalizado."""
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int = 0
    slen: int = 0
    positive: int = 0
    ppos: float = 0.0
    qcovhsp: float = 0.0
    stitle: str = ""
    qframe: int = 1
    sframe: int = 1
    query_seq: str = ""
    subject_seq: str = ""

    @property
    def q_start(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def q_end(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def s_start(self) -> int:
        return min(self.sstart, self.send)

    @property
    def s_end(self) -> int:
        return max(self.sstart, self.send)

    @property
    def s_strand(self) -> str:
        return '+' if self.send >= self.sstart else '-'

    @property
    def harmonic_score(self) -> float:
        cov = self.qcovhsp if self.qcovhsp > 0 else (self.length / self.qlen * 100.0 if self.qlen > 0 else 0.0)
        if self.pident <= 0 or cov <= 0:
            return 0.0
        return (2.0 * self.pident * cov) / (self.pident + cov)


# =============================================================================
# 1. PARSING DE DADOS DE ALINHAMENTO (TABULAR OUTFMT 6 E XML OUTFMT 5)
# =============================================================================

def parse_alignment_tabular(filepath: str) -> List[AlignmentHsp]:
    """
    Lê arquivo TSV/tabular de alinhamento (DIAMOND ou NCBI BLAST).
    Suporta formato estendido de 18 colunas, padrão de 12 colunas ou variações legadas.
    """
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return []

    hsps: List[AlignmentHsp] = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        first_line = True
        header_map = {}

        for row in reader:
            if not row or not any(row):
                continue
            first_cell = row[0].strip()
            if first_cell.startswith("#"):
                # Linha de comentário do BLAST
                continue

            if first_line:
                first_line = False
                # Checa se é cabeçalho nominal
                lower_row = [c.lower().strip() for c in row]
                if "qseqid" in lower_row and ("sseqid" in lower_row or "bitscore" in lower_row):
                    header_map = {name: idx for idx, name in enumerate(lower_row)}
                    continue

            # Mapeamento dinâmico se houver cabeçalho
            if header_map:
                try:
                    qseqid = row[header_map["qseqid"]].strip().split()[0]
                    sseqid = row[header_map.get("sseqid", header_map.get("subject", 1))].strip().split()[0]
                    pident = float(row[header_map.get("pident", 2)])
                    length = int(row[header_map.get("length", 3)])
                    mismatch = int(row[header_map.get("mismatch", 4)]) if "mismatch" in header_map else 0
                    gapopen = int(row[header_map.get("gapopen", 5)]) if "gapopen" in header_map else 0
                    qstart = int(row[header_map.get("qstart", 6)])
                    qend = int(row[header_map.get("qend", 7)])
                    sstart = int(row[header_map.get("sstart", 8)])
                    send = int(row[header_map.get("send", 9)])
                    evalue = float(row[header_map.get("evalue", 10)])
                    bitscore = float(row[header_map.get("bitscore", 11)])
                    qlen = int(row[header_map["qlen"]]) if "qlen" in header_map else 0
                    slen = int(row[header_map["slen"]]) if "slen" in header_map else 0
                    positive = int(row[header_map["positive"]]) if "positive" in header_map else int(length * (pident / 100.0))
                    ppos = float(row[header_map["ppos"]]) if "ppos" in header_map else pident
                    qcovhsp = float(row[header_map["qcovhsp"]]) if "qcovhsp" in header_map else 0.0
                    stitle = row[header_map["stitle"]].strip() if "stitle" in header_map else sseqid

                    hsps.append(AlignmentHsp(
                        qseqid=qseqid, sseqid=sseqid, pident=pident, length=length,
                        mismatch=mismatch, gapopen=gapopen, qstart=qstart, qend=qend,
                        sstart=sstart, send=send, evalue=evalue, bitscore=bitscore,
                        qlen=qlen, slen=slen, positive=positive, ppos=ppos,
                        qcovhsp=qcovhsp, stitle=stitle
                    ))
                except Exception:
                    continue
                continue

            # Mapeamento posicional baseado no número de colunas
            n_cols = len(row)
            try:
                if n_cols >= 18:
                    # Formato estendido recomendado:
                    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive ppos qcovhsp stitle
                    qseqid = row[0].strip().split()[0]
                    sseqid = row[1].strip().split()[0]
                    pident = float(row[2])
                    length = int(row[3])
                    mismatch = int(row[4])
                    gapopen = int(row[5])
                    qstart = int(row[6])
                    qend = int(row[7])
                    sstart = int(row[8])
                    send = int(row[9])
                    evalue = float(row[10])
                    bitscore = float(row[11])
                    qlen = int(row[12]) if row[12].isdigit() else 0
                    slen = int(row[13]) if row[13].isdigit() else 0
                    positive = int(float(row[14])) if row[14].replace('.', '', 1).isdigit() else int(length * (pident / 100.0))
                    ppos = float(row[15]) if row[15] else pident
                    qcovhsp = float(row[16]) if row[16] else 0.0
                    stitle = row[17].strip() if len(row) > 17 else sseqid

                elif n_cols >= 12:
                    # Formato blast 6 clássico
                    qseqid = row[0].strip().split()[0]
                    sseqid = row[1].strip().split()[0]
                    pident = float(row[2])
                    length = int(row[3])
                    mismatch = int(row[4])
                    gapopen = int(row[5])
                    qstart = int(row[6])
                    qend = int(row[7])
                    sstart = int(row[8])
                    send = int(row[9])
                    evalue = float(row[10])
                    bitscore = float(row[11])
                    qlen = 0
                    slen = 0
                    positive = int(length * (pident / 100.0))
                    ppos = pident
                    qcovhsp = 0.0
                    stitle = row[12].strip() if n_cols > 12 else sseqid

                elif n_cols == 5:
                    # Formato simplificado (qseqid, sseqid, bitscore, evalue, sseq)
                    qseqid = row[0].strip().split()[0]
                    sseqid = row[1].strip().split()[0]
                    bitscore = float(row[2])
                    evalue = float(row[3])
                    sseq = row[4].strip()
                    length = len(sseq)
                    pident = 50.0  # estimativa neutra
                    qstart, qend, sstart, send = 1, length * 3, 1, length
                    mismatch, gapopen = 0, 0
                    qlen, slen = qend, length
                    positive = length
                    ppos = 50.0
                    qcovhsp = 100.0
                    stitle = sseqid

                elif n_cols == 11 and ("ORF" in row[0] or "type:complete" in row[9] or "type:" in row[0]):
                    # Formato da busca primária (qseqid qlen sseqid slen stitle pident evalue bitscore qcovhsp qtitle qseq)
                    qseqid = row[0].strip().split()[0]
                    qlen = int(row[1]) if row[1].isdigit() else 0
                    sseqid = row[2].strip().split()[0]
                    slen = int(row[3]) if row[3].isdigit() else 0
                    stitle = row[4].strip()
                    pident = float(row[5])
                    evalue = float(row[6])
                    bitscore = float(row[7])
                    qcovhsp = float(row[8]) if row[8] else 0.0
                    length = int(qlen * (qcovhsp / 100.0)) if qcovhsp > 0 else 100
                    positive = int(length * (pident / 100.0))
                    ppos = pident
                    qstart = 1
                    qend = max(1, int(qlen * (qcovhsp / 100.0))) if qlen > 0 else 300
                    sstart = 1
                    send = max(1, length)
                    mismatch, gapopen = 0, 0

                else:
                    continue

                hsps.append(AlignmentHsp(
                    qseqid=qseqid, sseqid=sseqid, pident=pident, length=length,
                    mismatch=mismatch, gapopen=gapopen, qstart=qstart, qend=qend,
                    sstart=sstart, send=send, evalue=evalue, bitscore=bitscore,
                    qlen=qlen, slen=slen, positive=positive, ppos=ppos,
                    qcovhsp=qcovhsp, stitle=stitle
                ))
            except Exception:
                continue

    return hsps


def parse_alignment_xml(filepath: str) -> List[AlignmentHsp]:
    """
    Lê arquivo XML do BLAST (outfmt 5) via iterparse para eficiência de memória.
    Extrai coordenadas completas, identidade, positivos, cobertura e anotações.
    """
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return []

    hsps: List[AlignmentHsp] = []
    try:
        context = ET.iterparse(filepath, events=('end',))
        _, root = next(context)
        qseqid = ""
        qlen = 0

        for event, elem in context:
            if elem.tag == 'Iteration_query-def':
                # Captura ID principal da query (antes do espaço)
                text = elem.text or ""
                qseqid = text.strip().split()[0] if text else "Query"
            elif elem.tag == 'Iteration_query-len':
                try:
                    qlen = int(elem.text or 0)
                except ValueError:
                    qlen = 0
            elif elem.tag == 'Hit':
                hit_id_elem = elem.find('Hit_id')
                hit_def_elem = elem.find('Hit_def')
                hit_len_elem = elem.find('Hit_len')

                sseqid = (hit_id_elem.text if hit_id_elem is not None and hit_id_elem.text else "").strip().split()[0]
                stitle = (hit_def_elem.text if hit_def_elem is not None and hit_def_elem.text else sseqid).strip()
                try:
                    slen = int(hit_len_elem.text or 0) if hit_len_elem is not None else 0
                except ValueError:
                    slen = 0

                hsps_elem = elem.find('Hit_hsps')
                if hsps_elem is not None:
                    for hsp in hsps_elem.findall('Hsp'):
                        try:
                            bitscore = float(hsp.findtext('Hsp_bit-score', '0'))
                            evalue = float(hsp.findtext('Hsp_evalue', '1.0'))
                            qstart = int(hsp.findtext('Hsp_query-from', '1'))
                            qend = int(hsp.findtext('Hsp_query-to', '1'))
                            sstart = int(hsp.findtext('Hsp_hit-from', '1'))
                            send = int(hsp.findtext('Hsp_hit-to', '1'))
                            qframe = int(hsp.findtext('Hsp_query-frame', '1'))
                            sframe = int(hsp.findtext('Hsp_hit-frame', '1'))
                            nident = int(hsp.findtext('Hsp_identity', '0'))
                            pos = int(hsp.findtext('Hsp_positive', '0'))
                            gaps = int(hsp.findtext('Hsp_gaps', '0'))
                            align_len = int(hsp.findtext('Hsp_align-len', '0'))

                            pident = (nident / align_len * 100.0) if align_len > 0 else 0.0
                            ppos = (pos / align_len * 100.0) if align_len > 0 else pident
                            mismatch = max(0, align_len - nident - gaps)
                            gapopen = 1 if gaps > 0 else 0
                            hsp_span = abs(qend - qstart) + 1
                            qcov = (hsp_span / qlen * 100.0) if qlen > 0 else 0.0

                            hsps.append(AlignmentHsp(
                                qseqid=qseqid, sseqid=sseqid, pident=round(pident, 2),
                                length=align_len, mismatch=mismatch, gapopen=gapopen,
                                qstart=qstart, qend=qend, sstart=sstart, send=send,
                                evalue=evalue, bitscore=bitscore, qlen=qlen, slen=slen,
                                positive=pos, ppos=round(ppos, 2), qcovhsp=round(qcov, 2),
                                stitle=stitle, qframe=qframe, sframe=sframe
                            ))
                        except Exception:
                            continue

                # Libera memória de elementos já processados
                elem.clear()
                root.clear()
    except Exception as erro:
        print(f"[aviso] Falha ao processar arquivo XML {filepath}: {erro}", file=sys.stderr)

    return hsps


def parse_alignment_file(filepath: str) -> List[AlignmentHsp]:
    """Detecta automaticamente se o arquivo é XML ou Tabular e realiza o parsing."""
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return []

    # Lê os primeiros 200 bytes para detecção
    with open(filepath, 'rb') as f:
        head = f.read(200)

    if b"<?xml" in head or b"<BlastOutput" in head or filepath.endswith(".xml"):
        return parse_alignment_xml(filepath)
    else:
        return parse_alignment_tabular(filepath)


# =============================================================================
# 2. SELEÇÃO DE INTERVALOS NÃO-SOBREPOSTOS (WEIGHTED INTERVAL SCHEDULING)
# =============================================================================

def select_non_overlapping_hsps(
    hsps: List[AlignmentHsp],
    target_is_subject: bool = True,
    metric: str = "bitscore",
    max_overlap_bp: int = 15
) -> List[AlignmentHsp]:
    """
    Seleciona o subconjunto ótimo global de HSPs que NÃO se sobrepõem
    (tolerando no máximo 'max_overlap_bp' nucleotídeos de sobreposição),
    maximizando o somatório da métrica ('bitscore' ou 'harmonic_score').
    Usa Programação Dinâmica (Weighted Interval Scheduling) em O(n log n).
    """
    if not hsps:
        return []

    # Extrai tuplas: (start, end, weight, hsp)
    intervals = []
    for h in hsps:
        if target_is_subject:
            start, end = h.s_start, h.s_end
        else:
            start, end = h.q_start, h.q_end

        if metric == "harmonic_score":
            weight = h.harmonic_score
        else:
            weight = h.bitscore

        intervals.append((start, end, weight, h))

    # 1. Ordena os intervalos pelo ponto final
    intervals.sort(key=lambda x: x[1])
    n = len(intervals)

    # 2. Pré-computa p(j): índice do último intervalo compatível que termina antes do início de j
    ends = [inv[1] for inv in intervals]
    p = [0] * (n + 1)
    for j in range(1, n + 1):
        target_start = intervals[j - 1][0] + max_overlap_bp
        # Busca o maior índice cujo fim <= target_start
        idx = bisect.bisect_right(ends, target_start)
        p[j] = idx

    # 3. DP para calcular a pontuação máxima
    opt = [0.0] * (n + 1)
    for j in range(1, n + 1):
        weight = intervals[j - 1][2]
        opt[j] = max(opt[j - 1], weight + opt[p[j]])

    # 4. Backtracking para recuperar os intervalos escolhidos
    selected: List[AlignmentHsp] = []
    curr = n
    while curr > 0:
        weight = intervals[curr - 1][2]
        if weight + opt[p[curr]] >= opt[curr - 1]:
            selected.append(intervals[curr - 1][3])
            curr = p[curr]
        else:
            curr -= 1

    selected.reverse()
    # Ordena pelo início para apresentação limpa
    if target_is_subject:
        selected.sort(key=lambda x: x.s_start)
    else:
        selected.sort(key=lambda x: x.q_start)

    return selected


# =============================================================================
# 3. CLASSIFICAÇÃO TAXONÔMICA E IDENTIFICAÇÃO DE ORGANISMOS
# =============================================================================

VIRAL_PATTERNS = [
    r"virus", r"viral", r"viridae", r"virales", r"viroid", r"viricetes", r"phage",
    r"bacteriophage", r"polyhedrovirus", r"granulovirus", r"poxvir",
    r"parvovir", r"rhabdovir", r"flavivir", r"reovir", r"bunyavir", r"orthomyxovir",
    r"paramyxovir", r"coronavir", r"hepadnavir", r"nodavir", r"anellovir", r"circovir",
    r"partitivir", r"totivir", r"chuvir", r"jingmenvir", r"quaranjavir", r"nyavir",
    r"bornavir", r"filovir", r"picobirnavir", r"polyprotein", r"capsid", r"nucleocapsid",
    r"glycoprotein", r"envelope", r"large structural protein", r"rna-dependent rna polymerase",
    r"rdrp"
]
_VIRAL_REGEX = re.compile(r"(" + "|".join(VIRAL_PATTERNS) + r")", re.IGNORECASE)

CELLULAR_PATTERNS = [
    r"actin", r"tubulin", r"ribosomal", r"histone", r"ubiquitin", r"atp synthase",
    r"atpase", r"kinase", r"dehydrogenase", r"elongation factor", r"heat shock",
    r"chaperone", r"cytochrome", r"myosin", r"collagen", r"polymerase ii", r"peptidase",
    r"synthetase", r"transferase", r"synthase", r"isomerase", r"ligase", r"transporter",
    r"channel", r"zinc finger", r"dynein", r"kinesin", r"troponin", r"tropomyosin"
]
_CELLULAR_REGEX = re.compile(r"(" + "|".join(CELLULAR_PATTERNS) + r")", re.IGNORECASE)

ORTHOLOGOUS_EVE_PATTERNS = [
    r"uncharacterized", r"hypothetical", r"endogenous viral",
    r"envelope-like", r"polyprotein-like", r"capsid-like",
    r"non-structural protein-like", r"rdrp-like"
]
_ORTHOLOGOUS_REGEX = re.compile(r"(" + "|".join(ORTHOLOGOUS_EVE_PATTERNS) + r")", re.IGNORECASE)

HAM_PATTERNS = [
    r"ham1", r"ham1-like", r"ham-like", r"pyrophosphatase", r"fam003967", r"pf01725",
    r"cassava brown streak", r"maf-like"
]
_HAM_REGEX = re.compile(r"(" + "|".join(HAM_PATTERNS) + r")", re.IGNORECASE)

RETRO_PATTERNS = [
    r"retrovir", r"retroviridae", r"retrovirus", r"retroviral", r"endogenous retrovirus",
    r"\berv\b", r"spumavir", r"lentivir", r"alpharetrovir", r"betaretrovir",
    r"gammaretrovir", r"deltaretrovir", r"foamy virus", r"oncovirus",
    r"\bgag[- ]pol\b", r"\bgag protein\b", r"\bpol protein\b", r"\bretroviral gag\b", r"\bretroviral pol\b"
]
_RETRO_REGEX = re.compile(r"(" + "|".join(RETRO_PATTERNS) + r")", re.IGNORECASE)

RT_PATTERNS = [
    r"reverse transcriptase", r"transcriptase reversa", r"rna-dependent dna polymerase",
    r"\brvt_", r"\brvt\b", r"retrotransposon", r"retroposon", r"ty1/copia", r"ty3/gypsy",
    r"bel-pao", r"line-1", r"\bline\b", r"sine element", r"ltr retrotransposon",
    r"non-ltr retrotransposon", r"\bcopia\b", r"\bgypsy\b"
]
_RT_REGEX = re.compile(r"(" + "|".join(RT_PATTERNS) + r")", re.IGNORECASE)

# Blacklist de accessions RefSeq virais com falsos positivos conhecidos (bacteriófagos,
# genes celulares capturados por vírus como HAM1/GroEL/Calcineurina, contaminantes ou vírus de vertebrados espúrios)
FALSOS_POSITIVOS_REFSEQ_ACCESSION = {
    "NP_040937.1", "NP_620053.1", "NP_776266.1", "YP_001876447.1", "YP_003359177.1",
    "YP_003773387.1", "YP_004063681.1", "YP_004063983.1", "YP_004598981.2", "YP_004767305.1",
    "YP_005271224.1", "YP_005352845.1", "YP_005352853.1", "YP_007027011.1", "YP_007032446.1",
    "YP_007697629.1", "YP_007697636.1", "YP_007697642.1", "YP_009026378.1",
    "YP_009111319.1", "YP_009333345.1", "YP_009333466.1", "YP_009337442.1", "YP_009342322.1",
    "YP_009345058.1", "YP_009448185.1", "YP_009505587.1", "YP_009505589.1",
    "YP_009513020.1", "YP_009551696.1", "YP_009551698.1", "YP_010085116.1", "YP_010086809.1",
    "YP_010087319.1", "YP_010356086.1", "YP_010356144.1", "YP_010356448.1", "YP_010356503.1",
    "YP_010769290.1", "YP_010800375.1", "YP_010800379.1", "YP_010800854.1", "YP_010840026.1",
    "YP_460036.1", "YP_007714631.1", "YP_009506248.1", "YP_009333273.1", "YP_009305422.1"
}
FALSOS_POSITIVOS_REFSEQ_ACCESSION.update({a.split('.')[0] for a in list(FALSOS_POSITIVOS_REFSEQ_ACCESSION)})

REFSEQ_SPURIOUS_PATTERNS = [
    r"\bphage\b", r"\bbacteriophage\b", r"\bprophage\b", r"archaeal virus",
    r"\bgroel\b", r"\bchaperonin\b", r"\bcalcineurin\b", r"\bcdc48\b", r"\baaa-atpase\b",
    r"pestivirus", r"bovine viral diarrhea"
]
_REFSEQ_SPURIOUS_REGEX = re.compile(r"(" + "|".join(REFSEQ_SPURIOUS_PATTERNS) + r")", re.IGNORECASE)

# Lista de domínios HMM / conservados (RVDB FAM / eggNOG) a serem estritamente eliminados
DOMINIOS_EXCLUSAO_HMM = {
    "FAM006507", "FAM010466", "FAM005753", "4R0ZU", "FAM003478", "FAM003906", "FAM007500",
    "FAM020061", "FAM021751", "FAM033562", "4R0XZ", "FAM008720", "4QP6Y", "4R0UT", "4R13J",
    "FAM001232", "FAM001338", "FAM002630", "FAM024525", "FAM024750", "FAM023132", "FAM024758",
    "FAM025425", "FAM025723", "FAM030806", "FAM020950", "FAM022720", "FAM023030", "FAM010785",
    "FAM010808", "FAM019848", "FAM009125", "FAM010068", "FAM008389", "FAM008332", "FAM003101",
    "FAM005852", "FAM005924", "FAM006212", "FAM006283", "FAM006950", "FAM007012", "FAM007251",
    "FAM007270", "FAM008135", "4QP5T", "4R0V2", "4R10Q", "FAM003092", "FAM006504", "FAM006504",
    "FAM006619", "FAM011606","4R0Z1", "FAM008567", "FAM025671", "FAM003100", "FAM003400", "4R0R8",
    "FAM003075"
}
DOMINIOS_EXCLUSAO_HMM.update({d.split('.')[0].upper() for d in list(DOMINIOS_EXCLUSAO_HMM)})


def load_rvdb_domain_keywords(search_paths: Optional[List[str]] = None) -> Dict[str, str]:
    """
    Carrega o dicionário de anotações funcionais e palavras-chave de domínios HMM do RVDB:
    {HMM_ID: Protein_Keywords} (ex: FAM003967 -> '...ham1,inosine,xanthosine,pyrophosphatase...').
    """
    default_paths = [
        config.hmm_path("RVDB", "rvdb_base_annotations.tsv"),
        "rvdb_base_annotations.tsv",
        "../rvdb_base_annotations.tsv",
        "Data_base_hmm_eves/RVDB/rvdb_base_annotations.tsv"
    ]
    paths = (search_paths or []) + default_paths
    rvdb_dict = {}
    for p in paths:
        if p and os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                with open(p, 'r', encoding='utf-8', errors='replace') as f:
                    r = csv.reader(f, delimiter='\t')
                    for row in r:
                        if len(row) >= 3 and not row[0].startswith('#'):
                            rvdb_dict[row[0].strip()] = row[2].strip()
                if rvdb_dict:
                    break
            except Exception:
                continue
    return rvdb_dict


def check_candidate_exclusion(
    eve_id: str,
    top_virus_id: str,
    top_virus_title: str,
    hmm_domain: str,
    hmm_keywords: str,
    general_hits: List[Dict[str, Any]],
    consensus_summary: Dict[str, Any],
    all_domains: Optional[Set[str]] = None
) -> Tuple[bool, str, str]:
    """
    Avalia os critérios estritos de exclusão de candidatos a NIRVs (Non-Retroviral EVEs):
    1. Domínios HMM Proibidos (RVDB FAM / eggNOG especificados pelo usuário).
    2. Domínio HAM / HAM1 (pirofosfatase / HAM1-like CBSV / FAM003967 / PF01725).
    3. Retroviridae / Retrovirais / ERVs.
    4. Transcriptase Reversa (RT / RVT / Ty1-Copia / Ty3-Gypsy / Bel-Pao / LINE).
    5. Transferência Horizontal Pura (HGT) / Presença exclusiva de segmentos do hospedeiro na região sem corroboração viral.

    Retorna:
      (deve_excluir: bool, tipo_filtro: str, motivo_detalhado: str)
    """
    title_comb = f"{top_virus_title} {top_virus_id}".strip().lower()
    hmm_comb = f"{hmm_domain} {hmm_keywords}".strip().lower()

    # 0. Checagem Estrita de Domínios HMM Proibidos pelo Usuário (FAM / eggNOG)
    detected_prohibited_doms = set()
    if all_domains:
        for d in all_domains:
            d_norm = d.split('.')[0].strip().upper()
            if d_norm in DOMINIOS_EXCLUSAO_HMM or d.strip().upper() in DOMINIOS_EXCLUSAO_HMM:
                detected_prohibited_doms.add(d_norm)
    if hmm_domain and hmm_domain != "-":
        h_norm = hmm_domain.split('.')[0].strip().upper()
        if h_norm in DOMINIOS_EXCLUSAO_HMM or hmm_domain.strip().upper() in DOMINIOS_EXCLUSAO_HMM:
            detected_prohibited_doms.add(h_norm)
    for p_dom in DOMINIOS_EXCLUSAO_HMM:
        if p_dom in hmm_comb.upper() or p_dom in title_comb.upper():
            detected_prohibited_doms.add(p_dom)

    if detected_prohibited_doms:
        motivo = f"Presença de domínio(s) conservado(s) proibido(s) [{', '.join(sorted(detected_prohibited_doms))}]"
        return True, "DOMINIO_EXCLUIDO_HMM", motivo

    # SALVAGUARDA OBRIGATÓRIA: Densovírus (Parvoviridae) são EVEs artrópodes genuínas e JAMAIS devem ser excluídas
    if "densovirus" in title_comb or "densovir" in title_comb or "parvovir" in title_comb or "densovirus" in hmm_comb or "parvovir" in hmm_comb:
        return False, "VALIDO", ""

    # 0. Blacklist de Acessos Falsos Positivos do RefSeq Viral (e.g. bacteriófagos, fagos de Flavobacterium, Pestivirus, BVDV, GroEL, etc.)
    top_clean_id = top_virus_id.split('.')[0] if '.' in top_virus_id else top_virus_id
    if top_virus_id in FALSOS_POSITIVOS_REFSEQ_ACCESSION or top_clean_id in FALSOS_POSITIVOS_REFSEQ_ACCESSION:
        motivo = f"Acesso identificado na Blacklist de Falsos Positivos do RefSeq Viral [{top_virus_id}: {top_virus_title[:60] if top_virus_title else 'Falso positivo conhecido'}]"
        return True, "FALSO_POSITIVO_REFSEQ_VIRAL", motivo

    # 1. Bacteriófagos, vírus arqueais, chaperoninas, calcineurina ou vírus espúrios
    if _REFSEQ_SPURIOUS_REGEX.search(title_comb) or _REFSEQ_SPURIOUS_REGEX.search(hmm_comb):
        motivo = f"Bacteriófago / vírus arqueal ou gene celular capturado por vírus [{top_virus_id or hmm_domain}: {top_virus_title[:60] if top_virus_title else hmm_comb[:60]}]"
        return True, "FALSO_POSITIVO_REFSEQ_VIRAL", motivo

    # 2. Domínio HAM / HAM1
    if hmm_domain in ("FAM003967", "PF01725") or _HAM_REGEX.search(hmm_comb) or _HAM_REGEX.search(title_comb):
        motivo = f"Presença de Domínio HAM/HAM1 (pirofosfatase/HAM1-like) [{hmm_domain or top_virus_id}: {top_virus_title[:60] if top_virus_title else hmm_keywords[:60]}]"
        return True, "DOMINIO_HAM", motivo

    # 2. Retroviridae / Retrovirais
    if _RETRO_REGEX.search(title_comb) or _RETRO_REGEX.search(hmm_comb):
        motivo = f"Homologia com vírus da família Retroviridae / elemento retroviral [{top_virus_id or hmm_domain}: {top_virus_title[:60] if top_virus_title else hmm_comb[:60]}]"
        return True, "RETROVIRIDAE", motivo

    # 4. Filtro Estrito para 'Zinc finger' no banco geral
    for gh in general_hits:
        if "zinc finger" in gh.get("stitle", "").lower():
            motivo = f"Presença de domínio 'Zinc finger' detectada no banco geral [{gh['sseqid']}: {gh['stitle'][:60]}]"
            return True, "EXCLUSAO_ZINC_FINGER", motivo

    # 5. Genes de Transcriptase Reversa (RT / Retrotransposons)
    if _RT_REGEX.search(title_comb) or _RT_REGEX.search(hmm_comb):
        motivo = f"Homologia com gene de Transcriptase Reversa (RT/retrotransposon) [{top_virus_id or hmm_domain}: {top_virus_title[:60] if top_virus_title else hmm_comb[:60]}]"
        return True, "TRANSCRIPTASE_REVERSA", motivo

    # 4. HGT / Região Exclusivamente Celular Hospedeira sem Suporte Viral
    if general_hits:
        v_hits = [gh for gh in general_hits if gh["is_viral"]]
        cel_hits = [gh for gh in general_hits if gh["category"] == "SOBREPOSICAO_HOSPEDEIRO" or _CELLULAR_REGEX.search(gh.get("stitle", ""))]
        has_viral_hmm = bool(hmm_domain and hmm_domain != "-" and not _CELLULAR_REGEX.search(hmm_comb))
        if len(cel_hits) > 0 and len(v_hits) == 0 and not has_viral_hmm:
            top_cel = max(cel_hits, key=lambda x: x["bitscore"])
            motivo = f"Região apresenta exclusivamente segmentos do hospedeiro no BLAST [{top_cel['organism']} - {top_cel['stitle'][:50]}] sem corroboração viral (provável HGT / falso positivo celular)"
            return True, "HGT_INCONCLUSIVA_HOSPEDEIRO", motivo

    return False, "VALIDO", ""


def extract_organism_from_title(title: str) -> str:
    """Extrai o nome da espécie/organismo dos colchetes [Organism Name] no cabeçalho NCBI/RefSeq."""
    if not title:
        return "Desconhecido"

    matches = re.findall(r'\[([^\]]+)\]', title)
    if matches:
        # Pega a última ocorrência em colchetes, padrão comum no NCBI
        return matches[-1].strip()

    # Tentativa com pipes do RefSeq: |organism|
    pipe_parts = [p.strip() for p in title.split('|') if p.strip()]
    if len(pipe_parts) >= 3:
        for p in pipe_parts[1:]:
            if " " in p and not p.lower().startswith("gi") and not p.lower().startswith("ref"):
                return p

    return "Desconhecido"


def build_stitle_lookup_dict(search_paths: List[str]) -> Dict[str, str]:
    """
    Constrói um dicionário {accession: title_completo} varrendo arquivos de anotação
    ou FASTAs de referência existentes para resgatar descrições detalhadas.
    Prioriza descrições ricas (com título, organismo e anotação) sobre IDs curtos.
    """
    title_dict = {}
    for p in search_paths:
        if not p or not os.path.exists(p):
            continue
        try:
            if p.endswith((".fasta", ".faa", ".fa")):
                for rec in SeqIO.parse(p, "fasta"):
                    clean_id = rec.id.strip().split()[0]
                    desc = rec.description.strip()
                    if desc and clean_id:
                        if clean_id not in title_dict or len(desc) > len(title_dict[clean_id]):
                            title_dict[clean_id] = desc
            elif p.endswith((".tsv", ".tab")):
                with open(p, 'r', encoding='utf-8', errors='replace') as f:
                    for line in f:
                        if not line.strip() or line.startswith("#"):
                            continue
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) >= 5 and cols[4].strip():
                            # Se for triagem_viral_diamond_refseq_filtrado.tsv (cols[2]=sseqid, cols[4]=stitle)
                            clean_k = cols[2].strip()
                            t = cols[4].strip()
                            if clean_k not in title_dict or len(t) > len(title_dict[clean_k]):
                                title_dict[clean_k] = t
                        elif len(cols) >= 18 and cols[17].strip():
                            clean_k = cols[1].strip()
                            t = cols[17].strip()
                            if clean_k not in title_dict or len(t) > len(title_dict[clean_k]):
                                title_dict[clean_k] = t
        except Exception:
            continue
    return title_dict


def classify_hit_taxonomy(
    stitle: str,
    sseqid: str = "",
    known_viral_accessions: Optional[Set[str]] = None,
    title_dict: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """
    Classifica um hit do BLASTX em VIRAL ou NÃO-VIRAL, discriminando entre:
      - VIRAL: Hit genuíno em linhagem de vírus de RNA não-retroviral conhecida.
      - NÃO-VIRAL - PERSISTÊNCIA EM OUTRAS ESPÉCIES: Hit em hospedeiro eucarioto/artrópode
        anotado como proteína hipotética/não-caracterizada ou elemento endógeno (evidência
        de persistência e conservação ortóloga da EVE entre espécies).
      - NÃO-VIRAL - SOBREPOSIÇÃO HOSPEDEIRO: Hit em gene celular conservado (actina, tubulina, etc.),
        indicando risco de falso positivo caso o score celular supere o viral.
      - FILTRO ESPECÍFICO (HAM / RETRO / RT): Removido da linha de frente de NIRVs genuínos.
    """
    if title_dict and sseqid in title_dict:
        stitle = title_dict[sseqid]

    title_norm = stitle.lower()
    organism = extract_organism_from_title(stitle)
    org_norm = organism.lower()

    # 0. Blacklist de Falsos Positivos RefSeq (Bacteriófagos, Genes Celulares Capturados, Vírus Espúrios)
    # Salvaguarda: Densovírus / Parvoviridae são expressamente preservados
    is_denso = ("densovirus" in title_norm or "densovir" in title_norm or "parvovir" in title_norm or "densovirus" in org_norm or "parvovir" in org_norm)
    sseqid_base = sseqid.split('.')[0] if '.' in sseqid else sseqid
    if not is_denso and (sseqid in FALSOS_POSITIVOS_REFSEQ_ACCESSION or sseqid_base in FALSOS_POSITIVOS_REFSEQ_ACCESSION or _REFSEQ_SPURIOUS_REGEX.search(title_norm) or _REFSEQ_SPURIOUS_REGEX.search(org_norm)):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "EXCLUSAO_FALSO_POSITIVO",
            "label": "Falso Positivo RefSeq (Bacteriófago / Gene Capturado / Não-NIRV)",
            "description": f"Acesso ou termo na blacklist de falsos positivos do RefSeq [{sseqid}: {organism}]."
        }

    # 0.1 Domínios HMM Proibidos (RVDB FAM / eggNOG)
    if sseqid_base.upper() in DOMINIOS_EXCLUSAO_HMM or any(dom in title_norm.upper() for dom in DOMINIOS_EXCLUSAO_HMM):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "EXCLUSAO_DOMINIO_HMM",
            "label": "Domínio HMM Excluído",
            "description": f"Hit associado a domínio conservado proibido [{sseqid}: {organism}]."
        }

    # 1. Filtros de Exclusão Específica: HAM, Retroviridae, Transcriptase Reversa
    if _HAM_REGEX.search(title_norm) or _HAM_REGEX.search(sseqid):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "EXCLUSAO_HAM",
            "label": "Domínio HAM / HAM1 (Exclusão NIRV)",
            "description": f"Proteína ou domínio HAM1 (pirofosfatase/CBSV) [{organism}]."
        }

    if _RETRO_REGEX.search(title_norm) or _RETRO_REGEX.search(org_norm):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "EXCLUSAO_RETROVIRUS",
            "label": "Retroviridae / Retroviral (Exclusão NIRV)",
            "description": f"Vírus da família Retroviridae / elemento retroviral [{organism}]."
        }

    if _RT_REGEX.search(title_norm):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "EXCLUSAO_TRANSCRIPTASE_REVERSA",
            "label": "Transcriptase Reversa / Retrotransposon (Exclusão NIRV)",
            "description": f"Gene de transcriptase reversa ou retroelemento móvel [{organism}]."
        }

    # Se o sseqid pertence ao pool de referência viral e não foi excluído pelos filtros acima, é viral confirmado
    is_known_ref_virus = bool(known_viral_accessions and sseqid in known_viral_accessions)

    # 2. Checa se o organismo ou título indica vírus explícito
    is_viral_organism = bool(_VIRAL_REGEX.search(org_norm))
    is_viral_title = bool(_VIRAL_REGEX.search(title_norm))

    # Exceções onde "virus" ou "viral" está em contexto de "endogenous viral element"
    if "endogenous viral" in title_norm or "endogenous retrovirus" in title_norm:
        # Se está em um artrópode/hospedeiro, é persistência da EVE!
        if not is_viral_organism and organism != "Desconhecido":
            return {
                "is_viral": False,
                "organism": organism,
                "category": "PERSISTENCIA_OUTRAS_ESPECIES",
                "label": "EVE Ortóloga Persistente em Outro Genoma",
                "description": f"Elemento viral endógeno anotado no hospedeiro [{organism}]."
            }

    if is_known_ref_virus or is_viral_organism or (is_viral_title and not _CELLULAR_REGEX.search(title_norm)):
        # Hit viral genuíno
        viral_family = "Vírus Não Classificado"
        for fam in ["Rhabdoviridae", "Flaviviridae", r"Parvoviridae", "Bunyaviridae", "Phenuiviridae",
                    "Totiviridae", "Partitiviridae", "Chuviridae", "Orthomyxoviridae", "Poxviridae"]:
            if re.search(fam, title_norm, re.IGNORECASE):
                viral_family = fam
                break

        return {
            "is_viral": True,
            "organism": organism if organism != "Desconhecido" else "Viral RefSeq",
            "category": "VIRAL",
            "label": f"Vírus ({viral_family})",
            "description": f"Hit viral de referência [{organism}]."
        }

    # 3. Hit não-viral: avaliar se é gene celular ou persistência
    if _CELLULAR_REGEX.search(title_norm):
        return {
            "is_viral": False,
            "organism": organism,
            "category": "SOBREPOSICAO_HOSPEDEIRO",
            "label": "Gene Celular do Hospedeiro",
            "description": f"Homologia com maquinária celular do hospedeiro [{organism}]."
        }

    if _ORTHOLOGOUS_REGEX.search(title_norm) or organism != "Desconhecido":
        return {
            "is_viral": False,
            "organism": organism,
            "category": "PERSISTENCIA_OUTRAS_ESPECIES",
            "label": "Persistência em Outros Genomas",
            "description": f"Segmento homólogo conservado no genoma de [{organism}]."
        }

    return {
        "is_viral": False,
        "organism": organism,
        "category": "OUTRO_NAO_VIRAL",
        "label": "Não-Viral Diversificado",
        "description": f"Proteína não-viral [{organism}]."
    }


def load_hmm_domain_hits(search_paths: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Carrega resultados de anotação de domínios conservados (HMMER domtblout / triagem_profunda_hmm.tsv).
    Mapeia tanto pelo ID da ORF/alvo quanto por coordenadas genômicas [inicio-fim].
    """
    hmm_map = defaultdict(list)
    for p in search_paths:
        if not p or not os.path.exists(p):
            continue
        try:
            with open(p, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 8:
                        continue
                    target_name = parts[0].strip()
                    domain_name = parts[3].strip()
                    try:
                        evalue = float(parts[6])
                        score = float(parts[7])
                    except (ValueError, IndexError):
                        continue

                    desc = " ".join(parts[22:]) if len(parts) > 22 else domain_name
                    hit_info = {
                        "domain": domain_name,
                        "evalue": evalue,
                        "score": score,
                        "description": desc
                    }
                    hmm_map[target_name].append(hit_info)

                    # Extrai coordenada [inicio-fim](strand) e cria alias scaffold:inicio-fim(strand)
                    m = re.search(r'\[(\d+)-(\d+)\]\(([+-])\)', desc)
                    if m and "_" in target_name:
                        scaff = target_name.split("_ORF")[0].strip()
                        alias = f"{scaff}:{m.group(1)}-{m.group(2)}({m.group(3)})"
                        hmm_map[alias].append(hit_info)
        except Exception:
            continue

    # Ordena hits de cada alvo por e-value crescente
    for k in hmm_map:
        hmm_map[k].sort(key=lambda x: x["evalue"])
    return hmm_map


def load_blastn_validation_hits(search_paths: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Carrega validações de alinhamento a nível de nucleotídeo (BLASTN).
    Mapeia query ID para o melhor alinhamento nucleotídico viral encontrado.
    """
    blastn_map = {}
    for p in search_paths:
        if not p or not os.path.exists(p):
            continue
        try:
            with open(p, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.split("\t") if "\t" in line else line.split()
                    if len(parts) < 11:
                        continue
                    qid = parts[0].strip().split()[0]
                    sid = parts[1].strip().split()[0]
                    try:
                        pident = float(parts[2])
                        length = int(parts[3])
                        evalue = float(parts[10])
                        bitscore = float(parts[11]) if len(parts) > 11 else 0.0
                    except (ValueError, IndexError):
                        continue

                    curr = blastn_map.get(qid)
                    if not curr or bitscore > curr["bitscore"]:
                        blastn_map[qid] = {
                            "sseqid": sid,
                            "pident": pident,
                            "length": length,
                            "evalue": evalue,
                            "bitscore": bitscore
                        }
        except Exception:
            continue
    return blastn_map


def load_cluster_segmentation_map(search_paths: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Carrega classificação de grupos/clusters (flebo_cluster_v2_segmentos_classificacao.tsv).
    Mapeia segmentos aos seus respectivos clusters e membros parceiros.
    """
    cluster_map = {}
    group_members = defaultdict(list)
    for p in search_paths:
        if not p or not os.path.exists(p):
            continue
        try:
            with open(p, 'r', encoding='utf-8', errors='replace') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    seg = row.get("segmento", "").strip()
                    grp = row.get("grupo", "").strip()
                    if seg and grp:
                        group_members[grp].append(seg)
                        cluster_map[seg] = {
                            "grupo": grp,
                            "categoria": row.get("categoria", ""),
                            "status_acao": row.get("status_acao", ""),
                            "identidade_max": float(row.get("identidade_max", 0.0) or 0.0),
                            "cobertura_max": float(row.get("cobertura_max", 0.0) or 0.0),
                            "melhor_parceira": row.get("melhor_parceira", "")
                        }
        except Exception:
            continue

    # Anota parceiros do mesmo grupo
    for seg, info in cluster_map.items():
        grp = info["grupo"]
        info["membros_grupo"] = group_members.get(grp, [])
        info["n_membros_grupo"] = len(info["membros_grupo"])

    return cluster_map


def evaluate_viral_quality_tags(
    eve_id: str,
    consensus_summary: Dict[str, Any],
    general_hits: List[Dict[str, Any]],
    hmm_hits: Dict[str, List[Dict[str, Any]]],
    blastn_hits: Dict[str, Dict[str, Any]],
    cluster_map: Dict[str, Dict[str, Any]],
    all_virus_to_eves: Dict[str, List[str]]
) -> Dict[str, Any]:
    """
    Avalia minuciosamente a qualidade biológica e realidade viral do candidato a EVE,
    gerando as tags evolutivas, razão de positivos/identidade (BLOSUM62) e faixas de confiança:

    1. tag_ranking_viral:
       - BESTHIT_VIRAL (1º hit no banco geral é vírus)
       - VIRUS_TOP2_SOBREPOSTO (2º hit é vírus na mesma posição de hit não-viral, overlap >= 30%)
       - VIRUS_TOP2_ADJACENTE (2º hit é vírus em posição adjacente)
       - VIRUS_TOP3_PLUS (Rank N) (hit viral a partir do 3º lugar)
       - VIRAL_EXCLUSIVO_TRIAGEM (apenas hits virais, nenhum não-viral no banco geral)
       - SEM_HIT_VIRAL_GERAL (nenhum hit viral no banco geral)

    2. tag_dominio_conservado:
       - DOMINIO_HMM_CONFIRMADO: Pfam / vFam identificado via HMMER (RdRp, Capsídeo, etc.)
       - SEM_DOMINIO_HMM

    3. tag_blastn:
       - BLASTN_VIRAL_CONFIRMADO (alinhamento nucleotídico com genoma viral)
       - BLASTN_DIVERGENTE_ANCESTRAL (divergência nucleotídica antiga com saturação de bases)

    4. tag_divisao_proteina:
       - PROTEINA_VIRAL_FRAGMENTADA (N segmentos agrupados ou mapeando a mesma proteína viral de referência)
       - PROTEINA_INTEGRA (Locus Único)
       - SEM_MAPEAMENTO

    5. tag_deriva_positivos e razao_positivos_identidade:
       - R_pos = ppos / pident (índice de conservação de substituições conservativas)
       - DERIVA_CONSERVATIVA_ELEVADA: R_pos >= 1.30 com conservação estrutural de aminoácidos
       - HOMOLOGIA_ESTRITA_CONSERVADA: pident >= 55%
       - COMPETICAO_CELULAR_DESFAVORAVEL: hit celular hospedeiro supera o viral
       - DERIVA_MODERADA

    6. score_realidade_viral (0 - 100) e faixa_realidade_viral:
       - Faixa 1 (Tier 1): Alta Certeza (>= 70)
       - Faixa 2 (Tier 2): Ortóloga Persistente (50 - 69.9)
       - Faixa 3 (Tier 3): Ancestral Divergente / Fragmentada (20 - 49.9)
       - Faixa 4 (Tier 4): Alerta Celular / Incerteza (< 20)
    """
    c_sum = consensus_summary.get(eve_id, {})
    top_v_id = c_sum.get("top_virus_id", "Nenhum")
    v_cov = float(c_sum.get("viral_coverage_pct", 0.0) or 0.0)
    v_bs = float(c_sum.get("max_consensus_score", 0.0) or c_sum.get("top_virus_bitscore", 0.0) or 0.0)
    v_ident = float(c_sum.get("top_virus_ident", 0.0) or 0.0)
    v_ppos = float(c_sum.get("top_virus_ppos", 0.0) or 0.0)
    v_depth = int(c_sum.get("max_viral_depth", 0) or 0)

    # 1. Ranking Viral no Banco Geral (Top 100)
    tag_ranking = "SEM_HIT_VIRAL_GERAL"
    viral_rank = None
    primeiro_viral_hit = None
    melhor_nv_hit = None

    if general_hits:
        for idx, gh in enumerate(general_hits):
            if gh["is_viral"] and primeiro_viral_hit is None:
                primeiro_viral_hit = gh
                viral_rank = idx + 1
            elif not gh["is_viral"] and melhor_nv_hit is None:
                melhor_nv_hit = gh

        if viral_rank == 1:
            tag_ranking = "BESTHIT_VIRAL"
        elif viral_rank == 2:
            # Checa sobreposição na query entre Hit 1 (não-viral) e Hit 2 (viral)
            h1 = general_hits[0]
            h2 = general_hits[1]
            ov_start = max(h1["q_start"], h2["q_start"])
            ov_end = min(h1["q_end"], h2["q_end"])
            ov_len = max(0, ov_end - ov_start + 1)
            min_len = min(abs(h1["q_end"] - h1["q_start"]) + 1, abs(h2["q_end"] - h2["q_start"]) + 1)
            ov_pct = (ov_len / min_len * 100.0) if min_len > 0 else 0.0

            if ov_pct >= 30.0:
                tag_ranking = "VIRUS_TOP2_SOBREPOSTO"
            else:
                tag_ranking = "VIRUS_TOP2_ADJACENTE"
        elif viral_rank and viral_rank >= 3:
            tag_ranking = f"VIRUS_TOP3_PLUS (Rank {viral_rank})"
        elif viral_rank is None and melhor_nv_hit is not None:
            tag_ranking = "SEM_HIT_VIRAL_GERAL"
        elif viral_rank is None and v_cov > 0:
            tag_ranking = "VIRAL_EXCLUSIVO_TRIAGEM"
    else:
        if v_cov > 0 or top_v_id != "Nenhum":
            tag_ranking = "VIRAL_EXCLUSIVO_TRIAGEM"
        else:
            tag_ranking = "SEM_HOMOLOGIA_DETECTADA"

    # Se hit viral do banco geral existir, atualiza métricas virais de referência
    if primeiro_viral_hit:
        if v_ident <= 0:
            v_ident = float(primeiro_viral_hit["pident"])
        if v_ppos <= 0:
            v_ppos = float(primeiro_viral_hit["ppos"])
        if v_bs <= 0:
            v_bs = float(primeiro_viral_hit["bitscore"])

    # 2. Domínio Conservado (HMMER)
    h_hits = hmm_hits.get(eve_id, [])
    if not h_hits:
        for k, vlist in hmm_hits.items():
            if k == eve_id or (eve_id in k) or (k in eve_id):
                h_hits = vlist
                break

    if h_hits:
        best_hmm = h_hits[0]
        tag_hmm = f"DOMINIO_HMM_CONFIRMADO ({best_hmm['domain']})"
        hmm_domain_name = best_hmm["domain"]
        hmm_evalue = best_hmm["evalue"]
        hmm_score = best_hmm["score"]
    else:
        tag_hmm = "SEM_DOMINIO_HMM"
        hmm_domain_name = "-"
        hmm_evalue = 1.0
        hmm_score = 0.0

    # 3. BLASTN (Nível Nucleotídico)
    b_hit = blastn_hits.get(eve_id)
    if not b_hit:
        for k, bval in blastn_hits.items():
            if k == eve_id or (eve_id in k) or (k in eve_id):
                b_hit = bval
                break

    if b_hit and (b_hit["pident"] >= 55.0 or b_hit["evalue"] <= 1e-4):
        tag_blastn = f"BLASTN_VIRAL_CONFIRMADO ({b_hit['pident']:.1f}% id)"
        blastn_ident = b_hit["pident"]
        blastn_eval = b_hit["evalue"]
    else:
        tag_blastn = "BLASTN_DIVERGENTE_ANCESTRAL" if blastn_hits else "SEM_DADOS_BLASTN"
        blastn_ident = 0.0
        blastn_eval = 1.0

    # 4. Divisão da Proteína no Cluster / Segmentos
    mapped_eves = all_virus_to_eves.get(top_v_id, []) if top_v_id != "Nenhum" else []
    cl_info = cluster_map.get(eve_id, {})
    cl_members = cl_info.get("membros_grupo", [])

    n_frags = max(len(mapped_eves), len(cl_members))
    if n_frags >= 2:
        tag_divisao = f"PROTEINA_VIRAL_FRAGMENTADA ({n_frags} segmentos)"
    elif top_v_id != "Nenhum" or v_cov > 0:
        tag_divisao = "PROTEINA_INTEGRA (Locus Único)"
    else:
        tag_divisao = "SEM_MAPEAMENTO"

    # 5. Deriva Evolutiva & Razão Positivos/Identidade (BLOSUM62)
    nv_ident = float(melhor_nv_hit["pident"]) if melhor_nv_hit else 0.0
    nv_ppos = float(melhor_nv_hit["ppos"]) if melhor_nv_hit else 0.0
    nv_bs = float(melhor_nv_hit["bitscore"]) if melhor_nv_hit else 0.0
    nv_cat = melhor_nv_hit["category"] if melhor_nv_hit else "NENHUM"

    razao_pos_ident = round(v_ppos / max(1.0, v_ident), 2) if v_ident > 0 else 0.0
    delta_ident = round(v_ident - nv_ident, 1)
    delta_ppos = round(v_ppos - nv_ppos, 1)
    delta_bs = round(v_bs - nv_bs, 1)

    # Heurística para falsos positivos do Swiss-Prot
    if top_v_id != "Nenhum" and ("sp|" in top_v_id or "tr|" in top_v_id):
        nv_cat = "SOBREPOSICAO_HOSPEDEIRO"
        nv_bs = max(nv_bs, v_bs + 50)
        nv_ident = max(nv_ident, v_ident + 30)

    tem_alerta_celular_forte = (nv_cat == "SOBREPOSICAO_HOSPEDEIRO" and (nv_bs > v_bs + 25 or nv_ident > v_ident + 15))

    if tem_alerta_celular_forte:
        tag_deriva = "COMPETICAO_CELULAR_DESFAVORAVEL"
    elif nv_ident > v_ident + 20 and nv_bs > v_bs + 30:
        tag_deriva = "HOMOLOGIA_NAO_VIRAL_SUPERIOR"
    elif razao_pos_ident >= 1.30 and v_ppos >= 40.0:
        tag_deriva = "DERIVA_CONSERVATIVA_ELEVADA"
    elif v_ident >= 55.0:
        tag_deriva = "HOMOLOGIA_ESTRITA_CONSERVADA"
    elif v_ident > 0:
        tag_deriva = "DERIVA_MODERADA"
    else:
        tag_deriva = "SEM_DADOS_DERIVA"

    # 6. Score de Realidade Viral (0 - 100 pontos)
    score_realidade = 0.0

    # 6.1 Cobertura e Força Viral (0 a 25 pts)
    if v_cov >= 70.0: score_realidade += 15.0
    elif v_cov >= 35.0: score_realidade += 10.0
    elif v_cov >= 15.0: score_realidade += 5.0

    if v_bs >= 100.0: score_realidade += 10.0
    elif v_bs >= 50.0: score_realidade += 7.0
    elif v_bs >= 30.0: score_realidade += 4.0

    # 6.2 Ranking no Banco Geral (0 a 25 pts)
    if tag_ranking in ("BESTHIT_VIRAL", "VIRAL_EXCLUSIVO_TRIAGEM"):
        score_realidade += 25.0
    elif tag_ranking == "VIRUS_TOP2_SOBREPOSTO":
        score_realidade += 18.0
    elif tag_ranking == "VIRUS_TOP2_ADJACENTE":
        score_realidade += 14.0
    elif "VIRUS_TOP3_PLUS" in tag_ranking:
        score_realidade += 8.0

    # 6.3 Domínio HMM Conservado (0 a 20 pts)
    if "DOMINIO_HMM_CONFIRMADO" in tag_hmm:
        score_realidade += 20.0

    # 6.4 Assinatura de Deriva e Positivos (0 a 15 pts)
    if tag_deriva == "DERIVA_CONSERVATIVA_ELEVADA":
        score_realidade += 15.0
    elif tag_deriva == "HOMOLOGIA_ESTRITA_CONSERVADA":
        score_realidade += 15.0
    elif tag_deriva == "DERIVA_MODERADA":
        score_realidade += 8.0
    elif tag_deriva == "COMPETICAO_CELULAR_DESFAVORAVEL":
        score_realidade -= 25.0

    # 6.5 Consenso Posicional & Fragmentação (0 a 15 pts)
    if v_depth >= 3: score_realidade += 6.0
    elif v_depth >= 2: score_realidade += 3.0

    if "FRAGMENTADA" in tag_divisao:
        score_realidade += 5.0
    if "CONFIRMADO" in tag_blastn:
        score_realidade += 4.0

    score_realidade = max(0.0, min(100.0, round(score_realidade, 1)))

    # Faixas de Realidade Viral (Tiers 1 a 4)
    if score_realidade >= 70.0 and not tem_alerta_celular_forte:
        faixa_realidade = "Faixa 1 - Alta Certeza"
        faixa_rotulo = "Tier 1: Alta Certeza Viral"
    elif score_realidade >= 50.0 and not tem_alerta_celular_forte:
        faixa_realidade = "Faixa 2 - Ortóloga Persistente"
        faixa_rotulo = "Tier 2: Ortóloga Persistente Multi-Espécies"
    elif score_realidade >= 20.0 and not tem_alerta_celular_forte:
        faixa_realidade = "Faixa 3 - Ancestral Divergente"
        faixa_rotulo = "Tier 3: Ancestral Divergente / Fragmentada"
    else:
        faixa_realidade = "Faixa 4 - Baixa Certeza / Sobreposição Celular"
        faixa_rotulo = "Tier 4: Baixa Certeza / Sobreposição Celular"

    # Tags Resumo (para filtro rápido no DataTable)
    tags_ativas = []
    if "BESTHIT" in tag_ranking: tags_ativas.append("BESTHIT_VIRAL")
    elif "TOP2_SOBREPOSTO" in tag_ranking: tags_ativas.append("VIRUS_TOP2_SOBREPOSTO")
    elif "TOP2_ADJACENTE" in tag_ranking: tags_ativas.append("VIRUS_TOP2_ADJACENTE")
    elif "TOP3" in tag_ranking: tags_ativas.append("VIRUS_TOP3_PLUS")
    elif "EXCLUSIVO" in tag_ranking: tags_ativas.append("VIRAL_EXCLUSIVO")

    if "DOMINIO_HMM_CONFIRMADO" in tag_hmm: tags_ativas.append(f"HMM_{hmm_domain_name}")
    if "CONFIRMADO" in tag_blastn: tags_ativas.append("BLASTN_CONFIRMADO")
    if "FRAGMENTADA" in tag_divisao: tags_ativas.append("PROTEINA_FRAGMENTADA")
    if "ELEVADA" in tag_deriva: tags_ativas.append("DERIVA_CONSERVATIVA")
    if tem_alerta_celular_forte: tags_ativas.append("SOBREPOSICAO_CELULAR")

    return {
        "faixa_realidade_viral": faixa_realidade,
        "faixa_rotulo": faixa_rotulo,
        "score_realidade_viral": score_realidade,
        "tag_ranking_viral": tag_ranking,
        "tag_dominio_conservado": tag_hmm,
        "hmm_domain_name": hmm_domain_name,
        "hmm_evalue": hmm_evalue,
        "hmm_score": hmm_score,
        "tag_blastn": tag_blastn,
        "blastn_ident": blastn_ident,
        "blastn_eval": blastn_eval,
        "tag_divisao_proteina": tag_divisao,
        "tag_deriva_positivos": tag_deriva,
        "razao_positivos_identidade": razao_pos_ident,
        "delta_ident_v_vs_nv": delta_ident,
        "delta_ppos_v_vs_nv": delta_ppos,
        "delta_bitscore_v_vs_nv": delta_bs,
        "tags_resumo": "; ".join(tags_ativas) if tags_ativas else "SEM_TAGS"
    }


# =============================================================================
# 4. tBLASTn RECÍPROCO (PROTEÍNA VIRAL AA VS. SEGMENTO EVE NT)
# =============================================================================

def run_reciprocal_tblastn(
    fasta_candidatas_nt: str,
    fasta_virais_aa: str,
    diretorio_saida: str,
    threads: int = 4,
    evalue: float = 1e-3,
    export_xml: bool = False
) -> Dict[str, Any]:
    """
    Executa tBLASTn recíproco das proteínas virais de referência (AA)
    contra os segmentos nucleotídicos de putative EVEs (NT).
    Anota as posições de alinhamento que NÃO se sobrepõem via Weighted Interval Scheduling
    e agrupa todas as posições para pontuação de consenso.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    db_candidatas = os.path.join(diretorio_saida, "eves_candidatas_blastdb")
    tsv_saida = os.path.join(diretorio_saida, "tblastn_reciproco_virais_vs_eves.tsv")
    xml_saida = os.path.join(diretorio_saida, "tblastn_reciproco_virais_vs_eves.xml") if export_xml else None

    # 1. makeblastdb do banco de candidatas NT
    cmd_db = [
        "makeblastdb",
        "-in", fasta_candidatas_nt,
        "-dbtype", "nucl",
        "-out", db_candidatas
    ]
    try:
        subprocess.run(cmd_db, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as err:
        raise RuntimeError(f"Erro ao formatar makeblastdb de candidatas NT:\n{err.stderr}")

    # 2. Execução do tBLASTn (Query: Proteínas Virais, Subject: Segmentos de EVEs)
    # Formato com 18 campos completos para máxima granularidade
    outfmt_str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen positive ppos qcovhsp stitle"
    cmd_tblastn = [
        "tblastn",
        "-query", fasta_virais_aa,
        "-db", db_candidatas,
        "-out", tsv_saida,
        "-outfmt", outfmt_str,
        "-evalue", str(evalue),
        "-num_threads", str(threads)
    ]

    try:
        subprocess.run(cmd_tblastn, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as err:
        raise RuntimeError(f"Erro ao executar tBLASTn recíproco:\n{err.stderr}")

    # Execução opcional em XML (outfmt 5)
    if export_xml and xml_saida:
        cmd_xml = [
            "tblastn",
            "-query", fasta_virais_aa,
            "-db", db_candidatas,
            "-out", xml_saida,
            "-outfmt", "5",
            "-evalue", str(evalue),
            "-num_threads", str(threads)
        ]
        try:
            subprocess.run(cmd_xml, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except Exception as err:
            print(f"[aviso] Falha ao exportar XML do tBLASTn recíproco: {err}", file=sys.stderr)

    # 3. Parsing e seleção de HSPs não-sobrepostos por segmento de EVE (sseqid)
    all_hsps = parse_alignment_tabular(tsv_saida)
    hsps_por_eve = defaultdict(list)
    for h in all_hsps:
        hsps_por_eve[h.sseqid].append(h)

    non_overlapping_por_eve = {}
    tsv_nao_sobreposto = os.path.join(diretorio_saida, "validacao_posicional_tblastn_reciproco_nao_sobreposto.tsv")

    with open(tsv_nao_sobreposto, 'w', encoding='utf-8') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            "eve_id", "eve_nt_inicio", "eve_nt_fim", "eve_fita", "comprimento_alinhado_nt",
            "viral_proteina_id", "viral_titulo", "identidade_pct", "positivos_pct",
            "bitscore", "evalue", "cobertura_hsp_pct", "harmonic_score"
        ])

        for eve_id, lista in hsps_por_eve.items():
            # Seleciona HSPs não-sobrepostos nas coordenadas do subject (o segmento da EVE)
            chosen = select_non_overlapping_hsps(lista, target_is_subject=True, metric="bitscore")
            non_overlapping_por_eve[eve_id] = chosen

            for h in chosen:
                span_nt = abs(h.send - h.sstart) + 1
                writer.writerow([
                    eve_id, h.s_start, h.s_end, h.s_strand, span_nt,
                    h.qseqid, h.stitle or h.qseqid, h.pident, h.ppos,
                    h.bitscore, h.evalue, h.qcovhsp, round(h.harmonic_score, 2)
                ])

    print(f"[run] tBLASTn recíproco: {len(all_hsps)} hits brutos processados. "
          f"Anotações não-sobrepostas salvas em: {tsv_nao_sobreposto}")

    return {
        "tsv_bruto": tsv_saida,
        "xml_saida": xml_saida,
        "tsv_nao_sobreposto": tsv_nao_sobreposto,
        "all_hsps": all_hsps,
        "non_overlapping_por_eve": non_overlapping_por_eve
    }


# =============================================================================
# 5. CÁLCULO DO POSITIONAL VIRAL CONSENSUS SCORE
# =============================================================================

def compute_viral_positional_consensus(
    fasta_candidatas_nt: str,
    viral_hsps: List[AlignmentHsp],
    eve_lengths: Optional[Dict[str, int]] = None
) -> Dict[str, Any]:
    """
    Calcula a pontuação de consenso posicional dos vírus sobre cada segmento de EVE:
    - Para cada base ou janela nucleotídica:
        * Identifica quantos vírus distintos se alinharam naquela coordenada (viral_depth).
        * Aplica o Viral Consensus Score:
            Score(p) = (Soma dos Bitscores) * (1.0 + 0.25 * (viral_depth - 1))
        * Coleta a lista de linhagens/vírus corroboradores.
    - Condensa as posições em intervalos contíguos com resumo estatístico.
    """
    # Carrega comprimentos reais dos segmentos de EVE a partir do FASTA
    lengths = {}
    if eve_lengths:
        lengths.update(eve_lengths)
    if os.path.exists(fasta_candidatas_nt):
        for rec in SeqIO.parse(fasta_candidatas_nt, "fasta"):
            clean_id = rec.id.split()[0]
            lengths[clean_id] = len(rec.seq)

    # Agrupa HSPs virais por ID de EVE
    hsps_por_eve = defaultdict(list)
    for h in viral_hsps:
        # Se for hit de tBLASTn recíproco: subject é a EVE (sseqid), query é o vírus (qseqid)
        # Se for blastx inicial: query é a EVE (qseqid), subject é o vírus (sseqid)
        eve_id = h.sseqid if h.sseqid in lengths or "ORF" in h.sseqid else h.qseqid
        hsps_por_eve[eve_id].append(h)

    consensus_summary_por_eve = {}
    intervals_detalhados = []
    base_depth_records = []

    for eve_id, seq_len in lengths.items():
        hits = hsps_por_eve.get(eve_id, [])
        if not hits:
            consensus_summary_por_eve[eve_id] = {
                "eve_length": seq_len,
                "viral_covered_bases": 0,
                "viral_coverage_pct": 0.0,
                "max_viral_depth": 0,
                "mean_viral_depth": 0.0,
                "max_consensus_score": 0.0,
                "distinct_viruses_count": 0,
                "distinct_viruses": [],
                "top_virus_id": "Nenhum",
                "top_virus_bitscore": 0.0,
                "top_virus_ident": 0.0,
                "top_virus_ppos": 0.0
            }
            continue

        # Array base por base de suporte viral (1-indexed, tamanho seq_len + 1)
        depth_array = [0] * (seq_len + 1)
        bitscore_sum = [0.0] * (seq_len + 1)
        pident_sum = [0.0] * (seq_len + 1)
        ppos_sum = [0.0] * (seq_len + 1)
        viruses_per_pos = [set() for _ in range(seq_len + 1)]

        distinct_viruses_all = set()
        top_hit = max(hits, key=lambda h: h.bitscore)

        for h in hits:
            # Determina coordenadas e identificador do vírus
            if h.sseqid == eve_id:
                # Subject é a EVE
                start, end = h.s_start, h.s_end
                virus_id = h.qseqid
            else:
                # Query é a EVE
                start, end = h.q_start, h.q_end
                virus_id = h.sseqid

            distinct_viruses_all.add(virus_id)
            start_bounded = max(1, min(start, seq_len))
            end_bounded = max(1, min(end, seq_len))

            for pos in range(start_bounded, end_bounded + 1):
                depth_array[pos] += 1
                bitscore_sum[pos] += h.bitscore
                pident_sum[pos] += h.pident
                ppos_sum[pos] += h.ppos
                viruses_per_pos[pos].add(virus_id)

        # Segmentação em blocos contíguos de cobertura viral
        covered_bases = 0
        in_block = False
        block_start = 0
        block_scores = []
        block_depths = []
        block_pidents = []
        block_ppos = []
        block_viruses = set()

        for pos in range(1, seq_len + 1):
            d = len(viruses_per_pos[pos])
            if d > 0:
                covered_bases += 1
                bs = bitscore_sum[pos]
                # Fórmula do Consenso Posicional: amplifica regiões suportadas por múltiplos vírus independentes
                c_score = bs * (1.0 + 0.25 * max(0, d - 1))

                if not in_block:
                    in_block = True
                    block_start = pos
                    block_scores = [c_score]
                    block_depths = [d]
                    block_pidents = [pident_sum[pos] / d]
                    block_ppos = [ppos_sum[pos] / d]
                    block_viruses = set(viruses_per_pos[pos])
                else:
                    block_scores.append(c_score)
                    block_depths.append(d)
                    block_pidents.append(pident_sum[pos] / d)
                    block_ppos.append(ppos_sum[pos] / d)
                    block_viruses.update(viruses_per_pos[pos])
            else:
                if in_block:
                    in_block = False
                    block_end = pos - 1
                    intervals_detalhados.append({
                        "eve_id": eve_id,
                        "pos_inicio": block_start,
                        "pos_fim": block_end,
                        "comprimento_nt": block_end - block_start + 1,
                        "profundidade_viral_pico": max(block_depths),
                        "profundidade_viral_media": round(sum(block_depths) / len(block_depths), 1),
                        "score_consenso_pico": round(max(block_scores), 1),
                        "score_consenso_medio": round(sum(block_scores) / len(block_scores), 1),
                        "media_identidade": round(sum(block_pidents) / len(block_pidents), 1),
                        "media_positivos": round(sum(block_ppos) / len(block_ppos), 1),
                        "virus_presentes": "; ".join(sorted(block_viruses)[:8])
                    })

        if in_block:
            block_end = seq_len
            intervals_detalhados.append({
                "eve_id": eve_id,
                "pos_inicio": block_start,
                "pos_fim": block_end,
                "comprimento_nt": block_end - block_start + 1,
                "profundidade_viral_pico": max(block_depths),
                "profundidade_viral_media": round(sum(block_depths) / len(block_depths), 1),
                "score_consenso_pico": round(max(block_scores), 1),
                "score_consenso_medio": round(sum(block_scores) / len(block_scores), 1),
                "media_identidade": round(sum(block_pidents) / len(block_pidents), 1),
                "media_positivos": round(sum(block_ppos) / len(block_ppos), 1),
                "virus_presentes": "; ".join(sorted(block_viruses)[:8])
            })

        # Resumo por EVE
        max_d = max(depth_array) if depth_array else 0
        cov_pct = round((covered_bases / seq_len) * 100.0, 2) if seq_len > 0 else 0.0
        max_score = 0.0
        for pos in range(1, seq_len + 1):
            d = len(viruses_per_pos[pos])
            if d > 0:
                sc = bitscore_sum[pos] * (1.0 + 0.25 * max(0, d - 1))
                if sc > max_score:
                    max_score = sc

        consensus_summary_por_eve[eve_id] = {
            "eve_length": seq_len,
            "viral_covered_bases": covered_bases,
            "viral_coverage_pct": cov_pct,
            "max_viral_depth": max_d,
            "max_consensus_score": round(max_score, 1),
            "distinct_viruses_count": len(distinct_viruses_all),
            "distinct_viruses": sorted(list(distinct_viruses_all))[:15],
            "top_virus_id": top_hit.qseqid if top_hit.sseqid == eve_id else top_hit.sseqid,
            "top_virus_bitscore": top_hit.bitscore,
            "top_virus_ident": top_hit.pident,
            "top_virus_ppos": top_hit.ppos
        }

    return {
        "consensus_summary": consensus_summary_por_eve,
        "intervals": intervals_detalhados
    }


# =============================================================================
# 6. 2º BLASTX (BANCO GERAL - TOP 100 HITS) E CLASSIFICAÇÃO TAXONÔMICA
# =============================================================================

def analyze_general_blastx_hits(
    blastx_file: str,
    max_hits_per_query: int = 100,
    known_viral_accessions: Optional[Set[str]] = None,
    title_dict: Optional[Dict[str, str]] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Processa os hits do segundo BLASTX contra o banco geral (ex: nr).
    Limita aos 'max_hits_per_query' (padrão: 100) melhores hits por query,
    captura coordenadas na query (qstart, qend), métricas de similaridade e
    classifica cada hit como Viral ou Não-Viral (Persistência vs. Gene Celular).
    """
    raw_hsps = parse_alignment_file(blastx_file)
    hits_por_query = defaultdict(list)

    for h in raw_hsps:
        hits_por_query[h.qseqid].append(h)

    annotated_por_query = {}
    for qid, lista in hits_por_query.items():
        # Ordena por bitscore decrescente e seleciona até top 100
        lista.sort(key=lambda x: x.bitscore, reverse=True)
        top_hits = lista[:max_hits_per_query]

        annotated_list = []
        for h in top_hits:
            tax_info = classify_hit_taxonomy(
                h.stitle, h.sseqid,
                known_viral_accessions=known_viral_accessions,
                title_dict=title_dict
            )
            annotated_list.append({
                "hsp": h,
                "sseqid": h.sseqid,
                "stitle": h.stitle,
                "bitscore": h.bitscore,
                "evalue": h.evalue,
                "pident": h.pident,
                "ppos": h.ppos,
                "length": h.length,
                "qcovhsp": h.qcovhsp,
                "q_start": h.q_start,
                "q_end": h.q_end,
                "is_viral": tax_info["is_viral"],
                "organism": tax_info["organism"],
                "category": tax_info["category"],
                "category_label": tax_info["label"],
                "category_desc": tax_info["description"]
            })
        annotated_por_query[qid] = annotated_list

    return annotated_por_query


# =============================================================================
# 7. COMPARAÇÃO POSICIONAL MINUCIOSA (VIRAL VS. NÃO-VIRAL)
# =============================================================================

def compare_viral_and_nonviral_positions(
    consensus_data: Dict[str, Any],
    general_hits_annotated: Dict[str, List[Dict[str, Any]]],
    eve_lengths: Optional[Dict[str, int]] = None,
    hmm_hits: Optional[Dict[str, List[Dict[str, Any]]]] = None,
    blastn_hits: Optional[Dict[str, Dict[str, Any]]] = None,
    cluster_map: Optional[Dict[str, Dict[str, Any]]] = None,
    rvdb_anno: Optional[Dict[str, str]] = None,
    title_dict: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """
    Compara minuciosamente, intervalo por intervalo, as coordenadas alinhadas virais
    com as posições alinhadas celulares/não-virais obtidas no banco geral (top 100 hits).

    Determina:
      - Métricas quantitativas de diferença de alinhamento quando há sobreposição com genes celulares:
          * delta_bitscore_viral_vs_celular (Bitscore Viral - Bitscore Celular)
          * delta_cobertura_viral_vs_celular (% Cobertura Viral - % Cobertura Celular)
          * delta_identidade_viral_vs_celular (% Identidade Viral - % Identidade Celular)
          * razao_bitscore_viral_vs_celular (Razão Bitscore Viral / Bitscore Celular)
          * status_sobreposicao_celular (VIRAL_EXCLUSIVO, VIRAL_DOMINANTE, HOMOLOGIA_EQUILIBRADA, CELULAR_DOMINANTE, POSSIVEL_EVE_INCONCLUSIVA_HGT)
      - Filtro e segregação estrita de NIRVs genuínos vs. elementos excluídos:
          * Domínio HAM / HAM1 (pirofosfatase CBSV / FAM003967 / PF01725).
          * Retroviridae / Retrovirais / ERVs.
          * Genes de Transcriptase Reversa (RT / RVT / retrotransposons / Ty1-Copia / Ty3-Gypsy / Bel-Pao / LINE).
          * Transferência Horizontal Pura (HGT) onde a região apresenta apenas genes do hospedeiro sem corroboração viral.
          * Elementos removidos são salvos no relatório de auditoria 'validacao_posicional_removidas_retro_ham_rt.tsv'.
      - Avaliação de qualidade e estratificação em Faixas de Realidade Viral (Tiers 1 a 4).
    """
    consensus_summary = consensus_data["consensus_summary"]
    viral_intervals = consensus_data["intervals"]

    intervals_by_eve = defaultdict(list)
    for inv in viral_intervals:
        intervals_by_eve[inv["eve_id"]].append(inv)

    # Mapeamento reverso de proteínas virais para detecção de fragmentação
    all_virus_to_eves = defaultdict(list)
    for eid, cs in consensus_summary.items():
        vid = cs.get("top_virus_id")
        if vid and vid != "Nenhum":
            all_virus_to_eves[vid].append(eid)

    detalhes_comparacao = []
    sumario_eves = []
    removidas_retro_ham_rt = []

    for eve_id, c_sum in consensus_summary.items():
        gen_hits = general_hits_annotated.get(eve_id, [])
        v_invs = intervals_by_eve.get(eve_id, [])

        top_v_id = c_sum.get("top_virus_id", "Nenhum")
        top_v_title = title_dict.get(top_v_id, "") if title_dict else ""

        # Domínio HMM da EVE
        h_hits = hmm_hits.get(eve_id, []) if hmm_hits else []
        if not h_hits and hmm_hits:
            for k, vlist in hmm_hits.items():
                if k == eve_id or (eve_id in k) or (k in eve_id):
                    h_hits = vlist
                    break
        hmm_dom = h_hits[0]["domain"] if h_hits else "-"
        hmm_kw = rvdb_anno.get(hmm_dom, "") if (rvdb_anno and hmm_dom != "-") else ""
        candidate_all_doms = {h["domain"] for h in h_hits} if h_hits else set()
        if hmm_dom != "-":
            candidate_all_doms.add(hmm_dom)

        # Separa hits gerais
        nonviral_hits = [gh for gh in gen_hits if not gh["is_viral"]]
        viral_gen_hits = [gh for gh in gen_hits if gh["is_viral"]]

        # Rastreia espécies onde o elemento persiste e identifica genes celulares
        especies_persistencia = set()
        genes_celulares_detectados = []
        melhor_nv_hit = None
        melhor_v_hit = None

        if nonviral_hits:
            melhor_nv_hit = max(nonviral_hits, key=lambda x: x["bitscore"])
            for nv in nonviral_hits:
                if nv["category"] == "PERSISTENCIA_OUTRAS_ESPECIES" and nv["organism"] != "Desconhecido":
                    especies_persistencia.add(nv["organism"])
                elif nv["category"] == "SOBREPOSICAO_HOSPEDEIRO" or _CELLULAR_REGEX.search(nv.get("stitle", "")):
                    genes_celulares_detectados.append(nv)

        if viral_gen_hits:
            melhor_v_hit = max(viral_gen_hits, key=lambda x: x["bitscore"])

        # Identifica melhor hit celular que se sobreponha às coordenadas virais da query
        melhor_celular_sobreposto = None
        melhor_celular_bs = 0.0
        melhor_celular_cov = 0.0
        melhor_celular_ident = 0.0
        melhor_celular_ppos = 0.0
        melhor_celular_org = "Nenhum"
        melhor_celular_desc = "-"
        melhor_celular_faixa = "-"

        if genes_celulares_detectados and v_invs:
            overlapping_cel = []
            for inv in v_invs:
                v_st = inv["pos_inicio"]
                v_ed = inv["pos_fim"]
                for ch in genes_celulares_detectados:
                    c_st = ch["q_start"]
                    c_ed = ch["q_end"]
                    ov_st = max(v_st, c_st)
                    ov_ed = min(v_ed, c_ed)
                    if ov_ed >= ov_st:
                        ov_len = ov_ed - ov_st + 1
                        overlapping_cel.append((ov_len, ch))
            if overlapping_cel:
                overlapping_cel.sort(key=lambda x: (x[1]["bitscore"], x[0]), reverse=True)
                top_cel = overlapping_cel[0][1]
                melhor_celular_sobreposto = top_cel
                melhor_celular_bs = float(top_cel["bitscore"])
                melhor_celular_cov = float(top_cel["qcovhsp"])
                melhor_celular_ident = float(top_cel["pident"])
                melhor_celular_ppos = float(top_cel["ppos"])
                melhor_celular_org = top_cel["organism"]
                melhor_celular_desc = top_cel["stitle"][:80]
                melhor_celular_faixa = f"{top_cel['q_start']}-{top_cel['q_end']}"

        # Métricas Virais de Referência
        v_score_ref = float(c_sum.get("max_consensus_score", 0.0) or c_sum.get("top_virus_bitscore", 0.0) or 0.0)
        v_cov_ref = float(c_sum.get("viral_coverage_pct", 0.0) or 0.0)
        v_ident_ref = float(c_sum.get("top_virus_ident", 0.0) or 0.0)
        v_ppos_ref = float(c_sum.get("top_virus_ppos", 0.0) or 0.0)

        # Cálculo das Métricas Quantitativas de Comparação Posicional Viral vs. Celular
        tem_sobreposicao_celular = "SIM" if (melhor_celular_sobreposto is not None) else "NAO"

        if tem_sobreposicao_celular == "SIM":
            delta_bs_v_vs_c = round(v_score_ref - melhor_celular_bs, 1)
            delta_cov_v_vs_c = round(v_cov_ref - melhor_celular_cov, 1)
            delta_ident_v_vs_c = round(v_ident_ref - melhor_celular_ident, 1)
            delta_ppos_v_vs_c = round(v_ppos_ref - melhor_celular_ppos, 1)
            razao_bs_v_vs_c = round(v_score_ref / max(1.0, melhor_celular_bs), 2)
            if delta_bs_v_vs_c >= 15.0 and delta_cov_v_vs_c >= 0:
                status_sobreposicao = "VIRAL_DOMINANTE"
            elif abs(delta_bs_v_vs_c) < 15.0:
                status_sobreposicao = "HOMOLOGIA_EQUILIBRADA"
            else:
                status_sobreposicao = "CELULAR_DOMINANTE"
        else:
            delta_bs_v_vs_c = round(v_score_ref, 1)
            delta_cov_v_vs_c = round(v_cov_ref, 1)
            delta_ident_v_vs_c = round(v_ident_ref, 1)
            delta_ppos_v_vs_c = round(v_ppos_ref, 1)
            razao_bs_v_vs_c = round(v_score_ref, 2)
            status_sobreposicao = "VIRAL_EXCLUSIVO"

        # Checagem Estrita de Exclusão (HAM, Retroviridae, Transcriptase Reversa, HGT puro)
        deve_excluir, tipo_excl, motivo_excl = check_candidate_exclusion(
            eve_id=eve_id,
            top_virus_id=top_v_id,
            top_virus_title=top_v_title,
            hmm_domain=hmm_dom,
            hmm_keywords=hmm_kw,
            general_hits=gen_hits,
            consensus_summary=c_sum,
            all_domains=candidate_all_doms
        )

        if deve_excluir:
            doms_str = ", ".join(sorted(candidate_all_doms)) if candidate_all_doms else (hmm_dom if hmm_dom != "-" else "-")
            removidas_retro_ham_rt.append({
                "eve_id": eve_id,
                "eve_comprimento_nt": c_sum["eve_length"],
                "tipo_filtro": tipo_excl,
                "motivo_exclusao": motivo_excl,
                "top_virus_id": top_v_id,
                "top_virus_titulo": top_v_title[:80] if top_v_title else "-",
                "top_virus_bitscore": v_score_ref,
                "top_virus_ident": v_ident_ref,
                "top_virus_positivos": v_ppos_ref,
                "cobertura_viral_pct": v_cov_ref,
                "tag_dominio_conservado": f"DOMINIO_HMM_CONFIRMADO ({doms_str})" if doms_str != "-" else "SEM_DOMINIO_HMM",
                "hmm_keywords": hmm_kw[:80] if hmm_kw else "-",
                "melhor_hit_celular_organismo": melhor_celular_org,
                "melhor_hit_celular_titulo": melhor_celular_desc,
                "melhor_hit_celular_bitscore": melhor_celular_bs,
                "delta_bitscore_viral_vs_celular": delta_bs_v_vs_c,
                "status_final": f"REMOVIDO_{tipo_excl}"
            })
            continue

        # Comparação detalhada por intervalo viral
        tem_alerta_falso_positivo = (status_sobreposicao == "CELULAR_DOMINANTE")
        tem_persistencia_multiespecies = bool(especies_persistencia)

        for inv in v_invs:
            v_start = inv["pos_inicio"]
            v_end = inv["pos_fim"]
            v_score = inv["score_consenso_pico"]
            v_ident = inv["media_identidade"]
            v_ppos = inv["media_positivos"]

            overlapping_nv = []
            for nv in nonviral_hits:
                nv_start = nv["q_start"]
                nv_end = nv["q_end"]
                ov_start = max(v_start, nv_start)
                ov_end = min(v_end, nv_end)
                if ov_end >= ov_start:
                    ov_len = ov_end - ov_start + 1
                    overlapping_nv.append((ov_len, nv))

            if overlapping_nv:
                overlapping_nv.sort(key=lambda item: (item[1]["bitscore"], item[0]), reverse=True)
                _, top_nv = overlapping_nv[0]

                delta_bs_v_c = round(v_score - top_nv["bitscore"], 1)
                delta_cov_v_c = round(v_score - top_nv["qcovhsp"], 1)
                delta_ident_v_c = round(v_ident - top_nv["pident"], 1)
                razao_bs_v_c = round(v_score / max(1.0, top_nv["bitscore"]), 2)

                if top_nv["category"] == "SOBREPOSICAO_HOSPEDEIRO":
                    if delta_bs_v_c >= 15.0:
                        veredito_inv = "VIRAL_DOMINANTE_SOBRE_CELULAR"
                    elif delta_bs_v_c <= -15.0:
                        veredito_inv = "CELULAR_DOMINANTE_SOBRE_VIRAL"
                    else:
                        veredito_inv = "HOMOLOGIA_DUPLA_EQUILIBRADA"
                elif top_nv["category"] == "PERSISTENCIA_OUTRAS_ESPECIES":
                    veredito_inv = "PERSISTENCIA_OUTRAS_ESPECIES_HOMOLOGIA"
                else:
                    veredito_inv = "HOMOLOGIA_NAO_VIRAL_MODERADA"

                detalhes_comparacao.append({
                    "eve_id": eve_id,
                    "intervalo_inicio": v_start,
                    "intervalo_fim": v_end,
                    "comprimento_nt": inv["comprimento_nt"],
                    "tipo_regiao": "SOBREPOSICAO_VIRAL_E_NAO_VIRAL",
                    "viral_score_consenso": v_score,
                    "viral_profundidade_num_virus": inv["profundidade_viral_pico"],
                    "viral_virus_principais": inv["virus_presentes"],
                    "viral_media_identidade": v_ident,
                    "viral_media_positivos": v_ppos,
                    "celular_sobreposto_organismo": top_nv["organism"],
                    "celular_sobreposto_categoria": top_nv["category"],
                    "celular_sobreposto_titulo": top_nv["stitle"][:80],
                    "celular_sobreposto_identidade": top_nv["pident"],
                    "celular_sobreposto_positivos": top_nv["ppos"],
                    "celular_sobreposto_cobertura": top_nv["qcovhsp"],
                    "celular_sobreposto_bitscore": top_nv["bitscore"],
                    "delta_bitscore_viral_vs_celular": delta_bs_v_c,
                    "delta_cobertura_viral_vs_celular": delta_cov_v_c,
                    "delta_identidade_viral_vs_celular": delta_ident_v_c,
                    "razao_bitscore_viral_vs_celular": razao_bs_v_c,
                    "status_intervalo": veredito_inv
                })
            else:
                detalhes_comparacao.append({
                    "eve_id": eve_id,
                    "intervalo_inicio": v_start,
                    "intervalo_fim": v_end,
                    "comprimento_nt": inv["comprimento_nt"],
                    "tipo_regiao": "VIRAL_EXCLUSIVO",
                    "viral_score_consenso": v_score,
                    "viral_profundidade_num_virus": inv["profundidade_viral_pico"],
                    "viral_virus_principais": inv["virus_presentes"],
                    "viral_media_identidade": v_ident,
                    "viral_media_positivos": v_ppos,
                    "celular_sobreposto_organismo": "Nenhum",
                    "celular_sobreposto_categoria": "NENHUM_HIT",
                    "celular_sobreposto_titulo": "-",
                    "celular_sobreposto_identidade": 0.0,
                    "celular_sobreposto_positivos": 0.0,
                    "celular_sobreposto_cobertura": 0.0,
                    "celular_sobreposto_bitscore": 0.0,
                    "delta_bitscore_viral_vs_celular": v_score,
                    "delta_cobertura_viral_vs_celular": 100.0,
                    "delta_identidade_viral_vs_celular": v_ident,
                    "razao_bitscore_viral_vs_celular": v_score,
                    "status_intervalo": "VIRAL_ESPECIFICO_EXCLUSIVO"
                })

        # Veredito global da EVE candidata
        if tem_alerta_falso_positivo and not (c_sum["viral_coverage_pct"] > 60 and c_sum["max_viral_depth"] >= 3):
            veredito_global = "ALERTA_GENE_CELULAR_HOSPEDEIRO"
            grau_confianca = "Baixo (Suspeita de Falso Positivo Celular)"
        elif tem_persistencia_multiespecies and c_sum["viral_coverage_pct"] >= 20:
            veredito_global = "EVE_CONFIRMADA_PERSISTENTE_MULTIESPECIES"
            grau_confianca = "Alto (EVE Ortóloga Conservada em Múltiplas Espécies)"
        elif c_sum["viral_coverage_pct"] >= 30 and c_sum["max_consensus_score"] >= 40:
            veredito_global = "EVE_CONFIRMADA_VIRAL_ESPECIFICA"
            grau_confianca = "Alto (Específica Viral)"
        elif c_sum["viral_coverage_pct"] > 0:
            veredito_global = "EVE_CANDIDATA_PARCIAL_EM_VALIDACAO"
            grau_confianca = "Moderado"
        else:
            veredito_global = "SEM_HOMOLOGIA_VIRAL_DETECTADA"
            grau_confianca = "Indeterminado"

        # Avaliação de qualidade e faixas
        quality_eval = evaluate_viral_quality_tags(
            eve_id=eve_id,
            consensus_summary=consensus_summary,
            general_hits=gen_hits,
            hmm_hits=hmm_hits or {},
            blastn_hits=blastn_hits or {},
            cluster_map=cluster_map or {},
            all_virus_to_eves=all_virus_to_eves
        )

        sumario_eves.append({
            "eve_id": eve_id,
            "eve_comprimento_nt": c_sum["eve_length"],
            "cobertura_viral_nt": c_sum["viral_covered_bases"],
            "cobertura_viral_pct": c_sum["viral_coverage_pct"],
            "score_consenso_viral_max": c_sum["max_consensus_score"],
            "profundidade_viral_max": c_sum["max_viral_depth"],
            "num_virus_alinhados": c_sum["distinct_viruses_count"],
            "top_virus_id": top_v_id,
            "top_virus_ident": c_sum["top_virus_ident"],
            "top_virus_positivos": c_sum["top_virus_ppos"],
            "top_virus_bitscore": c_sum["top_virus_bitscore"],
            "num_hits_geral_top100": len(gen_hits),
            "num_hits_nao_virais": len(nonviral_hits),
            "persistencia_outras_especies": "SIM" if especies_persistencia else "NAO",
            "lista_outras_especies_persistencia": "; ".join(sorted(especies_persistencia)[:10]) if especies_persistencia else "-",
            "sobreposicao_com_celular": tem_sobreposicao_celular,
            "delta_bitscore_viral_vs_celular": delta_bs_v_vs_c,
            "delta_cobertura_viral_vs_celular": delta_cov_v_vs_c,
            "delta_identidade_viral_vs_celular": delta_ident_v_vs_c,
            "delta_positivos_viral_vs_celular": delta_ppos_v_vs_c,
            "razao_bitscore_viral_vs_celular": razao_bs_v_vs_c,
            "celular_sobreposto_organismo": melhor_celular_org,
            "celular_sobreposto_titulo": melhor_celular_desc,
            "celular_sobreposto_faixa": melhor_celular_faixa,
            "celular_sobreposto_bitscore": melhor_celular_bs,
            "celular_sobreposto_cobertura": melhor_celular_cov,
            "status_sobreposicao_celular": status_sobreposicao,
            "veredito_posicional": veredito_global,
            "grau_confianca": grau_confianca,
            # Tags biológicas de alinhamento e realidade viral
            "faixa_realidade_viral": quality_eval["faixa_realidade_viral"],
            "faixa_rotulo": quality_eval["faixa_rotulo"],
            "score_realidade_viral": quality_eval["score_realidade_viral"],
            "tag_ranking_viral": quality_eval["tag_ranking_viral"],
            "tag_dominio_conservado": quality_eval["tag_dominio_conservado"],
            "tag_blastn": quality_eval["tag_blastn"],
            "tag_divisao_proteina": quality_eval["tag_divisao_proteina"],
            "tag_deriva_positivos": quality_eval["tag_deriva_positivos"],
            "razao_positivos_identidade": quality_eval["razao_positivos_identidade"],
            "delta_ident_v_vs_nv": quality_eval["delta_ident_v_vs_nv"],
            "delta_ppos_v_vs_nv": quality_eval["delta_ppos_v_vs_nv"],
            "delta_bitscore_v_vs_nv": quality_eval["delta_bitscore_v_vs_nv"],
            "tags_resumo": quality_eval["tags_resumo"]
        })

    return {
        "sumario_eves": sumario_eves,
        "detalhes_comparacao": detalhes_comparacao,
        "removidas_retro_ham_rt": removidas_retro_ham_rt
    }


# =============================================================================
# 8. EXPORTAÇÃO DE RELATÓRIOS E JSON ESTRUTURADO
# =============================================================================

def export_positional_validation_reports(
    comparison_results: Dict[str, Any],
    consensus_data: Dict[str, Any],
    diretorio_saida: str,
    prefixo: str = "validacao_posicional"
) -> Dict[str, str]:
    """
    Exporta relatórios finais em TSV e JSON estruturado para o dashboard web:
      - validacao_posicional_eves_sumario.tsv (Candidatas Curadas NIRVs)
      - validacao_posicional_removidas_retro_ham_rt.tsv (Elementos Removidos por Filtro Retro/HAM/RT/HGT)
      - validacao_posicional_intervalos_detalhados.tsv
      - perfil_cobertura_consenso_posicional.tsv
      - validacao_posicional.json
    """
    os.makedirs(diretorio_saida, exist_ok=True)

    sumario_tsv = os.path.join(diretorio_saida, f"{prefixo}_eves_sumario.tsv")
    removidas_tsv = os.path.join(diretorio_saida, f"{prefixo}_removidas_retro_ham_rt.tsv")
    intervalos_tsv = os.path.join(diretorio_saida, f"{prefixo}_intervalos_detalhados.tsv")
    consenso_tsv = os.path.join(diretorio_saida, "perfil_cobertura_consenso_posicional.tsv")
    json_path = os.path.join(diretorio_saida, f"{prefixo}.json")

    # 1. Exporta sumário de EVEs curadas (NIRVs genuínas)
    sumario_eves = comparison_results["sumario_eves"]
    if sumario_eves:
        keys = list(sumario_eves[0].keys())
        with open(sumario_tsv, 'w', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            writer.writerows(sumario_eves)

    # 2. Exporta elementos excluídos pelo filtro Retro/HAM/RT/HGT
    removidas = comparison_results.get("removidas_retro_ham_rt", [])
    if removidas:
        keys_rem = list(removidas[0].keys())
        with open(removidas_tsv, 'w', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys_rem, delimiter='\t')
            writer.writeheader()
            writer.writerows(removidas)

    # 3. Exporta intervalos detalhados
    detalhes = comparison_results["detalhes_comparacao"]
    if detalhes:
        keys_det = list(detalhes[0].keys())
        with open(intervalos_tsv, 'w', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys_det, delimiter='\t')
            writer.writeheader()
            writer.writerows(detalhes)

    # 4. Exporta perfil de consenso posicional
    consenso_intervals = consensus_data.get("intervals", [])
    if consenso_intervals:
        keys_con = list(consenso_intervals[0].keys())
        with open(consenso_tsv, 'w', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys_con, delimiter='\t')
            writer.writeheader()
            writer.writerows(consenso_intervals)

    # 5. Exporta JSON estruturado para o Dashboard Web
    total_eves = len(sumario_eves)
    total_removidas = len(removidas)
    total_brutas = total_eves + total_removidas

    persistentes = sum(1 for e in sumario_eves if e["persistencia_outras_especies"] == "SIM")
    virais_puras = sum(1 for e in sumario_eves if e["veredito_posicional"] == "EVE_CONFIRMADA_VIRAL_ESPECIFICA")
    sobreposicoes_celulares = sum(1 for e in sumario_eves if e.get("sobreposicao_com_celular") == "SIM")

    # Contagem de exclusões por categoria
    removidas_por_tipo = defaultdict(int)
    for r in removidas:
        removidas_por_tipo[r.get("tipo_filtro", "OUTRO")] += 1

    # Contagem de status de sobreposição celular
    status_sobreposicao_counts = defaultdict(int)
    for e in sumario_eves:
        status_sobreposicao_counts[e.get("status_sobreposicao_celular", "VIRAL_EXCLUSIVO")] += 1

    # Coleta espécies hospedeiras externas identificadas
    especies_counts = defaultdict(int)
    for e in sumario_eves:
        sp_str = e.get("lista_outras_especies_persistencia", "-")
        if sp_str and sp_str != "-":
            for sp in sp_str.split("; "):
                sp_clean = sp.strip()
                if sp_clean:
                    especies_counts[sp_clean] += 1

    ranking_especies = [
        {"especie": sp, "num_eves_compartilhadas": count}
        for sp, count in sorted(especies_counts.items(), key=lambda x: x[1], reverse=True)
    ]

    # Distribuição das Faixas de Realidade Viral
    faixa_counts = defaultdict(int)
    for e in sumario_eves:
        fx = e.get("faixa_realidade_viral", "Faixa 4 - Baixa Certeza / Sobreposição Celular")
        faixa_counts[fx] += 1

    resumo_faixas = {
        "faixa_1_alta_certeza": {
            "contagem": faixa_counts["Faixa 1 - Alta Certeza"],
            "pct": round(faixa_counts["Faixa 1 - Alta Certeza"] / total_eves * 100.0, 1) if total_eves > 0 else 0.0,
            "descricao": "Alta certeza biológica com BestHit viral ou domínio conservado sem interferência celular"
        },
        "faixa_2_ortologa_persistente": {
            "contagem": faixa_counts["Faixa 2 - Ortóloga Persistente"],
            "pct": round(faixa_counts["Faixa 2 - Ortóloga Persistente"] / total_eves * 100.0, 1) if total_eves > 0 else 0.0,
            "descricao": "EVEs ortólogas com conservação inter-espécies e assinatura de substituição conservativa (Rpos >= 1.3)"
        },
        "faixa_3_ancestral_divergente": {
            "contagem": faixa_counts["Faixa 3 - Ancestral Divergente"],
            "pct": round(faixa_counts["Faixa 3 - Ancestral Divergente"] / total_eves * 100.0, 1) if total_eves > 0 else 0.0,
            "descricao": "Fósseis virais antigos erodidos, proteínas fragmentadas no cluster ou domínio HMM com divergência nt"
        },
        "faixa_4_alerta_incerteza": {
            "contagem": faixa_counts["Faixa 4 - Baixa Certeza / Sobreposição Celular"],
            "pct": round(faixa_counts["Faixa 4 - Baixa Certeza / Sobreposição Celular"] / total_eves * 100.0, 1) if total_eves > 0 else 0.0,
            "descricao": "Sobreposição com gene celular do hospedeiro dominante ou homologia viral incerta"
        }
    }

    # Distribuição das Tags de Evidência
    tag_counts = {
        "ranking_viral": defaultdict(int),
        "dominio_hmm": defaultdict(int),
        "blastn": defaultdict(int),
        "divisao_proteina": defaultdict(int),
        "deriva_positivos": defaultdict(int)
    }
    for e in sumario_eves:
        tag_counts["ranking_viral"][e.get("tag_ranking_viral", "OUTRO")] += 1
        tag_counts["dominio_hmm"]["DOMINIO_HMM_CONFIRMADO" if "CONFIRMADO" in e.get("tag_dominio_conservado", "") else "SEM_DOMINIO_HMM"] += 1
        tag_counts["blastn"]["BLASTN_VIRAL_CONFIRMADO" if "CONFIRMADO" in e.get("tag_blastn", "") else "BLASTN_DIVERGENTE_ANCESTRAL"] += 1
        tag_counts["divisao_proteina"]["PROTEINA_FRAGMENTADA" if "FRAGMENTADA" in e.get("tag_divisao_proteina", "") else ("PROTEINA_INTEGRA" if "INTEGRA" in e.get("tag_divisao_proteina", "") else "SEM_MAPEAMENTO")] += 1
        tag_counts["deriva_positivos"][e.get("tag_deriva_positivos", "OUTRO")] += 1

    distribuicao_tags = {k: dict(v) for k, v in tag_counts.items()}

    json_payload = {
        "resumo_geral": {
            "total_eves_analisadas_brutas": total_brutas,
            "total_eves_curadas_nirv": total_eves,
            "total_removidas_retro_ham_rt": total_removidas,
            "removidas_por_tipo": dict(removidas_por_tipo),
            "eves_com_persistencia_multiespecies": persistentes,
            "eves_confirmadas_virais_especificas": virais_puras,
            "eves_com_sobreposicao_celular": sobreposicoes_celulares,
            "pct_persistencia_evolutiva": round((persistentes / total_eves * 100.0), 1) if total_eves > 0 else 0.0
        },
        "resumo_faixas_realidade": resumo_faixas,
        "distribuicao_tags": distribuicao_tags,
        "distribuicao_status_sobreposicao": dict(status_sobreposicao_counts),
        "ranking_especies_persistencia": ranking_especies[:20],
        "sumario_eves": sumario_eves,
        "elementos_removidos_retro_ham_rt": removidas,
        "intervalos_detalhados": detalhes[:2000]
    }

    import json
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(json_payload, f, indent=2, ensure_ascii=False)

    print(f"[run] Relatórios de validação posicional gerados com sucesso:")
    print(f"       * Sumário EVEs Curadas (NIRVs)    : {sumario_tsv} ({total_eves} candidatos)")
    print(f"       * Elementos Excluídos Retro/HAM/RT: {removidas_tsv} ({total_removidas} removidos)")
    print(f"       * Intervalos Detalhados (Deltas)  : {intervalos_tsv}")
    print(f"       * Consenso Posicional             : {consenso_tsv}")
    print(f"       * JSON Dashboard                  : {json_path}")

    return {
        "sumario_tsv": sumario_tsv,
        "removidas_tsv": removidas_tsv,
        "intervalos_tsv": intervalos_tsv,
        "consenso_tsv": consenso_tsv,
        "json_path": json_path
    }


# =============================================================================
# 9. FLUXO CONSOLIDADO (PIPELINE CALL)
# =============================================================================

def executar_analise_posicional_completa(
    fasta_candidatas_nt: str,
    fasta_virais_ref: str,
    diretorio_saida: str,
    blastx_geral_tsv_ou_xml: Optional[str] = None,
    diamond_db_geral: Optional[str] = None,
    hmm_file: Optional[str] = None,
    blastn_file: Optional[str] = None,
    cluster_tsv: Optional[str] = None,
    threads: int = 4,
    evalue: float = 1e-3,
    hits_por_query: int = 100,
    export_xml: bool = False
) -> Dict[str, Any]:
    """
    Executa o fluxo completo solicitado pelo usuário:
    1. Executa ou lê o tBLASTn recíproco de proteínas virais contra o segmento EVE.
    2. Anota posições não-sobrepostas via Weighted Interval Scheduling.
    3. Calcula o Viral Positional Consensus Score para cada coordenada/intervalo.
    4. Executa ou lê o 2º BLASTX contra banco geral coletando os top 100 hits por query
       com coordenadas qstart, qend, identidade, positivos e cobertura.
    5. Classifica hits do 2º BLASTX em Virais e Não-Virais (Persistência vs. Gene Hospedeiro).
    6. Compara os locais virais e não-virais na query, calculando deltas quantitativos (bitscore e cobertura).
    7. Aplica filtro de remoção estrita para Domínio HAM, Retroviridae, Transcriptase Reversa e HGT celular.
    8. Avalia tags de qualidade (BestHit, Top2 sobreposto, Top3+, Domínio HMM, BLASTN, Deriva BLOSUM62 e Fragmentação).
    9. Exporta sumários tabulares, relatório de exclusões auditadas e JSON estruturado para o Dashboard Web.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    if os.path.basename(os.path.normpath(diretorio_saida)) == "validacao_posicional":
        pasta_posicional = diretorio_saida
    else:
        pasta_posicional = os.path.join(diretorio_saida, "validacao_posicional")
    os.makedirs(pasta_posicional, exist_ok=True)

    print(f"\n{'='*75}\n[Validação Posicional] Análise Posicional e Evolutiva das EVEs\n{'='*75}")
    # 1. Carregamos as HSPs virais diretamente do DIAMOND inicial, ignorando o tBLASTn recíproco
    # Procuramos o TSV de triagem viral
    import glob
    parent_dir = os.path.dirname(os.path.abspath(diretorio_saida))
    possiveis_tsv_viral = [
        os.path.join(parent_dir, "02_alinhamento_diamond", "triagem_viral_diamond.tsv"),
        os.path.join(parent_dir, "triagem_viral_diamond.tsv")
    ]
    viral_tsv = None
    for p in possiveis_tsv_viral:
        if os.path.exists(p):
            viral_tsv = p
            break
            
    all_viral_hsps = []
    if viral_tsv:
        print(f"[run] Carregando HSPs virais do arquivo: {viral_tsv}")
        all_viral_hsps = parse_alignment_tabular(viral_tsv)
    else:
        print("[aviso] triagem_viral_diamond.tsv não encontrado, o consenso posicional ficará vazio.")

    res_reciproco = {"tsv_nao_sobreposto": None, "xml": None, "all_hsps": all_viral_hsps}



    # 2. Viral Positional Consensus Score
    consensus_data = compute_viral_positional_consensus(
        fasta_candidatas_nt=fasta_candidatas_nt,
        viral_hsps=res_reciproco["all_hsps"]
    )

    parent_dir = os.path.dirname(os.path.abspath(diretorio_saida))

    # 3. 2º BLASTX contra Banco Geral (ou parsing de arquivo existente / autodetectado)
    if not blastx_geral_tsv_ou_xml or not os.path.exists(blastx_geral_tsv_ou_xml):
        candidatos_blastx = [
            os.path.join(pasta_posicional, "diamond_blastx_referencia_geral.tsv"),
            os.path.join(diretorio_saida, "diamond_blastx_referencia_geral.tsv"),
            os.path.join(parent_dir, "02_alinhamento_diamond", "diamond_blastx_referencia_geral.tsv"),
            os.path.join(parent_dir, "04_candidatas_eves", "diamond_blastx_referencia_geral.tsv"),
            os.path.join(parent_dir, "diamond_blastx_referencia_geral.tsv"),
            "resultados_organizados/02_alinhamento_diamond/diamond_blastx_referencia_geral.tsv",
            "teste_app/diamond_blastx_referencia_geral.tsv"
        ]
        for b_cand in candidatos_blastx:
            if b_cand and os.path.exists(b_cand) and os.path.getsize(b_cand) > 0:
                blastx_geral_tsv_ou_xml = b_cand
                print(f"[info] Arquivo do 2º BLASTX autodetectado: {b_cand}")
                break

    if not blastx_geral_tsv_ou_xml or not os.path.exists(blastx_geral_tsv_ou_xml):
        if diamond_db_geral and os.path.exists(diamond_db_geral):
            from evetool.discovery.sequence_searcher import buscar_proteinas_referencia_diamond
            print(f"[run] Executando 2º DIAMOND BLASTX estendido contra o banco geral {diamond_db_geral}...")
            blastx_geral_tsv_ou_xml = buscar_proteinas_referencia_diamond(
                fasta_query=fasta_candidatas_nt,
                banco_geral=diamond_db_geral,
                diretorio_saida=pasta_posicional,
                threads=threads,
                hits_por_query=hits_por_query,
                evalue=evalue
            )
        else:
            print("[aviso] Nenhum banco geral ou resultado de 2º BLASTX fornecido. "
                  "Validação posicional será computada com dados virais puros.")
            blastx_geral_tsv_ou_xml = None

    # Extrai conjunto de acessos virais conhecidos a partir do fasta de referências
    known_viral_accessions = set()
    if fasta_virais_ref and os.path.exists(fasta_virais_ref):
        for rec in SeqIO.parse(fasta_virais_ref, "fasta"):
            clean_id = rec.id.strip().split()[0]
            known_viral_accessions.add(clean_id)

    # Constrói mapa de títulos estendido varrendo arquivos no ambiente
    search_paths = [
        fasta_virais_ref,
        os.path.join(diretorio_saida, "triagem_viral_diamond_refseq_filtrado.tsv"),
        os.path.join(parent_dir, "triagem_viral_diamond_refseq_filtrado.tsv"),
        os.path.join(parent_dir, "02_alinhamento_diamond", "triagem_viral_diamond_refseq_filtrado.tsv"),
        os.path.join(diretorio_saida, "triagem_viral_diamond.tsv"),
        os.path.join(parent_dir, "triagem_viral_diamond.tsv"),
        os.path.join(parent_dir, "02_alinhamento_diamond", "triagem_viral_diamond.tsv"),
        os.path.join(parent_dir, "04_candidatas_eves", "proteinas_referencia_diamond_top1758.fasta"),
    ]
    title_dict = build_stitle_lookup_dict(search_paths)

    general_annotated = {}
    if blastx_geral_tsv_ou_xml and os.path.exists(blastx_geral_tsv_ou_xml):
        general_annotated = analyze_general_blastx_hits(
            blastx_file=blastx_geral_tsv_ou_xml,
            max_hits_per_query=hits_por_query,
            known_viral_accessions=known_viral_accessions,
            title_dict=title_dict
        )

    # Descoberta automática de HMM, BLASTN e Clusters
    hmm_search = [
        hmm_file,
        os.path.join(diretorio_saida, "triagem_profunda_hmm.tsv"),
        os.path.join(parent_dir, "03_triagem_hmm", "triagem_profunda_hmm.tsv"),
        os.path.join(parent_dir, "triagem_profunda_hmm.tsv"),
        "resultados_organizados/03_triagem_hmm/triagem_profunda_hmm.tsv"
    ]
    hmm_hits = load_hmm_domain_hits(hmm_search)

    # Carrega anotações do RVDB para checagem estrita de palavras-chave
    rvdb_search = [
        config.hmm_path("RVDB", "rvdb_base_annotations.tsv"),
        os.path.join(parent_dir, "rvdb_base_annotations.tsv"),
        "rvdb_base_annotations.tsv"
    ]
    rvdb_anno = load_rvdb_domain_keywords(rvdb_search)

    blastn_search = [
        blastn_file,
        os.path.join(diretorio_saida, "validacao_ncbi_blastn.tsv"),
        os.path.join(pasta_posicional, "validacao_ncbi_blastn.tsv"),
        os.path.join(parent_dir, "validacao_ncbi_blastn.tsv"),
        "validacao_ncbi_blastn.tsv"
    ]
    blastn_hits = load_blastn_validation_hits(blastn_search)

    cluster_search = [
        cluster_tsv,
        os.path.join(parent_dir, "08_clusters_e_tensores", "flebo_cluster_v2_segmentos_classificacao.tsv"),
        "resultados_organizados/08_clusters_e_tensores/flebo_cluster_v2_segmentos_classificacao.tsv"
    ]
    cluster_map = load_cluster_segmentation_map(cluster_search)

    # 4. Comparação Posicional (Viral vs. Não-Viral & Filtro Retro/HAM/RT/HGT)
    comparison_results = compare_viral_and_nonviral_positions(
        consensus_data=consensus_data,
        general_hits_annotated=general_annotated,
        hmm_hits=hmm_hits,
        blastn_hits=blastn_hits,
        cluster_map=cluster_map,
        rvdb_anno=rvdb_anno,
        title_dict=title_dict
    )

    # 5. Exportação dos Relatórios
    relatorios = export_positional_validation_reports(
        comparison_results=comparison_results,
        consensus_data=consensus_data,
        diretorio_saida=pasta_posicional
    )

    return {
        "pasta_posicional": pasta_posicional,
        "reciproco": res_reciproco,
        "consenso": consensus_data,
        "comparacao": comparison_results,
        "relatorios": relatorios
    }


def main():
    parser = argparse.ArgumentParser(
        description="Módulo de Validação Posicional e Comparação Evolutiva Viral vs. Não-Viral de EVEs."
    )
    parser.add_argument("--candidatas-nt", required=True, help="FASTA nucleotídico de EVEs candidatas/ORFs unificadas")
    parser.add_argument("--proteinas-virais", required=True, help="FASTA de proteínas virais de referência (AA)")
    parser.add_argument("--diamond-db-geral", help="Banco DIAMOND geral (.dmnd) para o 2º blastx")
    parser.add_argument("--blastx-geral-tsv", help="Arquivo TSV existente do 2º blastx (formato estendido ou padrão)")
    parser.add_argument("--blastx-geral-xml", help="Arquivo XML existente do 2º blastx (outfmt 5)")
    parser.add_argument("--hmm-file", help="Arquivo de resultados HMMER (triagem_profunda_hmm.tsv)")
    parser.add_argument("--blastn-file", help="Arquivo de validação BLASTN a nível de nucleotídeo")
    parser.add_argument("--cluster-tsv", help="Arquivo de classificação de clusters de segmentos")
    parser.add_argument("-o", "--outdir", default="validacao_posicional_saida", help="Diretório de saída")
    parser.add_argument("-t", "--threads", default=4, type=int, help="Número de threads")
    parser.add_argument("-e", "--evalue", default=1e-3, type=float, help="E-value de corte")
    parser.add_argument("--hits-per-query", default=100, type=int, help="Máximo de hits coletados por query (padrão: 100)")
    parser.add_argument("--xml", action="store_true", help="Gera também arquivos no formato BLAST XML (outfmt 5)")

    args = parser.parse_args()

    arquivo_geral = args.blastx_geral_xml if args.blastx_geral_xml else args.blastx_geral_tsv

    executar_analise_posicional_completa(
        fasta_candidatas_nt=args.candidatas_nt,
        fasta_virais_ref=args.proteinas_virais,
        diretorio_saida=args.outdir,
        blastx_geral_tsv_ou_xml=arquivo_geral,
        diamond_db_geral=args.diamond_db_geral,
        hmm_file=args.hmm_file,
        blastn_file=args.blastn_file,
        cluster_tsv=args.cluster_tsv,
        threads=args.threads,
        evalue=args.evalue,
        hits_por_query=args.hits_per_query,
        export_xml=args.xml
    )


if __name__ == "__main__":
    main()
