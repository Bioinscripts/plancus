import os
import re
import argparse

import sys
from Bio import SeqIO
from evetool.core.utils import verificar_dependencias, filtrar_taxonomia_refseq, unir_e_extrair_orfs
from evetool.discovery.orf_predictor import extrair_orfs
from evetool.discovery.cluster_filter import remover_redundancia, resgatar_nucleotideos_representantes, caminho_tabela_clusters
from evetool.discovery.sequence_searcher import preparar_e_buscar, buscar_proteinas_referencia_diamond
from evetool.discovery.hmm_scanner import resgatar_peptideos_nao_alinhados, escanear_hmm
from evetool.discovery.genome_mapper import preparar_banco_blast, rodar_tblastn, mesclar_tabelas_tblastn, gerar_gff_clusters, converter_gff_para_bed
from evetool.discovery.locus_extraction import agrupar_e_extrair_contexto
from evetool.discovery.online_blaster import rodar_blastx_api_em_lotes, extrair_melhores_hits_proteicos
from evetool.annotation.rnaseq_coverage import run_eve_coverage
from evetool.visualization.plot_genomic_context import gerar_figuras_contexto_lote
from evetool.comparative.cross_species_comparator import carregar_lista_especies, consolidar_eves_multiespecies, comparar_clusters_entre_especies
from evetool.core.output_layout import init_pipeline_dirs, generate_readme_structure
from evetool.validation.viral_positional_analyzer import executar_analise_posicional_completa
from evetool import config

# Code and results live apart: runs go to results_dir (config.py) unless -o says otherwise
PASTA_CODIGO = config.PASTA_CODIGO
RESULTADOS_PADRAO = config.RESULTS_DIR


def versao_banco_viral(hmm_db: str) -> dict:
    """Versions of the viral databases in use, from the build manifest next to the HMMs (if any)."""
    import json
    caminho = os.path.join(os.path.dirname(hmm_db or ""), "build_manifest.json")
    try:
        with open(caminho, encoding="utf-8") as fh:
            m = json.load(fh)
        f = m.get("sources", {})
        return {"build": m.get("build"), "date": m.get("date"), "recipe": m.get("recipe", {}).get("recipe_version"),
                "rvdb": f.get("rvdb", {}).get("version"), "eggnog": m.get("recipe", {}).get("eggnog", {}).get("version"),
                "refseq_viral": f.get("refseq", {}).get("date"), "swissprot": f.get("swissprot", {}).get("release")}
    except Exception:
        return {"build": "no build manifest (bundle v1: RVDB 32.0 + eggNOG 5.0, or a local database)"}


def _registrar_info_pipeline(outdir: str, genoma: str, species: str = None, gff: str = None,
                             flanco: int = 20000, hmm_db: str = None, diamond_db: str = None):
    """Record the genome and annotation paths so the curation step can compute GC and context."""
    import json
    try:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "pipeline_info.json"), "w", encoding="utf-8") as fh:
            json.dump({"genome": os.path.abspath(genoma) if genoma else None,
                       "hmm_db": hmm_db,
                       "viral_db": versao_banco_viral(hmm_db),
                       # the curation realigns candidates to these viral proteins (protein map)
                       "diamond_db": os.path.abspath(diamond_db) if diamond_db else None,
                       "species": species, "gff": os.path.abspath(gff) if gff else None,
                       "flank": flanco}, fh, ensure_ascii=False, indent=1)
    except Exception as erro:
        print(f"[warning] Could not write pipeline_info.json: {erro}")


def _cdsearch(args):
    """Optional extra layer: CD-Search over every candidate (local rpsblast+ on CDD by default)."""
    if not getattr(args, "cdsearch", False):
        return
    try:
        from evetool.annotation.cdsearch import anotar as anotar_cdsearch
        anotar_cdsearch(args.outdir, modo=args.cdsearch_mode, threads=args.threads)
    except Exception as erro:
        print(f"[warning] CD-Search annotation failed (it can be rerun from the page): {erro}")


def _pfam(args):
    """Optional extra layer: Pfam-A domains over every candidate (hmmscan --cut_ga)."""
    if not getattr(args, "pfam", False):
        return
    try:
        from evetool.annotation.pfam import anotar as anotar_pfam
        anotar_pfam(args.outdir, pfam_db=args.pfam_db, threads=args.threads)
    except Exception as erro:
        print(f"[warning] Pfam annotation failed (it can be rerun from the page): {erro}")


def _varredura_tblastn(args, db_diamond: str = None):
    """
    Last analysis: tblastn of the best viral hits of the screening over the genome, segments of
    nearby HSPs and characterisation of every ORF inside them (CD-Search and Pfam on the new ones).
    """
    if getattr(args, "no_tblastn_scan", False):
        return
    print(f"\n{'='*75}\n[tblastn scan] Viral references back on the genome\n{'='*75}")
    try:
        from evetool.discovery.tblastn_scan import varrer
        r = varrer(args.outdir, hits_por_orf=args.scan_hits_per_orf, evalue=args.scan_evalue,
                   gap=args.scan_gap, threads=args.threads, dmnd_db=db_diamond)
        print(f"[run] tblastn scan: {r['hsps']} HSPs in {r['segments']} segments, "
              f"{r['new_orfs']} ORFs not analysed before")
    except Exception as erro:
        print(f"[warning] tblastn scan failed (it can be rerun from the page): {erro}")


def _repeticoes(args):
    """Optional: repeats, low complexity and transposon domains around the EVEs (repeats.py)."""
    if not getattr(args, "repeats", False):
        return
    print(f"\n{'='*75}\n[repeats] Repeats and low complexity around the EVEs\n{'='*75}")
    try:
        from evetool.annotation.repeats import analisar
        r = analisar(args.outdir, flanco=args.repeats_flank, threads=args.threads)
        print(f"[run] repeats: {r['loci']} EVE loci · " + ", ".join(f"{v} {k}" for k, v in r["summary"].items()))
    except Exception as erro:
        print(f"[warning] repeat analysis failed (it can be rerun from the page): {erro}")


def _classificar_eves(outdir: str, multi: bool, ignorar_proteoma=()):
    """Classify the candidates by viral vs cellular score (see classify_eves.py)."""
    import subprocess
    cmd = [sys.executable, "-m", "evetool.validation.classify_eves", "-o", outdir]
    if multi:
        cmd.append("--multi")
    if ignorar_proteoma:
        cmd += ["--ignorar-proteoma", *ignorar_proteoma]
    try:
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            [PASTA_CODIGO] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]))
        subprocess.run(cmd, check=True, env=env)
    except Exception as erro:
        print(f"[warning] Score-based EVE classification failed: {erro}")


def _cobertura_rnaseq(outdir: str, args, genome: str, gff: str = None, r1: str = None,
                      r2: str = None, bams=None, species: str = None):
    """RNA-seq coverage over each EVE locus and its flanks (rnaseq_coverage.py), when reads or BAMs exist."""
    bams = [b for b in (bams or []) if b]
    if r1 and not os.path.exists(r1):
        print(f"[warning] RNA-seq reads not found: {r1}")
        r1 = r2 = None
    if not r1 and not bams:
        return
    print(f"\n{'='*75}\n[RNA-seq] EVE and flank coverage{f' - {species}' if species else ''}\n{'='*75}")
    try:
        run_eve_coverage(
            results_dir=outdir, genome=genome, gff=gff, flank=args.flank, r1=r1, r2=r2, bams=bams,
            species=species, threads=args.threads, library=args.rnaseq_library,
            min_mapq=args.rnaseq_min_mapq, read_type=args.rnaseq_read_type,
            reference=args.rnaseq_reference, plots=not args.no_plots)
    except Exception as erro:
        print(f"[warning] RNA-seq coverage failed: {erro}")


def _criar_links_conveniencia(diretorio_raiz: str, mapa_arquivos: dict):
    """Create relative symlinks at the root, for compatibility with legacy scripts."""
    for nome_link, alvo in mapa_arquivos.items():
        if alvo and os.path.exists(alvo):
            link_path = os.path.join(diretorio_raiz, nome_link)
            if not os.path.exists(link_path):
                try:
                    rel_target = os.path.relpath(alvo, diretorio_raiz)
                    os.symlink(rel_target, link_path)
                except Exception:
                    pass


def _coordenadas_do_cabecalho(descricao: str):
    """
    Extract (contig, start, end, strand) from an ORF header, in BED convention
    (0-based, end exclusive), which is what orfipy and bedtools getfasta use.

    Accepted layouts:
      1. 'contig_ORF.N [start-end](+) ...'     -> orfipy on the genome: exact contig coordinates.
      2. 'contig:start-end(+)_ORF.N [a-b](+)'  -> orfipy on pre-extracted segments: uses the
         segment interval in the genome (the ORF offset inside it is absorbed by the
         bedtools merge of loci).

    Returns None when the header carries no usable coordinates.
    """
    qseqid = descricao.split()[0]

    m = re.match(r"^(.+?):(\d+)-(\d+)\(([+-])\)", qseqid)
    if m:
        contig, start, end, strand = m.groups()
        return contig, int(start), int(end), strand

    if "_ORF." in qseqid:
        m = re.search(r"\[(\d+)-(\d+)\]\(([+-])\)", descricao)
        if m:
            start, end, strand = m.groups()
            return qseqid.rsplit("_ORF.", 1)[0], int(start), int(end), strand

    return None


def executar_triagem_especie(
    genoma_fasta: str,
    outdir: str,
    threads: int,
    args,
    db_diamond: str,
    is_fasta_diamond: bool,
    db_hmm: str,
    species: str = None
) -> dict:
    """
    Stage 1: ORF extraction, redundancy reduction and primary screening (DIAMOND + HMM).
    Files are written to the structured folders:
      - 01_orfs_cds/
      - 02_alinhamento_diamond/
      - 03_triagem_hmm/
      - 04_candidatas_eves/
    Returns the paths of the unified candidate FASTAs, in nucleotide and peptide form.
    """
    dir_map = init_pipeline_dirs(outdir)
    dir_cds = dir_map["01_orfs_cds"]
    dir_diamond = dir_map["02_alinhamento_diamond"]
    dir_hmm = dir_map["03_triagem_hmm"]
    dir_candidatas = dir_map["04_candidatas_eves"]

    tag_sp = f"[{species}] " if species else ""
    print(f"\n{'-'*60}\n{tag_sp}ORF mining and viral screening: {genoma_fasta}\n{'-'*60}")

    # 1. ORF extraction
    print(f"{tag_sp}[run] Extracting ORFs with orfipy...")
    orfs = extrair_orfs(genoma_fasta, dir_cds)

    # 2. Deduplication (90%)
    print(f"{tag_sp}[run] Removing redundancy with MMseqs2...")
    pep_nr = remover_redundancia(orfs["pep_fasta"], dir_cds, threads=threads)
    nt_nr = resgatar_nucleotideos_representantes(pep_nr, orfs["nt_fasta"], dir_cds)

    # 3. Primary screening (DIAMOND)
    print(f"{tag_sp}[run] Primary screening with DIAMOND blastx...")
    tsv_diamond, fasta_unaligned = preparar_e_buscar(
        db_diamond, is_fasta_diamond, nt_nr, dir_diamond, threads,
        hits_por_query=getattr(args, "viral_hits_per_query", 50)
    )
    tsv_diam_limpo = filtrar_taxonomia_refseq(
        tsv_diamond, dir_diamond, is_custom_db=is_fasta_diamond,
        manual_removal_list=getattr(args, "lista_remocao_manual", None),
        margem_ambiguo=getattr(args, "margem_ambiguo", 0.9)
    )

    # 4. Deep screening (HMMER)
    print(f"{tag_sp}[run] Rescuing deep viral scars with HMM...")
    pep_unaligned = resgatar_peptideos_nao_alinhados(fasta_unaligned, pep_nr, dir_hmm)
    tsv_hmm = escanear_hmm(
        pep_unaligned, dir_hmm, hmm_dir=args.hmm_dir, hmm_db=args.hmm_db, threads=threads
    )

    # 5. Unification (exports both nucleotide and peptide, for clean tblastn anchoring)
    candidatas = unir_e_extrair_orfs(
        tsv_diam_limpo, tsv_hmm, nt_nr, dir_candidatas, pep_fasta=pep_nr,
        # Recover the other copies of each validated cluster; the member sequences live
        # in the original FASTAs, before deduplication.
        tsv_clusters=caminho_tabela_clusters(pep_nr),
        nt_completo=orfs["nt_fasta"], pep_completo=orfs["pep_fasta"]
    )
    fasta_final_nt = candidatas["nt"]
    fasta_final_pep = candidatas["pep"]
    print(f"{tag_sp}[run] Candidates identified:")
    print(f"       Nucleotide: {fasta_final_nt}")
    print(f"       Peptide   : {fasta_final_pep}")

    # Convenience links at the root, for compatibility
    _criar_links_conveniencia(outdir, {
        "eve_candidatas_unificadas.fasta": fasta_final_nt,
        "eve_candidatas_unificadas_pep.fasta": fasta_final_pep
    })

    return {
        "candidatas_nt": fasta_final_nt,
        "candidatas_pep": fasta_final_pep,
        "diamond_tsv": tsv_diam_limpo,
        "hmm_tsv": tsv_hmm
    }



def executar_clusterizacao_hierarquica(fasta_exato: str, diretorio_saida: str, threads: int = 4):
    import os, subprocess
    if not fasta_exato or not os.path.exists(fasta_exato):
        return
    print(f"[run] Running the advanced hierarchical clustering on the EVEs...")
    os.makedirs(diretorio_saida, exist_ok=True)
    db_path = os.path.join(diretorio_saida, "eves_db")
    xml_out = os.path.join(diretorio_saida, "tblastx_eves.xml")
    
    try:
        subprocess.run(["makeblastdb", "-in", fasta_exato, "-dbtype", "nucl", "-out", db_path], check=True, stdout=subprocess.DEVNULL)
        cmd = [
            "tblastx", "-query", fasta_exato, "-db", db_path,
            "-out", xml_out, "-evalue", "1e-3", "-num_threads", str(threads),
            "-outfmt", "5"
        ]
        subprocess.run(cmd, check=True)
        
        from evetool.comparative.hierarchical_clustering import run_clustering
        run_clustering(fmt5=xml_out, out=os.path.join(diretorio_saida, "eve_clusters_hierarch"), no_r=False)
    except Exception as e:
        print(f"[warning] Hierarchical clustering failed: {e}")

def executar_ancoragem_e_contexto_especie(
    genoma_fasta: str,
    fasta_candidatas_pep: str,
    fasta_proteinas_ref: str,
    outdir: str,
    threads: int,
    gff_path: str = None,
    species: str = None,
    flanco: int = 20000,
    gerar_plots: bool = True
) -> dict:
    """
    Stage 3: tblastn anchoring on the genome (local ORFs plus the consolidated viral reference pool),
    GFF3 generation, 20 kb context extraction and context plots.
    Files are written to the structured folders:
      - 05_ancoragem_tblastn/
      - 06_anotacao_genomica/
      - 09_graficos_e_figuras/figuras_contexto_20kb/
    """
    dir_map = init_pipeline_dirs(outdir)
    dir_tblastn = dir_map["05_ancoragem_tblastn"]
    dir_anotacao = dir_map["06_anotacao_genomica"]
    dir_plots = os.path.join(dir_map["09_graficos_e_figuras"], "figuras_contexto_20kb")
    dir_rnaseq = dir_map["10_rnaseq_expressao"]

    tag_sp = f"[{species}] " if species else ""
    print(f"\n{'-'*60}\n{tag_sp}Genomic anchoring and context extraction (20 kb): {genoma_fasta}\n{'-'*60}")

    try:
        print(f"{tag_sp}[run] Skipping tblastn: taking coordinates straight from the EVE headers...")
        from Bio import SeqIO
        import subprocess
        from evetool.core.utils import is_hit_celular

        queries_proibidas = set()
        geral_tsv = os.path.join(os.path.dirname(fasta_candidatas_pep), "diamond_blastx_referencia_geral.tsv")
        if os.path.exists(geral_tsv):
            with open(geral_tsv, 'r') as fg:
                for line_tsv in fg:
                    if line_tsv.startswith('#') or not line_tsv.strip(): continue
                    parts = line_tsv.strip().split('\t')
                    if len(parts) < 2: continue
                    qseqid = parts[0]
                    stitle = parts[16] if len(parts) >= 18 else (parts[-1] if len(parts) >= 15 else "")
                    title_lower = stitle.lower()
                    sseqid = parts[1]
                    if is_hit_celular(sseqid) or "zinc finger" in title_lower or "retrovir" in title_lower or "reverse transcriptase" in title_lower or "pol gag" in title_lower:
                        queries_proibidas.add(qseqid)

        # Coordinates from the ORF headers; only valid when most of them are usable
        registros = list(SeqIO.parse(fasta_candidatas_pep, "fasta"))
        coordenadas = [(r.description, _coordenadas_do_cabecalho(r.description)) for r in registros]
        n_ok = sum(1 for _, c in coordenadas if c)
        tem_coordenadas = n_ok > 0 and n_ok >= 0.5 * len(coordenadas)

        if tem_coordenadas:
            print(f"{tag_sp}[run] Direct mapping: coordinates taken from the headers "
                  f"({n_ok}/{len(coordenadas)} ORFs).")
            bed_saida = os.path.join(dir_anotacao, "eve_hits_diretos.bed")
            escritas, sem_coord, descartadas = 0, [], 0
            with open(bed_saida, 'w') as fbed:
                for descricao, coord in coordenadas:
                    qseqid = descricao.split()[0]
                    if qseqid in queries_proibidas:
                        descartadas += 1
                        continue  # Falso positivo celular
                    if not coord:
                        # Without coordinates there is nothing to anchor: record it and move
                        # on, rather than writing a broken line that corrupts the whole BED.
                        sem_coord.append(qseqid)
                        continue
                    contig, start, end, strand = coord
                    fbed.write(f"{contig}\t{start}\t{end}\t{qseqid}\t1000\t{strand}\n")
                    escritas += 1

            if descartadas:
                print(f"{tag_sp}[run] {descartadas} ORFs dropped for a cellular hit in the general database.")
            if sem_coord:
                print(f"{tag_sp}[warning] {len(sem_coord)} ORFs without coordinates in the header were "
                      f"skipped (e.g. {', '.join(sem_coord[:3])}).")
            if escritas == 0:
                raise RuntimeError("No ORF with valid coordinates to anchor on the genome.")

            bed_ordenado = os.path.join(dir_anotacao, "eve_hits_diretos.sorted.bed")
            bed_mesclado = os.path.join(dir_anotacao, "eve_hits_diretos.merged.bed")
            subprocess.run(f"bedtools sort -i {bed_saida} > {bed_ordenado}", shell=True, check=True)
            # Strand becomes '.' when the cluster mixes both strands, so column 6 never
            # holds '+,-' and break whatever reads the BED downstream.
            subprocess.run(
                f'bedtools merge -i {bed_ordenado} -d 1000 -c 4,5,6 -o count,mean,distinct '
                f'| awk \'BEGIN{{OFS="\\t"}}{{$4="Cluster_"NR; if ($6 ~ /,/) $6="."; print $0}}\' '
                f'> {bed_mesclado}',
                shell=True, check=True)
            bed_final = bed_mesclado

            # Build the structured GFF3 (clusters -> ORFs) by distance
            from evetool.discovery.gff_grouping import cluster_bed_to_gff
            gff_final = os.path.join(dir_anotacao, "eve_clusters.gff3")
            cluster_bed_to_gff(bed_ordenado, gff_final, max_gap=1000)
            print(f"{tag_sp}[run] Structured cluster GFF3 written to: {gff_final}")
        else:
            print(f"{tag_sp}[run] Anchoring the candidates on the host genome with tblastn...")
            db_blast = preparar_banco_blast(genoma_fasta, dir_tblastn)
            tabela_tblastn = rodar_tblastn(fasta_candidatas_pep, db_blast, dir_tblastn, threads=threads)
            gff_final = os.path.join(dir_anotacao, "eve_clusters.gff3")
            gerar_gff_clusters(tabela_tblastn, gff_final)
            bed_final = os.path.join(dir_anotacao, "eve_hits_diretos.merged.bed")
            converter_gff_para_bed(gff_final, bed_final)

        from evetool.discovery.locus_extraction import agrupar_e_extrair_contexto
        contexto = agrupar_e_extrair_contexto(bed_final, genoma_fasta, dir_anotacao, flanco=flanco)



        print(f"{tag_sp}[run] Genomic context extracted to: {contexto['fasta_contexto']}")
        
        executar_clusterizacao_hierarquica(
            fasta_exato=contexto.get("fasta_exato"), 
            diretorio_saida=os.path.join(outdir, "08_clusters_e_tensores"), 
            threads=threads
        )


    except RuntimeError as erro:
        print(f"{tag_sp}[warning] Genomic anchoring produced no hit or failed: {erro}")
        return {"gff": None, "bed": None, "contexto": None}

    # RNA-seq coverage runs after the classification (see _cobertura_rnaseq), which
    # gives the loci; the context figures read its per-base table when it already exists.
    rnaseq_depth_tsv = os.path.join(dir_rnaseq, "rnaseq_depth_per_base.tsv")
    if not os.path.exists(rnaseq_depth_tsv):
        rnaseq_depth_tsv = None

    # Genomic context figures (with annotation and RNA-seq, when available)
    if gerar_plots and bed_final and contexto:
        try:
            print(f"{tag_sp}[run] Drawing the genomic context figures (20 kb flanks)...")
            gerar_figuras_contexto_lote(
                bed_original=bed_final,
                bed_flanked=contexto["bed_final"],
                gff_annotation=gff_path,
                depth_tsv=rnaseq_depth_tsv,
                output_dir=dir_plots,
                species=species,
                flank=flanco
            )
        except Exception as erro:
            print(f"{tag_sp}[warning] Could not draw the context figures: {erro}")

    # Convenience links at the species root, for legacy tools
    _criar_links_conveniencia(outdir, {
        "eve_tblastn_hits.tsv": os.path.join(dir_tblastn, "eve_tblastn_hits.tsv"),
        "eve_clusters.gff3": os.path.join(dir_anotacao, "eve_clusters.gff3"),
        "eve_clusters.bed": os.path.join(dir_anotacao, "eve_clusters.bed"),
        "eve_clusters_context.fasta": contexto["fasta_contexto"] if contexto else None,
        "figuras_contexto_20kb": dir_plots
    })

    return {
        "gff": gff_final,
        "bed": bed_final,
        "contexto": contexto
    }


def main():
    parser = argparse.ArgumentParser(
        description="PLANCUS (EVE tool): mining, validation, RNA-seq coverage and cross-species comparison."
    )

    # Modos de Entrada
    input_group = parser.add_argument_group("Input data")
    input_group.add_argument("-i", "--input", help="Host genome in FASTA format (single-species mode)")
    input_group.add_argument("--genomes-dir", help="Folder with FASTA genomes for the multi-species analysis")
    input_group.add_argument("--species-tsv", help="TSV listing the species and their files (genome, gff, rnaseq_r1, rnaseq_r2; any cell ending in .bam is used as an aligned RNA-seq BAM)")
    input_group.add_argument("--species", default=None, help="Host species name (single-species mode), used in the dashboard instead of guessing from the contig prefix")
    input_group.add_argument("-o", "--outdir", default=None,
                             help=f"Root output folder (default: {RESULTADOS_PADRAO}/results_<date>_<time>; "
                                  "a relative name is placed under that folder too)")
    input_group.add_argument("-t", "--threads", default=4, type=int, help="Number of threads for parallel execution")

    # Optional context and RNA-seq modules
    mod_group = parser.add_argument_group("Optional context and RNA-seq modules (single-species mode)")
    mod_group.add_argument("--gff", help="GFF3/GTF file with the host genome annotation")
    mod_group.add_argument("--rnaseq-r1", help="Raw RNA-seq FASTQ (R1 or single-end)")
    mod_group.add_argument("--rnaseq-r2", help="Raw RNA-seq FASTQ (R2 for paired-end)")
    mod_group.add_argument("--rnaseq-bam", nargs="+", default=[],
                           help="Aligned RNA-seq BAM(s), to the genome or to EVE context windows; used without re-aligning")
    mod_group.add_argument("--rnaseq-library", default="auto", choices=["auto", "fr", "rf", "unstranded"],
                           help="RNA-seq strandedness; auto infers it from the host exons of the GFF (default)")
    mod_group.add_argument("--rnaseq-read-type", default="sr", choices=["sr", "ont", "pb", "hifi"],
                           help="Read technology: sr = Illumina (minimap2 -ax sr), ont/pb/hifi = long reads (minimap2 splice)")
    mod_group.add_argument("--rnaseq-reference", default="context", choices=["context", "genome"],
                           help="Align reads to the EVE windows (fast, default) or to the whole genome (no paralogue recruitment)")
    mod_group.add_argument("--rnaseq-min-mapq", default=0, type=int,
                           help="Ignore RNA-seq alignments below this mapping quality (default 0)")
    mod_group.add_argument("--flank", default=20000, type=int, help="Upstream and downstream flank size (default: 20000 nt)")
    mod_group.add_argument("--no-plots", action="store_true", help="Turn off the genomic locus figures")
    mod_group.add_argument("--cdsearch", action="store_true",
                           help="Annotate every candidate with CD-Search (CDD domains; RNA and non-retroviral "
                                "ssDNA virus scopes)")
    mod_group.add_argument("--cdsearch-mode", choices=("local", "online"), default="local",
                           help="local: rpsblast+ on the local CDD (default, offline); "
                                "online: NCBI Batch CD-Search")
    mod_group.add_argument("--pfam", action="store_true",
                           help="Annotate every candidate with Pfam-A domains (hmmscan with the Pfam "
                                "gathering thresholds)")
    mod_group.add_argument("--pfam-db", default=config.PFAM_HMM,
                           help="Pfam-A.hmm pressed with hmmpress (default: %(default)s)")
    mod_group.add_argument("--repeats", action="store_true",
                           help="Also analyse repeats, low complexity and transposon domains around the EVEs")
    mod_group.add_argument("--repeats-flank", type=int, default=20000,
                           help="Flank read on each side of every EVE for the repeat analysis (default 20000 nt)")
    mod_group.add_argument("--no-tblastn-scan", action="store_true",
                           help="Skip the last analysis: tblastn of the best viral hits over the genome")
    mod_group.add_argument("--scan-hits-per-orf", type=int, default=1,
                           help="Viral hits of each screened ORF taken to the tblastn scan (default 1)")
    mod_group.add_argument("--scan-evalue", type=float, default=1e-5, help="tblastn scan E-value (default 1e-5)")
    mod_group.add_argument("--scan-gap", type=int, default=1000,
                           help="tblastn HSPs closer than this (nt) form one segment (default 1000)")

    # DIAMOND and HMM databases
    # Not given: the databases installed by prepare.py (paths in config.py)
    diamond_group = parser.add_mutually_exclusive_group()
    diamond_group.add_argument("--dmnd-db", help="Formatted viral DIAMOND database (.dmnd; default: mydb.dmnd "
                                                 "of the configured diamond_dir)")
    diamond_group.add_argument("--fasta-dir", help="Folder with viral FASTAs to build the DIAMOND database")

    hmm_group = parser.add_mutually_exclusive_group()
    hmm_group.add_argument("--hmm-db", help="Pre-indexed HMM database (default: the calibrated combined "
                                            "database of the configured hmm_dir)")
    hmm_group.add_argument("--hmm-dir", help="Folder with raw viral .hmm files")

    # Reference databases for the second anchoring
    ref_group = parser.add_argument_group("Second genomic anchoring (optional)")
    ref_group.add_argument("--diamond-db-geral", help="General DIAMOND database (.dmnd, e.g. nr) for local reference-protein collection")
    ref_group.add_argument("--diamond-hits-per-query", default=100, type=int,
                           help="Maximum best hits collected per query in the DIAMOND run against the general database (default: 100)")
    ref_group.add_argument("--diamond-top-hits", default=0, type=int,
                           help="Global cap on reference proteins (default: 0 = no cap, keep every unique top hit)")
    ref_group.add_argument("--validar-ncbi", action="store_true", help="Submit batches to the NCBI online blastx to fetch reference proteins")
    ref_group.add_argument("--ncbi-email", default=config.NCBI_EMAIL,
                           help="Contact e-mail sent with every NCBI Entrez request, required by the NCBI usage "
                                "policy (default: EVE_NCBI_EMAIL, then [ncbi] email in eve_config.ini)")
    ref_group.add_argument("--ncbi-lote-max-seqs", default=10, type=int, help="Maximum sequences per batch in the NCBI blastx")
    ref_group.add_argument("--ncbi-lote-max-nt", default=20000, type=int, help="Maximum nucleotides per batch in the NCBI blastx")
    ref_group.add_argument("--ncbi-top-hits", default=50, type=int, help="Number of best NCBI hits (default: 50)")
    ref_group.add_argument("--ncbi-timeout", default=1800, type=int, help="Per-batch timeout for the NCBI blastx, in seconds")
    ref_group.add_argument("--export-xml", action="store_true", help="Also write alignment reports in BLAST XML format (outfmt 5)")
    ref_group.add_argument("--no-validacao-posicional", action="store_true", help="Turn off the positional and evolutionary validation module")
    ref_group.add_argument("--viral-hits-per-query", default=50, type=int,
                           help="Hits per ORF in the primary RefSeq viral screening (default: 50). "
                                "The mix of genome types among them decides whether the ORF is an EVE or a host gene")
    ref_group.add_argument("--lista-remocao-manual", nargs="*", default=None,
                           help="One or more files with accessions/terms to remove (one per line). "
                                "They add to the built-in filter instead of replacing it: use them for the curation of "
                                "viral_screening_unknown_to_review.tsv and for accessions_suggested_removal.txt "
                                "written by audit_viral_db.py")
    ref_group.add_argument("--margem-ambiguo", default=0.9, type=float,
                           help="When the best dsDNA/RT hit reaches this fraction of the best RNA/ssDNA hit, "
                                "the ORF is removed as horizontal transfer (default: 0.9)")
    # Cross-species clustering and synteny
    multi_group = parser.add_argument_group("Multi-species clustering and synteny")
    multi_group.add_argument("--min-identity", default=30.0, type=float,
                             help="Minimum identity (%%) to join two EVEs in stage 1 (default: 30)")
    multi_group.add_argument("--min-coverage", default=30.0, type=float,
                             help="Minimum coverage (%%) to join two EVEs in stage 1 (default: 30)")
    multi_group.add_argument("--filter-mode", default="max", choices=["max", "mean", "both"],
                             help="How to apply the thresholds to a pair (default: max)")
    multi_group.add_argument("--max-clusters-heatmap", default=40, type=int,
                             help="Maximum shared clusters shown in the heatmap")
    multi_group.add_argument("--renomear-especies",
                             help="TSV 'current_name<TAB>correct_name' to fix the species labels")
    multi_group.add_argument("--sintenia", action="store_true",
                             help="Evaluate flank synteny for the clusters (off by default); "
                                  "it can also be requested per species, through the 'sintenia' column of --species-tsv")
    multi_group.add_argument("--proteoma-ref",
                             help="Reference proteome used to label the flanking genes in the synteny step")
    multi_group.add_argument("--flanco-sintenia", default=50000, type=int,
                             help="Flank size used by the synteny step (default: 50000 nt)")

    ref_group.add_argument("--ignorar-proteoma", nargs="*", default=[],
                           help="sseqid prefixes ignored by the classification (e.g. the accession of the species' "
                                "own proteome, so an ORF does not find itself)")

    args = parser.parse_args()

    # Validate the input arguments
    if not args.input and not args.species_tsv:
        parser.error("Provide either a single genome (-i / --input) OR a species list (--species-tsv).")

    # Databases left out: the installed ones (./eve-tool databases; ./eve-tool config shows them)
    if not args.dmnd_db and not args.fasta_dir:
        if not os.path.exists(config.VIRAL_DMND):
            parser.error(f"No viral DIAMOND database: {config.VIRAL_DMND} is missing. Install it with "
                         f"`./eve-tool databases` or pass --dmnd-db / --fasta-dir.")
        args.dmnd_db = config.VIRAL_DMND
    if not args.hmm_db and not args.hmm_dir:
        if not os.path.exists(config.HMM_DB):
            parser.error(f"No viral HMM database: {config.HMM_DB} is missing. Install it with "
                         f"`./eve-tool databases` or pass --hmm-db / --hmm-dir.")
        args.hmm_db = config.HMM_DB
    if not args.diamond_db_geral and os.path.exists(config.GERAL_DMND):
        args.diamond_db_geral = config.GERAL_DMND
    if args.pfam and not os.path.exists(args.pfam_db):
        parser.error(f"--pfam: {args.pfam_db} is missing (./eve-tool databases --pfam).")
    # Checked here, before hours of screening, rather than when the download is reached
    if (args.diamond_db_geral or args.validar_ncbi) and not args.ncbi_email:
        parser.error("The reference viral proteins are downloaded from the NCBI, which asks for a contact "
                     "e-mail. Pass --ncbi-email you@institution.org, export EVE_NCBI_EMAIL, or set "
                     "`email` under [ncbi] in eve_config.ini.")

    # Results never inside the code folder: a bare name goes under Resultados/
    if not args.outdir:
        import datetime
        args.outdir = os.path.join(RESULTADOS_PADRAO, datetime.datetime.now().strftime("results_%Y%m%d_%H%M"))
    elif not os.path.isabs(args.outdir) and os.sep not in args.outdir.rstrip(os.sep):
        args.outdir = os.path.join(RESULTADOS_PADRAO, args.outdir)
    args.outdir = os.path.abspath(args.outdir)
    if os.path.commonpath([args.outdir, PASTA_CODIGO]) == PASTA_CODIGO:
        parser.error(f"The output folder is inside the code folder ({PASTA_CODIGO}). "
                     f"Use a folder under {RESULTADOS_PADRAO}.")
    print(f"[info] Results: {args.outdir}")
    print("[info] Viral databases: " + " · ".join(f"{k} {v}" for k, v in versao_banco_viral(args.hmm_db).items() if v))
    os.makedirs(args.outdir, exist_ok=True)

    # 1. Check the external dependencies
    try:
        import pyhmmer
        tem_pyhmmer = True
    except ImportError:
        tem_pyhmmer = False

    dependencias_obrigatorias = [
        "orfipy", "mmseqs", "diamond",
        "makeblastdb", "tblastn", "tblastx", "bedtools", "samtools"
    ]
    if not tem_pyhmmer:
        dependencias_obrigatorias.extend(["hmmsearch", "hmmpress"])

    verificar_dependencias(dependencias_obrigatorias)

    # Database setup
    is_fasta_diamond = bool(args.fasta_dir)
    db_diamond = args.fasta_dir if is_fasta_diamond else args.dmnd_db
    db_hmm = args.hmm_dir if args.hmm_dir else args.hmm_db

    # =========================================================================
    # MODE 1: SYNCHRONISED TWO-STAGE MULTI-SPECIES ANALYSIS
    # =========================================================================
    if args.species_tsv:
        print(f"\n{'='*75}\n[Multi-species mode] Loading the species file: {args.species_tsv}\n{'='*75}")
        lista_especies = carregar_lista_especies(args.species_tsv, args.genomes_dir)

        # ---------------------------------------------------------------------
        # STAGE 1: PRIMARY SCREENING ACROSS EVERY SPECIES
        # ---------------------------------------------------------------------
        print(f"\n{'='*75}\n[STAGE 1/3] Mining and primary screening across all {len(lista_especies)} species\n{'='*75}")
        dados_triagem = {}

        for item in lista_especies:
            sp_nome = item["species"]
            genoma_arq = item["genome"]
            sp_outdir = os.path.join(args.outdir, "especies", sp_nome)

            _registrar_info_pipeline(sp_outdir, genoma_arq, sp_nome, item.get("gff"), args.flank, db_hmm, db_diamond)

            res_triagem = executar_triagem_especie(
                genoma_fasta=genoma_arq,
                outdir=sp_outdir,
                threads=args.threads,
                args=args,
                db_diamond=db_diamond,
                is_fasta_diamond=is_fasta_diamond,
                db_hmm=db_hmm,
                species=sp_nome
            )
            dados_triagem[sp_nome] = {
                "item": item,
                "outdir": sp_outdir,
                "screening": res_triagem
            }

        # ---------------------------------------------------------------------
        # FASE 2: SEGUNDO BLASTX (BANCO GERAL) COLETANDO OS TOP 100 HITS POR QUERY
        # ---------------------------------------------------------------------
        print(f"\n{'='*75}\n[STAGE 2/3] Consolidated collection of reference proteins (top {args.diamond_hits_per_query} hits/query)\n{'='*75}")
        pasta_referencias = os.path.join(args.outdir, "04_candidatas_eves")
        os.makedirs(pasta_referencias, exist_ok=True)

        # Merge every candidate EVE from every species
        fasta_todas_candidatas_nt = os.path.join(pasta_referencias, "todas_candidatas_todas_especies_nt.fasta")
        total_seqs = 0
        with open(fasta_todas_candidatas_nt, "w") as f_out:
            for sp_nome, info in dados_triagem.items():
                caminho_nt = info["screening"]["candidatas_nt"]
                if caminho_nt and os.path.exists(caminho_nt):
                    for record in SeqIO.parse(caminho_nt, "fasta"):
                        record.id = f"{sp_nome}|{record.id}"
                        record.description = ""
                        SeqIO.write(record, f_out, "fasta")
                        total_seqs += 1

        print(f"[run] {total_seqs} candidate EVEs consolidated from every species for the second blastx.")

        fasta_proteinas_ref_pool = None
        if args.diamond_db_geral:
            try:
                print(f"[run] Running the second DIAMOND blastx against the general database for every species...")
                buscar_proteinas_referencia_diamond(
                    fasta_query=fasta_todas_candidatas_nt,
                    banco_geral=args.diamond_db_geral,
                    diretorio_saida=pasta_referencias,
                    threads=args.threads,
                    hits_por_query=args.diamond_hits_per_query,
                    top_n_total=args.diamond_top_hits if args.diamond_top_hits > 0 else None,
                    export_xml=args.export_xml
                )
                
                # Consolidate viral_tsv paths
                # But wait, we don't have a single viral_tsv, we have one for each species!
                # Actually, triagem["diamond_tsv"] is in `dados_triagem[sp_nome]["triagem"]["diamond_tsv"]`.
                # Let's concatenate them!

                viral_tsvs = [info["screening"]["diamond_tsv"] for info in dados_triagem.values()]
                merged_viral_tsv = os.path.join(pasta_referencias, "todas_candidatas_diamond_viral.tsv")
                with open(merged_viral_tsv, "w") as f_out:
                    for v_tsv in viral_tsvs:
                        if os.path.exists(v_tsv):
                            with open(v_tsv, "r") as f_in:
                                f_out.write(f_in.read())

                geral_tsv = os.path.join(pasta_referencias, "diamond_blastx_referencia_geral.tsv")
                out_fasta = os.path.join(pasta_referencias, "referencias_virais_filtradas.fasta")
                from evetool.discovery.sequence_searcher import gerar_referencias_virais_filtradas
                fasta_proteinas_ref_pool = gerar_referencias_virais_filtradas(
                    merged_viral_tsv, geral_tsv, out_fasta, email_ncbi=args.ncbi_email
                )
                print(f"[run] Pool of filtered, genuine viral proteins ready at: {fasta_proteinas_ref_pool}")
            except Exception as erro:
                print(f"[warning] The consolidated DIAMOND search against the general database failed: {erro}")
        elif args.validar_ncbi and args.ncbi_email:
            try:
                print(f"[run] Submitting the unified candidates to the NCBI in batches...")
                tabelas_blastx = rodar_blastx_api_em_lotes(
                    fasta_todas_candidatas_nt, pasta_referencias,
                    max_seqs=args.ncbi_lote_max_seqs,
                    max_nt=args.ncbi_lote_max_nt,
                    max_espera_segundos=args.ncbi_timeout
                )
                fasta_proteinas_ref_pool = extrair_melhores_hits_proteicos(
                    tabelas_blastx, pasta_referencias,
                    top_n=args.ncbi_top_hits,
                    email_ncbi=args.ncbi_email
                )
            except Exception as erro:
                print(f"[warning] The NCBI online validation/collection failed: {erro}")

        # Consolidated positional and evolutionary validation
        if fasta_proteinas_ref_pool and not args.no_validacao_posicional:
            try:
                print(f"\n{'='*75}\n[Positional validation] Consolidated positional analysis and reciprocal tblastn\n{'='*75}")
                dir_posicional = os.path.join(args.outdir, "02_alinhamento_diamond")
                blastx_tsv = os.path.join(pasta_referencias, "diamond_blastx_referencia_geral.tsv")
                executar_analise_posicional_completa(
                    fasta_candidatas_nt=fasta_todas_candidatas_nt,
                    fasta_virais_ref=fasta_proteinas_ref_pool,
                    diretorio_saida=dir_posicional,
                    blastx_geral_tsv_ou_xml=blastx_tsv if os.path.exists(blastx_tsv) else None,
                    diamond_db_geral=args.diamond_db_geral,
                    threads=args.threads,
                    hits_por_query=args.diamond_hits_per_query,
                    export_xml=args.export_xml
                )
            except Exception as erro_pos:
                print(f"[warning] The consolidated positional validation failed: {erro_pos}")

        # ---------------------------------------------------------------------
        # STAGE 3: tblastn ANCHORING ON EVERY GENOME + 20 KB CONTEXT
        # ---------------------------------------------------------------------
        print(f"\n{'='*75}\n[STAGE 3/3] tblastn anchoring on the genomes and 20 kb context extraction\n{'='*75}")
        resultados_finais_especies = {}

        for sp_nome, info in dados_triagem.items():
            item = info["item"]
            sp_outdir = info["outdir"]
            screening = info["screening"]

            res_ancoragem = executar_ancoragem_e_contexto_especie(
                genoma_fasta=item["genome"],
                fasta_candidatas_pep=screening["candidatas_pep"],
                fasta_proteinas_ref=fasta_proteinas_ref_pool,  # reference proteins consolidated from every species
                outdir=sp_outdir,
                threads=args.threads,
                gff_path=item["gff"],
                species=sp_nome,
                flanco=args.flank,
                gerar_plots=not args.no_plots
            )
            if screening.get("candidatas_nt") and os.path.exists(screening["candidatas_nt"]):
                resultados_finais_especies[sp_nome] = screening["candidatas_nt"]

        # ---------------------------------------------------------------------
        # STAGE 4: CROSS-SPECIES COMPARISON (CLUSTERS, HEATMAP, PHYLOGENY)
        # ---------------------------------------------------------------------
        print(f"\n{'='*75}\n[STAGE 4/4] Cross-species comparison of the homologous clusters\n{'='*75}")
        pasta_comparativa = os.path.join(args.outdir, "11_sintenia_multiespecies")
        os.makedirs(pasta_comparativa, exist_ok=True)

        fasta_consolidado = consolidar_eves_multiespecies(
            resultados_finais_especies, pasta_comparativa
        )

        # Two-stage clustering (identity/coverage -> valleys), replacing the connected
        # components, which merged distinct families through single linkage.
        try:
            from evetool.comparative.multispecies_clustering import executar as clusterizar_multiespecies
            from argparse import Namespace
            sintenia_on = args.sintenia or any(i.get("sintenia") for i in lista_especies)
            locais, genomas_sin, origens_sin = {}, {}, {}
            if sintenia_on:
                # coordinates of every candidate and the genome of each species, for the flanks
                from evetool.curation import curation as _cur
                from evetool.comparative.flank_synteny import origens_de_segmentos, calibrar_origens
                for _sp, _pasta in _cur.pastas_especies(args.outdir).items():
                    for _c, _i, _f, _st, _orf in _cur._todos_os_loci(_pasta):
                        locais[f"{_sp}|{_orf}" if _sp else _orf] = (_sp, _c, _i, _f)
                for _sp, _g in _cur.genomas_para_sintenia(args.outdir).items():
                    genomas_sin[_sp] = _g["genome"]
                    if _g["segments"] and _g["own"]:
                        origens_sin[_sp] = calibrar_origens(origens_de_segmentos(_g["own"]), _g["own"], _g["genome"])
            cfg_cluster = Namespace(
                fasta=fasta_consolidado, xml=None, outdir=pasta_comparativa, threads=args.threads,
                min_identity=args.min_identity, min_coverage=args.min_coverage,
                filter_mode=args.filter_mode, valley_ratio=0.5, max_bandwidth=0.05,
                min_pairs=15, max_depth=6, w_cov=1.0, w_gap=1.0, max_overlap_nt=0,
                max_clusters_heatmap=args.max_clusters_heatmap, renomear=args.renomear_especies,
                # Synteny: enabled by the global flag or by the 'sintenia' column of the species TSV
                sintenia=sintenia_on,
                proteoma_ref=args.proteoma_ref, flanco_sintenia=args.flanco_sintenia,
                locais=locais, genomas_sintenia=genomas_sin, origens=origens_sin,
            )
            res_comparativo = clusterizar_multiespecies(cfg_cluster)
        except Exception as erro_cluster:
            print(f"[warning] The two-stage clustering failed ({erro_cluster}); "
                  f"falling back to the connected-components method.")
            try:
                res_comparativo = comparar_clusters_entre_especies(
                    fasta_consolidado=fasta_consolidado,
                    diretorio_saida=pasta_comparativa,
                    threads=args.threads
                )
            except Exception as erro_comp:
                # With a single species there may be nothing to compare; the classification
                # and the RNA-seq coverage below must still run.
                print(f"[warning] Cross-species comparison skipped: {erro_comp}")
                res_comparativo = None

        # Extra layer: HMM domains over every candidate of each species
        try:
            from evetool.annotation.hmm_domains import anotar as anotar_hmm_candidatas
            for sp_nome, info_sp in dados_triagem.items():
                anotar_hmm_candidatas(args.outdir, hmm_db=db_hmm, threads=args.threads,
                                      pasta_especie=info_sp["outdir"], forcar=True)
        except Exception as err_hmm:
            print(f"[warning] HMM annotation of the candidates failed: {err_hmm}")
        _cdsearch(args)
        _pfam(args)

        # Score-based classification of the candidates (viral vs cellular) and conservation tier
        _classificar_eves(args.outdir, multi=True, ignorar_proteoma=args.ignorar_proteoma)
        for sp_nome, info in dados_triagem.items():
            item = info["item"]
            _cobertura_rnaseq(info["outdir"], args, genome=item["genome"], gff=item["gff"],
                              r1=item["rnaseq_r1"], r2=item["rnaseq_r2"],
                              bams=item.get("rnaseq_bam"), species=sp_nome)

        _varredura_tblastn(args, db_diamond)
        _repeticoes(args)

        # Build the web dashboard automatically, when there is data to show
        try:
            from evetool.visualization.dashboard import gerar_dashboard
            gerar_dashboard(args.outdir,
                            genomas={sp: info["item"]["genome"] for sp, info in dados_triagem.items()},
                            hmm_db=args.hmm_db, hmm_dir=args.hmm_dir, threads=args.threads,
                            pastas_genomas=[p for p in (args.genomes_dir,) if p])
        except Exception as err_dash:
            print(f"[warning] Could not build the web dashboard: {err_dash}")

        # Convenience links at the root, for the positional validation
        _criar_links_conveniencia(args.outdir, {
            "validacao_posicional_eves_sumario.tsv": os.path.join(args.outdir, "02_alinhamento_diamond", "validacao_posicional", "validacao_posicional_eves_sumario.tsv"),
            "validacao_posicional.json": os.path.join(args.outdir, "02_alinhamento_diamond", "validacao_posicional", "validacao_posicional.json"),
            "perfil_cobertura_consenso_posicional.tsv": os.path.join(args.outdir, "02_alinhamento_diamond", "validacao_posicional", "perfil_cobertura_consenso_posicional.tsv"),
            "validacao_posicional_tblastn_reciproco_nao_sobreposto.tsv": os.path.join(args.outdir, "02_alinhamento_diamond", "validacao_posicional", "validacao_posicional_tblastn_reciproco_nao_sobreposto.tsv")
        })

        generate_readme_structure(args.outdir)

        print(f"\n{'='*75}\n[Success] Multi-species pipeline finished.\n{'='*75}")
        print(f"  * Per-species results      : {os.path.join(args.outdir, 'especies')}")
        print(f"  * Master viral protein pool: {pasta_referencias}")
        print(f"  * Cross-species comparison : {pasta_comparativa}")
        print(f"  * Folder layout documented : {os.path.join(args.outdir, 'README_STRUCTURE.md')}")

    # =========================================================================
    # MODE 2: SINGLE-SPECIES ANALYSIS
    # =========================================================================
    else:
        from evetool.discovery.ncbi_genome_downloader import e_acesso_ncbi, baixar_genoma_ncbi
        if e_acesso_ncbi(args.input):
            pasta_genomas = os.path.join(args.outdir, "genomas_ncbi")
            res_ncbi = baixar_genoma_ncbi(args.input, pasta_destino=pasta_genomas, baixar_gff=True)
            args.input = res_ncbi["genome"]
            if not args.gff and res_ncbi.get("gff"):
                args.gff = res_ncbi["gff"]
                print(f"[NCBI Datasets] GFF3 annotation attached automatically: {args.gff}")

        _registrar_info_pipeline(args.outdir, args.input, args.species, args.gff, args.flank, db_hmm, db_diamond)

        screening = executar_triagem_especie(
            genoma_fasta=args.input,
            outdir=args.outdir,
            threads=args.threads,
            args=args,
            db_diamond=db_diamond,
            is_fasta_diamond=is_fasta_diamond,
            db_hmm=db_hmm,
            species=None
        )

        # Viral reference collection
        dir_single = init_pipeline_dirs(args.outdir)
        pasta_candidatas_single = dir_single["04_candidatas_eves"]

        fasta_proteinas_ref = None
        if args.diamond_db_geral:
            try:
                # 1. Run blastx against the general database
                buscar_proteinas_referencia_diamond(
                    screening["candidatas_nt"],
                    args.diamond_db_geral,
                    pasta_candidatas_single,
                    threads=args.threads,
                    hits_por_query=args.diamond_hits_per_query,
                    top_n_total=args.diamond_top_hits if args.diamond_top_hits > 0 else None,
                    export_xml=args.export_xml
                )
                
                # 2. Use the viral and general tables to filter and download the right viral proteins

                from evetool.discovery.sequence_searcher import gerar_referencias_virais_filtradas
                viral_tsv = screening.get("diamond_tsv")
                geral_tsv = os.path.join(pasta_candidatas_single, "diamond_blastx_referencia_geral.tsv")
                out_fasta = os.path.join(pasta_candidatas_single, "referencias_virais_filtradas.fasta")
                
                fasta_proteinas_ref = gerar_referencias_virais_filtradas(
                    viral_tsv, geral_tsv, out_fasta, email_ncbi=args.ncbi_email
                )
            except Exception as erro:
                print(f"[warning] The DIAMOND search against the general database / NCBI download failed: {erro}")
        elif args.validar_ncbi and args.ncbi_email:
            try:
                tabelas_blastx = rodar_blastx_api_em_lotes(
                    screening["candidatas_nt"], pasta_candidatas_single,
                    max_seqs=args.ncbi_lote_max_seqs,
                    max_nt=args.ncbi_lote_max_nt,
                    max_espera_segundos=args.ncbi_timeout
                )
                fasta_proteinas_ref = extrair_melhores_hits_proteicos(
                    tabelas_blastx, pasta_candidatas_single,
                    top_n=args.ncbi_top_hits,
                    email_ncbi=args.ncbi_email
                )
            except Exception as erro:
                print(f"[warning] The NCBI online validation/collection failed: {erro}")

        # Positional and evolutionary validation
        if fasta_proteinas_ref and not args.no_validacao_posicional:
            try:
                print(f"\n{'='*75}\n[Positional validation] Positional analysis and reciprocal tblastn\n{'='*75}")
                dir_posicional_single = dir_single["02_alinhamento_diamond"]
                blastx_tsv_single = os.path.join(pasta_candidatas_single, "diamond_blastx_referencia_geral.tsv")
                executar_analise_posicional_completa(
                    fasta_candidatas_nt=screening["candidatas_nt"],
                    fasta_virais_ref=fasta_proteinas_ref,
                    diretorio_saida=dir_posicional_single,
                    blastx_geral_tsv_ou_xml=blastx_tsv_single if os.path.exists(blastx_tsv_single) else None,
                    diamond_db_geral=args.diamond_db_geral,
                    threads=args.threads,
                    hits_por_query=args.diamond_hits_per_query,
                    export_xml=args.export_xml
                )
            except Exception as erro_pos:
                print(f"[warning] The positional validation failed: {erro_pos}")

        executar_ancoragem_e_contexto_especie(
            genoma_fasta=args.input,
            fasta_candidatas_pep=screening["candidatas_pep"],
            fasta_proteinas_ref=fasta_proteinas_ref,
            outdir=args.outdir,
            threads=args.threads,
            gff_path=args.gff,
            species=None,
            flanco=args.flank,
            gerar_plots=not args.no_plots
        )

        # Extra layer: HMM domains over ALL candidates, not only those without a DIAMOND hit
        try:
            from evetool.annotation.hmm_domains import anotar as anotar_hmm_candidatas
            anotar_hmm_candidatas(args.outdir, hmm_db=db_hmm, threads=args.threads, forcar=True)
        except Exception as err_hmm:
            print(f"[warning] HMM annotation of the candidates failed: {err_hmm}")
        _cdsearch(args)
        _pfam(args)

        # Score-based classification of the candidates (viral vs cellular) and conservation tier
        _classificar_eves(args.outdir, multi=False, ignorar_proteoma=args.ignorar_proteoma)
        _cobertura_rnaseq(args.outdir, args, genome=args.input, gff=args.gff,
                          r1=args.rnaseq_r1, r2=args.rnaseq_r2, bams=args.rnaseq_bam)

        _varredura_tblastn(args, db_diamond)
        _repeticoes(args)

        # Build the web dashboard automatically, when there is data to show
        try:
            from evetool.visualization.dashboard import gerar_dashboard
            print(f"\n{'='*75}\n[Dashboard Web] Compilando Dashboard em: {os.path.join(args.outdir, '12_dashboard_web')}\n{'='*75}")
            gerar_dashboard(args.outdir, genomas={"": args.input},
                            nomes={"": args.species} if args.species else None,
                            hmm_db=args.hmm_db, hmm_dir=args.hmm_dir, threads=args.threads)
        except Exception as err_dash:
            print(f"[warning] Could not build the web dashboard: {err_dash}")
            import traceback
            traceback.print_exc()

        # Convenience links at the root, for the positional validation
        _criar_links_conveniencia(args.outdir, {
            "validacao_posicional_eves_sumario.tsv": os.path.join(dir_single["02_alinhamento_diamond"], "validacao_posicional", "validacao_posicional_eves_sumario.tsv"),
            "validacao_posicional.json": os.path.join(dir_single["02_alinhamento_diamond"], "validacao_posicional", "validacao_posicional.json"),
            "perfil_cobertura_consenso_posicional.tsv": os.path.join(dir_single["02_alinhamento_diamond"], "validacao_posicional", "perfil_cobertura_consenso_posicional.tsv"),
            "validacao_posicional_tblastn_reciproco_nao_sobreposto.tsv": os.path.join(dir_single["02_alinhamento_diamond"], "validacao_posicional", "validacao_posicional_tblastn_reciproco_nao_sobreposto.tsv")
        })

        generate_readme_structure(args.outdir)
        print(f"\n[run] Pipeline finished for the single genome. Results organised in: {args.outdir}")
        print(f"  * Folder layout documented in: {os.path.join(args.outdir, 'README_STRUCTURE.md')}")


if __name__ == "__main__":
    main()
