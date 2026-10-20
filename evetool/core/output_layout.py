#!/usr/bin/env python3
"""
output_layout.py
==========================
Governance and standardisation of the EVE_tool directory structure.
It isolates the outputs of each analysis phase in dedicated numbered folders:

  01_orfs_cds/             -> ORF and CDS prediction (orfipy, mmseqs2 deduplication)
  02_alinhamento_diamond/  -> Primary and secondary DIAMOND blastx screens
  03_triagem_hmm/          -> Deep HMMER screen for divergent viral scars
  04_candidatas_eves/      -> Consolidated EVE sequences before anchoring (NT/PEP)
  05_ancoragem_tblastn/    -> Genome anchoring with makeblastdb and tBLASTn
  06_anotacao_genomica/    -> Genomic coordinates, GFF3 annotation, BED and 20 kb context
  07_alinhamento_tblastx/  -> Reciprocal all-vs-all tBLASTx alignments (XML and tabular)
  08_clusters_e_tensores/  -> 3D CP/Tucker tensors, embeddings, linkage and clusters
  09_graficos_e_figuras/   -> Dendrograms, heatmaps and 20 kb genomic-context figures
  10_rnaseq_expressao/     -> RNA-seq mapping, BAMs, depth and quantification
  11_sintenia_multiespecies/ -> Multi-species comparison, synteny blocks and clinker
  12_dashboard_web/        -> Interactive web dashboard with dynamic tabs and JSON data
"""

import os
import sys
import re
import shutil
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple, List

# Standard folder structure of the pipeline
PIPELINE_DIRS = {
    "01_orfs_cds": {
        "alias": "cds",
        "title": "ORF/CDS prediction and extraction",
        "description": "ORF prediction with orfipy and redundancy removal with MMseqs2",
        "patterns": [
            r".*_orfs_nt\.fasta$",
            r".*_orfs_pep\.fasta$",
            r".*_orfs_pep_nr_.*\.fasta$",
            r".*_orfs_pep_nr_cluster\.tsv$",
            r".*tmp_mmseqs.*",
            r".*orfipy.*\.log$"
        ]
    },
    "02_alinhamento_diamond": {
        "alias": "diamond",
        "title": "DIAMOND screening and alignments",
        "description": "blastx alignments against reference viral databases and comparative positional validation",
        "patterns": [
            r".*diamond.*\.tsv$",
            r".*diamond.*\.xml$",
            r".*triagem_viral_diamond.*\.tsv$",
            r".*validacao_posicional.*",
            r".*perfil_cobertura_consenso.*",
            r".*tblastn_reciproco.*",
            r".*banco_customizado.*"
        ]
    },
    "03_triagem_hmm": {
        "alias": "hmm",
        "title": "Deep HMMER screen",
        "description": "Search for divergent viral scars with HMM profiles (hmmsearch)",
        "patterns": [
            r".*triagem_profunda_hmm\.tsv$",
            r".*orfs_sem_hit_viral.*\.fasta$",
            r".*banco_viral_customizado\.hmm.*"
        ]
    },
    "04_candidatas_eves": {
        "alias": "candidates",
        "title": "Candidate EVE sequences",
        "description": "Consolidated FASTAs of candidate EVEs (nucleotides and peptides) and viral references",
        "patterns": [
            r".*eve_candidatas_unificadas.*\.fasta$",
            r".*todas_candidatas_.*\.fasta$",
            r".*proteinas_referencia.*\.fasta$"
        ]
    },
    "05_ancoragem_tblastn": {
        "alias": "tblastn",
        "title": "tBLASTn genome anchoring",
        "description": "Host BLAST databases and tBLASTn hit-anchoring tables",
        "patterns": [
            r".*_blastdb\.(ndb|nhr|nin|not|nsq|ntf|nto|pdb|phr|pin|pot|psq|ptf|pto)$",
            r".*eve_tblastn_hits.*\.tsv$"
        ]
    },
    "06_anotacao_genomica": {
        "alias": "annotation",
        "title": "Genomic annotation and context",
        "description": "GFF3 and BED annotation files, EVE loci and sequences with 20 kb flanks",
        "patterns": [
            r".*eve_clusters.*\.gff3?$",
            r".*eve_clusters.*\.bed$",
            r".*eve_hits.*\.bed$",
            r".*eve_clusters_context\.fasta$",
            r".*eve_loci\.tsv$",
            r".*cromossomos\.sizes$"
        ]
    },
    "07_alinhamento_tblastx": {
        "alias": "tblastx",
        "title": "All-vs-all tBLASTx alignment",
        "description": "Reciprocal tBLASTx alignments (XML and tabular), non-overlapping HSPs and metric tables",
        "patterns": [
            r".*tblastx.*\.xml$",
            r".*tblastx.*\.tsv$",
            r".*_pares_metricas\.tsv$",
            r".*_pares\.tsv$",
            r".*_tabela_metricas\.tsv$",
            r".*db_flebo\..*",
            r".*tblastx_db\..*"
        ]
    },
    "08_clusters_e_tensores": {
        "alias": "clusters",
        "title": "Cluster and tensor results",
        "description": "3D CP/Tucker tensors, embeddings, metric loadings, linkage and cluster assignment",
        "patterns": [
            r".*_clusters\.tsv$",
            r".*_tensor_embedding\.tsv$",
            r".*_tensor_carga_metricas\.tsv$",
            r".*_segmentos_classificacao\.tsv$",
            r".*_componentes_resumo\.tsv$",
            r".*_linkage\.tsv$",
            r".*_tabela_ponderada\.tsv$",
            r".*_dissimilaridade\.tsv$",
            r".*_clusters_tolerancia\.tsv$",
            r".*_vizinhos_tolerancia\.tsv$",
            r".*_varredura_tolerancia\.tsv$",
            r".*_mcl_.*\.tsv$",
            r".*_limiar_.*\.tsv$",
            r".*_comp_feats_.*\.tsv$",
            r".*bolhas.*\.tsv$",
            r".*cone.*\.tsv$"
        ]
    },
    "09_graficos_e_figuras": {
        "alias": "plots",
        "title": "Charts and figures",
        "description": "Dendrograms, heatmaps, tolerance sweep and genomic-context figures (PNG/PDF)",
        "patterns": [
            r".*_dendrograma\.(png|pdf|svg)$",
            r".*_heatmap\.(png|pdf|svg)$",
            r".*_varredura_tolerancia\.(png|pdf|svg)$",
            r".*_limiar_grafico\.(png|pdf|svg)$",
            r".*_contexto_20kb\.(png|pdf|svg)$",
            r".*\.(png|pdf|svg)$"
        ]
    },
    "10_rnaseq_expressao": {
        "alias": "rnaseq",
        "title": "RNA-seq expression and coverage",
        "description": "BAM/BAI alignments, samtools depth profiles and RPKM quantification",
        "patterns": [
            r".*\.bam$",
            r".*\.bai$",
            r".*_depth\.tsv$",
            r".*cobertura.*\.tsv$",
            r".*rnaseq_analise.*"
        ]
    },
    "11_sintenia_multiespecies": {
        "alias": "synteny",
        "title": "Synteny and multi-species comparison",
        "description": "Homologous clusters across species, clinker reports, GenBank (.gbk) and phylogenetic trees",
        "patterns": [
            r".*comparacao_multiespecies.*",
            r".*clinker.*",
            r".*sintenia.*",
            r".*\.gbk$",
            r".*alinhamento.*\.aln$"
        ]
    },
    "12_dashboard_web": {
        "alias": "dashboard",
        "title": "Interactive web dashboard",
        "description": "Interactive web dashboard with dynamic tabs, JSON data and a local web server",
        "patterns": [
            r".*index\.html$",
            r".*server\.py$",
            r".*webpage.*"
        ]
    }
}


def init_pipeline_dirs(base_dir: str) -> Dict[str, str]:
    """
    Cria fisicamente todas as pastas estruturadas do pipeline dentro de base_dir
    and returns a dictionary mapping both the folder name (e.g. '01_orfs_cds') 
    quanto o alias curto (ex: 'cds') para o caminho absoluto correspondente.
    """
    base_dir = os.path.abspath(base_dir)
    os.makedirs(base_dir, exist_ok=True)
    
    dir_map = {"root": base_dir}
    for folder_name, info in PIPELINE_DIRS.items():
        folder_path = os.path.join(base_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)
        dir_map[folder_name] = folder_path
        dir_map[info["alias"]] = folder_path
        
    return dir_map


def classify_file(filename: str) -> Optional[str]:
    """
    Classifies a file by its name/pattern and returns the matching destination folder
    (ex: '01_orfs_cds', '08_clusters_e_tensores', etc.).
    """
    fname = os.path.basename(filename)
    
    # Priority order for disambiguation
    priority_order = [
        "01_orfs_cds",
        "03_triagem_hmm",
        "04_candidatas_eves",
        "02_alinhamento_diamond",
        "05_ancoragem_tblastn",
        "06_anotacao_genomica",
        "07_alinhamento_tblastx",
        "08_clusters_e_tensores",
        "09_graficos_e_figuras",
        "10_rnaseq_expressao",
        "11_sintenia_multiespecies",
        "12_dashboard_web"
    ]
    
    for folder_name in priority_order:
        patterns = PIPELINE_DIRS[folder_name]["patterns"]
        for pat in patterns:
            if re.match(pat, fname, re.IGNORECASE):
                return folder_name
                
    return None


def organize_directory(
    source_dir: str, 
    target_dir: Optional[str] = None, 
    mode: str = "link",
    verbose: bool = True
) -> Dict[str, List[str]]:
    """
    Scans the files in source_dir and sorts them into the pipeline folders.
    
    Parameters:
      source_dir: input directory with raw / unsorted files.
      target_dir: structured destination. When None, the folders are created inside source_dir.
      mode: 'link' (relative symlinks, saving space), 'copy' or 'move'.
      verbose: print each step.
      
    Retorna:
      Dictionary with the list of files allocated to each folder.
    """
    source_dir = os.path.abspath(source_dir)
    target_dir = os.path.abspath(target_dir) if target_dir else source_dir
    
    dir_map = init_pipeline_dirs(target_dir)
    allocated: Dict[str, List[str]] = {folder: [] for folder in PIPELINE_DIRS}
    allocated["unclassified"] = []
    
    items = sorted(os.listdir(source_dir))
    for item in items:
        # Skip the structured folders themselves when they already exist
        if item in PIPELINE_DIRS or item in ("webpage", "Web page"):
            continue
            
        src_path = os.path.join(source_dir, item)
        
        # Special case for the context-figures directory
        if item == "figuras_contexto_20kb" and os.path.isdir(src_path):
            dst_folder = os.path.join(dir_map["09_graficos_e_figuras"], "figuras_contexto_20kb")
            if not os.path.exists(dst_folder):
                if mode == "move":
                    shutil.move(src_path, dst_folder)
                elif mode == "copy":
                    shutil.copytree(src_path, dst_folder)
                else: # link
                    os.symlink(os.path.relpath(src_path, os.path.dirname(dst_folder)), dst_folder)
            allocated["09_graficos_e_figuras"].append(item)
            continue

        # Special case for the positional-validation folder
        if item == "validacao_posicional" and os.path.isdir(src_path):
            dst_folder = os.path.join(dir_map["02_alinhamento_diamond"], "validacao_posicional")
            if not os.path.exists(dst_folder):
                if mode == "move":
                    shutil.move(src_path, dst_folder)
                elif mode == "copy":
                    shutil.copytree(src_path, dst_folder)
                else: # link
                    os.symlink(os.path.relpath(src_path, os.path.dirname(dst_folder)), dst_folder)
            allocated["02_alinhamento_diamond"].append(item)
            continue
            
        # Caso especial para pasta de rnaseq_analise
        if item == "rnaseq_analise" and os.path.isdir(src_path):
            dst_folder = dir_map["10_rnaseq_expressao"]
            # Copy or link the content
            for sub in os.listdir(src_path):
                s_sub = os.path.join(src_path, sub)
                d_sub = os.path.join(dst_folder, sub)
                if not os.path.exists(d_sub):
                    if mode == "move":
                        shutil.move(s_sub, d_sub)
                    elif mode == "copy":
                        if os.path.isdir(s_sub):
                            shutil.copytree(s_sub, d_sub)
                        else:
                            shutil.copy2(s_sub, d_sub)
                    else:
                        os.symlink(os.path.relpath(s_sub, os.path.dirname(d_sub)), d_sub)
            allocated["10_rnaseq_expressao"].append(item)
            continue
            
        # Caso especial para tmp_mmseqs
        if item == "tmp_mmseqs" and os.path.isdir(src_path):
            dst_folder = os.path.join(dir_map["01_orfs_cds"], "tmp_mmseqs")
            if not os.path.exists(dst_folder):
                if mode == "move":
                    shutil.move(src_path, dst_folder)
                elif mode == "copy":
                    shutil.copytree(src_path, dst_folder)
                else:
                    os.symlink(os.path.relpath(src_path, os.path.dirname(dst_folder)), dst_folder)
            allocated["01_orfs_cds"].append(item)
            continue
            
        if os.path.isfile(src_path):
            folder_target = classify_file(item)
            if folder_target:
                dst_path = os.path.join(dir_map[folder_target], item)
                if not os.path.exists(dst_path):
                    if mode == "move":
                        shutil.move(src_path, dst_path)
                    elif mode == "copy":
                        shutil.copy2(src_path, dst_path)
                    else: # link
                        os.symlink(os.path.relpath(src_path, os.path.dirname(dst_path)), dst_path)
                allocated[folder_target].append(item)
            else:
                allocated["unclassified"].append(item)

    # Write the explanatory README_STRUCTURE.md
    generate_readme_structure(target_dir, allocated)
    
    if verbose:
        print(f"\n{'='*75}\n[EVE_tool file organiser] Directory organised at: {target_dir}\n{'='*75}")
        for folder_name, files in allocated.items():
            if files:
                title = PIPELINE_DIRS.get(folder_name, {}).get("title", folder_name)
                print(f"  * {folder_name}/ ({len(files)} arquivos) -> {title}")
        if allocated["unclassified"]:
            print(f"  [!] {len(allocated['unclassified'])} unclassified files kept at the root.")
        print(f"{'='*75}\n")
        
    return allocated


def generate_readme_structure(target_dir: str, allocated: Optional[Dict[str, List[str]]] = None):
    """Writes a README_STRUCTURE.md describing the purpose of each folder and its files."""
    readme_path = os.path.join(target_dir, "README_STRUCTURE.md")
    
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write("# EVE_tool file and directory structure\n\n")
        f.write("This directory holds the results of the endogenous viral element (EVE) mining and characterisation pipeline, organised by analysis step.\n\n")
        f.write("## Folder summary\n\n")
        
        for folder_name, info in PIPELINE_DIRS.items():
            f.write(f"### `{folder_name}/` — {info['title']}\n")
            f.write(f"**Description:** {info['description']}\n\n")
            
            # Check what is actually on disk
            folder_real_path = os.path.join(target_dir, folder_name)
            files_in_folder = []
            if os.path.exists(folder_real_path) and os.path.isdir(folder_real_path):
                files_in_folder = sorted(os.listdir(folder_real_path))
                
            if files_in_folder:
                f.write(f"**Files present ({len(files_in_folder)} items):**\n")
                for fname in files_in_folder[:12]:
                    f.write(f"- `{fname}`\n")
                if len(files_in_folder) > 12:
                    f.write(f"- *... e mais {len(files_in_folder) - 12} arquivos*\n")
            else:
                f.write("*No file produced at this step in this run.*\n")
            f.write("\n")
            
        f.write("---\n")
        f.write("*Gerado automaticamente pelo `output_layout.py`.*\n")


def main():
    parser = argparse.ArgumentParser(
        description="EVE_tool folder organiser: sorts files by analysis type."
    )
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory holding the files")
    parser.add_argument("-o", "--output-dir", help="Destination directory (optional; organises in place when omitted)")
    parser.add_argument("-m", "--mode", choices=["link", "copy", "move"], default="link",
                        help="Organisation mode: 'link' (symlinks, no extra space), 'copy' or 'move'")
    parser.add_argument("--init-only", action="store_true", help="Apenas inicializa as pastas sem mover arquivos")
    
    args = parser.parse_args()
    
    if args.init_only:
        out_dir = args.output_dir or args.input_dir
        dirs = init_pipeline_dirs(out_dir)
        print(f"[Sucesso] 12 pastas estruturadas inicializadas em: {out_dir}")
        for k in sorted([d for d in dirs.keys() if d.startswith("0") or d.startswith("1")]):
            print(f"  - {dirs[k]}")
    else:
        organize_directory(args.input_dir, args.output_dir, mode=args.mode)


if __name__ == "__main__":
    main()
