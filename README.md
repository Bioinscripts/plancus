<p align="center"><img src="docs/logo.jpeg" alt="PLANCUS logo" width="460"></p>

# EVE tool

Discovery, validation, curation and cross-species comparison of **endogenous viral elements (EVEs)** in eukaryotic genomes.

The pipeline goes from a genome (a FASTA file or an NCBI `GCA_`/`GCF_` accession) to a curated EVE catalogue:

1. **ORF screening.** ORFs are predicted with orfipy, deduplicated with MMseqs2 and screened with DIAMOND blastx against viral proteins.
2. **HMM rescue.** ORFs with no viral hit are searched against calibrated RVDB + eggNOG viral HMMs (pyhmmer).
3. **Validation.** Each candidate is compared against Swiss-Prot. Viral and cellular evidence are weighed against each other, and the viral family comes from the NCBI taxonomy.
4. **Domains.** Optional CD-Search (local CDD) and Pfam-A annotation.
5. **Genomic context.** Candidates are anchored in the genome with tblastn and grouped into loci with flanks. A last tblastn scan recovers fragments that have no ORF.
6. **Repeats.** Optional: low complexity (DUST/SEG), soft-masked repeats and transposon domains around each EVE.
7. **Cross-species.** Two-stage tblastx clustering, a heatmap of the supported sharing and flank synteny.
8. **Curation.** A local web page to inspect, approve or remove candidates.

## Installation

The tool needs Linux (or macOS) and conda, mamba or micromamba.

```bash
git clone https://github.com/<user>/<repo>.git
cd <repo>
conda env create -f environment.yml
conda activate eve_env
```

### Databases

By default the databases go to `~/eve_data`. To use another folder, copy `eve_config.example.ini` to `eve_config.ini` and set `data_dir`, or export `EVE_DATA_DIR=/path`.

```bash
./eve-tool databases                  # essentials, ~1.5 GB, ~10 min
./eve-tool databases --pfam --cdd     # optional: Pfam-A (~1.5 GB) and the local CDD (~4.5 GB)
./eve-tool check                      # checks tools, libraries and databases
```

| Part | Source | Used for |
|---|---|---|
| viral bundle | this project ([Zenodo DOI]) | viral proteins (DIAMOND), calibrated viral HMMs, annotations |
| Swiss-Prot | UniProt | cellular reference for the validation |
| taxdump | NCBI Taxonomy | viral family of each hit |
| proteome | NCBI (default *P. papatasi* GCF_024763615.1; `--proteome` picks another) | labels of the flanking genes in the synteny step |
| Pfam-A | EBI (optional) | domains and transposon detection |
| CDD | NCBI (optional) | CD-Search |

### Updating the viral databases (optional)

The viral bundle is a fixed release: cite its version, and runs made with it are reproducible. To follow newer RVDB releases, rebuild the HMMs with the same criteria (`evetool/databases/recipe_viral_db.json`):

```bash
./eve-tool build-databases -t 8          # latest RVDB + eggNOG 5.0; ~1 day (calibration)
./eve-tool databases --bundle <build folder>/eve_tool_db_<date>.tar --force
```

The recipe keeps:

- the RVDB models of viral families without dsDNA, dsDNA-RT and ssRNA-RT, plus the models with no family whose seed sequences are RNA viruses;
- the eggNOG RNA (and ssDNA) virus groups;
- the same HMM calibration (Swiss-Prot and six dipteran proteomes as negatives).

The RefSeq viral proteins and Swiss-Prot are not downloaded by the build. It uses the installed files, or the ones given with `--refseq-fasta` and `--swissprot`. The RefSeq file needs the NCBI Virus defline `accession |title|molecule type|species|genus|family|length`; a file without it is refused. [docs/refseq_viral_proteins.md](docs/refseq_viral_proteins.md) shows, step by step, how to download it from NCBI Virus.

The build is resumable: an interrupted run continues where it stopped. `build_report.md` compares the new database with the installed one: models added, removed or re-calibrated, and proteins per genome type. Every run records the database version in `pipeline.log` and `pipeline_info.json`. A newer database finds a somewhat different set of EVEs.

## Quick test

The `examples/` folder has three 50–60 kb windows of the *Lutzomyia longipalpis* genome (RefSeq GCF_024334085.1) that carry rhabdovirus and phasmavirus EVEs. The run takes about 3 minutes with 8 threads.

```bash
./run_pipeline.sh example -i examples/llongipalpis_3_windows.fna -t 8
./open_curation.sh ~/eve_data/Resultados/example
```

## Usage

```bash
# One genome (FASTA, with an optional GFF3 annotation)
./run_pipeline.sh my_run -i genome.fna --gff genome.gff3

# Several species, with the cross-species comparison
./run_pipeline.sh sandflies --species-tsv species.tsv --sintenia --repeats
```

`run_pipeline.sh` writes each run to `<results_dir>/<run_name>/`, keeps a `pipeline.log` and turns on CD-Search and Pfam when they are installed. `./eve-tool run --help` lists every option. The databases default to the installed ones and can be replaced with `--dmnd-db`, `--hmm-db` and `--diamond-db-geral`.

`species.tsv` (tab-separated). `genome` takes a local FASTA or a `GCA_`/`GCF_` accession, which is downloaded together with its GFF3:

```
species	genome	gff	rnaseq_r1	rnaseq_r2
Phlebotomus_papatasi	GCF_024763615.1	None	None	None
Lutzomyia_longipalpis	L_longipalpis.fna	L_longipalpis.gff3	R1.fq.gz	R2.fq.gz
```

## Curation page

```bash
./open_curation.sh <results folder>
```

The page runs on a local Python server (standard library only) and opens in the browser. It has tabs for the pipeline, an overview, per-EVE inspection (domains, ORFs and repeats in context), clusters, tblastn segments and the cross-species comparison. Curation decisions are saved in `<results>/13_curadoria/`.

When a run has RNA-seq (`--rnaseq-r1/--rnaseq-r2` or BAMs in the species TSV), the inspection view, the expanded context clusters and the tblastn segments show the **read coverage** under their maps: depth per bin, split by strand for stranded libraries, with the exonic and intron-spanning reads of the EVE. Other BAMs (aligned to the genome or to the pipeline's context windows) can be added on the page with **add BAM**. In runs on extracted segments, the reads are read on the original contig.

## Repository layout

```
eve-tool                command line: ./eve-tool run | curate | databases | check | config
run_pipeline.sh         runs the pipeline into <results_dir>/<run_name>/ with a log
open_curation.sh        opens the curation page for a results folder
evetool/
  main.py               pipeline orchestration (./eve-tool run)
  config.py             database and results paths
  check.py              installation check (./eve-tool check)
  core/                 shared helpers and the results folder layout
  discovery/            ORFs, MMseqs2, DIAMOND and HMM screening, tblastn anchoring, loci, tblastn scan
  validation/           viral vs cellular classification, positional validation, viral protein map
  annotation/           CD-Search, viral HMMs, Pfam, repeats and low complexity, RNA-seq coverage
  visualization/        locus figures and the static dashboard
  comparative/          cross-species clustering, hierarchical clustering, flank synteny
  curation/             curation server and web page (web/)
  databases/            database download, packing, audit and HMM calibration
examples/               small public test genome
docs/                   guides (downloading the RefSeq viral proteins from NCBI Virus)
```

Each module can also be run on its own, for example `python -m evetool.annotation.repeats -h` (with the repository folder on `PYTHONPATH`, or after `pip install -e .`, which also puts `eve-tool` on the PATH).

## Configuration

Each path is resolved in this order:

1. an environment variable (`EVE_DATA_DIR`, `EVE_HMM_DIR`, `EVE_DIAMOND_DIR`, `EVE_RESULTS_DIR`, `EVE_PFAM_HMM`);
2. `eve_config.ini`;
3. the defaults under `data_dir`.

`./eve-tool config` prints the paths in use.

## Citation

If you use this tool, please cite: _in preparation_. Please also cite the underlying tools and databases: DIAMOND, MMseqs2, HMMER/pyhmmer, BLAST+, orfipy, RVDB, eggNOG, Pfam, CDD and UniProt.

## License

The code is released under the [MIT License](LICENSE). The viral database bundle on Zenodo has its own license and derives from RVDB, eggNOG and NCBI RefSeq; cite them too.
