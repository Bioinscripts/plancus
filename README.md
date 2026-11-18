<p align="center"><img src="docs/logo.jpeg" alt="PLANCUS logo" width="460"></p>

# PLANCUS

Discovery, validation, curation and cross-species comparison of **endogenous viral elements (EVEs)** in eukaryotic genomes. Named after *Caracara plancus*, the Brazilian carcará — a scavenger with a sharp eye for what others leave behind.

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

The run downloads the reference viral proteins from the NCBI, whose usage policy asks for a contact e-mail. Set yours once, in `eve_config.ini`:

```ini
[ncbi]
email = you@institution.org
```

or with `export EVE_NCBI_EMAIL=you@institution.org`, or per run with `--ncbi-email`. A run without it stops at the start with a message.

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

The page runs on a local Python server (standard library only) and opens in the browser. It has tabs for the pipeline, an overview, per-EVE inspection (domains, ORFs and repeats in context), clusters, tblastn segments, the cross-species comparison and the composition of the EVEs next to circulating viruses. Curation decisions are saved in `<results>/13_curadoria/`.

The **Circulating viruses** tab compares the composition of the approved EVEs with coding sequences (CDS) of circulating viruses that you add on the page (FASTA; the family is read from an NCBI Virus header, a `family=` tag, or assigned with DIAMOND), optionally with host CDS as a reference. Each EVE is read in the frame of its viral protein (the stretch DIAMOND aligned, in the frame that matches the reference), not in its ORF frame. The tab shows %GC and GC1/2/3/GC3s, codon usage (RSCU, ENC, ENC and neutrality plots), dinucleotide odds ratios over the whole sequence and at codon positions 1-2, 2-3 and 3-1, and a PCA and t-SNE on the features that best separate the groups (Kruskal-Wallis, redundant ones skipped). With the genome annotation (the run's GFF, a GFF next to the genome, a third column `gff` in `synteny_genomes.tsv`, or one set on the page), the CDS of the host genes are read too: the genes within ±50 kb of an EVE are its neighbours, a sample of the rest is the host genome. Each EVE is then placed on the axis from the circulating viruses of its family (0) to its neighbouring host genes or the host genome (1), and shown EVE by EVE against them. A host gene model that covers an EVE is left out and reported: it is probably the EVE itself, annotated as a host gene. Points are coloured by viral family, drawn by origin (EVE, circulating virus, host gene near an EVE, host genome gene) and filled when they carry a viral HMM domain; clicking an EVE links it to its neighbouring genes in the PCA and t-SNE. The same analysis runs from the command line: `python -m evetool.comparative.codon_composition -o <results> --add cds.fasta --name <set>`, `--gff <species> <file.gff3>` and `--run`.

**Verdict.** Every candidate gets two labels, kept apart on purpose (`evetool/validation/verdict.py`). *Detection* says how it was found: which of the three independent lines of viral evidence it has (alignment to a viral protein, viral HMM, RNA/ssDNA-virus conserved domain). *Verdict* says what it most likely is, from explicit rules shown with their reasons for and against: **EVE (high confidence)**, **EVE (probable)**, **Ambiguous**, **Possible ancient insertion**, **No viral evidence**, **Host gene** (a label only, see *Host genes* below), and the two explicitly non-viral ones, **Transposon** and **Out of scope**. The scope is non-retroviral RNA viruses and ssDNA viruses of eukaryotes: retroelements, pararetroviruses, dsDNA viruses and bacteriophages are out. Mapped reads add to the verdict when a species has them (`evetool/validation/read_evidence.py`): DNA reads of the assembly that confirm the integration, and piRNA or siRNA production, lift a clean *probable* to *high confidence* (both together lift an *Ambiguous* to *probable*), while assembly reads that do not cover the locus, pile up far above its flanks or cross no junction count against it; reads never remove a candidate. Only the explicitly non-viral loci are moved to the Removed box by themselves, once, and restoring one keeps it; everything else is the researcher's call. A candidate whose cellular hit wins is *not* removed when that hit is an uncharacterised protein: an old EVE annotated as a host gene looks exactly like that.

**Genomic DNA.** The *Genomic DNA* tab maps DNA reads on the whole genome (`minimap2 -ax sr`, or the long-read presets) and keeps only the alignments over the EVE loci ± 5 kb, so a whole-genome BAM is never written (`evetool/annotation/dna_reads.py`). Each library is declared as the **reads of the assembly** or a **population sample**, and as an **individual** or a **pool**. Assembly reads test the integration: an integrated EVE is covered like its flanks and the genome, with reads across both junctions; a large excess of depth points to an exogenous virus or a contaminant. Population samples give a loci × samples matrix: present / partial / absent for an individual, an estimated frequency (EVE depth over flank depth) for a pool, and "no data" when the flanks themselves are not covered. Only measurements are stored (depths, breadth, junction reads, and reads clipped at the junction, which come from chromosomes without the insertion); the calls are made on the page from thresholds the researcher changes while looking at the matrix. From the command line: `python -m evetool.annotation.dna_reads -o <results> --species <name> --r1 R1.fq.gz --r2 R2.fq.gz --role population --sample-type pool --group <population>`. Many libraries can be mapped and measured elsewhere (e.g. on a cluster, with the same `_medir_locus`) and brought back as measurements only, without reads or BAMs: `dna_reads.importar_medidas(outdir, species, measurements.tsv, libraries.json)` registers them (marked *precomputed*); a precomputed value is used only while the locus keeps the coordinates it was measured on.

**Speed.** The page redraws only the charts of the visible tab (the others are redrawn when shown), builds the option lists of the row selects when one is opened, does not paint the saved snapshot when a server is answering, and an inspection tab loads neither the snapshot nor Plotly.

**Curation record in the inspection.** The inspection tab shows the row of the candidate table for the EVE it opened (species, state, group, conservation, viral family, protein class, scores, CDD), and viral family and protein class can be changed there: the change is applied in the curation tab, like the Qualify / Reject decision, and is kept when that tab is saved.

**Contig ends.** The inspection marks where the contig ends: a thick bar on the neighbourhood track (hatched beyond it; the track zooms to a contig shorter than itself), on the GC profile and on the cluster context, and the header gives the distance from the EVE to each end, with the tag *at the contig end* when it is 1 kb or less. A run made on extracted segments is placed back on the whole genome first (`synteny_genomes.tsv`), since the end of a segment is not the end of a contig; where the segment stops and the contig goes on, the drawing says so.

**Simple sequence.** Poly-A, poly-anything (a run of one base of 10 nt or more) and short tandem repeats (a unit of 2–6 nt repeated at least 4 times, 16 nt or more) inside a candidate segment are flagged in the table (`poly-N` / `repeat`), and listed against it in the verdict. From 5 % of the segment, or a single homopolymer of 20 nt, they hold *high confidence* back to *probable*; a segment that is 30 % or more simple sequence is held back to *Ambiguous* (`evetool/validation/low_complexity.py`, measured on the candidate sequence itself). When the repeats step ran, 20 % or more low complexity (DUST) and 50 % or more soft-masked repeat count against it too. Nothing is removed for this.

**Host genes.** A host gene is not always a non-EVE (an old insertion is often annotated as a gene of the host), so the page never removes one by itself. The *Host genes* tab searches the peptide of every candidate against the installed host proteomes (DIAMOND blastp, best hit per ORF in each proteome, isoforms collapsed; `evetool/validation/host_genes.py`, also `python -m evetool.validation.host_genes -o <results>`) and lists the genes that were hit with their candidates and three flags: *uncharacterised* (no described function: possibly the EVE itself), *same species* (the candidate's own genome annotation contains the ORF) and *host ≥ viral*. Removing a gene rejects every candidate ORF that hit it, with its whole locus; restoring the gene brings them back. **Removed host genes (TSV)** is the list for the methods (gene, description, organism, proteome, the ORFs that hit it and the ORFs rejected with it), also written to `13_curadoria/host_genes/removed_host_genes.tsv` on save. In the verdict, *Host gene* is now only a label; a characterised host hit as strong as the viral one holds *high confidence* back to *probable*.

**Qualifying and neighbours.** Every candidate starts approved; an EVE the researcher has inspected and kept is **qualified** (the `✓` button of the table, or *Qualify* / key `A` in the inspection tab) and turns light green in the table and on the maps. The state is saved with the selection (`qualified` in `curation_selection.json`, a `qualified` column in the exported tables) and is dropped when the EVE is rejected. The inspection tab shows the other EVEs within 25 kb as a small map zoomed to them, with the context cluster underlined and each block coloured by its state (qualified, not decided yet, rejected): hovering a block or its table row shows the EVE, a click (or `←` / `→`) opens it in the same tab, and the tab stays open after a decision while neighbours remain.

**Removing.** Rejecting a candidate takes its whole locus (every overlapping reading frame) out of the table and into a folded **Removed** box under it, so the tables hold only what is still being analysed. The Clusters and tblastn tabs have a **remove** button on each cluster, segment and tblastn segment: it rejects the EVEs in the candidate table too, and **restore** in the Removed box of that tab brings them back together. Nothing is lost until the curation is saved, and a removal can always be restored.

**GFF3.** Four exports, in the coordinates of the contigs of the run: the tblastn segments (tblastn tab: one `match` per segment, its HSPs as `match_part` with the viral protein as `Target`), the genomic-context clusters (Clusters tab: one `region` per cluster, its segments as children), the unique EVEs (Inspection tab: one `viral_sequence` per approved locus, with family, protein class, best viral hit and cluster) and **every ORF, nested left out** (Inspection tab, and *Export selected* in the Sequences tab of `/v2`). That last one holds one `ORF` feature per approved ORF, under one `viral_sequence` per locus. An ORF lying entirely inside another ORF of the list on the same contig, in any frame or strand, is left out, so the overlapping frames of one insertion do not pile up. When two ORFs have the same span, the one with the best viral hit stays.

**Distance to the nearest EVE.** The Clusters tab answers *how far apart are the EVEs in the genome, and is there a pattern?* (`evetool/comparative/eve_spacing.py`).
- **Measure.** For each approved locus, it takes the gap to the nearest other approved locus of the same genome. The gap is measured on the original chromosome; a run on extracted segments is placed back on the whole genome first.
- **Comparison.** The same loci are dropped at random on the assembly 200 times: a contig is drawn in proportion to its length and a position uniformly on it.
- **Chart.** A ridgeline with one species per row: the density of the gaps on a log scale, with one tick per EVE and the random placement dashed.
- **Statistics.** It reports the median gap observed and at random, the share within 10 kb, p (closer together / more spread than random) and a Kolmogorov–Smirnov test.
- **Merge within N bp.** This option first joins loci closer than N bp into one insertion (the ORFs of one insertion broken by stop codons), so the gaps are between insertions.
- **Species left out.** A species of a segment run without its whole genome in `synteny_genomes.tsv` is reported as not measurable.
- **Command line.** `python -m evetool.comparative.eve_spacing -o <results> [--merge 1000] --tsv`.

**Primers.** The inspection view suggests primers with Primer3 (`primer3-py`): an **internal** pair, both primers inside the inserted sequence, and a **junction** pair for each side, one primer inside the EVE and one in the host flank, which gives a product only if the viral sequence is integrated at that point of the genome. *Primers of the approved EVEs (TSV)* designs them for every approved locus; from the command line, `python -m evetool.curation.primers -o <results> --all`. Each pair is checked against its own template only, so confirm the specificity on the whole genome before ordering.

**Reads.** The *Read mapping* card of the Inspection tab maps reads on the genome of a species, in the background, and the inspection view, the expanded context clusters and the tblastn segments then show them under their maps. In runs on extracted segments, the reads are mapped and read on the original contigs.

- **Small RNA** follows [sRNAProfiler](https://github.com/LymF/srnaprofiler). Give the reads with the adapters already removed. They are mapped with Bowtie1 allowing at most one mismatch (`-v 1 -k 20 --best --strata`), optionally after dropping the reads that match a Bowtie index of structural RNA, and a read placed at *N* positions counts 1/*N* at each. For every span, the reads with at least half their length inside it give the size profile (18-35 nt), strand, 5' U, A at position 10, the 10-nt ping-pong overlap (fraction and z-score), phasing and the class of sRNAProfiler's rule set (siRNA or piRNA, with its confidence and `review` flag). *Small RNA features of the approved EVEs (TSV)* writes them for every approved locus; from the command line, `python -m evetool.annotation.srna -o <results> --species Sp --reads trimmed.fq.gz` and `--table`.
- **Long and other RNA** follow [ViralQuest](https://github.com/gabrielvpina/viralquest): `minimap2 -ax <preset> | samtools sort` (`sr` for Illumina, the splice presets for Nanopore and PacBio), on the whole genome or on the EVE windows, and its coverage track: bar height is depth, bar colour is the mean base quality (green Q>=30, amber Q20-29, red Q<20), + strand above and - strand below for stranded libraries, with mean depth, breadth, CV and mean quality, and the internal stretches below 20% of the mean depth underlined. `minimap2 -ax sr` does not open introns, so a short read crossing one is clipped at the junction. The pipeline's own RNA-seq step (`--rnaseq-r1/--rnaseq-r2`, or BAMs in the species TSV) maps the same way.

BAMs made elsewhere can be added as they are, as small RNA or as long RNA.

**Report layout (`/v2`).** The same server also serves the page with the layout of the ViralQuest report, at `http://127.0.0.1:<port>/v2` (the classic page stays at `/`; both read and save the same curation). It adds two tabs:
- **Sequences:** one collapsible card per EVE locus. The header shows family, protein, length, verdict, viral bits and the read evidence of every library: DNA presence, RNA-seq transcription, and small RNA reads with the piRNA signature. The body, read when the card opens, shows the verdict, the genomic context with every read track, the domain architecture, the reads library by library, the tblastn segments, the clusters and the sequence. Cards can be searched, filtered, sorted, selected and exported (FASTA, TSV).
- **Libraries:** three heatmaps of locus × library. They show DNA presence (the calls of the *Genomic DNA* tab, columns grouped by population), RNA-seq transcription (reads, mean depth or share covered) and small RNA abundance with the piRNA signature.

For the piRNA signature, the piRNA-sized reads (24–31 nt) are read strand by strand. The signature is called *primary* when a strand carries a 5′ U bias, *secondary* when a strand carries an A at position 10, and *ping-pong* when the z-score of the 10-nt 5′–5′ overlap is significant. The thresholds can be changed on the page, and a table lists every locus with enough piRNA-sized reads.

The page is generated from `curation.html` by `evetool/curation/web/v2_src/build_v2.py`. Run the script again after changing the classic page.

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
  validation/           viral vs cellular classification, verdict rules, positional validation, viral protein map
  annotation/           CD-Search, viral HMMs, Pfam, repeats and low complexity, read coverage (long RNA),
                        small RNA profiles (srna.py), genomic DNA reads (dna_reads.py)
  visualization/        locus figures and the static dashboard
  comparative/          cross-species clustering, hierarchical clustering, flank synteny,
                        codon usage and dinucleotide composition vs circulating viruses
  curation/             curation server, web page (web/) and primer design (primers.py)
  databases/            database download, packing, audit and HMM calibration
examples/               small public test genome
docs/                   guides (downloading the RefSeq viral proteins from NCBI Virus)
```

Each module can also be run on its own, for example `python -m evetool.annotation.repeats -h` (with the repository folder on `PYTHONPATH`, or after `pip install -e .`, which also puts `eve-tool` on the PATH).

## Configuration

Each path is resolved in this order:

1. an environment variable (`EVE_DATA_DIR`, `EVE_HMM_DIR`, `EVE_DIAMOND_DIR`, `EVE_RESULTS_DIR`, `EVE_PFAM_HMM`; `EVE_NCBI_EMAIL` for the NCBI contact e-mail);
2. `eve_config.ini`;
3. the defaults under `data_dir`.

`./eve-tool config` prints the paths in use.

## Citation

If you use this tool, please cite: _in preparation_. Please also cite the underlying tools and databases: DIAMOND, MMseqs2, HMMER/pyhmmer, BLAST+, orfipy, RVDB, eggNOG, Pfam, CDD and UniProt.

## License

The code is released under the [MIT License](LICENSE). The viral database bundle on Zenodo has its own license and derives from RVDB, eggNOG and NCBI RefSeq; cite them too.
