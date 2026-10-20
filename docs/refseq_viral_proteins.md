# Downloading the RefSeq viral proteins

The viral DIAMOND database of the EVE tool (`DIAMOND/sequences.fasta`) holds every RefSeq protein of every virus, downloaded from NCBI Virus. The tool never downloads this file by itself. The published bundle ships one copy. To use a newer one, download it yourself as described here and give it to the database build:

```bash
./eve-tool build-databases --refseq-fasta refseq_viral_proteins.fasta
```

The screening reads the genome type, species, genus and family of each hit from the FASTA definition line. The file must therefore have the definition line described below. `build-databases` checks it and refuses a file that does not have it.

## 1. Search

1. Open NCBI Virus: <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/>.
2. In **Search by virus**, type `Viruses` and choose **Viruses, taxid:10239**. This includes all viruses, not only RNA viruses.
3. In the results, set **Sequence Type** to **Protein**.
4. Under **Nucleotide Sequence Type** (or **RefSeq Genome Completeness**, depending on the page version), keep only **RefSeq**.

Leave every other filter empty: no host, no genome type, no completeness, no date.

The screening needs the DNA and reverse-transcribing viruses too. They let it tell an EVE from a host gene, so do not filter by genome type here. The pipeline does that itself.

## 2. Download

1. Click **Download**.
2. Choose **Protein** and **FASTA**, then **Next**.
3. Choose **Download All Records**, then **Next**.
4. In **FASTA definition line**, choose **Build custom** and add these fields **in this order**:

   | # | Field |
   |---|---|
   | 1 | Accession (with version) |
   | 2 | GenBank Title |
   | 3 | Molecule type (genomic molecule type) |
   | 4 | Species |
   | 5 | Genus |
   | 6 | Family |
   | 7 | Length |

   The separator is `|`.
5. Click **Download**. The file is about 270 MB.

## 3. Check the file

Each header must look like this (accession, a space, then six `|`-separated fields):

```
>YP_013025916.1 |N protein [Parahenipavirus daeryongense]|ssRNA(-)|Parahenipavirus daeryongense|Parahenipavirus|Paramyxoviridae|574
>YP_013031981.1 |polyprotein [Potyvirus thevetiae]|ssRNA(+)|Potyvirus thevetiae|Potyvirus|Potyviridae|3182
>YP_011426526.1 |replication-associated protein [Pacific flying fox faeces associated circular DNA virus-2]|ssDNA|Protegovirus tongis|Protegovirus|Geplanaviridae|363
```

Quick checks:

```bash
grep -c '>' refseq_viral_proteins.fasta                       # number of proteins (~720,000 in 2026)
grep -m3 '>' refseq_viral_proteins.fasta                      # headers as above
grep '>' refseq_viral_proteins.fasta | awk -F'|' 'NF<7' | head  # should print nothing
```

The third field is the genome type (`ssRNA(+)`, `ssRNA(-)`, `dsRNA`, `ssDNA`, `dsDNA`, `ssRNA-RT`…). Count how many proteins each type has:

```bash
grep '>' refseq_viral_proteins.fasta | cut -d'|' -f3 | sort | uniq -c | sort -nr | head
```

Most proteins are dsDNA (phages), followed by ssRNA(+) and ssRNA(-). If almost everything is `unknown` or empty, the **Molecule type** field was missed.

## 4. Use it

```bash
./eve-tool build-databases --refseq-fasta refseq_viral_proteins.fasta
```

The build copies the file, rebuilds the DIAMOND database from it and calibrates the HMMs with it. It also records the file's path, date, number of proteins and md5 in `build_manifest.json`. Keep the downloaded file, and write down the download date: it is the version of this database.
