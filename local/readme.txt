# Ferramenta de análise de EVEs

Um pipeline de bioinformática escalável, modular e universal para a mineração, validação, genômica funcional e análise evolutiva de Elementos Virais Endógenos (EVEs) em genomas eucariotos hospedeiros (artrópodes, vertebrados, plantas e outros táxons). Este conjunto de ferramentas integra predição de ORFs, desduplicação de sequências, alinhamento de alta sensibilidade (DIAMOND + HMMER3), decomposição tensorial, análise posicional de consenso e visualização de sintenia genômica.

## Arquitetura do Sistema

O pipeline é modular e orquestrado por scripts em Python e R, permitindo execução paralela e personalização de bancos de dados.

* **`main.py`**: Orquestrador central que gerencia os argumentos de linha de comando, valida as dependências do sistema e aciona os módulos na sequência correta.


* **`orf_predictor.py`**: Extrai transcritos (ORFs) em nucleotídeos e peptídeos utilizando o `orfipy`, capturando sequências parciais para não perder cicatrizes virais antigas.


* **`cluster_filter.py`**: Remove redundâncias de sequências (limiar padrão de 90%) executando o algoritmo de agrupamento do `MMseqs2`.


* **`sequence_searcher.py`**: Executa a triagem primária em alta velocidade via `DIAMOND` (modo blastx) contra bancos de dados virais customizados ou RefSeq.


* **`hmm_scanner.py`**: Resgata EVEs altamente divergentes submetendo as sequências não alinhadas pelo DIAMOND a perfis HMM (ex: RVDB, vFam) usando o `hmmsearch`.


* **`online_blaster.py`**: Submete automaticamente as candidatas consolidadas à Common URL API do NCBI para validação contra falsos positivos do hospedeiro.


* **`genome_mapper.py`**: Ancara as EVEs de volta no genoma do hospedeiro usando `tBLASTn` para extrair as coordenadas exatas de integração.


* **`data_extraction_cluster_eves.py`**: Utiliza o `bedtools` para aglutinar segmentos próximos e extrair os flancos genômicos (ex: 1.000 nt de borda) para análise do contexto hospedeiro.


* **`utils.py`**: Funções genéricas de apoio, governança de diretórios e filtros taxonômicos.


* **`mapa_proteinas.py`**: Realinha cada EVE candidata (Smith-Waterman, BLOSUM62, nas 6 fases do trecho alinhado pelo DIAMOND) contra a proteína viral de referência e mostra onde ela cai na proteína (painel "EVEs on the viral proteins" da curadoria; `python mapa_proteinas.py -o resultados --tsv mapa.tsv`). Marca as EVEs cujo sinal viral está na fita oposta ou em outra fase do ORF.



### Módulos de Visualização e Pós-Processamento (R & Bash)

* **`plot_eve_gene_context.R`**: Gera mapas visuais interativos (HTML) e estáticos da vizinhança genômica, mapeando as EVEs entre íntrons e éxons do hospedeiro utilizando `gggenes` e `plotly`.


* **`pipeline_sintenia_eves.R`**: Automatiza a descoberta de blocos ortólogos via `DECIPHER` e desenha gráficos de conservação sintênica ancestral usando o `geneviewer`.


* **`smallrna_viral_pipeline.py`**: Mapeia bibliotecas de pequenos RNAs (siRNAs/piRNAs) contra as EVEs usando `bowtie`, verificando atividade imunológica ativa do hospedeiro.



---

## Dependências e Pré-requisitos

O pipeline interage fisicamente com ferramentas binárias otimizadas instaladas no sistema operacional.

**Softwares de Linha de Comando (PATH):**

* `orfipy`

* `mmseqs` (MMseqs2)


* `diamond`

* `hmmer` (`hmmsearch`, `hmmpress`)


* `blast+` (`makeblastdb`, `tblastn`)


* `bedtools` e `samtools`

* `bowtie` (apenas para análise de pequenos RNAs)


* `emboss` (para o comando `cons` gerador de sequências consenso)



**Bibliotecas Python 3:**

* `biopython`, `pandas`, `requests`.



**Bibliotecas R (Análises Visuais e Evolutivas):**

* `DECIPHER`, `geneviewer`, `gggenes`, `ggplot2`, `rtracklayer`, `Rsamtools`, `Biostrings`, `plotly`, `htmlwidgets`.



---

## Configuração dos Caminhos

Os bancos de dados e a pasta de resultados não ficam no código. Copie `eve_config.example.ini` para `eve_config.ini` (na mesma pasta do `main.py`) e ajuste os caminhos:

```ini
[paths]
data_dir = /caminho/eve_data          # raiz; as demais pastas são subpastas dela por padrão
# hmm_dir     = <data_dir>/Data_base_hmm_eves   (HMMs, RVDB, CDD, taxdump)
# diamond_dir = <data_dir>/DIAMOND              (mydb.dmnd, geral/uniprot_sprot.dmnd)
# results_dir = <data_dir>/Resultados
# pfam_hmm    = <data_dir>/Pfam/Pfam-A.hmm
```

Cada chave também pode ser definida por variável de ambiente (`EVE_DATA_DIR`, `EVE_HMM_DIR`, `EVE_DIAMOND_DIR`, `EVE_RESULTS_DIR`, `EVE_PFAM_HMM`), e `EVE_CONFIG` aponta para outro arquivo de configuração. Para conferir os caminhos em uso: `python config.py`.

## Como Usar (Uso Básico)

A execução base exige o genoma do hospedeiro e ao menos um banco de dados de referência (formatado ou um diretório de arquivos brutos para formatação automática).

**Exemplo com bancos de dados customizados formatados na hora:**

```bash
python main.py -i genoma_artropode.fasta \
               --fasta-dir /caminho/banco_viral_pep/ \
               --hmm-dir /caminho/modelos_rvdb/ \
               -o resultados_eve \
               -t 8

```


**Exemplo utilizando bancos já pré-indexados (Espécie Única):**

```bash
python main.py -i genoma_artropode.fasta \
               --dmnd-db nr_viral.dmnd \
               --hmm-db vfam.hmm \
               --gff anotacao_hospedeiro.gff3 \
               --rnaseq-r1 rna_R1.fastq.gz \
               --rnaseq-r2 rna_R2.fastq.gz \
               --flank 20000 \
               -o resultados_eve \
               -t 12
```

**Exemplo com Modo Multi-Espécies e Comparação Cruzada de Clusters:**

```bash
python main.py --genomes-dir /caminho/pasta_genomas/ \
               --species-tsv lista_especies.tsv \
               --dmnd-db nr_viral.dmnd \
               --hmm-db vfam.hmm \
               --flank 20000 \
               -o analise_multiespecies \
               -t 16
```

O arquivo `lista_especies.tsv` deve conter as seguintes colunas (separadas por tabulação). A coluna `genome` aceita tanto arquivos FASTA locais quanto números de acesso **GCA_** (GenBank) ou **GCF_** (RefSeq) do NCBI, que serão baixados automaticamente com suas respectivas anotações GFF3:
```tsv
species	genome	gff	rnaseq_r1	rnaseq_r2
Phlebotomus_sergenti	GCA_058968885.1	None	None	None
Phlebotomus_papatasi	GCF_000466445.1	None	None	None
Lutzomyia_longipalpis	L_longipalpis.fasta	L_longipalpis.gff3	rnaseq_Llong_R1.fq.gz	rnaseq_Llong_R2.fq.gz
```
*(Se for passado um accession GCA/GCF, o pipeline baixa o genoma FASTA e busca automaticamente o GFF3 correspondente no NCBI via `ncbi_genome_downloader.py`).*

**Download avulso de genomas via NCBI Datasets:**
```bash
python ncbi_genome_downloader.py -a GCA_058968885.1 -o ./pasta_genomas/
```

## Fluxo de Trabalho (Workflow Genômico)

1. **Descoberta e Triagem:** Transcritos putativos são extraídos com `orfipy` e colapsados matematicamente com `MMseqs2` para remover redundâncias.
2. **Dupla Camada de Mapeamento:** Segmentos passam por triagem primária via `DIAMOND`. O que não alinha é submetido a modelos de Markov (`HMMER/hmmsearch`) para captura de cicatrizes e regiões hiperdivergentes.
3. **Ancoragem:** As EVEs candidatas (peptídeos) são ancoradas nos scaffolds hospedeiros via `tBLASTn`. Os flancos de 20.000 nucleotídeos são exportados junto com a sequência viral para caracterização do contexto genômico.
4. **Verificação de Atividade Transcricional (RNA-seq):** O módulo `rnaseq_analyzer.py` mapeia leituras brutas de RNA-seq (Bowtie2) contra os loci e flancos de 20 kb, gerando perfis de profundidade e determinando o status de expressão (ativo, moderado, apenas flanco ou silenciado).
5. **Visualização do Contexto Genômico:** O módulo `plot_genomic_context.py` renderiza imagens em alta resolução (PNG e PDF) da inserção viral cercada pelos flancos de 20 kb, genes/éxons do hospedeiro (se anotados) e a curva de cobertura de RNA-seq.
6. **Comparação Cruzada Multi-Espécies:** O módulo `cross_species_comparator.py` consolida todas as EVEs encontradas, executa busca all-vs-all (tBLASTx), constrói redes de homologia com NetworkX e gera:
   - Resumo dos MacroClusters com status evolutivo (Ancestral / Compartilhado vs Específico de Espécie).
   - Matriz de presença/ausência e número de cópias por espécie.
   - Heatmap visual de compartilhamento entre espécies hospedeiras.
