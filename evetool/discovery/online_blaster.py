import os
import time
import requests
from pathlib import Path
from Bio import SeqIO, Entrez

def rodar_blastx_api(
    fasta_query: str,
    diretorio_saida: str,
    nome_arquivo_saida: str = "validacao_ncbi_blastx.tsv",
    max_espera_segundos: int = 1800
) -> str:
    """
    Submete ORFs unificadas à Common URL API do NCBI usando os comandos Put/Get.
    Gera uma tabela de resultados similar ao output do DIAMOND.
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    arquivo_saida = os.path.join(diretorio_saida, nome_arquivo_saida)
    url_base = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    
    # 1. Carrega as sequências e formata como string única
    sequencias = list(SeqIO.parse(fasta_query, "fasta"))
    if not sequencias:
        raise ValueError("O arquivo consolidado está vazio.")
    fasta_string = "".join([seq.format("fasta") for seq in sequencias])
    
    # 2. CMD = Put (Submissão da Busca)
    params_put = {
        "CMD": "Put",
        "PROGRAM": "blastx",
        "DATABASE": "ClusteredNR",
        "QUERY": fasta_string,
        "EXPECT": "1e-3",
        "HITLIST_SIZE": 5
    }
    
    print(f"[run] Submetendo {len(sequencias)} sequências via NCBI Common URL API...")
    resposta_put = requests.post(url_base, data=params_put, timeout=60)
    resposta_put.raise_for_status()
    
    # Extrai o Request ID (RID) e o tempo estimado (RTOE) do HTML retornado
    rid = None
    rtoe = 30
    for linha in resposta_put.text.split('\n'):
        if linha.startswith("    RID = "):
            rid = linha.split("=")[1].strip()
        elif linha.startswith("    RTOE = "):
            rtoe = int(linha.split("=")[1].strip())
            
    if not rid:
        raise RuntimeError("Falha ao obter o RID do NCBI. A submissão falhou.")
        
    print(f"[run] Busca iniciada com sucesso. RID: {rid}. Espera estimada: {rtoe} segundos.")
    
    # 3. CMD = Get (Checagem de Status via FORMAT_OBJECT=SearchInfo)
    params_status = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo"
    }
    
    inicio_espera = time.time()
    while True:
        tempo_decorrido = time.time() - inicio_espera
        if tempo_decorrido > max_espera_segundos:
            raise TimeoutError(
                f"Tempo limite de {max_espera_segundos}s excedido aguardando o RID {rid} no servidor do NCBI."
            )

        time.sleep(rtoe if rtoe > 0 else 60)
        rtoe = 60  # Política do NCBI: nunca checar o mesmo RID mais de 1x por minuto
        
        resposta_status = requests.get(url_base, params=params_status, timeout=60)
        if "Status=WAITING" in resposta_status.text:
            print("[run] Status: WAITING. O servidor ainda está processando...")
            continue
        elif "Status=FAILED" in resposta_status.text or "Status=UNKNOWN" in resposta_status.text:
            raise RuntimeError(f"O processamento falhou no servidor do NCBI para o RID {rid}.")
        elif "Status=READY" in resposta_status.text:
            if "ThereAreHits=yes" in resposta_status.text:
                print("[run] Status: READY. Fazendo download dos resultados...")
                break
            else:
                print("[run] Status: READY. Nenhum hit encontrado para as sequências.")
                with open(arquivo_saida, "w") as f:
                    f.write("Nenhum hit encontrado no banco nr do NCBI.\n")
                return arquivo_saida

    # 4. CMD = Get (Recuperação dos Dados em formato Tabular)
    params_get_data = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Text",
        "ALIGNMENT_VIEW": "Tabular"
    }
    
    resposta_dados = requests.get(url_base, params=params_get_data, timeout=120)
    resposta_dados.raise_for_status()
    
    with open(arquivo_saida, "w") as out_handle:
        out_handle.write(resposta_dados.text)
        
    print(f"[run] Tabela de validação salva em: {arquivo_saida}")
    return arquivo_saida


def dividir_fasta_para_blastx(fasta_path: str, diretorio_saida: str, max_seqs: int = 10, max_nt: int = 20000) -> list:
    """
    Divide um FASTA em lotes menores para submissão ao blastx online do NCBI,
    respeitando dois limites simultâneos por lote: no máximo `max_seqs`
    sequências E no máximo `max_nt` nucleotídeos somados nesse lote. Lotes
    menores reduzem drasticamente o tempo de fila estimado (RTOE) do NCBI
    em comparação a submeter tudo de uma vez.

    Uma sequência sozinha maior que max_nt vai para seu próprio lote (não dá
    para dividir uma ORF ao meio sem invalidá-la).
    """
    os.makedirs(diretorio_saida, exist_ok=True)
    nome_base = Path(fasta_path).stem

    sequencias = list(SeqIO.parse(fasta_path, "fasta"))
    if not sequencias:
        raise ValueError(f"O arquivo '{fasta_path}' está vazio, nada para dividir.")

    lotes = []
    lote_atual = []
    nt_atual = 0

    for seq in sequencias:
        tamanho = len(seq.seq)
        excederia_qtd = len(lote_atual) + 1 > max_seqs
        excederia_nt = bool(lote_atual) and (nt_atual + tamanho > max_nt)

        if lote_atual and (excederia_qtd or excederia_nt):
            lotes.append(lote_atual)
            lote_atual = []
            nt_atual = 0

        lote_atual.append(seq)
        nt_atual += tamanho

    if lote_atual:
        lotes.append(lote_atual)

    caminhos_lotes = []
    for i, lote in enumerate(lotes, start=1):
        caminho_lote = os.path.join(diretorio_saida, f"{nome_base}_lote{i:02d}.fasta")
        with open(caminho_lote, 'w') as f_out:
            for seq in lote:
                f_out.write(seq.format("fasta"))
        caminhos_lotes.append(caminho_lote)

    total_nt = sum(len(s.seq) for s in sequencias)
    print(f"[run] {len(sequencias)} sequências ({total_nt} nt) divididas em {len(caminhos_lotes)} lote(s) "
          f"(máx. {max_seqs} seqs / {max_nt} nt por lote).")

    return caminhos_lotes


def rodar_blastx_api_em_lotes(
    fasta_query: str,
    diretorio_saida: str,
    max_seqs: int = 10,
    max_nt: int = 20000,
    pausa_entre_lotes: int = 10,
    max_espera_segundos: int = 1800
) -> list:
    """
    Divide o FASTA de entrada em lotes menores (ver dividir_fasta_para_blastx)
    e submete cada lote separadamente ao blastx online do NCBI. Retorna a
    lista de tabelas de resultado geradas, uma por lote. Se um lote falhar
    (rede, timeout, fila do NCBI), os demais continuam sendo processados.
    """
    lotes = dividir_fasta_para_blastx(fasta_query, diretorio_saida, max_seqs=max_seqs, max_nt=max_nt)

    tabelas_resultado = []
    for i, lote in enumerate(lotes, start=1):
        print(f"[run] Submetendo lote {i}/{len(lotes)} ao blastx do NCBI ({os.path.basename(lote)})...")
        nome_saida = f"validacao_ncbi_blastx_lote{i:02d}.tsv"
        try:
            tabela = rodar_blastx_api(
                lote,
                diretorio_saida,
                nome_arquivo_saida=nome_saida,
                max_espera_segundos=max_espera_segundos
            )
            tabelas_resultado.append(tabela)
        except Exception as erro:
            print(f"[aviso] Lote {i}/{len(lotes)} falhou, seguindo para o próximo: {erro}")

        # Cortesia com o servidor do NCBI entre submissões de lotes distintos
        if i < len(lotes):
            time.sleep(pausa_entre_lotes)

    if not tabelas_resultado:
        raise RuntimeError("Nenhum lote foi processado com sucesso pelo blastx do NCBI.")

    return tabelas_resultado


def extrair_melhores_hits_proteicos(tabelas_blastx: list, diretorio_saida: str, top_n: int, email_ncbi: str) -> str:
    """
    Lê as tabelas de resultado do blastx (NCBI, um ou mais lotes), extrai os
    `top_n` melhores hits (por bitscore, deduplicados por accession) somando
    todos os lotes, busca as sequências de proteína correspondentes via NCBI
    Efetch, e salva um FASTA consolidado — pronto para ser usado como query
    no tblastn de ancoragem genômica.
    """
    os.makedirs(diretorio_saida, exist_ok=True)

    hits = []  # (accession, bitscore)
    for tabela in tabelas_blastx:
        if not tabela or not os.path.exists(tabela) or os.path.getsize(tabela) == 0:
            continue
        with open(tabela, 'r') as f_in:
            for linha in f_in:
                if linha.startswith("#") or not linha.strip():
                    continue
                campos = linha.strip().split("\t")
                if len(campos) < 12:
                    continue
                subject_id = campos[1]
                try:
                    bitscore = float(campos[11])
                except ValueError:
                    continue
                # subject_id pode vir em formatos diferentes:
                # "gi|123|ref|NP_001.1|" (antigo, com gi), "ref|NP_001.1|"
                # (atual, sem gi) ou apenas "NP_001.1" (accession puro).
                # Em todos os casos, o accession é o último token não-vazio
                # ao separar por "|".
                partes = [p for p in subject_id.split("|") if p]
                accession = partes[-1] if partes else subject_id
                hits.append((accession, bitscore))

    if not hits:
        raise RuntimeError("Nenhum hit encontrado nas tabelas de blastx fornecidas; "
                            "não há proteínas de referência para buscar.")

    # Mantém apenas o melhor bitscore por accession, depois pega os top_n
    melhor_por_accession = {}
    for accession, bitscore in hits:
        if accession not in melhor_por_accession or bitscore > melhor_por_accession[accession]:
            melhor_por_accession[accession] = bitscore

    top_accessions = sorted(melhor_por_accession.items(), key=lambda x: x[1], reverse=True)[:top_n]
    ids_para_buscar = [acc for acc, _ in top_accessions]

    print(f"[run] Buscando {len(ids_para_buscar)} proteínas de referência via NCBI Efetch...")

    Entrez.email = email_ncbi
    handle = Entrez.efetch(db="protein", id=",".join(ids_para_buscar), rettype="fasta", retmode="text")
    fasta_proteinas_texto = handle.read()
    handle.close()

    arquivo_saida = os.path.join(diretorio_saida, f"proteinas_referencia_ncbi_top{len(ids_para_buscar)}.fasta")
    with open(arquivo_saida, 'w') as f_out:
        f_out.write(fasta_proteinas_texto)

    print(f"[run] {len(ids_para_buscar)} proteínas de referência salvas em: {arquivo_saida}")
    return arquivo_saida
