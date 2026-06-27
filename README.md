# Detecção de Domínios Funcionais em Proteínas por Redes Complexas

Projeto desenvolvido como parte de uma disciplina de Redes Complexas. O objetivo é modelar estruturas proteicas como grafos e aplicar análise de redes para identificar domínios funcionais e comunidades estruturais, comparando os resultados com as anotações biológicas do RCSB PDB.

O pipeline completo é executado em dois estágios: primeiro em uma **proteína de teste** (1A3N) para validação metodológica, depois na **proteína-alvo** (6B1T).

---

## Estrutura do projeto

```
Projeto_Redes_Complexas_Proteina_Teste_1A3N_e_Proteina_Alvo_6B1T.ipynb
resultados/
├── 1A3N/
│   ├── 1A3N_residuos_ca.csv
│   ├── 1A3N_arestas.csv
│   ├── 1A3N_metricas_rede.csv
│   ├── 1A3N_distribuicao_graus.csv
│   ├── 1A3N_histograma_graus.png
│   ├── 1A3N_graus_loglog.png
│   ├── 1A3N_centralidades.csv
│   ├── 1A3N_top20_*.csv
│   ├── 1A3N_histogramas_centralidades.png
│   ├── 1A3N_comunidades_louvain.csv
│   ├── 1A3N_resumo_comunidades.csv
│   ├── 1A3N_comunidade_x_cadeia.csv
│   ├── 1A3N_heatmap_comunidade_x_cadeia.png
│   ├── 1A3N_metricas_validacao.csv
│   ├── 1A3N_3d_comunidades.html
│   └── resumo_execucao_1A3N.txt
└── 6B1T/
    ├── 6B1T_residuos_ca.csv
    ├── 6B1T_arestas.csv
    ├── 6B1T_metricas_rede.csv
    ├── 6B1T_distribuicao_graus.csv
    ├── 6B1T_histograma_graus.png
    ├── 6B1T_graus_loglog.png
    ├── 6B1T_centralidades.csv
    ├── 6B1T_top20_*.csv
    ├── 6B1T_histogramas_centralidades.png
    ├── 6B1T_comunidades_louvain.csv
    ├── 6B1T_resumo_comunidades.csv
    ├── 6B1T_comunidade_x_cadeia.csv
    ├── 6B1T_heatmap_comunidade_x_cadeia.png
    ├── 6B1T_comunidades_com_familias.csv
    ├── 6B1T_contagem_familias_rcsb.csv
    ├── 6B1T_tabela_comunidade_familia_rcsb.csv
    ├── 6B1T_heatmap_comunidade_familia_rcsb.png
    ├── 6B1T_comunidade_x_familia_rcsb.csv
    ├── 6B1T_heatmap_comunidade_x_familia_rcsb.png
    ├── 6B1T_comunidade_x_hexon_penton.csv
    ├── 6B1T_heatmap_comunidade_x_hexon_penton.png
    ├── 6B1T_metricas_validacao_biologica.csv
    ├── 6B1T_3d_comunidades.html
    └── resumo_execucao_6B1T.txt
```

---

## Proteínas analisadas

### 1A3N — Proteína de teste

Hemoglobina humana (desoxigenada), estrutura com múltiplas cadeias (A, B, C, D), utilizada para validar o pipeline antes da execução na proteína-alvo. Seu tamanho reduzido permite cálculos exatos de todas as centralidades.

- **Fonte:** [RCSB PDB — 1A3N](https://www.rcsb.org/structure/1A3N)
- **Betweenness:** cálculo exato

### 6B1T — Proteína-alvo

Adenovírus humano tipo 5 (HAdV-C5), estrutura de capsídeo viral resolvida por cryo-EM. Composta por dezenas de cadeias distribuídas em diferentes famílias proteicas estruturais. Por ser uma estrutura de grande porte, o arquivo PDB é distribuído pelo RCSB em dois bundles separados.

- **Fonte:** [RCSB PDB — 6B1T](https://www.rcsb.org/structure/6B1T)
- **Arquivos necessários:** `6b1t-pdb-bundle1.pdb` e `6b1t-pdb-bundle2.pdb`
- **Betweenness:** cálculo aproximado (`k=300`)

Famílias estruturais mapeadas:

| Cadeias | Família estrutural |
|---|---|
| A–L | Hexon protein |
| M | Penton protein |
| N | Pre-hexon-linking protein IIIa |
| O, P | Pre-hexon-linking protein VIII |
| Q–T | Hexon-interlacing protein |
| U, V, X, Y | Pre-protein VI |
| W | Pre-histone-like nucleoprotein |

---

## Metodologia

### Modelagem como grafo

Cada resíduo de aminoácido é representado pelo seu átomo Cα. Dois resíduos são conectados por uma aresta se a distância euclidiana entre seus Cα for menor ou igual ao limiar definido (`DISTANCE_THRESHOLD = 8.0 Å`). O grafo resultante é não-direcionado e não ponderado (peso = 1.0 para todas as arestas).

A busca de pares é feita com `cKDTree` (scipy), que torna o processo eficiente mesmo para proteínas grandes.

### Métricas da rede

Para cada proteína são calculadas:

- Número de nós e arestas
- Densidade
- Grau médio e máximo
- Coeficiente de agrupamento médio
- Número de componentes conectados e tamanho do maior componente

### Distribuição de graus

Histograma linear e gráfico log-log da distribuição de graus, com a média destacada.

### Centralidades

Calculadas para todos os resíduos e salvas com ranking Top 20 por métrica:

| Centralidade | Método |
|---|---|
| Degree | `nx.degree_centrality` |
| Closeness | `nx.closeness_centrality` |
| Eigenvector | `nx.eigenvector_centrality` (max_iter=1000) |
| PageRank | `nx.pagerank` (α=0.85) |
| Betweenness | Exato (1A3N) / Aproximado k=300 (6B1T) |

### Detecção de comunidades — Louvain

O algoritmo de Louvain (`python-louvain`) é aplicado para identificar comunidades, com resolução configurável separadamente para cada proteína (`LOUVAIN_RESOLUTION = 1.0` para ambas). A modularidade é calculada ao final.

### Validação biológica

**1A3N:** as comunidades detectadas são comparadas com as cadeias PDB via tabela cruzada, ARI e NMI.

**6B1T:** a validação é feita em três níveis:

1. Comunidade × cadeia PDB
2. Comunidade × família estrutural RCSB (7 famílias)
3. Comunidade × grupo mínimo (Hexon / Penton / Outras)

Métricas quantitativas calculadas: Modularidade, Purity, ARI, NMI, Homogeneity e Completeness.

---

## Parâmetros configuráveis

| Parâmetro | Valor padrão | Descrição |
|---|---|---|
| `DISTANCE_THRESHOLD` | `8.0` | Limiar de distância Cα-Cα (Å) |
| `LOUVAIN_RESOLUTION_TESTE` | `1.0` | Resolução Louvain para 1A3N |
| `LOUVAIN_RESOLUTION_6B1T` | `1.0` | Resolução Louvain para 6B1T |
| `RANDOM_STATE` | `42` | Semente para reprodutibilidade |
| `BETWEENNESS_K_6B1T` | `300` | Amostras para betweenness aproximada |
| `CALCULAR_CLOSENESS` | `True` | Ativa/desativa Closeness |
| `CALCULAR_EIGENVECTOR` | `True` | Ativa/desativa Eigenvector |
| `CALCULAR_PAGERANK` | `True` | Ativa/desativa PageRank |

---

## Como executar

### Pré-requisitos

O notebook instala automaticamente todas as dependências na primeira célula. As bibliotecas utilizadas são:

- `biopython` — leitura de arquivos PDB
- `networkx` — construção e análise de grafos
- `python-louvain` — detecção de comunidades
- `plotly` — visualizações 3D interativas
- `pandas`, `numpy`, `scipy` — manipulação de dados e KDTree
- `matplotlib` — gráficos estáticos
- `scikit-learn` — métricas de validação (ARI, NMI, homogeneity, completeness)
- `tqdm` — barras de progresso
- `kaleido` — exportação de figuras Plotly para PNG

### Execução no Google Colab

1. Abra o notebook no Google Colab.
2. Execute todas as células em ordem (**Runtime → Run all**).
3. Quando solicitado na célula 13, faça upload dos dois arquivos da 6B1T:
   - `6b1t-pdb-bundle1.pdb`
   - `6b1t-pdb-bundle2.pdb`

   Os arquivos estão disponíveis no bundle `tar.gz` da página do RCSB: [https://www.rcsb.org/structure/6B1T](https://www.rcsb.org/structure/6B1T) → *Download Files → PDB Format (gz bundle)*.

4. Ao final, a célula 26 compacta todos os resultados em `resultados_redes_complexas_1A3N_6B1T.zip` e inicia o download automaticamente.

---

## Saídas geradas

Todos os arquivos são salvos em `resultados/1A3N/` e `resultados/6B1T/`:

| Tipo | Arquivos |
|---|---|
| Dados tabulares | `.csv` com resíduos, arestas, centralidades, comunidades, validação |
| Figuras estáticas | `.png` — histogramas de graus, log-log, centralidades, heatmaps |
| Visualização 3D | `.html` interativo (Plotly) com resíduos coloridos por comunidade |
| Resumo textual | `.txt` com métricas e lista de arquivos gerados |

---

## Dependências

```
biopython
networkx
python-louvain
plotly
pandas
numpy
scipy
matplotlib
scikit-learn
tqdm
kaleido
```

---

## References

- PDB 1A3N: Deoxyhemoglobin — [https://www.rcsb.org/structure/1A3N](https://www.rcsb.org/structure/1A3N)
- PDB 6B1T: Human Adenovirus Type 5 — [https://www.rcsb.org/structure/6B1T](https://www.rcsb.org/structure/6B1T)
- Blondel, V. D. et al. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics*, 2008(10), P10008.
- Newman, M. E. J. (2004). Analysis of weighted networks. *Physical Review E*, 70(5), 056131.
