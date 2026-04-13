# 🧬 Projeto Banco Genoma

Projeto de construção de um banco de dados de sequências genômicas para micro, meso e macroorganismos do solo, com armazenamento, filtragem e disponibilização de dados biológicos.

---

## 🔬 Ferramenta `primer_blast_local`

Ferramenta de validação de primers *in silico*, desenvolvida e validada contra o NCBI Primer-BLAST para 12 organismos (fungos e protozoários).

### ✨ Funcionalidades

* ✅ Suporte nativo a degenerações IUPAC (S, Y, R, W, K, M, B, D, H, V)
* ✅ Interface web simples (Flask + TailwindCSS)
* ✅ Linha de comando
* ✅ Extração automática de sequências de amplicons
* ✅ Download de resultados em CSV, FASTA e HTML
* ✅ Configuração avançada de parâmetros (E-value, mismatches 3', cobertura, etc.)
* ✅ Validação comprovada contra o NCBI (identidade > 99,9%)

---

## 🚀 Como Usar

### 🌐 Opção 1: Interface Web

```bash
# 1. Instale as dependências
pip install flask biopython pandas

# 2. Navegue até a pasta Web
cd Web

# 3. Execute o servidor
python app.py

# 4. Acesse no navegador
http://localhost:5000
```

---

### 💻 Opção 2: Linha de Comando

#### Comando básico

```bash
python primer_blast_local.py -g genoma.fna -p primers.fasta -o resultados
```

#### Comando completo

```bash
python primer_blast_local.py \
    -g genoma.fna \
    -p primers.fasta \
    -o resultados \
    -e 30000 \
    --min_size 4000 \
    --max_size 6000 \
    -t 4 \
    --max_target_seqs 100 \
    --max_3prime_mismatches 2 \
    --qcov_hsp_perc 80 \
    --amp_seq \
    -m 0
```

---

## ⚙️ Explicação dos Parâmetros

| Parâmetro               | Descrição                    | Valor recomendado |
| ----------------------- | ---------------------------- | ----------------- |
| -g                      | Arquivo FASTA do genoma      | genoma.fna        |
| -p                      | Arquivo FASTA com os primers | primers.fasta     |
| -o                      | Prefixo de saída             | resultados        |
| -e                      | E-value do BLAST             | 30000             |
| --min_size              | Tamanho mínimo do amplicon   | 4000              |
| --max_size              | Tamanho máximo do amplicon   | 6000              |
| -t                      | Número de threads            | 4                 |
| --max_target_seqs       | Limite de hits               | 100               |
| --max_3prime_mismatches | Mismatches no 3'             | 2                 |
| --qcov_hsp_perc         | Cobertura mínima (%)         | 80                |
| --amp_seq               | Extrai amplicons             | flag              |
| -m                      | Limiar de Tm                 | 0                 |

---

## ⚙️ Parâmetros Avançados (Opcionais)

| Parâmetro         | Descrição              | Padrão |
| ----------------- | ---------------------- | ------ |
| --no_blast        | Pula execução do BLAST | Off    |
| --use_existing_db | Reutiliza banco BLAST  | Off    |
| --na              | Na⁺ (mM)               | 50     |
| -k / --pot        | K⁺ (mM)                | 0      |
| --tris            | Tris (mM)              | 0      |
| --mg              | Mg²⁺ (mM)              | 0      |
| --dntps           | dNTPs (mM)             | 0      |
| --saltcorr        | Correção de sal        | 5      |

---

## 🧪 Exemplo completo

```bash
python primer_blast_local.py \
    -g genoma_saccharomyces_cerevisiae.fna \
    -p primers.fasta \
    -o resultados_sc \
    -e 30000 \
    --min_size 4000 \
    --max_size 6000 \
    -t 4 \
    --max_target_seqs 100 \
    --max_3prime_mismatches 2 \
    --qcov_hsp_perc 80 \
    --amp_seq \
    -m 0
```

---

## 📁 Estrutura do Projeto

```
primer_blast_local/
├── primer_blast_local.py
├── reformat.py
├── run_parse_blastn.py
├── primers.fasta
├── extrair_amplicons.py
├── comparar_amplicons.py
├── gerar_relatorio_amplicons.py
└── Web/
    ├── app.py
    ├── templates/
    │   └── index.html
    ├── uploads/
    └── results/
```

---

## 📜 Scripts Auxiliares

### extrair_amplicons.py

```bash
# Arquivos separados
python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo S_cerevisiae

# Arquivo único
python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo S_cerevisiae --unico
```

---

### comparar_amplicons.py

```bash
# Comparação única
python comparar_amplicons.py --ncbi ncbi.fasta --local local.fasta --org "S. cerevisiae"

# Em lote
python comparar_amplicons.py --batch comparacoes.csv --relatorio resultados.txt
```

---

### gerar_relatorio_amplicons.py

```bash
python gerar_relatorio_amplicons.py \
    --ncbi amplicon_ncbi.fasta \
    --local amplicon_local.fasta \
    --org "Saccharomyces cerevisiae" \
    --cromossomo "NC_001144.5"
```

---

## 📊 Validação

A ferramenta foi validada contra o NCBI Primer-BLAST:

| #  | Organismo                 | Tipo        | Produtos | Tamanho        |
| -- | ------------------------- | ----------- | -------- | -------------- |
| 1  | Saccharomyces cerevisiae  | Fungo       | 2        | 4.817 bp       |
| 2  | Phaeodactylum tricornutum | Protozoário | 1        | 4.886 bp       |
| 3  | Rhizophagus irregularis   | Fungo       | 10       | 4.675–4.704 bp |
| 4  | Dictyostelium discoideum  | Protozoário | 1        | 5.530 bp       |
| 5  | Fusarium oxysporum        | Fungo       | 3        | 4.480 bp       |
| 6  | Neurospora crassa         | Fungo       | 1        | 4.522 bp       |
| 7  | Phytophthora infestans    | Protozoário | 33       | 4.023–5.471 bp |
| 8  | Candida albicans          | Fungo       | 1        | 4.474 bp       |
| 9  | Aspergillus fumigatus     | Fungo       | 1        | 4.553 bp       |
| 10 | Blastocystis hominis      | Protozoário | 21       | 4.345–4.492 bp |
| 11 | Babesia microti           | Protozoário | 2        | 4.765 bp       |
| 12 | Sclerotinia sclerotiorum  | Fungo       | 4        | 5.113-5.122 bp |

### ✅ Resultados

* Identidade média: > 99,9%
* Diferença média: 1–2 bp (0,02–0,04%)
* 100% de concordância nos cromossomos

---

## 📚 Referência

Latz, M. A. C., et al. (2022).
*Short- and long-read metabarcoding of the eukaryotic rRNA operon.*
Molecular Ecology Resources, 22, 2304–2318.
https://doi.org/10.1111/1755-0998.13623
