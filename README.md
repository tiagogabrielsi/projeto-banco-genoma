# 🧬 Projeto Banco Genoma

Projeto de construção de um banco de dados de sequências genômicas para micro, meso e macroorganismos do solo, com armazenamento, filtragem e disponibilização de dados biológicos.

---

## 🔬 Ferramenta `primer_blast_local`

Ferramenta de validação de primers _in silico_, desenvolvida e validada contra o NCBI Primer-BLAST para 12 organismos (fungos e protozoários).

### ✨ Funcionalidades

- ✅ Suporte nativo a degenerações IUPAC (S, Y, R, W, K, M, B, D, H, V)
- ✅ Interface web moderna (Flask + TailwindCSS)
- ✅ Linha de comando para automação e integração com pipelines
- ✅ Extração automática de sequências de amplicons
- ✅ Download de resultados em CSV, FASTA e HTML
- ✅ Configuração avançada de parâmetros (E-value, mismatches 3', cobertura, etc.)
- ✅ Validação comprovada contra o NCBI (identidade > 99,9%)

---

## 🚀 Como Usar

### Opção 1: Interface Web

cmd
# 1. Instale as dependências
pip install flask biopython pandas

# 2. Navegue até a pasta Web
cd Web

# 3. Execute o servidor
python app.py

# 4. Acesse no navegador
http://localhost:5000


Opção 2: Linha de Comando (Terminal)

# Comando Básico
cmd
python primer_blast_local.py -g genoma.fna -p primers.fasta -o resultados

# Comando Completo (Parâmetros Recomendados)
cmd
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
    
Explicação dos Parâmetros
Parâmetro	Descrição	Valor Recomendado
-g	Arquivo FASTA do genoma	genoma.fna
-p	Arquivo FASTA com os primers	primers.fasta
-o	Prefixo para arquivos de saída	resultados
-e	E-value do BLAST	30000
--min_size	Tamanho mínimo do amplicon (bp)	4000
--max_size	Tamanho máximo do amplicon (bp)	6000
-t	Número de threads	4
--max_target_seqs	Limite de hits do BLAST	100
--max_3prime_mismatches	Mismatches permitidos na extremidade 3'	2
--qcov_hsp_perc	Cobertura mínima do primer (%)	80
--amp_seq	Extrai sequências dos amplicons	(flag)
-m	Limiar de Tm (°C)	0 (desativado)
Parâmetros Avançados (Opcionais)
Parâmetro	Descrição	Padrão
--no_blast	Pula execução do BLAST (usa resultados existentes)	Desligado
--use_existing_db	Reutiliza banco BLAST existente	Desligado
--na	Concentração de Na⁺ (mM) para cálculo de Tm	50
-k / --pot	Concentração de K⁺ (mM) para cálculo de Tm	0
--tris	Concentração de Tris (mM) para cálculo de Tm	0
--mg	Concentração de Mg²⁺ (mM) para cálculo de Tm	0
--dntps	Concentração de dNTPs (mM) para cálculo de Tm	0
--saltcorr	Método de correção de sal (1, 2, ou 5)	5
Exemplo Completo
cmd
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
📁 Estrutura do Projeto
text
primer_blast_local/
├── primer_blast_local.py          # Script principal (linha de comando)
├── reformat.py                    # Formatação de arquivos
├── run_parse_blastn.py            # Execução do BLAST e parsing
├── primers.fasta                  # Arquivo de primers (exemplo)
├── extrair_amplicons.py           # Extrai sequências dos resultados CSV
├── comparar_amplicons.py          # Compara múltiplos amplicons (NCBI vs Local)
├── gerar_relatorio_amplicons.py   # Gera relatórios HTML individuais
└── Web/                           # Aplicação Web (Flask)
    ├── app.py                     # Servidor Flask
    ├── templates/
    │   └── index.html             # Interface web
    ├── uploads/                   # Uploads de genomas
    └── results/                   # Resultados das análises
📜 Scripts Auxiliares
extrair_amplicons.py
Extrai sequências de amplicons do CSV e gera arquivos FASTA.

cmd
# Extrair todos os amplicons em arquivos separados
python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo S_cerevisiae

# Extrair todos em um único arquivo
python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo S_cerevisiae --unico
comparar_amplicons.py
Compara múltiplos pares de sequências (NCBI vs Local) e gera relatório.

cmd
# Comparação única
python comparar_amplicons.py --ncbi ncbi.fasta --local local.fasta --org "S. cerevisiae"

# Múltiplas via arquivo CSV
python comparar_amplicons.py --batch comparacoes.csv --relatorio resultados.txt
gerar_relatorio_amplicons.py
Gera relatórios HTML individuais para cada organismo.

cmd
python gerar_relatorio_amplicons.py \
    --ncbi amplicon_ncbi.fasta \
    --local amplicon_local.fasta \
    --org "Saccharomyces cerevisiae" \
    --cromossomo "NC_001144.5"
📊 Validação
A ferramenta foi validada contra o NCBI Primer-BLAST para 12 organismos:

#	Organismo	Tipo	Produtos	Tamanho (Local)
1	Saccharomyces cerevisiae	Fungo	2	4.817 bp
2	Phaeodactylum tricornutum	Protozoário	1	4.886 bp
3	Rhizophagus irregularis	Fungo	10	4.675-4.704 bp
4	Dictyostelium discoideum	Protozoário	1	5.530 bp
5	Fusarium oxysporum	Fungo	3	4.480 bp
6	Neurospora crassa	Fungo	1	4.522 bp
7	Phytophthora infestans	Protozoário	33	4.023-5.471 bp
8	Candida albicans	Fungo	1	4.474 bp
9	Aspergillus fumigatus	Fungo	1	4.553 bp
10	Blastocystis hominis	Protozoário	21	4.345-4.492 bp
11	Babesia microti	Protozoário	2	4.765 bp
12	Sclerotinia sclerotiorum	Fungo	-	-
Resultado da Validação:

✅ Identidade média das sequências: > 99,9%

✅ Diferença média de tamanho: 1-2 bp (0,02-0,04%)

✅ 100% de concordância nos cromossomos/scaffolds identificados

📚 Referência
Latz, M. A. C., Grujic, V., Brugel, S., Lycken, J., John, U., Karlson, B., Andersson, A., & Andersson, A. F. (2022). Short‐ and long‐read metabarcoding of the eukaryotic rRNA operon: Evaluation of primers and comparison to shotgun metagenomics sequencing. Molecular Ecology Resources, 22, 2304–2318. https://doi.org/10.1111/1755-0998.13623
