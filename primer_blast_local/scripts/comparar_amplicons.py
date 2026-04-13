#!/usr/bin/env python
"""
Script universal para comparar sequências de amplicons do NCBI e da ferramenta local.
Suporta múltiplas comparações de uma só vez e gera um relatório detalhado.

Uso:
    # Comparação única
    python comparar_amplicons.py --ncbi amplicon_ncbi.fasta --local amplicon_local.fasta --org "S. cerevisiae"
    
    # Múltiplas comparações via arquivo CSV de configuração
    python comparar_amplicons.py --batch comparacoes.csv
    
    # Múltiplas comparações via argumentos
    python comparar_amplicons.py --pares "ncbi1.fasta,local1.fasta,Org1" "ncbi2.fasta,local2.fasta,Org2"

Autor: Tiago
Data: 2026
"""

import argparse
import csv
import os
import sys
from datetime import datetime
from Bio import SeqIO
from Bio import Align

def ler_fasta(arquivo):
    """Lê um arquivo FASTA e retorna a sequência em maiúsculas."""
    if not os.path.exists(arquivo):
        raise FileNotFoundError(f"Arquivo não encontrado: {arquivo}")
    return str(SeqIO.read(arquivo, "fasta").seq).upper()

def comparar_sequencias(seq_ncbi, seq_local, organismo="Desconhecido"):
    """
    Compara duas sequências e retorna um dicionário com os resultados.
    """
    resultado = {
        "organismo": organismo,
        "tamanho_ncbi": len(seq_ncbi),
        "tamanho_local": len(seq_local),
        "diferenca_tamanho": abs(len(seq_ncbi) - len(seq_local)),
        "identidade": 0.0,
        "matches": 0,
        "alinhado": 0,
        "diferencas": [],
        "contida": False,
        "direcao_contida": "",
        "status": ""
    }
    
    # Verifica se uma sequência está contida na outra
    if seq_local in seq_ncbi:
        resultado["contida"] = True
        resultado["direcao_contida"] = "Local contida na NCBI"
        pos = seq_ncbi.find(seq_local)
        resultado["diferencas"] = [{
            "tipo": "contida",
            "posicao_ncbi": pos + 1,
            "tamanho_local": len(seq_local),
            "extra_inicio": seq_ncbi[:pos],
            "extra_fim": seq_ncbi[pos + len(seq_local):]
        }]
        resultado["identidade"] = 100.0
        resultado["matches"] = len(seq_local)
        resultado["alinhado"] = len(seq_local)
        resultado["status"] = "✅ 100% CONTIDA"
    elif seq_ncbi in seq_local:
        resultado["contida"] = True
        resultado["direcao_contida"] = "NCBI contida na Local"
        pos = seq_local.find(seq_ncbi)
        resultado["diferencas"] = [{
            "tipo": "contida",
            "posicao_local": pos + 1,
            "tamanho_ncbi": len(seq_ncbi),
            "extra_inicio": seq_local[:pos],
            "extra_fim": seq_local[pos + len(seq_ncbi):]
        }]
        resultado["identidade"] = 100.0
        resultado["matches"] = len(seq_ncbi)
        resultado["alinhado"] = len(seq_ncbi)
        resultado["status"] = "✅ 100% CONTIDA"
    else:
        # Faz alinhamento local para encontrar a melhor correspondência
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        alignments = aligner.align(seq_ncbi, seq_local)
        
        if alignments:
            best = alignments[0]
            aligned_ncbi = str(best[0])
            aligned_local = str(best[1])
            
            # Conta matches (excluindo gaps)
            matches = 0
            diffs = []
            for i, (a, b) in enumerate(zip(aligned_ncbi, aligned_local)):
                if a == b and a != '-' and b != '-':
                    matches += 1
                elif a != '-' or b != '-':
                    diffs.append({
                        "posicao": i + 1,
                        "ncbi": a,
                        "local": b
                    })
            
            aligned_length = sum(1 for a, b in zip(aligned_ncbi, aligned_local) if a != '-' or b != '-')
            identidade = (matches / aligned_length * 100) if aligned_length > 0 else 0.0
            
            resultado["identidade"] = identidade
            resultado["matches"] = matches
            resultado["alinhado"] = aligned_length
            resultado["diferencas"] = diffs[:50]  # Limita a 50 diferenças para não sobrecarregar
            
            if identidade >= 99.99:
                resultado["status"] = f"✅ {identidade:.2f}% idêntico"
            elif identidade >= 99.0:
                resultado["status"] = f"⚠️ {identidade:.2f}% similar"
            else:
                resultado["status"] = f"❌ {identidade:.2f}% divergente"
    
    return resultado

def gerar_relatorio(resultados, arquivo_saida=None):
    """Gera um relatório formatado com os resultados das comparações."""
    
    linhas = []
    linhas.append("=" * 80)
    linhas.append("🔬 RELATÓRIO DE COMPARAÇÃO DE AMPLICONS")
    linhas.append("=" * 80)
    linhas.append(f"Data/Hora: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    linhas.append(f"Total de comparações: {len(resultados)}")
    linhas.append("=" * 80)
    linhas.append("")
    
    # Estatísticas gerais
    contidas = sum(1 for r in resultados if r["contida"])
    identicas = sum(1 for r in resultados if r["identidade"] >= 99.99)
    similares = sum(1 for r in resultados if 99.0 <= r["identidade"] < 99.99)
    divergentes = sum(1 for r in resultados if r["identidade"] < 99.0)
    
    linhas.append("📊 ESTATÍSTICAS GERAIS")
    linhas.append("-" * 40)
    linhas.append(f"  ✅ Sequências idênticas (100% contidas): {contidas}")
    linhas.append(f"  ✅ Identidade >= 99.99%: {identicas}")
    linhas.append(f"  ⚠️ Similaridade 99.0% - 99.99%: {similares}")
    linhas.append(f"  ❌ Divergentes (< 99.0%): {divergentes}")
    linhas.append("")
    
    # Resultados por organismo
    for i, r in enumerate(resultados, 1):
        linhas.append("=" * 80)
        linhas.append(f"🔍 COMPARAÇÃO {i}: {r['organismo']}")
        linhas.append("=" * 80)
        linhas.append(f"  Tamanho NCBI:  {r['tamanho_ncbi']} bp")
        linhas.append(f"  Tamanho Local: {r['tamanho_local']} bp")
        linhas.append(f"  Diferença:     {r['diferenca_tamanho']} bp")
        linhas.append(f"  Status:        {r['status']}")
        
        if r["contida"]:
            linhas.append(f"  Tipo:          {r['direcao_contida']}")
            if r["diferencas"]:
                d = r["diferencas"][0]
                if d.get("extra_inicio"):
                    linhas.append(f"  Extra NCBI (início): {d['extra_inicio'][:50]}{'...' if len(d['extra_inicio']) > 50 else ''}")
                if d.get("extra_fim"):
                    linhas.append(f"  Extra NCBI (fim):    {d['extra_fim'][:50]}{'...' if len(d['extra_fim']) > 50 else ''}")
        else:
            linhas.append(f"  Identidade:    {r['identidade']:.2f}%")
            linhas.append(f"  Matches:       {r['matches']} / {r['alinhado']}")
            
            if r["diferencas"]:
                linhas.append("")
                linhas.append(f"  📍 Diferenças encontradas ({len(r['diferencas'])}):")
                for d in r["diferencas"][:10]:
                    linhas.append(f"     Posição {d['posicao']}: NCBI={d['ncbi']} | Local={d['local']}")
                if len(r["diferencas"]) > 10:
                    linhas.append(f"     ... e mais {len(r['diferencas']) - 10} diferenças")
        linhas.append("")
    
    linhas.append("=" * 80)
    linhas.append("✅ RELATÓRIO CONCLUÍDO")
    linhas.append("=" * 80)
    
    # Salva ou imprime
    conteudo = "\n".join(linhas)
    
    if arquivo_saida:
        with open(arquivo_saida, "w", encoding="utf-8") as f:
            f.write(conteudo)
        print(f"\n✅ Relatório salvo em: {arquivo_saida}")
    else:
        print(conteudo)
    
    return conteudo

def carregar_pares_csv(arquivo_csv):
    """Carrega pares de arquivos de um CSV."""
    pares = []
    with open(arquivo_csv, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2 and not row[0].startswith("#"):
                ncbi = row[0].strip()
                local = row[1].strip()
                org = row[2].strip() if len(row) > 2 else "Desconhecido"
                pares.append((ncbi, local, org))
    return pares

def main():
    parser = argparse.ArgumentParser(
        description="Compara sequências de amplicons do NCBI e da ferramenta local.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  # Comparação única
  python comparar_amplicons.py --ncbi S_cerevisiae_ncbi.fasta --local S_cerevisiae_local.fasta --org "S. cerevisiae"

  # Múltiplas comparações via argumentos
  python comparar_amplicons.py --pares "ncbi1.fasta,local1.fasta,Org1" "ncbi2.fasta,local2.fasta,Org2"

  # Múltiplas comparações via arquivo CSV
  python comparar_amplicons.py --batch comparacoes.csv

  # Salvar relatório em arquivo
  python comparar_amplicons.py --batch comparacoes.csv --relatorio resultados.txt

Formato do arquivo CSV (comparacoes.csv):
  ncbi_S_cerevisiae.fasta,local_S_cerevisiae.fasta,S. cerevisiae
  ncbi_P_tricornutum.fasta,local_P_tricornutum.fasta,P. tricornutum
  ncbi_R_irregularis.fasta,local_R_irregularis.fasta,R. irregularis
        """
    )
    
    # Modos de entrada
    grupo = parser.add_mutually_exclusive_group(required=True)
    grupo.add_argument("--ncbi", help="Arquivo FASTA da sequência do NCBI")
    grupo.add_argument("--pares", nargs="+", help="Pares de arquivos no formato 'ncbi,local,organismo'")
    grupo.add_argument("--batch", help="Arquivo CSV com lista de pares (formato: ncbi,local,organismo)")
    
    parser.add_argument("--local", help="Arquivo FASTA da sequência da ferramenta local (obrigatório com --ncbi)")
    parser.add_argument("--org", default="Desconhecido", help="Nome do organismo (usado com --ncbi)")
    parser.add_argument("--relatorio", "-r", help="Arquivo para salvar o relatório (se não especificado, imprime na tela)")
    
    args = parser.parse_args()
    
    pares = []
    
    if args.ncbi:
        if not args.local:
            parser.error("--local é obrigatório quando --ncbi é usado")
        pares = [(args.ncbi, args.local, args.org)]
    
    elif args.pares:
        for p in args.pares:
            partes = p.split(",")
            if len(partes) >= 2:
                ncbi = partes[0].strip()
                local = partes[1].strip()
                org = partes[2].strip() if len(partes) > 2 else "Desconhecido"
                pares.append((ncbi, local, org))
    
    elif args.batch:
        if not os.path.exists(args.batch):
            print(f"❌ ERRO: Arquivo de batch não encontrado: {args.batch}")
            sys.exit(1)
        pares = carregar_pares_csv(args.batch)
    
    if not pares:
        print("❌ Nenhum par de arquivos para comparar.")
        sys.exit(1)
    
    print(f"\n🔬 Iniciando comparação de {len(pares)} organismo(s)...\n")
    
    resultados = []
    for ncbi_file, local_file, org in pares:
        print(f"  Comparando: {org}...")
        try:
            seq_ncbi = ler_fasta(ncbi_file)
            seq_local = ler_fasta(local_file)
            resultado = comparar_sequencias(seq_ncbi, seq_local, org)
            resultados.append(resultado)
            print(f"    ✅ {resultado['status']}")
        except Exception as e:
            print(f"    ❌ ERRO: {e}")
            resultados.append({
                "organismo": org,
                "tamanho_ncbi": 0,
                "tamanho_local": 0,
                "diferenca_tamanho": 0,
                "identidade": 0.0,
                "matches": 0,
                "alinhado": 0,
                "diferencas": [],
                "contida": False,
                "direcao_contida": "",
                "status": f"❌ ERRO: {str(e)[:50]}"
            })
    
    print()
    gerar_relatorio(resultados, args.relatorio)


if __name__ == "__main__":
    main()