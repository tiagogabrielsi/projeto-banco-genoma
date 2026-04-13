#!/usr/bin/env python
"""
Script universal para extrair sequências de amplicons de arquivos CSV gerados pelo primer_blast_local.
Gera um arquivo FASTA individual para cada produto encontrado.

Uso:
    python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo organismo
    python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo organismo --pasta amplicons
    python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo organismo --unico

Autor: Tiago
Data: 2026
"""

import csv
import argparse
import os
import sys

def extrair_amplicons(arquivo_csv, prefixo="amplicon", pasta_saida=".", unico_arquivo=False):
    """
    Extrai sequências de amplicons de um arquivo CSV e salva em arquivos FASTA.
    
    Parâmetros:
    - arquivo_csv: Caminho do arquivo .pass.csv
    - prefixo: Prefixo para os nomes dos arquivos gerados
    - pasta_saida: Pasta onde os arquivos serão salvos
    - unico_arquivo: Se True, salva todos os amplicons em um único arquivo FASTA
    
    Retorna:
    - Número de amplicons extraídos
    """
    
    # Verifica se o arquivo existe
    if not os.path.exists(arquivo_csv):
        print(f"❌ ERRO: Arquivo não encontrado: {arquivo_csv}")
        return 0
    
    # Cria a pasta de saída se não existir
    if not os.path.exists(pasta_saida):
        os.makedirs(pasta_saida)
        print(f"📁 Pasta criada: {pasta_saida}")
    
    amplicons = []
    
    # Lê o arquivo CSV
    with open(arquivo_csv, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)  # Pula o cabeçalho
        
        for i, row in enumerate(reader, start=1):
            if len(row) >= 10:
                assay = row[0]
                cromossomo = row[3]
                tamanho = row[6]
                inicio = row[7]
                fim = row[8]
                sequencia = row[9]
                tm_amp = row[10] if len(row) > 10 else "N/A"
                
                amplicons.append({
                    "indice": i,
                    "assay": assay,
                    "cromossomo": cromossomo,
                    "tamanho": tamanho,
                    "inicio": inicio,
                    "fim": fim,
                    "tm_amp": tm_amp,
                    "sequencia": sequencia
                })
    
    if not amplicons:
        print("❌ Nenhum amplicon encontrado no arquivo CSV.")
        return 0
    
    # Salva os amplicons
    if unico_arquivo:
        # Salva todos em um único arquivo
        nome_arquivo = os.path.join(pasta_saida, f"{prefixo}_todos_amplicons.fasta")
        with open(nome_arquivo, "w", encoding="utf-8") as out:
            for amp in amplicons:
                cabecalho = f">{prefixo}_{amp['cromossomo']}_produto{amp['indice']} | tamanho={amp['tamanho']}bp | posicao={amp['inicio']}-{amp['fim']} | Tm={amp['tm_amp']}C"
                out.write(f"{cabecalho}\n")
                out.write(f"{amp['sequencia']}\n")
        print(f"✅ {len(amplicons)} amplicon(s) salvos em: {nome_arquivo}")
    else:
        # Salva cada amplicon em um arquivo separado
        for amp in amplicons:
            # Sanitiza o nome do cromossomo para ser válido como nome de arquivo
            cromossomo_limpo = amp['cromossomo'].replace("/", "_").replace("\\", "_").replace(":", "_")
            nome_arquivo = os.path.join(pasta_saida, f"{prefixo}_{cromossomo_limpo}_produto{amp['indice']}.fasta")
            
            with open(nome_arquivo, "w", encoding="utf-8") as out:
                cabecalho = f">{prefixo}_{amp['cromossomo']}_produto{amp['indice']} | tamanho={amp['tamanho']}bp | posicao={amp['inicio']}-{amp['fim']} | Tm={amp['tm_amp']}C"
                out.write(f"{cabecalho}\n")
                out.write(f"{amp['sequencia']}\n")
            
            print(f"  ✅ Produto {amp['indice']}: {amp['cromossomo']} ({amp['tamanho']} bp) → {nome_arquivo}")
    
    return len(amplicons)


def main():
    parser = argparse.ArgumentParser(
        description="Extrai sequências de amplicons de arquivos CSV do primer_blast_local e gera arquivos FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo S_cerevisiae
  python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo Babesia --pasta ./amplicons
  python extrair_amplicons.py --csv resultados__results.pass.csv --prefixo Rhizophagus --unico
        """
    )
    
    parser.add_argument("--csv", required=True, help="Arquivo CSV de resultados (ex: resultados__results.pass.csv)")
    parser.add_argument("--prefixo", required=True, help="Prefixo para os nomes dos arquivos gerados (ex: S_cerevisiae)")
    parser.add_argument("--pasta", default=".", help="Pasta de saída para os arquivos (padrão: diretório atual)")
    parser.add_argument("--unico", action="store_true", help="Salva todos os amplicons em um único arquivo FASTA")
    
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print("🔬 EXTRAÇÃO DE AMPLICONS - primer_blast_local")
    print("="*60)
    print(f"📄 Arquivo CSV: {args.csv}")
    print(f"🏷️  Prefixo: {args.prefixo}")
    print(f"📁 Pasta de saída: {args.pasta}")
    print(f"📦 Modo: {'Único arquivo' if args.unico else 'Arquivos separados'}")
    print("="*60 + "\n")
    
    num_extraidos = extrair_amplicons(
        arquivo_csv=args.csv,
        prefixo=args.prefixo,
        pasta_saida=args.pasta,
        unico_arquivo=args.unico
    )
    
    print("\n" + "="*60)
    print(f"✅ EXTRAÇÃO CONCLUÍDA: {num_extraidos} amplicon(s) processado(s)")
    print("="*60)


if __name__ == "__main__":
    main()