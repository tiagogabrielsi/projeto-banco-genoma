#!/usr/bin/env python
"""
Script para gerar relatórios individuais de comparação de amplicons (NCBI vs Local).
Cada organismo terá um arquivo HTML e TXT com a comparação lado a lado.

Uso:
    python gerar_relatorio_amplicons.py --ncbi amplicon_ncbi.fasta --local amplicon_local.fasta --org "Babesia microti" --cromossomo "NC_027207.2" --pasta relatorios

Autor: Tiago
Data: 2026
"""

import argparse
import os
from datetime import datetime
from Bio import SeqIO

def ler_fasta(arquivo):
    """Lê um arquivo FASTA e retorna a sequência em maiúsculas."""
    if not os.path.exists(arquivo):
        raise FileNotFoundError(f"Arquivo não encontrado: {arquivo}")
    record = SeqIO.read(arquivo, "fasta")
    return str(record.seq).upper(), record.description

def formatar_sequencia(seq, tamanho_linha=80):
    """Formata uma sequência para exibição com linhas de tamanho fixo."""
    return "\n".join([seq[i:i+tamanho_linha] for i in range(0, len(seq), tamanho_linha)])

def encontrar_diferencas(seq1, seq2, nome1="NCBI", nome2="Local"):
    """Encontra e marca as diferenças entre duas sequências."""
    diffs = []
    min_len = min(len(seq1), len(seq2))
    
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            diffs.append({
                "posicao": i + 1,
                nome1: seq1[i],
                nome2: seq2[i]
            })
    
    # Verifica diferenças de tamanho
    if len(seq1) > len(seq2):
        for i in range(min_len, len(seq1)):
            diffs.append({
                "posicao": i + 1,
                nome1: seq1[i],
                nome2: "-"
            })
    elif len(seq2) > len(seq1):
        for i in range(min_len, len(seq2)):
            diffs.append({
                "posicao": i + 1,
                nome1: "-",
                nome2: seq2[i]
            })
    
    return diffs

def gerar_relatorio_txt(org, cromossomo, seq_ncbi, seq_local, desc_ncbi, desc_local, arquivo_saida):
    """Gera um relatório em formato TXT."""
    
    tamanho_ncbi = len(seq_ncbi)
    tamanho_local = len(seq_local)
    diferenca = abs(tamanho_ncbi - tamanho_local)
    
    # Verifica contenção
    contida = ""
    if seq_local in seq_ncbi:
        contida = "✅ A sequência LOCAL está 100% CONTIDA na sequência NCBI"
        pos = seq_ncbi.find(seq_local)
        extra_inicio = seq_ncbi[:pos]
        extra_fim = seq_ncbi[pos + len(seq_local):]
    elif seq_ncbi in seq_local:
        contida = "✅ A sequência NCBI está 100% CONTIDA na sequência LOCAL"
        pos = seq_local.find(seq_ncbi)
        extra_inicio = seq_local[:pos]
        extra_fim = seq_local[pos + len(seq_ncbi):]
    else:
        contida = "⚠️ As sequências NÃO estão contidas uma na outra"
        extra_inicio = ""
        extra_fim = ""
        diffs = encontrar_diferencas(seq_ncbi, seq_local)
    
    with open(arquivo_saida, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("🔬 RELATÓRIO DE COMPARAÇÃO DE AMPLICONS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Organismo:        {org}\n")
        f.write(f"Cromossomo:       {cromossomo}\n")
        f.write(f"Data/Hora:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("📊 RESUMO\n")
        f.write("-" * 40 + "\n")
        f.write(f"Tamanho NCBI:     {tamanho_ncbi} bp\n")
        f.write(f"Tamanho Local:    {tamanho_local} bp\n")
        f.write(f"Diferença:        {diferenca} bp ({diferenca/max(tamanho_ncbi, tamanho_local)*100:.2f}%)\n")
        f.write(f"Status:           {contida}\n")
        f.write("\n")
        
        if "CONTIDA" in contida:
            f.write("🔍 DETALHES DA CONTENÇÃO\n")
            f.write("-" * 40 + "\n")
            if extra_inicio:
                f.write(f"Bases extras no INÍCIO da NCBI: {extra_inicio}\n")
            if extra_fim:
                f.write(f"Bases extras no FIM da NCBI:    {extra_fim}\n")
            f.write("\n")
        else:
            f.write(f"🔍 DIFERENÇAS ENCONTRADAS ({len(diffs)})\n")
            f.write("-" * 40 + "\n")
            for d in diffs[:20]:
                f.write(f"  Posição {d['posicao']:6d}: NCBI={d['NCBI']} | Local={d['Local']}\n")
            if len(diffs) > 20:
                f.write(f"  ... e mais {len(diffs) - 20} diferenças\n")
            f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("🧬 SEQUÊNCIA NCBI\n")
        f.write("=" * 80 + "\n")
        f.write(f"> {desc_ncbi}\n")
        f.write(formatar_sequencia(seq_ncbi, 80) + "\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("🧬 SEQUÊNCIA LOCAL\n")
        f.write("=" * 80 + "\n")
        f.write(f"> {desc_local}\n")
        f.write(formatar_sequencia(seq_local, 80) + "\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("✅ FIM DO RELATÓRIO\n")
        f.write("=" * 80 + "\n")
    
    return arquivo_saida

def gerar_relatorio_html(org, cromossomo, seq_ncbi, seq_local, desc_ncbi, desc_local, arquivo_saida):
    """Gera um relatório em formato HTML."""
    
    tamanho_ncbi = len(seq_ncbi)
    tamanho_local = len(seq_local)
    diferenca = abs(tamanho_ncbi - tamanho_local)
    
    # Verifica contenção
    if seq_local in seq_ncbi:
        status = "✅ A sequência LOCAL está 100% CONTIDA na sequência NCBI"
        status_class = "success"
        pos = seq_ncbi.find(seq_local)
        extra_inicio = seq_ncbi[:pos]
        extra_fim = seq_ncbi[pos + len(seq_local):]
        contida = True
    elif seq_ncbi in seq_local:
        status = "✅ A sequência NCBI está 100% CONTIDA na sequência LOCAL"
        status_class = "success"
        pos = seq_local.find(seq_ncbi)
        extra_inicio = seq_local[:pos]
        extra_fim = seq_local[pos + len(seq_ncbi):]
        contida = True
    else:
        status = "⚠️ As sequências NÃO estão contidas uma na outra"
        status_class = "warning"
        contida = False
        diffs = encontrar_diferencas(seq_ncbi, seq_local)
    
    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Relatório: {org}</title>
    <style>
        body {{ font-family: 'Courier New', monospace; max-width: 1200px; margin: 0 auto; padding: 20px; background-color: #f5f5f5; }}
        .container {{ background-color: white; border-radius: 10px; padding: 30px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .summary {{ background-color: #ecf0f1; padding: 20px; border-radius: 8px; margin: 20px 0; }}
        .summary p {{ margin: 5px 0; }}
        .success {{ color: #27ae60; font-weight: bold; }}
        .warning {{ color: #e67e22; font-weight: bold; }}
        .sequence {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; font-family: 'Courier New', monospace; font-size: 14px; overflow-x: auto; white-space: pre-wrap; word-break: break-all; border: 1px solid #dee2e6; max-height: 400px; overflow-y: auto; }}
        .info {{ color: #7f8c8d; font-size: 0.9em; }}
        .badge {{ display: inline-block; padding: 3px 8px; border-radius: 4px; font-size: 0.8em; font-weight: bold; }}
        .badge-ncbi {{ background-color: #3498db; color: white; }}
        .badge-local {{ background-color: #2ecc71; color: white; }}
        .diff-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        .diff-table th, .diff-table td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
        .diff-table th {{ background-color: #34495e; color: white; }}
        .diff-table tr:hover {{ background-color: #f1f1f1; }}
        .footer {{ margin-top: 40px; text-align: center; color: #95a5a6; font-size: 0.9em; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 Relatório de Comparação de Amplicons</h1>
        
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0;">
            <div>
                <p><strong>Organismo:</strong> {org}</p>
                <p><strong>Cromossomo:</strong> {cromossomo}</p>
            </div>
            <div style="text-align: right;">
                <p class="info">Data/Hora: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
        </div>
        
        <h2>📊 Resumo</h2>
        <div class="summary">
            <p><span class="badge badge-ncbi">NCBI</span> <strong>Tamanho:</strong> {tamanho_ncbi} bp</p>
            <p><span class="badge badge-local">Local</span> <strong>Tamanho:</strong> {tamanho_local} bp</p>
            <p><strong>Diferença:</strong> {diferenca} bp ({diferenca/max(tamanho_ncbi, tamanho_local)*100:.2f}%)</p>
            <p class="{status_class}"><strong>Status:</strong> {status}</p>
"""
    
    if contida:
        html += f"""
            <p><strong>Detalhes da contenção:</strong></p>
            <ul>
                <li>Bases extras no INÍCIO da NCBI: <code>{extra_inicio if extra_inicio else '(nenhuma)'}</code></li>
                <li>Bases extras no FIM da NCBI: <code>{extra_fim if extra_fim else '(nenhuma)'}</code></li>
            </ul>
        """
    else:
        html += f"""
            <p><strong>Diferenças encontradas:</strong> {len(diffs)}</p>
        </div>
        
        <h2>🔍 Diferenças Detectadas</h2>
        <table class="diff-table">
            <tr>
                <th>Posição</th>
                <th>NCBI</th>
                <th>Local</th>
            </tr>
"""
        for d in diffs[:50]:
            html += f"""
            <tr>
                <td>{d['posicao']}</td>
                <td>{d['NCBI']}</td>
                <td>{d['Local']}</td>
            </tr>"""
        if len(diffs) > 50:
            html += f"""
            <tr>
                <td colspan="3" style="text-align: center; color: #7f8c8d;">... e mais {len(diffs) - 50} diferenças</td>
            </tr>"""
        html += """
        </table>
"""
    
    html += f"""
        <h2>🧬 Sequência NCBI</h2>
        <p class="info">{desc_ncbi}</p>
        <div class="sequence">{formatar_sequencia(seq_ncbi, 80)}</div>
        
        <h2>🧬 Sequência Local</h2>
        <p class="info">{desc_local}</p>
        <div class="sequence">{formatar_sequencia(seq_local, 80)}</div>
        
        <div class="footer">
            <p>✅ Relatório gerado automaticamente pelo script comparar_amplicons.py</p>
        </div>
    </div>
</body>
</html>"""
    
    with open(arquivo_saida, "w", encoding="utf-8") as f:
        f.write(html)
    
    return arquivo_saida

def main():
    parser = argparse.ArgumentParser(
        description="Gera relatórios individuais de comparação de amplicons (NCBI vs Local).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplo:
  python gerar_relatorio_amplicons.py --ncbi Babesia_ncbi.fasta --local Babesia_local.fasta --org "Babesia microti" --cromossomo "NC_027207.2"
        """
    )
    
    parser.add_argument("--ncbi", required=True, help="Arquivo FASTA da sequência do NCBI")
    parser.add_argument("--local", required=True, help="Arquivo FASTA da sequência da ferramenta local")
    parser.add_argument("--org", required=True, help="Nome do organismo")
    parser.add_argument("--cromossomo", default="Desconhecido", help="Identificador do cromossomo/scaffold")
    parser.add_argument("--pasta", default="relatorios", help="Pasta de saída para os relatórios")
    parser.add_argument("--formato", choices=["txt", "html", "ambos"], default="ambos", help="Formato do relatório")
    
    args = parser.parse_args()
    
    # Cria a pasta de saída
    if not os.path.exists(args.pasta):
        os.makedirs(args.pasta)
        print(f"📁 Pasta criada: {args.pasta}")
    
    # Sanitiza o nome do organismo para usar no nome do arquivo
    org_limpo = args.org.replace(" ", "_").replace(".", "").replace("/", "_")
    
    print(f"\n🔬 Gerando relatório para: {args.org}")
    print(f"   Cromossomo: {args.cromossomo}")
    print()
    
    try:
        seq_ncbi, desc_ncbi = ler_fasta(args.ncbi)
        seq_local, desc_local = ler_fasta(args.local)
        
        print(f"   ✅ NCBI:  {len(seq_ncbi)} bp")
        print(f"   ✅ Local: {len(seq_local)} bp")
        
        if args.formato in ["txt", "ambos"]:
            arquivo_txt = os.path.join(args.pasta, f"{org_limpo}_{args.cromossomo.replace(':', '_')}.txt")
            gerar_relatorio_txt(args.org, args.cromossomo, seq_ncbi, seq_local, desc_ncbi, desc_local, arquivo_txt)
            print(f"   📄 Relatório TXT: {arquivo_txt}")
        
        if args.formato in ["html", "ambos"]:
            arquivo_html = os.path.join(args.pasta, f"{org_limpo}_{args.cromossomo.replace(':', '_')}.html")
            gerar_relatorio_html(args.org, args.cromossomo, seq_ncbi, seq_local, desc_ncbi, desc_local, arquivo_html)
            print(f"   🌐 Relatório HTML: {arquivo_html}")
        
        print(f"\n✅ Relatório gerado com sucesso!")
        
    except Exception as e:
        print(f"\n❌ ERRO: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())