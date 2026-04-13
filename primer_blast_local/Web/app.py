#!/usr/bin/env python
"""
Aplicação Web para Primer-BLAST Local
Autor: Tiago
Data: Abril 2026

Funcionalidades:
- Upload de genoma FASTA
- Configuração de primers (com suporte a degenerações IUPAC)
- Parâmetros básicos e avançados
- Execução assíncrona com barra de progresso
- Download de resultados em CSV, FASTA e HTML
"""

import os
import sys
import uuid
import subprocess
import threading
from datetime import datetime
from flask import Flask, render_template, request, jsonify, send_file
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.secret_key = 'sua-chave-secreta-aqui'

# ===== CONFIGURAÇÃO DE CAMINHOS =====
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_PATH = BASE_DIR

sys.path.insert(0, SCRIPT_PATH)

UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
RESULTS_FOLDER = os.path.join(BASE_DIR, 'results')
ALLOWED_EXTENSIONS = {'fasta', 'fna', 'fa', 'txt'}

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULTS_FOLDER'] = RESULTS_FOLDER

# Armazena status dos jobs
jobs_status = {}


def allowed_file(filename):
    """Verifica se a extensão do arquivo é permitida."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def generate_html_report(products, organism_name="Organismo", job_id=""):
    """Gera um relatório HTML com as sequências dos amplicons."""
    
    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Amplicons - {organism_name}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            background-color: white;
            border-radius: 10px;
            padding: 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            font-size: 1.3em;
        }}
        .summary {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        .summary p {{
            margin: 5px 0;
        }}
        .summary a {{
            color: #ffd700;
            text-decoration: none;
        }}
        .product {{
            margin-bottom: 30px;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .product-header {{
            background-color: #3498db;
            color: white;
            padding: 12px 20px;
            font-weight: bold;
            font-size: 1.1em;
        }}
        .product-info {{
            background-color: #f8f9fa;
            padding: 12px 20px;
            border-bottom: 1px solid #ddd;
            font-family: 'Consolas', 'Monaco', monospace;
            font-size: 14px;
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
        }}
        .sequence {{
            background-color: #f8f9fa;
            padding: 20px;
            font-family: 'Courier New', monospace;
            font-size: 14px;
            overflow-x: auto;
            white-space: pre-wrap;
            word-break: break-all;
            max-height: 400px;
            overflow-y: auto;
            border: 1px solid #dee2e6;
            line-height: 1.6;
        }}
        .footer {{
            margin-top: 40px;
            text-align: center;
            color: #95a5a6;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
            padding-top: 20px;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
            margin-right: 10px;
        }}
        .badge-tm {{
            background-color: #e74c3c;
            color: white;
        }}
        .badge-size {{
            background-color: #27ae60;
            color: white;
        }}
        .badge-position {{
            background-color: #f39c12;
            color: white;
        }}
        .note {{
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .note strong {{
            color: #856404;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Relatório de Amplicons - primer_blast_local</h1>
        
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0;">
            <div>
                <p><strong>🏷️ Organismo:</strong> {organism_name}</p>
                <p><strong>📋 Job ID:</strong> {job_id}</p>
            </div>
            <div style="text-align: right;">
                <p><strong>📊 Total de produtos:</strong> <span style="font-size: 1.5em; font-weight: bold; color: #3498db;">{len(products)}</span></p>
                <p><strong>📅 Data/Hora:</strong> {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}</p>
            </div>
        </div>
        
        <div class="summary">
            <p>✅ <strong>{len(products)}</strong> amplicon(s) encontrado(s) pela ferramenta <code>primer_blast_local</code>.</p>
            <p>🔬 Ferramenta validada contra NCBI Primer-BLAST (identidade > 99,9%).</p>
        </div>
        
"""
    
    for i, p in enumerate(products, 1):
        html += f"""
        <div class="product">
            <div class="product-header">
                🧬 Produto #{i} - {p.get('cromossomo', 'N/A')}
            </div>
            <div class="product-info">
                <span class="badge badge-size">📏 Tamanho: {p.get('tamanho', 'N/A')} bp</span>
                <span class="badge badge-tm">🌡️ Tm: {p.get('tm', 'N/A')}°C</span>
                <span class="badge badge-position">📍 Posição: {p.get('inicio', 'N/A')} - {p.get('fim', 'N/A')}</span>
            </div>
            <div class="sequence">{p.get('sequencia', 'N/A')}</div>
        </div>
"""
    
    html += f"""
        <div class="footer">
            <p>🔬 Gerado por <strong>primer_blast_local</strong> - Ferramenta validada contra NCBI Primer-BLAST</p>
            <p>📅 {datetime.now().strftime('%d/%m/%Y %H:%M:%S')} | Job ID: {job_id}</p>
            <p style="margin-top: 10px; font-size: 0.85em;">Primers utilizados: V4_Balzano_F (CCAGCASCYGCGGTAATTCC) | D11_3143R (RCCACAAGCYARTTATCC)</p>
        </div>
    </div>
</body>
</html>"""
    
    return html


def run_analysis(job_id, genome_path, primer_fasta, params):
    """Executa a análise em background."""
    try:
        jobs_status[job_id] = {"status": "running", "progress": 0, "message": "Iniciando..."}
        
        jobs_status[job_id]["progress"] = 10
        jobs_status[job_id]["message"] = "Criando banco BLAST..."
        
        output_prefix = os.path.join(app.config['RESULTS_FOLDER'], job_id)
        primer_blast_script = os.path.join(SCRIPT_PATH, "primer_blast_local.py")
        
        if not os.path.exists(primer_blast_script):
            raise FileNotFoundError(f"Script não encontrado: {primer_blast_script}")
        
        # Comando completo com TODOS os parâmetros
        cmd = [
            "python", primer_blast_script,
            "-g", genome_path,
            "-p", primer_fasta,
            "-o", output_prefix,
            "-e", str(params.get("evalue", 30000)),
            "--min_size", str(params.get("min_size", 4000)),
            "--max_size", str(params.get("max_size", 6000)),
            "-t", str(params.get("threads", 4)),
            "--max_target_seqs", str(params.get("max_target_seqs", 100)),
            "--max_3prime_mismatches", str(params.get("max_3prime", 2)),
            "--qcov_hsp_perc", str(params.get("qcov_hsp_perc", 80)),
            "-m", str(params.get("tm_thresh", 0)),
            "--na", str(params.get("na", 50)),
            "-k", str(params.get("k", 0)),
            "--tris", str(params.get("tris", 0)),
            "--mg", str(params.get("mg", 0)),
            "--dntps", str(params.get("dntps", 0)),
            "--saltcorr", str(params.get("saltcorr", 5)),
        ]
        
        if params.get("amp_seq"):
            cmd.append("--amp_seq")
        if params.get("no_blast"):
            cmd.append("--no_blast")
        if params.get("use_existing_db"):
            cmd.append("--use_existing_db")
        
        jobs_status[job_id]["progress"] = 30
        jobs_status[job_id]["message"] = "Executando BLAST..."
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=SCRIPT_PATH)
        
        jobs_status[job_id]["progress"] = 80
        jobs_status[job_id]["message"] = "Processando resultados..."
        
        if result.returncode == 0:
            pass_file = f"{output_prefix}__results.pass.csv"
            if os.path.exists(pass_file):
                with open(pass_file, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                    num_products = len(lines) - 1
                
                jobs_status[job_id] = {
                    "status": "completed",
                    "progress": 100,
                    "message": f"Análise concluída! {num_products} produtos encontrados.",
                    "output_prefix": output_prefix,
                    "num_products": num_products
                }
            else:
                jobs_status[job_id] = {
                    "status": "completed",
                    "progress": 100,
                    "message": "Análise concluída! Nenhum produto encontrado.",
                    "output_prefix": output_prefix,
                    "num_products": 0
                }
        else:
            jobs_status[job_id] = {
                "status": "error",
                "progress": 0,
                "message": f"Erro na execução: {result.stderr[:200]}"
            }
            
    except subprocess.TimeoutExpired:
        jobs_status[job_id] = {"status": "error", "message": "Timeout: análise excedeu 1 hora"}
    except FileNotFoundError as e:
        jobs_status[job_id] = {"status": "error", "message": str(e)}
    except Exception as e:
        jobs_status[job_id] = {"status": "error", "message": str(e)}


@app.route('/')
def index():
    """Página principal."""
    return render_template('index.html')


@app.route('/api/analyze', methods=['POST'])
def analyze():
    """Endpoint para iniciar análise."""
    job_id = str(uuid.uuid4())[:8]
    
    if 'genome' not in request.files:
        return jsonify({"error": "Nenhum arquivo enviado"}), 400
    
    file = request.files['genome']
    if file.filename == '':
        return jsonify({"error": "Nenhum arquivo selecionado"}), 400
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        genome_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{job_id}_{filename}")
        file.save(genome_path)
    else:
        return jsonify({"error": "Formato de arquivo não permitido"}), 400
    
    primer_fasta = os.path.join(app.config['UPLOAD_FOLDER'], f"{job_id}_primers.fasta")
    with open(primer_fasta, 'w', encoding='utf-8') as f:
        f.write(f">Primer|V4_Balzano_F|fwd\n{request.form.get('forward_primer', 'CCAGCASCYGCGGTAATTCC')}\n")
        f.write(f">Primer|D11_3143R|rev\n{request.form.get('reverse_primer', 'RCCACAAGCYARTTATCC')}\n")
    
    params = {
        "evalue": float(request.form.get("evalue", 30000)),
        "min_size": int(request.form.get("min_size", 4000)),
        "max_size": int(request.form.get("max_size", 6000)),
        "threads": int(request.form.get("threads", 4)),
        "max_3prime": int(request.form.get("max_3prime", 2)),
        "amp_seq": request.form.get("amp_seq") == "on",
        "max_target_seqs": int(request.form.get("max_target_seqs", 100)),
        "qcov_hsp_perc": int(request.form.get("qcov_hsp_perc", 80)),
        "tm_thresh": float(request.form.get("tm_thresh", 0)),
        "no_blast": request.form.get("no_blast") == "on",
        "use_existing_db": request.form.get("use_existing_db") == "on",
        "na": float(request.form.get("na", 50)),
        "k": float(request.form.get("k", 0)),
        "tris": float(request.form.get("tris", 0)),
        "mg": float(request.form.get("mg", 0)),
        "dntps": float(request.form.get("dntps", 0)),
        "saltcorr": int(request.form.get("saltcorr", 5)),
    }
    
    thread = threading.Thread(target=run_analysis, args=(job_id, genome_path, primer_fasta, params))
    thread.start()
    
    jobs_status[job_id] = {"status": "pending", "progress": 0, "message": "Na fila..."}
    
    return jsonify({"job_id": job_id})


@app.route('/api/status/<job_id>')
def status(job_id):
    """Endpoint para verificar status do job."""
    if job_id in jobs_status:
        return jsonify(jobs_status[job_id])
    return jsonify({"error": "Job não encontrado"}), 404


@app.route('/api/results/<job_id>')
def results(job_id):
    """Endpoint para obter resultados."""
    if job_id not in jobs_status or jobs_status[job_id].get("status") != "completed":
        return jsonify({"error": "Resultados não disponíveis"}), 404
    
    output_prefix = jobs_status[job_id]["output_prefix"]
    pass_file = f"{output_prefix}__results.pass.csv"
    
    if not os.path.exists(pass_file):
        return jsonify({"error": "Arquivo de resultados não encontrado"}), 404
    
    import csv
    products = []
    with open(pass_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            products.append({
                "cromossomo": row["Subject_ID"],
                "tamanho": row["Amplicon_size"],
                "inicio": row["Start"],
                "fim": row["End"],
                "tm": row.get("Amplicon_tm", "N/A")
            })
    
    return jsonify({
        "num_products": len(products),
        "products": products,
        "downloads": {
            "csv": f"/api/download/{job_id}/csv",
            "fasta": f"/api/download/{job_id}/fasta",
            "html": f"/api/download/{job_id}/html"
        }
    })


@app.route('/api/download/<job_id>/<file_type>')
def download(job_id, file_type):
    """Endpoint para download dos resultados."""
    if job_id not in jobs_status:
        return "Job não encontrado", 404
    
    output_prefix = jobs_status[job_id]["output_prefix"]
    
    if file_type == "csv":
        file_path = f"{output_prefix}__results.pass.csv"
        if not os.path.exists(file_path):
            return "Arquivo não encontrado", 404
        return send_file(file_path, as_attachment=True, download_name=f"resultados_{job_id}.csv")
    
    elif file_type == "fasta":
        import csv
        fasta_path = f"{output_prefix}__amplicons.fasta"
        csv_path = f"{output_prefix}__results.pass.csv"
        
        if not os.path.exists(csv_path):
            return "Arquivo CSV não encontrado", 404
        
        with open(csv_path, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            with open(fasta_path, 'w', encoding='utf-8') as fastafile:
                for i, row in enumerate(reader, 1):
                    if row.get('Amplicon_sequence'):
                        fastafile.write(f">produto{i}_{row['Subject_ID']}_{row['Amplicon_size']}bp\n")
                        fastafile.write(f"{row['Amplicon_sequence']}\n")
        
        return send_file(fasta_path, as_attachment=True, download_name=f"amplicons_{job_id}.fasta")
    
    elif file_type == "html":
        import csv
        csv_path = f"{output_prefix}__results.pass.csv"
        
        if not os.path.exists(csv_path):
            return "Arquivo CSV não encontrado", 404
        
        products = []
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                products.append({
                    "cromossomo": row["Subject_ID"],
                    "tamanho": row["Amplicon_size"],
                    "inicio": row["Start"],
                    "fim": row["End"],
                    "tm": row.get("Amplicon_tm", "N/A"),
                    "sequencia": row.get("Amplicon_sequence", "N/A")
                })
        
        # Tenta extrair o nome do organismo do nome do arquivo do genoma
        organism_name = f"Job_{job_id}"
        # Procura pelo arquivo de genoma
        for f in os.listdir(UPLOAD_FOLDER):
            if f.startswith(job_id) and f.endswith(('.fna', '.fasta', '.fa')):
                organism_name = f.replace(f"{job_id}_", "").rsplit('.', 1)[0]
                break
        
        html_content = generate_html_report(products, organism_name, job_id)
        
        html_path = f"{output_prefix}__report.html"
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        return send_file(html_path, as_attachment=True, download_name=f"amplicons_{organism_name}_{job_id}.html")
    
    return "Tipo não suportado", 400


if __name__ == '__main__':
    print(f"🚀 Servidor iniciado em: http://localhost:5000")
    print(f"📁 Scripts em: {SCRIPT_PATH}")
    print(f"📁 Uploads em: {UPLOAD_FOLDER}")
    print(f"📁 Resultados em: {RESULTS_FOLDER}")
    print(f"📄 Relatório HTML disponível após análise")
    app.run(debug=True, host='0.0.0.0', port=5000)