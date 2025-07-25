from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def identidade_score(a, b):
    return 0 if '-' in (a, b) else int(a == b)

substituicao = {
    'A': {'A':  1, 'T': -1, 'G': -1, 'C': -1, '-': 0},
    'T': {'A': -1, 'T':  1, 'G': -1, 'C': -1, '-': 0},
    'G': {'A': -1, 'T': -1, 'G':  1, 'C': -1, '-': 0},
    'C': {'A': -1, 'T': -1, 'G': -1, 'C':  1, '-': 0},
    '-': {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
}

def ler_fasta_duplo(caminho):
    with open(r"C:\Users\juhdi\OneDrive\Área de Trabalho\Python\programa 8\Seq1Seq2_corrigido.txt") as f:
        linhas = [l.strip() for l in f if l.strip() and not l.startswith(">")]
    meio = len(linhas) // 2
    return ''.join(linhas[:meio]), ''.join(linhas[meio:])

def calcular_scores(seq1, seq2):
    id_total = sub_total = 0
    for a, b in zip(seq1, seq2):
        id_total += identidade_score(a, b)
        sub_total += substituicao[a][b]
    comprimento = len(seq1)
    return {
        "identidade_total": id_total,
        "identidade_pct": id_total / comprimento * 100,
        "substituicao_total": sub_total,
        "substituicao_media": sub_total / comprimento
    }

def alinhamento_pairwise(seq1, seq2):
    seq1_limpa = seq1.replace("-", "")
    seq2_limpa = seq2.replace("-", "")
    alinhamentos = pairwise2.align.globalms(seq1_limpa, seq2_limpa, 1, -1, -0.5, -0.1)
    return alinhamentos[0] 
# Carregando arquivo
seq1, seq2 = ler_fasta_duplo("/mnt/data/Seq1Seq2_corrigido.txt")

res = calcular_scores(seq1, seq2)

print("🔍 RESULTADOS COM MATRIZES DADAS:")
print(f"- Identidade total: {res['identidade_total']} nucleotídeos")
print(f"- Identidade média: {res['identidade_pct']:.2f}%")
print(f"- Escore total (matriz substituição): {res['substituicao_total']}")
print(f"- Escore médio por base: {res['substituicao_media']:.2f}")

alinhamento = alinhamento_pairwise(seq1, seq2)
print("\n📌 ALINHAMENTO COM Bio.pairwise2:")
print(format_alignment(*alinhamento))
