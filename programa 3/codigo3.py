from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import math


def calcular_gc_proporcao(seq):
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    if total == 0:
        return 0
    return (g + c) / total  


def calcular_tm(seq, gc_proporcao, na_conc=0.1):
    comprimento = len(seq)
    if comprimento == 0:
        return 0
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * (gc_proporcao * 100) - (675 / comprimento)
    return tm


def processar_fasta(arquivo_fasta):
    dados_sequencias = []
    dados_gc = []
    dados_tm_gc = []

    for record in SeqIO.parse(arquivo_fasta, "fasta"):
        seq_id = record.id
        seq = str(record.seq).upper()

        a = seq.count('A')
        t = seq.count('T')
        c = seq.count('C')
        g = seq.count('G')
        total = len(seq)

        gc_proporcao = calcular_gc_proporcao(seq)
        tm = calcular_tm(seq, gc_proporcao, na_conc=0.1)  

        dados_sequencias.append({
            "Sequência": seq_id,
            "A": a,
            "T": t,
            "C": c,
            "G": g,
            "Total": total
        })

        dados_gc.append({
            "Sequência": seq_id,
            "GC_proporcao": f"{gc_proporcao:.3f}"
        })

        dados_tm_gc.append({
            "Tm_C": f"{tm:.3f}",
            "GC_percent": f"{gc_proporcao * 100:.2f}"
        })

    return dados_sequencias, dados_gc, dados_tm_gc

def salvar_csv_dados(dados, nome_arquivo):
    df = pd.DataFrame(dados)
    df.to_csv(nome_arquivo, index=False)
    print(f"Arquivo '{nome_arquivo}' salvo com sucesso.")


def plotar_gc_vs_tm(dados_tm_gc):
    df = pd.DataFrame(dados_tm_gc)
    df["Tm_C"] = df["Tm_C"].astype(float)
    df["GC_percent"] = df["GC_percent"].astype(float)

    plt.figure(figsize=(10, 6))
    plt.scatter(df['Tm_C'], df['GC_percent'], color='blue', alpha=0.7)
    plt.xlabel("Temperatura de Melting (°C)")
    plt.ylabel("Conteúdo GC (%)")
    plt.title("Conteúdo GC vs Temperatura de Melting", fontsize=14, y=1)
    plt.suptitle("Júlia Cristina", fontsize=10, y=0.92)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("grafico_conteudoGC_vs_Tm.png")
    plt.show()


if __name__ == "__main__":
    arquivo_fasta = r"C:\Users\juhdi\Downloads\Ecoli_Sakai_cds_from_genomic.fna"

    dados_seqs, dados_gc, dados_tm_gc = processar_fasta(arquivo_fasta)
    salvar_csv_dados(dados_seqs, "Dados das sequencias.csv")
    salvar_csv_dados(dados_gc, "Conteudo_GC.csv")
    salvar_csv_dados(dados_tm_gc, "Temperatura _x_GC.csv")
    plotar_gc_vs_tm(dados_tm_gc)
