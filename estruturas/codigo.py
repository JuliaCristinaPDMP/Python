from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import math


def calcular_gc(seq):
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    if total == 0:
        return 0
    return ((g + c) / total) * 100


def calcular_tm(seq, gc_percent, na_conc=0.1):
    comprimento = len(seq)
    if comprimento == 0:
        return 0
    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (500 / comprimento)
    return tm


def processar_fasta(arquivo_fasta):
    dados = []

    for record in SeqIO.parse(arquivo_fasta, "fasta"):
        seq_id = record.id
        seq = str(record.seq).upper()

        comprimento = len(seq)
        a_count = seq.count('A')
        t_count = seq.count('T')
        c_count = seq.count('C')
        g_count = seq.count('G')

        gc_percent = calcular_gc(seq)
        tm = calcular_tm(seq, gc_percent)

        dados.append({
            "ID": seq_id,
            "Comprimento": comprimento,
            "A": a_count,
            "T": t_count,
            "C": c_count,
            "G": g_count,
            "GC_%": gc_percent,
            "Tm_C": tm
        })

    return pd.DataFrame(dados)


def salvar_csv(df, nome_saida="resultado.csv"):
    df.to_csv(nome_saida, index=False)
    print(f"Arquivo {nome_saida} salvo com sucesso!")


def plotar_gc_vs_tm(df):
    plt.figure(figsize=(10, 6))
    plt.scatter(df['Tm_C'], df['GC_%'], alpha=0.7)
    plt.xlabel("Temperatura de Anelamento (Tm) [°C]")
    plt.ylabel("Conteúdo GC (%)")
    plt.title("Conteúdo GC vs Temperatura de Anelamento")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("grafico_gc_tm.png")
    plt.show()


if __name__ == "__main__":
    arquivo_fasta = r"C:\Users\juhdi\Downloads\Ecoli_Sakai_cds_from_genomic.fna"  
    df = processar_fasta(arquivo_fasta)
    salvar_csv(df, "dados_saida.csv")
    plotar_gc_vs_tm(df)
