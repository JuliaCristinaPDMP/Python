from Bio import SeqIO
import csv

files = {
    "Genoma_Completo": "Ecoli Sakai sequence.fasta",
    "Plasmidio_1": "Ecoli Sakai Plasmid 1 sequence (1).fasta",
    "Plasmidio_2": "Ecoli Sakai plasmid 2 sequence (2).fasta"
}

def contar_nucleotideos(file_path):
    contagem = {"A": 0, "C": 0, "G": 0, "T": 0}
    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper()
        for nucleotideo in contagem.keys():
            contagem[nucleotideo] += seq.count(nucleotideo)
    return contagem

for nome, caminho in files.items():
    contagem = contar_nucleotideos(caminho)
    total = sum(contagem.values())
    output_file = f"{nome}_contagem.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([nome])
        writer.writerow(["A", "C", "G", "T", "Total"])
        writer.writerow([contagem["A"], contagem["C"], contagem["G"], contagem["T"], total])

print("Contagem concluída! Arquivos CSV gerados.")


