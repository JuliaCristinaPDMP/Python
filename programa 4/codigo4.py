from Bio import SeqIO
from Bio.Seq import Seq

codon_table = {
    'AUG':'M', 'UGG':'W',
    'UUU':'F', 'UUC':'F',
    'UUA':'L', 'UUG':'L', 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'AUU':'I', 'AUC':'I', 'AUA':'I',
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
    'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S', 'AGU':'S', 'AGC':'S',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'UAU':'Y', 'UAC':'Y',
    'CAU':'H', 'CAC':'H',
    'CAA':'Q', 'CAG':'Q',
    'AAU':'N', 'AAC':'N',
    'AAA':'K', 'AAG':'K',
    'GAU':'D', 'GAC':'D',
    'GAA':'E', 'GAG':'E',
    'UGU':'C', 'UGC':'C',
    'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'UAA':'*', 'UAG':'*', 'UGA':'*'
}

def transcrever_dna_para_rna(seq):
    return str(seq).replace('T', 'U')

def traduzir_rna_para_proteina(rna_seq, frame):
    proteina = ''
    for i in range(frame, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        aminoacido = codon_table.get(codon, 'X') 
        proteina += aminoacido
    return proteina

dna_records = list(SeqIO.parse(r"C:\Users\juhdi\OneDrive\Área de Trabalho\Python\programa 4\Ecoli_Sakai_cds_from_genomic.fna", "fasta"))
rna_records = []
frames_proteinas = [[] for _ in range(6)]

for record in dna_records:
    seq_dna = str(record.seq)
    rna_seq = transcrever_dna_para_rna(seq_dna)
    rna_records.append(f">{record.id}_RNA\n{rna_seq}")

    for i in range(3):
        prot_seq = traduzir_rna_para_proteina(rna_seq, i)
        frames_proteinas[i].append(f">{record.id}_Proteína\n{prot_seq}")

    seq_dna_rev = str(Seq(seq_dna).reverse_complement())
    rna_seq_rev = transcrever_dna_para_rna(seq_dna_rev)
    for i in range(3):
        prot_seq = traduzir_rna_para_proteina(rna_seq_rev, i)
        frames_proteinas[i+3].append(f">{record.id}_Proteína\n{prot_seq}")

with open("Sakai_RNA.fasta", "w") as f:
    for r in rna_records:
        f.write(r + "\n")

for i in range(6):
    with open(f"Frame{i+1}.fasta", "w") as f:
        for p in frames_proteinas[i]:
            f.write(p + "\n")
