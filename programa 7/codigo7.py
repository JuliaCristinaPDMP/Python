from Bio import SeqIO
from Bio.Seq import Seq

def read_fasta(reads4):
    return [str(record.seq) for record in SeqIO.parse(r"C:\Users\juhdi\OneDrive\Área de Trabalho\Python\programa 7\reads4.fasta", "fasta")]

def find_overlap(a, b, min_length=1):
    max_len = 0
    for i in range(min_length, min(len(a), len(b)) + 1):
        if a[-i:] == b[:i]:
            max_len = i
    return max_len

def greedy_assembly(reads):
    while len(reads) > 1:
        max_overlap = -1
        best_pair = (0, 0)
        best_merged = ""
        
        for i in range(len(reads)):
            for j in range(len(reads)):
                if i != j:
                    overlap_len = find_overlap(reads[i], reads[j])
                    if overlap_len > max_overlap:
                        max_overlap = overlap_len
                        best_pair = (i, j)
                        best_merged = reads[i] + reads[j][overlap_len:]

        i, j = best_pair
        reads.pop(j)
        reads.pop(i)
        reads.append(best_merged)
    
    return reads[0]

def write_fasta(sequence, output_file, ra="170348"):
    with open(output_file, "w") as f:
        f.write(f">Contig montado utilizando Algoritmo Guloso - {ra}\n")
        f.write(sequence + "\n")

if __name__ == "__main__":
    reads = read_fasta("reads4.fasta")
    contig = greedy_assembly(reads)
    write_fasta(contig, "contig.fasta", ra="170348")