from Bio import SeqIO
from Bio.Seq import Seq

# Definições
start_codon = "ATG"
stop_codons = {"TAA", "TAG", "TGA"}
ra_number = "170348"

def find_orfs(seq, frame_offset, reverse=False):
    orfs = []
    seq_len = len(seq)
    frame = seq[frame_offset:]
    protein_counter = 1
    i = 0
    while i < len(frame) - 3:
        codon = frame[i:i+3]
        if codon == start_codon:
            for j in range(i+3, len(frame)-2, 3):
                stop_codon = frame[j:j+3]
                if stop_codon in stop_codons:
                    orf_dna = frame[i:j+3]
                    protein_seq = str(Seq(orf_dna).translate(to_stop=True))
                    if len(protein_seq) > 0:
                        start_pos = frame_offset + i + 1
                        end_pos = frame_offset + j + 3
                        if reverse:
                            start_pos, end_pos = seq_len - end_pos + 1, seq_len - start_pos + 1
                        orfs.append({
                            "protein": protein_seq,
                            "start": start_pos,
                            "end": end_pos,
                            "num": protein_counter
                        })
                        protein_counter += 1
                    i = j
                    break
        i += 3
    return orfs

def process_sequence(file_path, label_prefix, genbank_code):
    results = {}
    seq_record = SeqIO.read(file_path, "fasta")
    sequence = str(seq_record.seq).upper()

    # Frames diretos
    for frame in range(3):
        orfs = find_orfs(sequence, frame)
        results[f"frame{frame+1}"] = orfs

    # Frames reversos
    rev_seq = str(Seq(sequence).reverse_complement())
    for frame in range(3):
        orfs = find_orfs(rev_seq, frame, reverse=True)
        results[f"frame{frame+4}"] = orfs

    # Salvar arquivos
    for i in range(1, 7):
        frame_key = f"frame{i}"
        orf_list = results[frame_key]
        output_path = f"{label_prefix}_frame{i}_ativ5_{ra_number}.fasta"
        with open(output_path, "w") as f:
            for orf in orf_list:
                header = f">{genbank_code}, Frame {i}, proteína {orf['num']}, [location={orf['start']}..{orf['end']}], {ra_number}"
                f.write(header + "\n")
                f.write(orf["protein"] + "\n")

# ⚠️ Caminho absoluto completo para seus arquivos:
base_path = r"C:\Users\juhdi\OneDrive\Área de Trabalho\Python\programa 5"

process_sequence(f"{base_path}\\Ecoli Sakai sequence (1).fasta", "gc", "BA000007.3")
process_sequence(f"{base_path}\\Ecoli Sakai Plasmid 1 sequence.fasta", "p1", "AB011549.2")
process_sequence(f"{base_path}\\Ecoli Sakai plasmid 2 sequence.fasta", "p2", "AB011548.2")
