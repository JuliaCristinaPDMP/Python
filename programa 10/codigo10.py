def smith_waterman_full(seq1, seq2, match_score=5, mismatch_score=-3, gap_penalty=-4):
    m, n = len(seq1), len(seq2)

    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    traceback_matrix = [['' for _ in range(m + 1)] for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq2[i - 1] == seq1[j - 1]:
                score = match_score
            else:
                score = mismatch_score

            diag = score_matrix[i - 1][j - 1] + score
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0, diag, up, left)

            if score_matrix[i][j] == diag and score_matrix[i][j] != 0:
                traceback_matrix[i][j] = '↖'
            elif score_matrix[i][j] == up and score_matrix[i][j] != 0:
                traceback_matrix[i][j] = '↑'
            elif score_matrix[i][j] == left and score_matrix[i][j] != 0:
                traceback_matrix[i][j] = '←'

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    print("(C)Matriz de escore:")
    print("     " + "  ".join(seq1))
    for i in range(n + 1):
        if i == 0:
            row_label = " "
        else:
            row_label = seq2[i - 1]
        print(row_label + " " + " ".join(f"{score_matrix[i][j]:>2}" for j in range(m + 1)))

    aligned_seq1 = []
    aligned_seq2 = []
    alignment_score = 0
    i, j = max_pos
    matches = mismatches = gaps = 0

    while score_matrix[i][j] != 0:
        direction = traceback_matrix[i][j]
        if direction == '↖':
            aligned_seq1.append(seq1[j - 1])
            aligned_seq2.append(seq2[i - 1])
            if seq1[j - 1] == seq2[i - 1]:
                alignment_score += match_score
                matches += 1
            else:
                alignment_score += mismatch_score
                mismatches += 1
            i -= 1
            j -= 1
        elif direction == '↑':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[i - 1])
            alignment_score += gap_penalty
            gaps += 1
            i -= 1
        elif direction == '←':
            aligned_seq1.append(seq1[j - 1])
            aligned_seq2.append('-')
            alignment_score += gap_penalty
            gaps += 1
            j -= 1

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    print("\n(D)Alinhamento ótimo:")
    print("Seq1:", "".join(aligned_seq1))
    print("      ", "".join("|" if a == b else " " for a, b in zip(aligned_seq1, aligned_seq2)))
    print("Seq2:", "".join(aligned_seq2))

    print("\n(E)Detalhamento do escore:")
    print(f"Matches:    {matches} × {match_score} = {matches * match_score}")
    print(f"Mismatches: {mismatches} × {mismatch_score} = {mismatches * mismatch_score}")
    print(f"Gaps:       {gaps} × {gap_penalty} = {gaps * gap_penalty}")
    print(f"\n→ Escore total: {alignment_score}")

    return alignment_score


# Exemplo dado pelo prof Shida na descrição da atividade
seq1 = "GAATTCAGTTA"
seq2 = "GGATCGA"

smith_waterman_full(seq1, seq2)
