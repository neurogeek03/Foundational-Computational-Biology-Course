# Assignment 2 
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# --- EDIT PATHS ---
print("""
Edit the paths to match your file structure here: """) #TODO delete this path if not needed 
fasta_path = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/Foundational-Computational-Biology-Course/assignment_2/sample_input.fasta"

print("""
Part 1: Modifying the Smith-Waterman algorithm for 
""")

print("We have to make the following modifications to the Smith-Waterman algorithm")

print("""\
Align RNA1 (5’→3’) with the reverse complement of RNA2 (3’→5’).

Use the custom scoring scheme:
- G:C → +3
- A:U → +2
- G:U → +2
- Other mismatches → -1
- Gaps → -1 (only single-length gaps allowed)
""")

print("First, we have to import the fasta file. We do this ensuring that no comments such as (#Test case 1) are included (using fasta-blast):")
print("Following the assignment instructions, seq2 each case of the fasta file has to be reversed to stimulate a 3' -> 5' strand.")

# Defining a custom fasta parser as a function

def parse_fasta_pairs(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta-blast"))
    print(f"Found {len(records)} sequences in file.")
    
    # Verifying that the fasta sequences are an even number in total
    if len(records) % 2 != 0:
        raise ValueError("Expected an even number of sequences (seq1 + seq2 pairs).")
    
    # Defining the sequence pairs based on the order the file was read 
    pairs = []
    for i in range(0, len(records), 2):
        seq1 = str(records[i].seq)
        seq2 = str(records[i+1].seq)[::-1]  # Reverse for 3′→5′ hybridization
        pairs.append((seq1, seq2))
    
    return pairs

pairs = parse_fasta_pairs("test_input.fasta")

print("This is how the fasta sequences are stored in the 'pairs' variable:")
for idx, (s1, s2) in enumerate(pairs, 1):
    print(f"Test Case {idx}")
    print("Seq1 (5'–3') :", s1)
    print("Seq2 (3'–5') :", s2)
    print()

print("Defining a function that will return specific scores for different pairs")

def rna_score(b1, b2):
    if (b1, b2) in [('G', 'C'), ('C', 'G')]:
        return 3
    elif (b1, b2) in [('A', 'U'), ('U', 'A')]:
        return 2
    elif (b1, b2) in [('G', 'U'), ('U', 'G')]:
        return 2
    else:
        return -1

print("""Defining a function that will compute a 2d score matrix for a pair of RNA sequences using 
the scoring system defined in the function above. """)

# Essentially, one dimension will be the length of 
# sequence 1 (m) with an extra position included (m+1), which allows to account for the zero-based
# indexing of the alignment process. Similarly, the other dimension will be the lnegth of sequence 2 (n)
# plus an extra position (n+1). Each cell in the matrix will contain a score.

# The function fills in the matrix by computing the following: 

# MATCH
# The score of a cell. It takes into account the score of the previous diagonal cell in the matrix
# and computes the score of the current diagonal cell. 

# DELETE
# While calculating the cumulative score of the cells in the diagonal, the algorithm also considers the
# possibility of introducing a gap in seq2 (aligning a nucleotide from se1 to a gap). This value is 
# calculated by subtracting 1 (-1 penalty for gaps) to the score of the cell above (in the matrix diagonal)

# INSERT
# While calculating the cumulative score of the cells in the diagonal, the algorithm also considers the
# possibility of introducing a gap in seq1 (aligning a nucleotide from se1 to a gap). his value is 
# calculated by subtracting 1 (-1 penalty for gaps) to the score of the cell to the left (in the matrix diagonal)

# MAXIMUM SCORE
# This function tracks the highest score and its position(s) in the matrix. 
# It updates the max score and its position if the current calculated score is higher
# If the current score computed equals the max score, the value and its position are also recorded.

def local_align_rna(seq1, seq2):
    m, n = len(seq1), len(seq2)
    H = [[0]*(n+1) for _ in range(m+1)]  # the score matrix
    max_score = 0
    max_pos = []

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = H[i-1][j-1] + rna_score(seq1[i-1], seq2[j-1])
            delete = H[i-1][j] - 1
            insert = H[i][j-1] - 1

            H[i][j] = max(0, match, delete, insert)

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = [(i, j)]
            elif H[i][j] == max_score:
                max_pos.append((i, j))

    return H, max_pos, max_score

# The matrix H helps identify the best local alignment between two RNA sequences, considering matches, mismatches, 
# and gaps. The local alignment focuses on finding the highest-scoring subsequences rather than aligning the entire 
# sequences.

print("""When aligning 2 sequences, we can consider each sequence a row and the nucleotides 
that are laid on top of each other a column. The goal is to allow for consecutive gapa in a row 
but not gaps aligned to each other (in the same column). 
""")

print("""Using the matrix computed with the function above, we define the traceback function which 
performs the alignment based on the high-scoring subsequences""")

# The traceback function uses the max scores calculated in the matrix as the starting points 

def traceback_all(seq1, seq2, matrix, starts):
    alignments = []
    for start in starts:
        stack = [("", "", start[0], start[1], None)]  # (aligned1, aligned2, i, j, prev_move)

        while stack:
            a1, a2, i, j, prev_move = stack.pop()

            if matrix[i][j] == 0:
                alignments.append((a1[::-1], a2[::-1]))
                continue

            score = matrix[i][j]

            # Diagonal move: match/mismatch
            if i > 0 and j > 0:
                diag_score = matrix[i-1][j-1] + rna_score(seq1[i-1], seq2[j-1])
                if diag_score == score:
                    stack.append((a1 + seq1[i-1], a2 + seq2[j-1], i-1, j-1, 'diag'))

            # Left move: gap in seq2
            if j > 0 and prev_move != 'left':  # disallow repeated left if you want no extensions
                left_score = matrix[i][j-1] - 1
                if left_score == score:
                    if prev_move != 'up':  # prevent gap vs gap
                        stack.append((a1 + seq1[i-1] if i > 0 else "-", a2 + "-", i, j-1, 'left'))

            # Up move: gap in seq1
            if i > 0 and prev_move != 'up':
                up_score = matrix[i-1][j] - 1
                if up_score == score:
                    if prev_move != 'left':  # prevent gap vs gap
                        stack.append((a1 + "-", a2 + seq2[j-1] if j > 0 else "-", i-1, j, 'up'))

    return alignments











print("===== Section 3 =====")

print(f"SECTION 2 - ANSWER: The computed amino acid frequencies for the yeast proteome are")