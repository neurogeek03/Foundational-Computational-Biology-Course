# Assignment 2 
import sys
from collections import defaultdict
from collections import Counter 
from Bio import SeqIO
import pandas as pd
from io import StringIO
import zipfile
import gzip
import os
import shutil
import seaborn as sns 
import matplotlib.pyplot as plt

# --- EDIT PATHS ---
print("""
Edit the name of your fasta file here: """) 
fasta_file_name = "test_input.fasta"
hap_map_path = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/Foundational-Computational-Biology-Course/HapMap3.zip"

print("""
Part 1: Modifying the Smith-Waterman algorithm for RNA-RNA alignment, where the first sequence is at a 3'->5 direction and the 
second sequence at 5'->3' direction
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

print("===== Part 1: Defining Functions =====")
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
        seq1 = str(records[i].seq)[::-1]  # Reverse for 3′→5′ hybridization
        seq2 = str(records[i+1].seq)
        pairs.append((seq1, seq2))
    
    return pairs

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
#TODO INCLUDE SEQUENCE INDEXING INFO

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
# All optimal alignments are stored in the alignments list
# A stack is used because there may be more than 1 optimal path, and we wish to explore all optimal paths 

def traceback_all(seq1, seq2, matrix, starts):
    alignments = []
    for start in starts:
        stack = [("", "", start[0], start[1], None)]  # Initial stack: (aligned1, aligned2, i, j, prev_move)

        while stack:
            a1, a2, i, j, prev_move = stack.pop() # Takes the last thing we added to the stack & unpacks into variables

            # This line is only executed once we have reached the end of the alignment, so it is mostly skipped 
            if matrix[i][j] == 0: # if the score of a matrix cell is zero, then we are at the end of a local alignment
                alignments.append((a1[::-1], a2[::-1])) # The alignments that were read backwards, are reversed and saved
                continue

            score = matrix[i][j] 

            # Diagonal move: match/mismatch
            # Checks if we can move across the matrix diagonal to align 
            # If diag_score == score is true, it means that aligning seq1[i-1] with seq2[j-1] is part of the optimal 
            # alignment and results in the same score as the current matrix cell. In other words, this is a valid step 
            # in the traceback process, so the algorithm continues to move diagonally and adds the aligned 
            # characters to the growing sequence alignments (a1 and a2).
            if i > 0 and j > 0:
                diag_score = matrix[i-1][j-1] + rna_score(seq1[i-1], seq2[j-1])
                if diag_score == score:
                    stack.append((a1 + seq1[i-1], a2 + seq2[j-1], i-1, j-1, 'diag'))

            # Left move: gap in seq2
            if j > 0:
            # and prev_move != 'left':  # disallow repeated left if you want no extensions
                left_score = matrix[i][j-1] - 1
                if left_score == score:
                    if prev_move != 'up':  # prevent gap vs gap
                        stack.append((a1 + seq1[i-1] if i > 0 else "-", a2 + "-", i, j-1, 'left'))
                        # the algorithm will allow consecutive gaps in seq2 if the score matches.

            # Up move: gap in seq1
            if i > 0:
            # and prev_move != 'up':
                up_score = matrix[i-1][j] - 1
                if up_score == score:
                    if prev_move != 'left':  # prevent gap vs gap
                        stack.append((a1 + "-", a2 + seq2[j-1] if j > 0 else "-", i-1, j, 'up'))

    return alignments

# This function outputs the alignment in the same format as the sample output file 
def print_alignment(aln1, aln2, file = None):
    # reverse both sequences for biological convention: 3' to 5' on top, 5' to 3' on bottom
    top = aln1[::-1]
    bottom = aln2[::-1]

    # Match bars
    match_line = ""
    for a, b in zip(top, bottom):
        if (a, b) in [('G', 'C'), ('C', 'G'), ('A', 'U'), ('U', 'A'), ('G', 'U'), ('U', 'G')]:
            match_line += '|'
        else:
            match_line += ' '

    print("3' " + top + " 5'", file=file)
    print("   " + match_line, file=file)
    print("5' " + bottom + " 3'\n", file=file)


# Function that will save the output of the printing function above in a text file:
def write_alignment_output(case_id, score, alignments, file):
    print(f"# Test Case {case_id}", file=file)
    print(f"Score: {score}\n", file=file)
    for a1, a2 in alignments:
        print_alignment(a1, a2, file=file)

print("===== Part 1: Running the Program =====")

pairs = parse_fasta_pairs(fasta_file_name)

print("This is how the fasta sequences are stored in the 'pairs' variable:")
for idx, (s1, s2) in enumerate(pairs, 1):
    print(f"Test Case {idx}")
    print("Seq1 (5'–3') :", s1)
    print("Seq2 (3'–5') :", s2)
    print()

print("Iterating through the pairs and saving the sequences as a list for later use")
seq_pairs = []
for case_id, (seq1, seq2) in enumerate(pairs, 1):
    # Store the pair of sequences as a tuple (seq1, seq2)
    seq_pairs.append((seq1, seq2))

# Iterating through the sequence pairs list created above to run the program 
with open("output.txt", "w") as out_file:
    for case_id, (seq1, seq2) in enumerate(seq_pairs, 1):
        print(f"# Test Case {case_id}")
        
        # Running the local alignment function
        H, max_pos, max_score = local_align_rna(seq1, seq2)  # Passing seq1 and seq2 to the function
        
        # Running the traceback function to get the alignments
        alignments = traceback_all(seq1, seq2, H, max_pos)
        
        # Print or process the results (not saving to a file as you requested)
        print(f"Score: {max_score}")
        print(" ")

        for a1, a2 in alignments:
            print_alignment(a1, a2)

        write_alignment_output(case_id, max_score, alignments, file=out_file)


print("""
Part 2: Interpreting SNPs from the HapMap3 dataset
""")

print("===== Part 2: Task 1 - Unzipping data =====")

# Unzipping and loading the data 

print("Here's a preview of the HapMap3 dataset:")

#TODO uncomment before submission

# #Unzipping the HapMap3.zip file into the disk
# with zipfile.ZipFile(hap_map_path, 'r') as zip_ref:
#     zip_ref.extractall('hapmap3_data')

unzipped_files = os.listdir('hapmap3_data')
print(unzipped_files) # ['__MACOSX', 'HapMap3']

# Unzipping the files contained in the HapMap3 folder of the file we unzipped above 
input_dir = 'hapmap3_data/HapMap3/'
output_dir = 'hapmap3_unzipped'

os.makedirs(output_dir, exist_ok=True)

#TODO uncomment before submission

# for filename in os.listdir(input_dir):
#     if filename.endswith('.hmap.gz'):
#         input_path = os.path.join(input_dir, filename)
#         output_filename = filename[:-3]  # Remove '.gz'
#         output_path = os.path.join(output_dir, output_filename)

#         with gzip.open(input_path, 'rb') as f_in:
#             with open(output_path, 'wb') as f_out:
#                 shutil.copyfileobj(f_in, f_out)

#         print(f"Unzipped: {filename} → {output_filename}")


print("Previewing one of the .hapmap files:")
files = [f for f in os.listdir(output_dir) if f.endswith('.hmap')]

#TODO uncomment before submission

# # Pick the first one (just for preview)
# preview_file = os.path.join(output_dir, files[0])
# print(f"Previewing: {preview_file}")
# # Read and display first few rows
# df = pd.read_csv(preview_file, sep='\t')
# print(df.head())

print("===== Part 2: Task 1 - Functions =====")

# Defining necessary Functions

# This fucntion allows us to select for the SNPs of interest
def extract_snp_info(file_path, snps_of_interest):
    df = pd.read_csv(file_path, delim_whitespace=True, dtype=str, low_memory=False)
    print("Columns in the file:", df.columns.tolist())
    df_filtered = df[df['rs#'].isin(snps_of_interest)]
    return df_filtered

# Function to compute allele frequencies
def compute_allele_frequencies(snp_row):
    genotypes = snp_row.iloc[11:]  # skipping the metadata columns
    allele_counts = Counter()

    for gt in genotypes:
        if pd.isna(gt) or gt in ['NN', '--']:
            continue
        allele_counts[gt[0]] += 1
        allele_counts[gt[1]] += 1  # since it's diploid

    total = sum(allele_counts.values())
    frequencies = {allele: round(count / total, 4) for allele, count in allele_counts.items()}
    return frequencies

populations = ['ASW', 'CEU', 'CHB', 'LWK', 'MKK', 'YRI']
print(f"We will look for the SNPs in the following populations {populations}")

# TODO UNCOMMENT

# # print("===== Part 2: Task 1 - Running the program =====")

# # Defining the input directory 
# input_dir = 'hapmap3_unzipped'
# hapmap_files = [
#     os.path.join(input_dir, fname)
#     for fname in os.listdir(input_dir)
#     if any(pop in fname for pop in populations)
# ]

# print("Selected HapMap3 files:")
# for f in hapmap_files:
#     print(" →", os.path.basename(f))


# # Defining the SNPs of interest and the populations we wish to search for them 
# snps_of_interest = ['rs683', 'rs910']
# results = {}

# for pop in populations:
#     file_path = f'hapmap3_unzipped/{pop}.hmap'
#     snp_data = extract_snp_info(file_path, snps_of_interest)
#     print(snp_data)
#     results[pop] = {}

#     for _, row in snp_data.iterrows():
#         rsid = row['rs#']
#         freqs = compute_allele_frequencies(row)
#         is_fixed = len(freqs) == 1
#         fixed_allele = list(freqs.keys())[0] if is_fixed else None
#         results[pop][rsid] = {
#             'frequencies': freqs,
#             'is_fixed': is_fixed,
#             'fixed_allele': fixed_allele
#         }

# # Displaying the results per population 
# for pop in results:
#     print(f'\nPopulation: {pop}')
#     for snp in snps_of_interest:
#         data = results[pop].get(snp, {})
#         if not data:
#             print(f"  {snp}: Not found")
#         else:
#             print(f"  {snp}: Fixed? {data['is_fixed']}, Allele Frequencies: {data['frequencies']}, Fixed Allele: {data['fixed_allele']}")

# # Creating a figure of the output for easier visualization & interpretation
# records = []
# for pop in results:
#     for snp in snps_of_interest:
#         data = results[pop].get(snp, {})
#         if data and 'frequencies' in data:
#             for allele, freq in data['frequencies'].items():
#                 records.append({
#                     'Population': pop,
#                     'SNP': snp,
#                     'Allele': allele,
#                     'Frequency': freq
#                 })

# df = pd.DataFrame(records)

# # Optional: Sort populations consistently
# df['Population'] = pd.Categorical(df['Population'], categories=['ASW', 'CEU', 'CHB', 'LWK', 'MKK', 'YRI'], ordered=True)

# # Creating the plot using seaborn
# sns.set(style="whitegrid")
# g = sns.catplot(
#     data=df,
#     x='Population', y='Frequency', hue='Allele',
#     col='SNP', kind='bar', height=5, aspect=1,
#     palette='pastel', ci=None
# )
# g.set_titles('SNP: {col_name}')
# g.set_axis_labels('Population', 'Allele Frequency')
# g.add_legend(title='Allele')
# plt.tight_layout()

# # Saving the plot in my current working directory
# plot_path = os.path.join(os.getcwd(), "allele_frequencies_by_population.png")
# g.savefig(plot_path, dpi=300)

# print(f"Plot saved to: {plot_path}")

print("end of Part 2, Task 1")

print("===== Part 2: Task 2 - Obtaining the allele frequencies for the given subpopulations =====")
populations = ['ASW', 'CHD', 'GIH', 'LWK', 'MEX', 'MKK', 'TSI']
print(f"We will look for the SNPs in the following populations {populations}")

snps_of_interest = ['rs683', 'rs910']
results = {}

for pop in populations:
    file_path = f'hapmap3_unzipped/{pop}.hmap'
    snp_data = extract_snp_info(file_path, snps_of_interest)
    print(snp_data)
    results[pop] = {}

    for _, row in snp_data.iterrows():
        rsid = row['rs#']
        freqs = compute_allele_frequencies(row)
        is_fixed = len(freqs) == 1
        fixed_allele = list(freqs.keys())[0] if is_fixed else None
        results[pop][rsid] = {
            'frequencies': freqs,
            'is_fixed': is_fixed,
            'fixed_allele': fixed_allele
        }

# Displaying the results per population 
for pop in results:
    print(f'\nPopulation: {pop}')
    for snp in snps_of_interest:
        data = results[pop].get(snp, {})
        if not data:
            print(f"  {snp}: Not found")
        else:
            print(f"  {snp}: Fixed? {data['is_fixed']}, Allele Frequencies: {data['frequencies']}, Fixed Allele: {data['fixed_allele']}")

# Defining a simple heterozygosity function
def heterozygosity(pA):
    return 1 - (pA**2 + (1 - pA)**2)

# Looping through the existing results and computing heterozygosity
for pop in results:
    for snp in results[pop]:
        freqs = results[pop][snp]['frequencies']

        if len(freqs) == 2:
            # Use frequency of 'A' if available, otherwise we pick the first allele available
            pA = freqs.get('A', list(freqs.values())[0])
            hs = heterozygosity(pA)

        # Store the computed heterozygosity in results
        results[pop][snp]['heterozygosity'] = hs

# Printing the results as a table: 
for snp in snps_of_interest:
    print(f"\n### Results for {snp}\n")
    print("| Population | Allele Frequencies | Heterozygosity (H_S) |")
    print("|------------|--------------------|------------------------|")
    
    for pop in results:
        snp_data = results[pop].get(snp)
        if snp_data:
            freqs = snp_data['frequencies']
            hs = snp_data['heterozygosity']
            freqs_str = ', '.join([f"{allele}: {freq:.4f}" for allele, freq in freqs.items()])
            print(f"| {pop} | {freqs_str} | {hs:.4f} |")
        else:
            print(f"| {pop} | Not found | - |")

# Calculating and displaying the mean allele frequencies

# Initializing a dictionary that will hold the mean frequencies 
mean_frequencies_per_snp = {}

for snp in snps_of_interest:
    allele_totals = defaultdict(float)
    allele_counts = defaultdict(int)

    for pop in results:
        snp_data = results[pop].get(snp)
        if snp_data:
            freqs = snp_data['frequencies']
            for allele, freq in freqs.items():
                allele_totals[allele] += freq
                allele_counts[allele] += 1  # tracking how many populations had this allele

    # Computing mean frequencies 
    mean_freqs = {
        allele: allele_totals[allele] / allele_counts[allele]
        for allele in allele_totals
    }

    # Saving the frequency values as a dictionary
    mean_freq_a = mean_freqs.get('A', 0.0)
    mean_freq_c = mean_freqs.get('C', 0.0)

    # Optionally save to dictionary for later access
    mean_frequencies_per_snp[snp] = {
        'A': mean_freq_a,
        'C': mean_freq_c
    }

    # Printing the results as a table: 
    print(f"\n### Mean Allele Frequencies for {snp} (across populations)\n")
    print("| Allele | Mean Frequency |")
    print("|--------|----------------|")
    for allele in sorted(mean_freqs.keys()):
        print(f"| {allele} | {mean_freqs[allele]:.4f} |")

# Calculating the total heterozygosity and the mean subpopulation-specific heterozygosity:
for snp in snps_of_interest:
    pA_mean = mean_frequencies_per_snp[snp].get('A', 0.0)
    
    # Calculating H_T
    ht = heterozygosity(pA_mean)

    # 3. Calculate each population's H_S for the SNP
    hs_values = []
    for pop in results:
        snp_data = results[pop].get(snp)
        if snp_data:
            freqs = snp_data['frequencies']
            if len(freqs) == 2:  # only compute for biallelic SNPs
                pA = freqs.get('A', list(freqs.values())[0])
                hs_values.append(heterozygosity(pA))
            else:
                hs_values.append(0.0)

    # 4. Compute mean H_S
    hs_mean = sum(hs_values) / len(hs_values) if hs_values else 0.0

    # 5. Compute F_ST
    fst = (ht - hs_mean) / ht if ht > 0 else 0.0

    # 6. Print results
    print(f"\n### Fixation Index Summary for {snp}")
    print(f"H_T       : {ht:.4f}")
    print(f"H_S(mean) : {hs_mean:.4f}")
    print(f"F_ST      : {fst:.4f}")



print("===== Part 2: Task 2 - Running the program =====")