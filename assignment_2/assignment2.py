# Assignment 2 
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# --- EDIT PATHS ---
print("""
Edit the paths to match your file structure here: """)
fasta_path = "/Users/marlenfaf/Desktop/UofT - PhD/MMG1344H/Foundational-Computational-Biology-Course/assignment_2/sample_input.fasta"

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

print("First, we have to import the fasta file. We do this ensuring that no comments are included (using fasta-blast):")

records = list(SeqIO.parse(fasta_path, "fasta-blast"))

# Extract ID and sequence into a list of tuples
data = [(record.id, str(record.seq)) for record in records]

# Converting to a pandas dataframe
fasta_df = pd.DataFrame(data, columns=["ID", "Sequence"])

print("This is a preview of the .fasta file:")
print(fasta_df.head(5))




print("===== Section 3 =====")

print(f"SECTION 2 - ANSWER: The computed amino acid frequencies for the yeast proteome are")