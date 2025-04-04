from Bio import SeqIO # to read the .fasta file
import pandas as pd # to process dataframes
import re # to perform regex searches 
from collections import Counter # to count amino acids in the proteome
import pdfplumber # to read the pdf file from assignment 1 
import seaborn as sns # to create and a the matrix as an image
import matplotlib.pyplot as plt # to create and a the matrix as an image

print("""
Edit the paths to match your file structure here""")
fasta_path = "/Users/marlenfaf/Desktop/UofT - PhD/MMG1344H/Foundational-Computational-Biology-Course/orf_trans.fasta"
assignment_pdf_path = "/Users/marlenfaf/Desktop/UofT - PhD/MMG1344H/FCBI_Assignment_1_Instructions_2025-1.pdf"

print("""
Part 1: Probability Analysis
""")

print("===== Section 1: Defining Amino Acid Probabilities  =====")
# Based on frequency data from S. cerevisiae

P_R = 0.06   # Probability of Arginine (R)
P_S = 0.093  # Probability of Serine (S)
P_T = 0.05   # Probability of Threonine (T)
P_ST = P_S + P_T  # Probability of either Serine or Threonine

# Section 2: Compute Probability of One 6-Mer Matching the Pattern
P_6mer = P_R * P_ST * P_ST

# Output the result
print(f"Probability of a random 6-mer matching the pattern: {P_6mer:.6f}")

# The protein of interest is made of 150 amino acids, so we need to find 
# how many 6-mers are there in the protein 

def calculate_k_mers(sequence_length, k):
    """Calculates and returns the number of k-mers that can exist in a 
    protein sequence of given length."""
    if sequence_length < k:
        return 0  # Not enough amino acids to form even a single k-mer
    return sequence_length - k + 1

number_of_6mers = calculate_k_mers(150,6)

print(f"The number of 6-mers in the Yfp1 150 residue-long sequence is {number_of_6mers}")

# Calculate the probability of at least one of the 145 6-mers to match the query sequence. 
P_match_one = P_6mer * 145 

print(f"The probability of at least one of the 145 6-mers matching the query sequence is {P_match_one:.6f}")


print("===== Section 2: Posterior probability that Yfp1 is phosphorylated by Cmk2  =====")

total_proteins = 4000  # Total number of proteins on the array of Ptacek et al
cmk2_targets = 9  # Number of proteins phosphorylated by Cmk2 according to Ptacek et al

# Prior probability that a protein is phosphorylated by Cmk2
P_phosphorylated_by_Cmk2 = cmk2_targets / total_proteins

# Given probability that Yfp1 matches the R-x(2)-[ST]-x-[ST] pattern (specifically for Yfp1)
P_pattern_Yfp1 = P_match_one  # This number was identified in Section 1 

# Likelihood that a protein matching the pattern is phosphorylated by Cmk2
P_pattern_given_phosphorylated_by_Cmk2 = 1  # Assuming that all Cmk2 substrates match the pattern

# Apply Bayes' Theorem to calculate the posterior probability
posterior_probability = (P_pattern_given_phosphorylated_by_Cmk2 * P_phosphorylated_by_Cmk2) / P_pattern_Yfp1

# Print the result
print(f"The posterior probability that Yfp1 is phosphorylated by Cmk2, given that it matches the pattern, is: {posterior_probability:.4f}")

print("""
Part 2: Probability Analysis
""")
print("===== Section 1 =====")

# Specifying the location of the fasta file
fasta_file = fasta_path 

# Initializing a list to hold the data
data = []

# Read and parse the FASTA file
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Append the sequence ID and sequence to the data list
        data.append([record.id, str(record.seq)])

# Storing the data list that was created above in a pandas dataframe 
fasta_file = pd.DataFrame(data, columns=["ID", "Sequence"])

print("This is a preview of the .fasta file")
print(fasta_file.head(5))

# Defining the pattern we wish to look for 
print("The pattern we are searching for is:")
pattern = r'R.{2}[ST].{1}[ST]'
print(pattern)

print("A function to find the pattern in each sequence in the fasta file is defined")
# Defining a function that will search for the pattern
def contains_pattern(sequence):
    """ Returns the number of times the specified pattern appears in a sequence 

    Args:
        sequence
    """
    return len(re.findall(pattern, sequence))

print("""The function is applied to read the 'Sequence' column in the .fasta file and note the 
number of matches in the 'Match_count' column""")

fasta_file["Match_Count"] = fasta_file["Sequence"].apply(contains_pattern)

print("This is how the dataframe looks like after we count the pattern matches:")
print(fasta_file.head(5))

# Counting how many matches were found 
print("Counting how many proteins match the pattern at least once...")
proteins_with_matches = (fasta_file["Match_Count"] > 0).sum()
print(f"The total number of proteins with at least one match to the pattern is {proteins_with_matches}")

total_proteins = len(fasta_file)
print(f"The total number of proteins in the genome is {total_proteins}")

fraction_with_match = proteins_with_matches / total_proteins
print(f"SECTION 1 - ANSWER: The fraction of proteins that match the pattern is {fraction_with_match:.4f}")

print("===== Section 2 =====")
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
print (f"First we define the 20 amino acids in the genome as follows: {amino_acids}")

# Initialize a Counter to hold the amino acid counts
amino_acid_count = Counter()

print("We use the counter to find how many appearances of each amino acid are there in each protein sequence.")
# Iteratig over each protein sequence to count the number of amino acids 
for seq in fasta_file["Sequence"]:
    amino_acid_count.update(seq)

print("We sum the counts for each amino acid across sequences")
# Summing the number of amino acids in all sequences to get the count for the entire proteome
total_amino_acids = sum(amino_acid_count.values())
print(f"The total appearances of each amino acid are as follows: {total_amino_acids}")

print("We iterate over the string of amino acids to compute the frequency for each one and saving it all in a dictionary")
# Calculate the frequency of each amino acid
amino_acid_frequencies = {aa: amino_acid_count[aa] / total_amino_acids for aa in amino_acids}

print("Converting the dictionary into a vector")
# Convert to a vector (list) for output
frequency_vector = [(aa, amino_acid_frequencies[aa]) for aa in amino_acids]

print(f"SECTION 2 - ANSWER: The computed amino acid frequencies for the yeast proteome are {frequency_vector}")

print("===== Section 3 =====")
# The pdf path that was edited in the beginning of the script is saved into a new variable
pdf_path = assignment_pdf_path

print("Reading the assigment instructions pdf file where the amino acid frequency table is found")
print("Extracting the text from page 2 of the pdf...")
# Extracting the text from page 2 
with pdfplumber.open(pdf_path) as pdf:
    page = pdf.pages[1]  # to adhere to the zero based indexing that pdfplumber is using 
    text_page_2 = page.extract_text()

print("Here's how the text looks like:")
print(text_page_2)

print("Splitting the text by line and spaces, while removing the '%' sign. Then the data is saved into a list")
lines = text_page_2.split("\n")  # Split text by line
data = []

for line in lines:
    parts = line.split()  # Split line by spaces
    if len(parts) == 2:  # Ensure valid row format
        amino_acid, frequency = parts[0], parts[1].strip("%")  # Remove '%' sign
        data.append([amino_acid, float(frequency)])

print("Converting the table into a pandas dataframe...")
part1_table = pd.DataFrame(data, columns=["Amino Acid", "Frequency (%)"])

print("This is how the table looks like:")
print(part1_table)

print("To compare the two frequency tables we will place them in the same pandas dataframe")

print("Converting the frequency vector from Section 2 into a pandas dataframe")
section2_frequencies = pd.DataFrame(frequency_vector, columns=['Amino Acid', 'Frequency'])

print("Here is how the frequencies from Section 2 look like:")
print(section2_frequencies)
print("""This table contains frequencies that are computed as a fraction of the total amino acids and not as percentage. 
It also uses the single letter nomenclature rather than the 3-letter convention. Let's change that.
""")

section2_frequencies['Frequency'] = (section2_frequencies['Frequency'] * 100).round(1)
amino_acid_map = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His", "I": "Ile",
    "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser",
    "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"
}

# Replace the 1-letter codes with the 3-letter codes in the 'Amino Acid' column
section2_frequencies['Amino Acid'] = section2_frequencies['Amino Acid'].map(amino_acid_map)

print("Here's how the refined table looks like:")
print(section2_frequencies)

print("Merging the 2 dataframes...")
tables_merged = pd.merge(part1_table, section2_frequencies, left_on='Amino Acid', right_on='Amino Acid', how='left')
tables_merged.columns = ['Amino Acid', 'Part 1 Table', 'Section 2']

print("Here's how the merged dataframe looks like:")
print(tables_merged)

print("Computing the difference in frequencies in a new column of the merged dataframe...")
tables_merged['Difference'] = tables_merged['Part 1 Table'] - tables_merged['Section 2']

print("Here's the dataframe with the difference between frequencies:")
print(tables_merged)

print("""SECTION 3 - ANSWER: The two organisms differ in their amino acid frequencies especially in terma of Alanine,
Arginine, Phenyalanine, Threonine and Tryptophan. After looking at the column "Differences" in the pandas dataset, we can 
observe that there are both positive and negative values. The positive values mean in S. cevirisae, these amino acids are 
more abundant in its proteins compared to the S288C yeast strain. The negative values mean that these amino acids are more 
abundant in the proteins of the S288C yeast strain. In general, different values can reflect specific adaptations of each 
organism, while similar values (Difference is close to zero) can be interpreted as having some similarity in the proteome
and hence sharing certain pathways such as metabolic pathways. There are several biological implications of these 
differences, such as different adaptations that each organism needs to make to adjust to its environment. This could dictate
the need for a different protein repertoire. Further, the organisms might not be sharing the same codon frequencies or 
might not use the same codons (due to differences in open reading frames, for example), which would lead to differences in
amino acid usage. 
""")

print("===== Section 4 =====")
print("A frequency matrix will be created from the vector of 20 frequencies calculated in Step 2")
print(f"We will iterate over this vector {frequency_vector}")

print("Separating the ferquency values from the amino acis keys in the vector.")
# Converting the vector to a dictionary
aa_freq = dict(frequency_vector)
print(f"The amino acid frequencies are: {aa_freq}")

# Getting a list of amino acids (for rows and columns)
aa_list = list(aa_freq.keys())
print(f"The amino acid codes are: {aa_list}")

print("Creating an empty dataframe where the column represents the first amino acid and the row represents the second one")
freq_matrix = pd.DataFrame(index=aa_list, columns=aa_list)

print(f"Here's how it looks like: {freq_matrix}")

print("Now we will compute the expected frequencies assuming independece, following this formula: P(AA1, AA2) = P(AA1) * P(AA2)")
for col in aa_list:  
    for row in aa_list:  
        freq_matrix.loc[row, col] = aa_freq[col] * aa_freq[row]

# Converting the values to floats to ensure all values are of the same data type
freq_matrix = freq_matrix.astype(float)

print(f"SECTION 4 - ANSWER: This is how the frequency matrix looks like: {freq_matrix}")

print("I would like to save this as a figure, coloring the matrix cells based on their probability value")

plt.figure(figsize=(12, 10))
sns.heatmap(freq_matrix, cmap='viridis')

plt.title("Section 4: Expected Frequencies of Di-Amino Acid Words", fontsize=14)
plt.xlabel("First Amino Acid")
plt.ylabel("Second Amino Acid")

# Saving as PNG
plt.tight_layout()
plt.savefig("section4_dipeptide_freq_matrix.png", dpi=300)

print("===== Section 5 =====")
print("Here we will count how often the amino acid pairs actually occur in the proteome.")

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

# Initiating the matrix where the di-amino acid combinations will be stored 
di_counts = pd.DataFrame(0, index=amino_acids, columns=amino_acids)

print(f"Here's how the emmpty matrix looks like: {di_counts}")

total_pairs = 0

print("We will iterate over each protein sequence in the .fasta file")
print("This process takes a while ...")
# Parse FASTA file
for record in SeqIO.parse(fasta_path, "fasta"):
    seq = str(record.seq) # converting the a.a sequence to a string 
    for i in range(len(seq) - 1): #Ensuring that we do not get out of the string bounds as we are computing a.a pairs 
        a1 = seq[i] # 1st a.a
        a2 = seq[i+1] # 2nd a.a
        if a1 in amino_acids and a2 in amino_acids: # validatig that the amino acids exist in our list
            di_counts.at[a1, a2] += 1 # increasing the count of the speicifc pair by 1
            total_pairs += 1 # increasing the total number of amino acid pairs by 1 

print(f"Here's how the matrix containing the count looks like: {di_counts}")

print("Dividing each di-amino acid count by the total counts to get its frequency...")
di_freqs = di_counts / total_pairs

print(f"SECTION 5 - ANSWER: Here's how the 20x20 matrix of the empirial frequencies looks like:{di_freqs}")

print("""I would like to save this as a figure next to the matrix from section 4, 
coloring the matrix cells based on their probability value""")

# Create a figure with 1 row and 2 columns of subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 7))  # width x height in inches

# Plot from section 4
sns.heatmap(freq_matrix, cmap='viridis', ax = axes[0])
axes[0].set_title("Section 4: Expected Frequencies of Di-Amino Acid Words")
axes[0].set_xlabel("1st Amino Acid")
axes[0].set_ylabel("2nd Amino Acid")

# Plot from section 5
sns.heatmap(di_freqs, cmap='viridis', ax = axes[1])
axes[1].set_title("Section 5: Empirical Frequencies of Di-Amino Acid Words (yeast genome)")
axes[1].set_xlabel("1st Amino Acid")
axes[1].set_ylabel("2nd Amino Acid")

# Saving the figure from both plots 
plt.savefig("section4_5_side_by_side.png", dpi=300)

print("===== Section 6 =====")

print("""To compute the conditional amino acid frequencies given that a certain amino acid 
is preceded by Isoleucine (I) or Glutamine (Q), a relevant fucntion is defined, 
where the probability of the amino acid of interest (x) being preceded by a defined amino acid (y)
is divided with the total occurrences of y. This is explained in detail in the .doc file 
of this assignment. """)

def calculate_conditional_frequencies(preceding_aa, di_counts, total_pairs):
    # Get the total counts of the preceding amino acid
    total_preceding_count = di_counts.loc[preceding_aa].sum()
    
    # Calculate the conditional frequency for each amino acid in the next position
    conditional_frequencies = {aa: di_counts.loc[preceding_aa, aa] / total_preceding_count
                               for aa in amino_acids}
    
    return conditional_frequencies

# Step 2: Specify the amino acids for conditioning
isoleucine_as_first = "I"  # Isoleucine
glutamine_as_first = "Q"  # Glutamine

print("After computing the frequencies, I will organize them in a pandas dataframe")

# Step 4: Calculate conditional frequencies for Isoleucine (I) and Glutamine (Q)
conditional_frequencies_I = calculate_conditional_frequencies(isoleucine_as_first, di_counts, total_pairs)
conditional_frequencies_Q = calculate_conditional_frequencies(glutamine_as_first, di_counts, total_pairs)

# Create a list of dictionaries with amino acid and frequency
data_I = [{'Preceding AA': 'Isoleucine', 'Following AA': aa, 'Frequency': freq} 
          for aa, freq in conditional_frequencies_I.items()]

data_Q = [{'Preceding AA': 'Glutamine', 'Following AA': aa, 'Frequency': freq} 
          for aa, freq in conditional_frequencies_Q.items()]

# Combine data for both Isoleucine and Glutamine
combined_data = data_I + data_Q

# Create a DataFrame from the combined data
section6_conditional_freqs = pd.DataFrame(combined_data)

print(f"SECTION 6 - ANSWER: The conditional probabilities for amino acids preceded by Isoleusine (I) and Glutamine (Q) are: {section6_conditional_freqs}")

print("===== Section 7 =====")