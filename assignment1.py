from Bio import SeqIO # to read the .fasta file
import pandas as pd # to process dataframes
import re # to perform regex searches 
from collections import Counter # to count amino acids in the proteome

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
fasta_file = "/Users/marlenfaf/Desktop/UofT - PhD/MMG1344H/Foundational-Computational-Biology-Course/orf_trans.fasta"

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
print("===== Section 4 =====")
print("===== Section 5 =====")
print("===== Section 6 =====")
print("===== Section 7 =====")