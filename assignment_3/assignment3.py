import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from scipy.stats import expon
from statsmodels.stats.multitest import multipletests

print("""
Part 1: Constructing and implementing the k-means algorithm.
""")

print("===== Question 1: Importing & Standardizing data, Defining the k-means function =====")

# Importing the dataset 
project_path = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/ass_3"
data_path = os.path.join(project_path, "Part1_Biase_2014.csv")
data = pd.read_csv(data_path, index_col=0)

# Previewing the dataset
print(data.shape) # this is a matrix with genes as rows, where each sample is a row
#(6812 rows , 40 columns)
print(data.head()) 

# Transposing the data to keep in line with the required input of standard clustering algorithms
sc_data = data.T.values # making this an array, meaning it loses some dataframe properties (like head)
print(sc_data.shape) #(40 rows, 6812 columns)

# Manually standardizing the data 
# Standardizing each gene (column) - Z-score normalization
means = np.mean(sc_data, axis=0) # compute across rows
stds = np.std(sc_data, axis=0) # compute across rows
# Implementing it on the data matrix 
standardized_data = (sc_data - means) / stds
# Verifying that it works
st_mean = np.mean(standardized_data, axis=0)  # close to 0
print(f"The mean of the data is now {st_mean}")
st_stdev = np.std(standardized_data, axis=0)  # 1
print(f"The standard deviation of the data is now {st_stdev}")

# Creating a manual function for k-means clustering 

# Defining Eucledian distance 
def eucledian_distance(a,b):
    return np.linalg.norm(a - b)

def kmeans(X, k=4, max_iters=100, tol=1e-4):
    """
    Performs K-Means clustering on the input data using Euclidean distance.

    Parameters:
    -----------
    X : np.ndarray
        The input data of shape (n_samples, n_features), where each row is a sample and each column is a feature (gene).
        For scRNA-seq data, rows are blastomeres and columns are gene expression levels.
    k : int, optional
        The number of clusters to form (default is 4).
    max_iters : int, optional
        The maximum number of iterations to run the algorithm (default is 100).
    tol : float, optional
        The tolerance for convergence. If the centroids move less than `tol`, the algorithm stops (default is 1e-4).

    Returns:
    --------
    labels : np.ndarray
        An array of shape (n_samples,) indicating the cluster assignment for each sample.
    centroids : np.ndarray
        An array of shape (k, n_features) representing the final centroid positions.

    Notes:
    ------
    This function implements the K-Means algorithm from scratch without using pre-built libraries such as scikit-learn.
    It uses Euclidean distance and standard iterative updates until convergence or max iterations.
    """

    n_samples, n_features = X.shape # samples = rows, genes/features = columns
    centroids = X[np.random.choice(n_samples, k, replace = False)] # randomly pick k initial steroids, avoiding duplicates 

    for iteration in range(max_iters): # steps below will be repeated accoding to the max_iters value
        clusters = [[]for _ in range(k)] # initializing k empty lists to hold point indices per cluster 
        # Step 1: cluster assignment 
        for i, x in enumerate(X):
            distances = [eucledian_distance(x, c) for c in centroids] # compute distance to all centroids 
            cluster_idx = np.argmin(distances) # finding the closest centroid by getting the index of the smallest value in the distances list 
            clusters[cluster_idx].append(i) # appending the index determined above to the list of clusters  
        # Step 2: Recomputing centroids 
        new_centroids = np.zeros_like(centroids) # initializing a zero-array that will store updated centroids 
        for idx, cluster in enumerate(clusters):
            if cluster:  # checking that it's not empty
                new_centroids[idx] = np.mean(X[cluster], axis=0) # mean of all samples assigned to a specific cluster index
        # Step 3: Check convergence 
        if np.allclose(centroids, new_centroids, atol=tol): # if the centroids moved less than the tolerance, stop
            break
        centroids = new_centroids # updating the centroids for the next iteration (if any)

    # Create final labels
    labels = np.zeros(n_samples) # zero-array of cluster labels
    for cluster_idx, cluster in enumerate(clusters): 
        for sample_idx in cluster:
            labels[sample_idx] = cluster_idx # filling in the array by through each cluster's assigned samples
    return labels, centroids # output: a cluster label for each sample, final centroid positions

# Running the algorithm on the standardized data 
np.random.seed(42)

# Run K-means
labels, centroids = kmeans(standardized_data, k=4)

# Initialize dictionary to hold sample names per cluster
clusters = {i: [] for i in range(4)}

# Map each sample to its cluster
for idx, label in enumerate(labels.astype(int)):
    clusters[label].append(data.columns[idx])  # assuming sample names are in data.index

# Print results
print("K-means clustering results (seed=42):")
for cid in range(4):
    print(f"\nCluster {cid} ({len(clusters[cid])} samples):")
    for sample in clusters[cid]:
        print(f"  - {sample}")

print("===== Question 2: Running K-means 10x using different seeds =====")
# Defining a function that implements the within-cluster sum of squares to evaluate clustering quality
def compute_wcss(X, labels, centroids, k):
    """
    Compute the Within-Cluster Sum of Squares (WCSS) for one cluster.

    Parameters:
    -----------
    - cluster_data: ndarray of shape (n_samples_in_cluster, n_features)

    Returns:
    --------
    - wcss: float
    """
    total_wcss = 0
    for cluster_id in range(k):
        # Select data points in this cluster
        cluster_points = X[labels == cluster_id]
        centroid = centroids[cluster_id]

        # Sum squared distances to the centroid
        squared_distances = np.sum((cluster_points - centroid) ** 2, axis=1)
        total_wcss += np.sum(squared_distances)
    return total_wcss

# We can automate this process by looping over using the function and recording the results for 10 times ß

cluster_matrix = np.zeros((40, 10), dtype=int) # initializing empty matrix to hold the results 
results = []

for seed in range(10):
    np.random.seed(seed) # setting ramdom seed

    labels, centroids = kmeans(standardized_data, k=4) # running the function

    # Ensuring that the seed is different every time and the centroids are also different: 
    print(f"Seed {seed}, First centroid: {centroids[0][:5]}")

    cluster_matrix[:, seed] = labels.astype(int)
    
    # computing wcss
    wcss = compute_wcss(standardized_data, labels, centroids, k=4)

    results.append({
        'seed': seed,
        'labels': labels,
        'wcss': wcss
    })

# Create DataFrame for labeling
cluster_df = pd.DataFrame(cluster_matrix, index=data.columns, columns=[f"Seed {s}" for s in range(10)])

# ============================== Visualization of iterations (heatmap) ==============================
# Setting color preferences: 
colors = ["#1f77b4", "#ff7f0e", "#d62728", "#2ca02c"]
cmap = ListedColormap(colors)
bounds = [0, 1, 2, 3, 4] 
norm = BoundaryNorm(bounds, cmap.N)

plt.figure(figsize=(12, 10))
ax = sns.heatmap(cluster_df, cmap=cmap, cbar = True, norm=norm, linewidths=0.1, linecolor='gray')

colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.5, 1.5, 2.5, 3.5])
colorbar.set_ticklabels(['1', '2', '3', '4'])
colorbar.set_label('Cluster ID')

plt.title("Cluster Assignment per Sample Across Seeds")
plt.xlabel("Random Seed")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig("part1_q2_iterations_heatmap.png", dpi=300)

print("Figure, question 2 saved!")

print("===== Question 3  =====")

# Printing the Within-Cluster Sum of Squares (WCSS) for all runs:
for run in results:
    print(f"Seed: {run['seed']}  |  WCSS: {run['wcss']:.2f}")

best_run = min(results, key=lambda x: x['wcss'])
print(f"\nBest run: Seed {best_run['seed']} with WCSS = {best_run['wcss']:.2f}")

# Saving the labels from the best run: 
best_labels = best_run['labels']

# =========== Function ===========
# Defining a function to detect the enriched genes per cluster 
def find_enriched_genes_per_cluster(data, labels, gene_names=None, alpha=0.05):
    """
    Perform Mann-Whitney U tests to find enriched/depleted genes per cluster.
    Apply Bonferroni correction for multiple testing.
    
    Parameters:
    -----------
        data: np.ndarray of shape (samples, genes), original unstandardized FPKM
        labels: np.ndarray of shape (samples,), cluster assignment for each sample
        alpha: significance threshold after correction

    Returns:
    --------
        results: dictionary with cluster ID → list of (gene_index, corrected_p)
    """
    n_genes = data.shape[1]
    print(n_genes)
    unique_clusters = np.unique(labels)
    print(unique_clusters)
    results = {}

    for cluster_id in unique_clusters:
        in_cluster = data[labels == cluster_id]
        out_cluster = data[labels != cluster_id]

        sig_genes = []

        for gene_idx in range(n_genes):
            gene_in = in_cluster[:, gene_idx]
            gene_out = out_cluster[:, gene_idx]

            # Mann-Whitney U test (two-sided)
            try:
                stat, p = mannwhitneyu(gene_in, gene_out, alternative='two-sided')
            except ValueError:
                # If the MannWhitneyU test fails we continue to the next gene 
                continue

            # Bonferroni correction
            n_tests = n_genes * len(unique_clusters) # computing the total tests that were performed 
            corrected_p = p * n_tests

            if corrected_p < alpha:
                    if gene_names is not None:
                        sig_genes.append((gene_names[gene_idx], corrected_p))
                    else:
                        sig_genes.append((gene_idx, corrected_p))

        # Sort by corrected p-value
        sig_genes.sort(key=lambda x: x[1])
        results[cluster_id] = sig_genes

    return results


# Implementing the function
# Ensuring that the raw data is transposed and in the format of a numpy array
gene_names = data.index.tolist()  # Keep gene names in order
significant_genes = find_enriched_genes_per_cluster(data.T.values, best_labels, gene_names=gene_names)

# Print example output
for cluster, genes in significant_genes.items():
    print(f"\nCluster {cluster} - {len(genes)} significant genes")
    for gene_name, pval in genes[:5]:  # print top 5
        print(f"  Gene {gene_name}, corrected p = {pval:.3e}")

# Saving significant genes per cluster to a .txt file
with open("significant_genes_per_cluster.txt", "w") as f:
    for cluster, genes in significant_genes.items():
        f.write(f"Cluster {cluster} - {len(genes)} significant genes\n")
        for gene_name, pval in genes:
            f.write(f"  Gene {gene_name}, corrected p = {pval:.3e}\n")
        f.write("\n") 

print("===== Question 4  =====")
print("This question is answered analytically in teh document")

print("""
Part 2: Constructing and implementing the k-means algorithm.
""")

print("===== Question 1: Variant Allele Frequency (VAF) for controls =====")

assignment_dir = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/ass_3"

# Loading the mutations CSV (variants of interest detected at diagnosis)
mutations_path = os.path.join(assignment_dir, 'input/mutations.csv')
mutations_df = pd.read_csv(mutations_path)

# Listing all control samples
control_dir = os.path.join(assignment_dir, "input", "controls")
control_files = [f for f in os.listdir(control_dir)]

# Creating an empty DataFrame for the result
vaf_table = pd.DataFrame(columns=['chr', 'pos', 'nucleotide'] + [f'{f.split(".")[0]}_VAF' for f in control_files])

# Iterating through each mutation (variant detected at diagnosis)
for index, mutation in mutations_df.iterrows():
    chr_ = mutation['chr']
    pos = mutation['pos']
    alt = mutation['alt']
    
    # Create a list of VAFs for each control sample (set to 0 if variant not found)
    vaf_row = [chr_, pos, alt]
    
    for control_file in control_files:
        # Loading the control file
        control_df = pd.read_csv(os.path.join(assignment_dir, control_dir, control_file))
        
        # Check if the mutation exists in the control file
        control_variant = control_df[(control_df['chr'] == chr_) & 
                                     (control_df['pos'] == pos) & 
                                     (control_df['nucleotide'] == alt)]
        
        if not control_variant.empty:
            # If variant exists, get the VAF (assuming it's the 'VAF' column)
            vaf_row.append(control_variant['VAF'].values[0])
        else:
            # If not found, append 0
            vaf_row.append(0)
    
    # Add the row to the result table
    vaf_table.loc[index] = vaf_row

# Saving table
vaf_table.to_csv('control_vaf_table.csv', index=False)

print("Control VAF table has been generated and saved to 'control_vaf_table.csv'.")


print("===== Question 2: Fitting exponential distributions =====")

# Creating a list of control VAF columns
control_columns = [col for col in vaf_table.columns if col.endswith('_VAF')]

# Initiating a list to store result rows
results = []

for _, row in vaf_table.iterrows():
    vaf_values = row[control_columns].astype(float).values # obtaining the VAF values 
    
    # Checking if at least one VAF is non-zero
    if np.any(vaf_values > 0):
        # Fitting exponential distribution (mu=0) where lamda = 1 / mean
        loc_fixed = 0
        _, scale = expon.fit(vaf_values[vaf_values > 0], floc=loc_fixed)  # fit only non-zero VAFs
        rate = 1 / scale
        
        # Creating result row
        variant_id = f"{row['chr']} {int(row['pos'])} {row['nucleotide']}"
        results.append([variant_id, rate] + list(vaf_values))
        
# saving as dataframe & then .csv
output_columns = ['variant', 'rate'] + [f'control{i+1}' for i in range(len(control_columns))]
output_df = pd.DataFrame(results, columns=output_columns)
output_df.to_csv("variant_rate_fitting.csv", index=False)

print("Output was saved to 'variant_rate_fitting.csv'")

print("===== Question 3: Comparing patient-variant pairs at diagnosis vs. after treatment =====")

# Making a table of VAFs for patients in the same way as we did for the controls 
patient_dir = os.path.join(assignment_dir, 'input/patients')
patient_files = [f for f in os.listdir(patient_dir) if f.endswith('.csv')]

patient_vaf_table = pd.DataFrame(columns=['chr', 'pos', 'nucleotide'] + [f'{f.split(".")[0]}_VAF' for f in patient_files])

for index, mutation in mutations_df.iterrows():
    chr_ = mutation['chr']
    pos = mutation['pos']
    alt = mutation['alt']
    
    vaf_row = [chr_, pos, alt]
    
    for patient_file in patient_files:
        patient_df = pd.read_csv(os.path.join(patient_dir, patient_file))
        
        patient_variant = patient_df[(patient_df['chr'] == chr_) & 
                                     (patient_df['pos'] == pos) & 
                                     (patient_df['nucleotide'] == alt)]
        
        if not patient_variant.empty:
            vaf_row.append(patient_variant['VAF'].values[0])
        else:
            vaf_row.append(0) 
    
    patient_vaf_table.loc[index] = vaf_row

patient_vaf_table.to_csv('patient_vaf_table.csv', index=False)

# We will also use the table with the exponential modelling we did on control data 
background_table = output_df 

# Inspecting table structure
print(background_table.head())
print(patient_vaf_table.head())

# Ensuring the information appears in the same format betweeen the two tables 
background_table['chr'], background_table['pos'], background_table['nucleotide'] = zip(*background_table['variant'].str.split(" "))
background_table['pos'] = background_table['pos'].astype(int)

# Joining the two datasets based on variant 
merged = pd.merge(background_table[['chr', 'pos', 'nucleotide', 'rate']], 
                  patient_vaf_table, 
                  on=['chr', 'pos', 'nucleotide'])

# Preparing a list of the patient columns
vaf_columns = [col for col in merged.columns if col.endswith('_VAF')]

results = []

for _, row in merged.iterrows():
    variant_id = f"{row['chr']} {row['pos']} {row['nucleotide']}"
    rate = row['rate']  # exponential rate for this variant
    
    for patient_col in vaf_columns:
        patient_id = patient_col.replace('_VAF', '')
        vaf = row[patient_col]
        
        # p-value under exponential null
        p_val = np.exp(-rate * vaf)
        results.append([variant_id, patient_id, vaf, p_val])

# Saving 
results_df = pd.DataFrame(results, columns=['variant', 'patient', 'VAF', 'p_value'])
results_df.to_csv('variant_patient_nominal_pvalues.csv', index=False)
print(results_df.head())


print("===== Question 4: Multiple testing correction =====")

# We will use the dataframe with the p-values from the previous question

print(results_df.head())
adjusted_results = multipletests(results_df['p_value'], method='fdr_bh')

results_df['FDR'] = adjusted_results[1]

# Defining the FDR thresholds required 
fdr_001 = results_df[results_df['FDR'] <= 0.01]
fdr_005 = results_df[results_df['FDR'] <= 0.05]

output_df = fdr_005.copy()
output_df.rename(columns={'patient': 'sample', 'p_value': 'p'}, inplace=True)

output_df = output_df[['variant', 'sample', 'p', 'FDR']]

output_df.to_csv('significant_variants_fdr05.csv', index=False)
print("The significant variants at FDR ≤ 0.05 have been saved to 'significant_variants_fdr05.csv'")

# Summary of results 
print(f"Significant at FDR ≤ 0.01: {len(fdr_001)}")
print(f"Significant at FDR ≤ 0.05: {len(fdr_005)}")

print("===== Question 5: Fisher's exact test =====")

# determining the dataframes that I will use for this analysis 
sig_variants = output_df

print(mutations_df.head())
print(sig_variants.head())

# Building "variant" column in mutations_df to match sig_variants
mutations_df['variant'] = mutations_df['chr'] + ' ' + mutations_df['pos'].astype(str) + ' ' + mutations_df['alt']

# Conveting df with all patient variants to long format
print(patient_vaf_table.head())  # assuming this is the wide-format VAF table

long_df = patient_vaf_table.melt(
    id_vars=['chr', 'pos', 'nucleotide'],
    var_name='sample',
    value_name='VAF'
)

# Clean up column values and create matching 'variant' column
long_df['sample'] = long_df['sample'].str.replace('_VAF', '', regex=False)
long_df['variant'] = long_df['chr'] + ' ' + long_df['pos'].astype(str) + ' ' + long_df['nucleotide']

# Join with variant info
long_df = long_df.merge(
    sig_variants[['variant', 'sample', 'FDR']],
    on=['variant', 'sample'],
    how='left'
)

# Adding 'significant' column
long_df['significant'] = long_df['FDR'] <= 0.05
long_df['significant'] = long_df['significant'].fillna(False)

# Adding 'at_diagnosis' flag 
diagnosis_pairs = set(zip(mutations_df['variant'], mutations_df['sample']))
long_df['at_diagnosis'] = long_df.apply(
    lambda row: (row['variant'], row['sample']) in diagnosis_pairs,
    axis=1
)

# Optional: check result
print("Now the full long-format variant dataframe looks like this:")
print(long_df.head())

# Contignency table
A = ((long_df['significant']) & (long_df['at_diagnosis'])).sum()
B = ((~long_df['significant']) & (long_df['at_diagnosis'])).sum()
C = ((long_df['significant']) & (~long_df['at_diagnosis'])).sum()
D = ((~long_df['significant']) & (~long_df['at_diagnosis'])).sum()

contingency = pd.DataFrame(
    [[A, B], [C, D]],
    index=['Observed at diagnosis', 'Not observed at diagnosis'],
    columns=['Significant', 'Not significant']
)

print("\nContingency Table:")
print(contingency)

# Running Fisher's exact test 
from scipy.stats import fisher_exact

oddsratio, p_value = fisher_exact([[A, B], [C, D]], alternative='greater')
print(f"\nFisher's exact test p-value (one-tailed, 'greater'): {p_value}")

print("===== Question 6: Fraction of diagnosed variants at different FDR thresholds =====")

# Loading or define nominal p-values DataFrame
nominal_values_path = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/Foundational-Computational-Biology-Course/assignment_3/variant_patient_nominal_pvalues.csv"
nominal_df = pd.read_csv(nominal_values_path)

# Filtering to diagnostic variants
diagnostic_variants = set(mutations_df['variant'] + '_' + mutations_df['sample'])  # or however you defined it
nominal_df['variant_sample'] = nominal_df['variant'] + '_' + nominal_df['patient']
diagnostic_df = nominal_df[nominal_df['variant_sample'].isin(diagnostic_variants)].copy()

# Benjamini-Hochberg FDR correction
diagnostic_df['FDR'] = multipletests(diagnostic_df['p_value'], method='fdr_bh')[1]

# Output
diagnostic_df.rename(columns={
    'patient': 'sample',
    'p_value': 'p',
    'VAF': 'rate'  # assuming 'rate' is your observed VAF here
}, inplace=True)

output_df = diagnostic_df[['variant', 'sample', 'rate', 'p', 'FDR']]
output_df.to_csv("significant_diagnostic_variants.csv", index=False)

# Summary
total = len(diagnostic_df)
n_fdr_01 = (diagnostic_df['FDR'] <= 0.01).sum()
n_fdr_05 = (diagnostic_df['FDR'] <= 0.05).sum()

print(f"Fraction no longer significant at FDR ≤ 0.01: {1 - n_fdr_01 / total:.2%}")
print(f"Fraction no longer significant at FDR ≤ 0.05: {1 - n_fdr_05 / total:.2%}")

print("===== Question 7: Fraction of diagnosed variants at different FDR thresholds =====")

print(mutations_df.head()) # diagnosed variants 
signif_variants_q6_path = "/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H/Foundational-Computational-Biology-Course/assignment_3/to_submit_Part2_q6_significant_diagnostic_variants.csv"
signif_variants_q6 = pd.read_csv(signif_variants_q6_path)

print(output_df.head()) # diagnosed variants which remained significant after treatment

# Step 1: Create a unique ID for each patient-variant pair
mutations_df['variant_sample'] = mutations_df['variant'] + '_' + mutations_df['sample']
signif_variants_q6['variant_sample'] = signif_variants_q6['variant'] + '_' + signif_variants_q6['sample']

print(mutations_df.head()) # diagnosed variants 
print(signif_variants_q6.head())

all_diagnosed = set(mutations_df['variant_sample'])

# Identifying all still-significant pairs after treatment
still_significant = set(signif_variants_q6['variant_sample'])

# pairs that lost significance
no_longer_significant = all_diagnosed - still_significant

# which diagnosed variants are no longer significant
mutations_df['responded'] = mutations_df['variant_sample'].isin(no_longer_significant)

# Group by patient and find those with at least one non-significant variant
responding_patients = (
    mutations_df[mutations_df['responded']]
    .groupby('sample')
    .size()
    .index
    .tolist()
)

# Compute the fraction of responding patients
total_patients = mutations_df['sample'].nunique()
fraction_responded = len(responding_patients) / total_patients

# Output
print("Patients who responded to treatment:")
print(responding_patients)
print(f"\nFraction of patients who responded: {fraction_responded:.2%}")

# Details 
print("Total diagnosed variant_sample entries:", len(mutations_df['variant_sample'].unique()))
print("Total significant variant_sample entries:", len(signif_variants_q6['variant_sample'].unique()))
print("Overlap:", len(set(mutations_df['variant_sample']) & set(signif_variants_q6['variant_sample'])))

diagnosed_set = set(mutations_df['variant_sample'])
significant_set = set(signif_variants_q6['variant_sample'])

# Check a few samples
print("Diagnosed examples:", list(diagnosed_set)[:5])
print("Significant after treatment examples:", list(significant_set)[:5])

no_longer_significant = diagnosed_set - significant_set
print("Variants no longer significant:", list(no_longer_significant)[:5])
