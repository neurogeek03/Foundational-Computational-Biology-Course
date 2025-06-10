import os 
import scanpy as sc 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import tensorflow as tf
from sklearn.decomposition import NMF
from sklearn.manifold import TSNE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report 
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout, ReLU
from tensorflow.keras import layers, models, regularizers, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras.utils import to_categorical
import time

# np.random.seed(723)

# Essential paths 
project_path = '/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H'
data_path = os.path.join(project_path, 'pt2_ass2/resources')
output_base = os.path.join(project_path, 'Foundational-Computational-Biology-Course/pt2_assignment_2')

# Data paths
# part 1
exp_mtx_path = os.path.join(data_path, 'expression_matrix.txt')
cell_labels_path = os.path.join(data_path, 'cell_labels.txt')

# part 2 
train_path = os.path.join(data_path, 'mnist_train.csv')
test_path = os.path.join(data_path, 'mnist_test.csv')

# ================ PART 1 ================
# ================ 0. DATA PREPARATION ================

# Loading data
expr_df = pd.read_csv(exp_mtx_path, sep="\t", index_col=0)  # rows = genes, columns = cells
labels_df = pd.read_csv(cell_labels_path, sep="\t")
labels_df = labels_df.set_index("NAME")  # Make cell names the index

# Data Preview 
print(expr_df.head())
print(labels_df.head())

# Ensuring columns in expression data match order in labels
expr_df = expr_df.loc[:, labels_df.index] # now, each column in the expression matrix corresponds to 1 row in the cell labels

# Transposing expression matrix to match AnnData format: observations = cells, variables = genes
adata = sc.AnnData(expr_df.T)

# Combining the 2 datasets: Add cell type labels as metadata
adata.obs['CellType'] = labels_df.loc[adata.obs_names, 'CellType'].values

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Creating a color dictionary, such that each cell type is consistently the same color: 
color_dict_hex = {
    'DC1': '#1f77b4',
    'DC5': '#ff7f0e',
    'DC2': '#2ca02c',
    'Mono1': '#d62728',
    'Mono4': '#9467bd',
    'DC3': '#8c564b',
    'DC4': '#e377c2',
    'Mono2': '#7f7f7f',
    'Mono3': '#bcbd22',
    'DC6': '#17becf'
}

# ================ 1. PRINCIPAL COMPONENT ANALYSIS ================
# ====== PCA PLOT
sc.tl.pca(adata, zero_center=False, svd_solver='arpack')  

# Extract variance explained
pc1_var = adata.uns['pca']['variance_ratio'][0] * 100
pc2_var = adata.uns['pca']['variance_ratio'][1] * 100

# Passing the variance explained to labels for both axes 
xlabel = f"PC1 ({pc1_var:.1f}%)"
ylabel = f"PC2 ({pc2_var:.1f}%)"

# Get PCA coordinates
pca_df = pd.DataFrame(adata.obsm['X_pca'][:, :2], columns=['PC1', 'PC2'], index=adata.obs_names)
pca_df['CellType'] = adata.obs['CellType']
pca_df['color'] = pca_df['CellType'].map(color_dict)

# Plot & save
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='CellType', s=30)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title('PCA: Human Blood scRNA-seq')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'pca_plot_all_genes.png'), dpi=300)
plt.close() 

# ====== SCREE PLOT
# Extracting variance explained from Scanpy's PCA result
pc_variance = adata.uns['pca']['variance_ratio']  
pc_variance_percent = pc_variance * 100  # Convert to %

# first 10 PCs
top10_var = pc_variance_percent[:10]

# Plot scree plot
plt.figure(figsize=(7, 4))
plt.bar(range(1, 11), top10_var, color='skyblue', edgecolor='black')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained (%)')
plt.title('Scree Plot of First 10 Principal Components')
plt.xticks(range(1, 11))
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'scree_plot_all_genes.png'), dpi=300)
plt.close() 

# How many PCs explain more than 1% of the variance? 
num_above_1 = np.sum(pc_variance_percent > 1)
print(f"{num_above_1} PCs explain more than 1% of the variance.")

# How much variance do the first 10 PCs explain?
total_top10 = np.sum(top10_var)
print(f"The first 10 PCs explain {total_top10:.2f}% of the total variance.")

# ====== REPEATING WITH HIGHLY VARIABLE GENES ONLY 
# Making a copy of the object 
adata_hvg = adata.copy()

# get highly variable genes
sc.pp.highly_variable_genes(adata_hvg, n_top_genes=500, flavor='seurat')

# subseting the genes (variable) to keep only the 500 most variable ones
adata_hvg = adata_hvg[:, adata_hvg.var['highly_variable']]

# ====== PCA PLOT (HVGs)
sc.tl.pca(adata_hvg, zero_center=False, svd_solver='arpack')  

pc1_var = adata_hvg.uns['pca']['variance_ratio'][0] * 100
pc2_var = adata_hvg.uns['pca']['variance_ratio'][1] * 100

xlabel = f"PC1 ({pc1_var:.1f}%)"
ylabel = f"PC2 ({pc2_var:.1f}%)"

pca_df = pd.DataFrame(adata_hvg.obsm['X_pca'][:, :2], columns=['PC1', 'PC2'], index=adata_hvg.obs_names)
pca_df['CellType'] = adata_hvg.obs['CellType']

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='CellType', s=30)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title('HVGs - PCA: Human Blood scRNA-seq')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'pca_plot_HVGs.png'), dpi=300)
plt.close() 


# ====== SCREE PLOT (HVGs)
pc_variance = adata_hvg.uns['pca']['variance_ratio']  
pc_variance_percent = pc_variance * 100  # Convert to %

top10_var = pc_variance_percent[:10]

plt.figure(figsize=(7, 4))
plt.bar(range(1, 11), top10_var, color='skyblue', edgecolor='black')
plt.xlabel('Principal Component')
plt.ylabel('Variance Explained (%)')
plt.title('HVGs - Scree Plot of First 10 Principal Components')
plt.xticks(range(1, 11))
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'scree_plot_HVGs.png'), dpi=300)
plt.close() 

# How many PCs explain more than 1% of the variance? 
num_above_1 = np.sum(pc_variance_percent > 1)
print(f"{num_above_1} PCs explain more than 1% of the variance.")

# How much variance do the first 10 PCs explain?
total_top10 = np.sum(top10_var)
print(f"The first 10 PCs explain {total_top10:.2f}% of the total variance.")

# ================ 2. NON-NEGATIVE MATRIX FACTORIZATION (NMF) ================
# Converting the expression matrix into a dense numpy array
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

# ====== Running NMF
# Initializing the model 
nmf = NMF(n_components=2, init='nndsvda', random_state=42, max_iter=1000)

# Fitting the model and transforming the data
X_nmf = nmf.fit_transform(X)

# Converting the nmf result into a dataframe, using the cell names as the row index.
nmf_df = pd.DataFrame(X_nmf, columns=['NMF1', 'NMF2'], index=adata.obs_names)
nmf_df['CellType'] = adata.obs['CellType'].values

# Plot & save
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=nmf_df,
    x='NMF1', y='NMF2',
    hue='CellType', 
    palette='tab10',  
    s=30, alpha=0.8, linewidth=0
)
plt.title('NMF: Human Blood scRNA-seq')
plt.xlabel('NMF1')
plt.ylabel('NMF2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Cell Type')
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'NMF_plot.png'), dpi=300)
plt.close() 

# # ====== DIFFERENT RANDOM INITIALIZATIONS
# list to store variance explained 
variance_explained_list = []

# 3 different random seeds
seeds = [0, 1, 2]

for seed in seeds:
    nmf = NMF(n_components=10, init='random', random_state=seed, max_iter=2000)
    W = nmf.fit_transform(X)
    H = nmf.components_
    X_hat = np.dot(W, H)
    
    # Computing how much variance is explained
    ss_total = np.sum(X ** 2)
    ss_residual = np.sum((X - X_hat) ** 2)
    variance_explained = 1 - (ss_residual / ss_total)

    
    variance_explained_list.append(variance_explained)
print(variance_explained_list)

# bar plot
plt.figure(figsize=(6, 4))
plt.bar(['Run 1', 'Run 2', 'Run 3'], variance_explained_list, color='skyblue')
plt.ylabel('Variance Explained')
plt.title('Variance Explained by NMF (3 initializations)')
for i, v in enumerate(variance_explained_list):
    plt.text(i, v + 0.001, f"{v:.6f}", ha='center', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'VAR_RANDOM_initializations_NMF.png'), dpi=300)
plt.close() 

# ====== HEATMAP OF THE H MATRIX

# running NMF 
nmf = NMF(n_components=10, init='nndsvda', random_state=42, max_iter=1000)
W = nmf.fit_transform(X)
H = nmf.components_  

# creating dataframe 
H_df = pd.DataFrame(H, 
                    index=[f'NMF{i+1}' for i in range(H.shape[0])],
                    columns=adata.var_names)


# Selecting top 10 genes per component
top_genes = set()
for i in range(H.shape[0]):
    top_indices = H[i].argsort()[-10:]  # top 10 genes
    top_genes.update(adata.var_names[top_indices])

# Subset the DataFrame
H_df_subset = H_df.loc[:, list(top_genes)]

plt.figure(figsize=(12, 6))
sns.heatmap(H_df_subset, cmap='viridis', xticklabels=True)
plt.title('NMF H Matrix: Top Genes per Component')
plt.xlabel('Genes')
plt.ylabel('NMF Components')
plt.tight_layout()
plt.savefig(os.path.join(output_base, 'NMF_Hmatrix_heatmap.png'), dpi=300)
plt.close() 

# ================ 3. t-DISTRIBUTED STOCHASTIC NEIGHBOR EMBEDDING (t-SNE) ================
# ====== Running t-SNE
# expression matrix (dense array)
X_full = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

# Run t-SNE
tsne = TSNE(n_components=2, perplexity=4, random_state=42)
X_tsne = tsne.fit_transform(X_full)

# saving tsne embeddings to anndata object
adata.obsm['X_tsne'] = X_tsne

# creating a dataframe for plotting 
tsne_df = pd.DataFrame(X_tsne, columns=['TSNE1', 'TSNE2'], index=adata.obs_names)
tsne_df['CellType'] = adata.obs['CellType']

# print(tsne_df.head())

plt.figure(figsize=(8,6))
sns.scatterplot(data=tsne_df, x='TSNE1', y='TSNE2', hue='CellType', palette = color_dict_hex, s=30)
plt.title('t-SNE plot (perplexity=4)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(output_base, 't-SNE_plot.png'), dpi=300)
plt.close() 

# ====== Trying different perplexity values 
# perplexities = [2, 12]
# tsne_results = {}

X_full = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
cell_types = adata.obs['CellType']

for p in perplexities:
    tsne = TSNE(n_components=2, perplexity=p, random_state=42)
    X_tsne = tsne.fit_transform(X_full)
    
    tsne_df = pd.DataFrame(X_tsne, columns=['TSNE1', 'TSNE2'], index=adata.obs_names)
    tsne_df['CellType'] = cell_types.values
    tsne_results[p] = tsne_df

# plot and save
for p, df in tsne_results.items():
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=df, x='TSNE1', y='TSNE2', hue='CellType', palette=color_dict_hex, s=30)
    plt.title(f't-SNE with Perplexity = {p}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_base, f"tsne_perplexity_{p}.png"), dpi=300)
    plt.close()

# ====== Trying different random initializations
seeds = [42, 7, 123]
tsne_random_results = {}
perplexity = 12

for seed in seeds:
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=seed)
    X_tsne = tsne.fit_transform(X_full)

    tsne_df = pd.DataFrame(X_tsne, columns=['TSNE1', 'TSNE2'], index=adata.obs_names)
    tsne_df['CellType'] = cell_types.values
    tsne_random_results[seed] = tsne_df

for seed, df in tsne_random_results.items():
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=df, x='TSNE1', y='TSNE2', hue='CellType', palette=color_dict_hex, s=30)
    plt.title(f't-SNE with Perplexity = {perplexity}, Random Seed = {seed}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'tsne_seed_{seed}.png', dpi=300)
    plt.close()

#sanity check
seed1, seed2 = seeds[0], seeds[1]
df1 = tsne_random_results[seed1]
df2 = tsne_random_results[seed2]

merged = pd.DataFrame({
    'TSNE1_seed1': df1['TSNE1'],
    'TSNE2_seed1': df1['TSNE2'],
    'TSNE1_seed2': df2['TSNE1'],
    'TSNE2_seed2': df2['TSNE2']
}, index=adata.obs_names)

plt.figure(figsize=(6, 6))
for i in range(len(merged)):
    plt.plot([merged.iloc[i, 0], merged.iloc[i, 2]],
             [merged.iloc[i, 1], merged.iloc[i, 3]],
             color='gray', linewidth=0.3)

plt.scatter(merged['TSNE1_seed1'], merged['TSNE2_seed1'], color='blue', label=f'Seed {seed1}', alpha=0.5)
plt.scatter(merged['TSNE1_seed2'], merged['TSNE2_seed2'], color='red', label=f'Seed {seed2}', alpha=0.5, s=10)
plt.legend()
plt.title(f"t-SNE Comparison: Seed {seed1} vs Seed {seed2}")
plt.tight_layout()
plt.savefig("tsne_sanity_comparison.png", dpi=300)

# ================ PART 2 ================
# Load data
train_df = pd.read_csv(train_path)
test_df = pd.read_csv(test_path)

# ===== Q1
# selecting the training images (x) and their labels (y)
X_train = train_df.iloc[:, 1:].values 
y_train = train_df.iloc[:, 0].values   

# selecting the test images (x) and their labels (y)
X_test = test_df.iloc[:, 1:].values
y_test = test_df.iloc[:, 0].values

# Defining regularization types
penalties = {
    "No Regularization": None,
    "L1 Regularization": "l1",
    "L2 Regularization": "l2"
}

results = {}

for name, penalty in penalties.items():
    print(f"\nTraining with {name}...")

    if penalty is None:
        model = LogisticRegression(penalty=None, solver='saga', multi_class='multinomial', max_iter=1000)
    else:
        model = LogisticRegression(penalty=penalty, solver='saga', multi_class='multinomial', max_iter=1000)

    start_time = time.time()
    model.fit(X_train, y_train)
    train_time = time.time() - start_time

    y_pred = model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)

    results[name] = {
        "accuracy": acc,
        "confusion_matrix": cm,
        "training_time": train_time
    }

# output results
for name, metrics in results.items():
    print(f"\n=== {name} ===")
    print(f"Accuracy: {metrics['accuracy']:.4f}")
    print(f"Training time: {metrics['training_time']:.2f} seconds")
    print("Confusion matrix:")
    print(metrics["confusion_matrix"])

    # Plot heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(metrics["confusion_matrix"], 
                annot=True, 
                fmt='d', 
                cmap='Blues', 
                xticklabels=range(10), 
                yticklabels=range(10))
    plt.xlabel("Predicted Label")
    plt.ylabel("True Label")
    plt.title(f"Confusion Matrix Heatmap: {name}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_base, f'CM_{name}.png'), dpi=300)


# ===== Q2
# Linear function kernel 
print("Training SVM with Linear Kernel...")
start_time = time.time()
linear_svm = SVC(kernel='linear')  # creating model 
linear_svm.fit(X_train, y_train)   # training 
end_time = time.time()
print(f"Training time (Linear Kernel): {end_time - start_time:.4f} seconds\n")

# Predict labels and evaluate
y_pred_linear = linear_svm.predict(X_test)
print("Linear Kernel Performance:\n")
print(classification_report(y_test, y_pred_linear))

# Radial basis function kernel (RBF)
print("Training SVM with RBF Kernel...")
start_time = time.time()
rbf_svm = SVC(kernel='rbf')
rbf_svm.fit(X_train, y_train)
end_time = time.time()
print(f"Training time (RBF Kernel): {end_time - start_time:.4f} seconds\n")

# Predict labels and evaluate
y_pred_rbf = rbf_svm.predict(X_test)
print("RBF Kernel Performance:\n")
print(classification_report(y_test, y_pred_rbf))


# ===== Q3
tf.random.set_seed(42)
tf.keras.backend.clear_session()

# Normalizing pixel values
X_train = X_train.astype('float32') / 255.0
X_test = X_test.astype('float32') / 255.0

# One-hot encoding of the labels
y_train = to_categorical(y_train, 10)
y_test = to_categorical(y_test, 10)

# Defining hyperparameters
learning_rates = [0.001, 0.01, 0.1]
regularizations = {
    'none': None,
    'l1': regularizers.l1(0.001),
    'l2': regularizers.l2(0.001)
}

results = []

# Looping over all combinations
for reg_name, reg in regularizations.items():
    for lr in learning_rates:
        print(f"\nTraining with {reg_name.upper()} regularization and learning rate = {lr}")

        # define model
        model = models.Sequential([
            Input(shape=(784,)),  
            layers.Dense(800, activation='relu', kernel_regularizer=reg),
            layers.Dense(10, activation='softmax')
        ])

        optimizer = tf.keras.optimizers.Adam(learning_rate=lr)

        model.compile(optimizer=optimizer,
                      loss='categorical_crossentropy',
                      metrics=['accuracy'])

        # training and timing 
        start_time = time.time()
        history = model.fit(X_train, y_train,
                            epochs=10,
                            batch_size=128,
                            validation_split=0.2,
                            verbose=0)
        elapsed_time = time.time() - start_time

        # Evaluation
        test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)

        # recording the result
        results.append({
            'reg': reg_name,
            'lr': lr,
            'accuracy': test_acc,
            'train_time': elapsed_time
        })

# displaying results
results = sorted(results, key=lambda x: x['accuracy'], reverse=True)

print("\nTop results (sorted by test accuracy):")
for res in results:
    print(f"Reg: {res['reg']:>4}, LR: {res['lr']}, Accuracy: {res['accuracy']:.4f}, Time: {res['train_time']:.2f}s")


# ===== Q4
# Normalizing pixel values
X_train = X_train.astype('float32') / 255.0
X_test = X_test.astype('float32') / 255.0

# Add channel dimension for grayscale images
X_train = X_train.reshape(-1, 28, 28, 1)
X_test = X_test.reshape(-1, 28, 28, 1)

# One-hot encoding of the labels
y_train = to_categorical(y_train, 10)
y_test = to_categorical(y_test, 10)

model = Sequential([
    # 32 3x3 filters will "scan" the image of 28x28 pixels, generating 32 new images, one per filter
    Conv2D(32, (3, 3), padding='same', input_shape=(28, 28, 1)),
    ReLU(),
    MaxPooling2D(pool_size=(2, 2)),

    Conv2D(64, (3, 3), padding='same'),
    ReLU(),
    MaxPooling2D(pool_size=(2, 2)),

    Flatten(),
    Dense(128, kernel_regularizer=l2(0.0005)),
    ReLU(),
    Dropout(0.5),
    Dense(10, activation='softmax')  # Softmax here because we'll use categorical_crossentropy
])

# compiling model
model.compile(
    optimizer=Adam(learning_rate=0.001),
    loss='categorical_crossentropy',
    metrics=['accuracy']
)

# training model
start = time.time()
history = model.fit(X_train, y_train, epochs=10, batch_size=64, validation_split=0.1)
end = time.time()
print(f"Training time: {end - start:.2f} seconds")

# Predict class probabilities
y_pred_probs = model.predict(X_test)
y_pred = np.argmax(y_pred_probs, axis=1)
y_true = np.argmax(y_test, axis=1)

# Accuracy and reports
print("Classification Report:\n")
print(classification_report(y_true, y_pred))

print("Confusion Matrix:\n")
cm=confusion_matrix(y_true, y_pred)

plt.figure(figsize=(8, 6))
sns.heatmap(cm, 
                annot=True, 
                fmt='d', 
                cmap='Blues', 
                xticklabels=range(10), 
                yticklabels=range(10))
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title(f"Confusion Matrix Heatmap: My CNN")
plt.tight_layout()
plt.savefig(os.path.join(output_base, f'CM_myCNN.png'), dpi=300)