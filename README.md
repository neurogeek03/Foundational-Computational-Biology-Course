# FCBII Assignment 2 - Dimensionality Reduction & Machine Learning 

## Overview

This project is comprised of two parts: 

- **Part 1**: Dimensionality Reduction using PCA, NMF, t-SNE in a subset of sc-RNAseq data. 
- **Part 2**: Machine learning approaches usign the MNIST dataset. Methods used: Multinomial logistic regression, Support Vector Machine (SVM) with linear or radial basis function kernel, Fully connected Neural Network, Convolutional Neural Network (CNN).

## Structure

- `pt2_assignment_2.py`: Main script for running the analysis
- `Assignment_2_FCBII.docx`: Document with the rationale about each question, icluding performance metrics, results and figures. 
- `MMG1344H_env.yml`: Conda environment created to perform the computations for this assignment. Python 3.12.9 is used. 
along with basic libraries such as pandas and matplotlib. Machine learning libraries used: scikit-learn, Keras, Tensorflow - all are embedded within the environment.

## How to Run

1. Unzip the .zip assignment file 

```bash
unzip mariaeleni_fafouti.zip
```

2. After adding the contents of the .zip file in a folder, the structure should look like this: 

```bash
.
├── Assignment_2_FCBII.docx
├── MMG1344H_env.yml
├── pt2_assignment_2.py
└── README.md

4 files
```

3. Create the conda environment 
```bash
conda env create -f MMG1344H_env.yml
```

4. Make sure that the input files are installed and unzipped. The input folder structure should look like this: 

```bash
.
└── resources
    ├── cell_labels.txt
    ├── expression_matrix.txt
    ├── mnist_test.csv
    ├── mnist_train.csv
    ├── PCam_test_10K.csv
    ├── PCam_train_40K.csv
    └── sample_pca.png

1 directory, 7 files
```
5. In the `pt2_assignment_2.py` file, change the following paths to match your file structrue. 

```python
project_path = "/path/to/your/ass_1_FCB2/folder"
data_path = "path/to/your/ass_2_FCB2_input/resources"
output_base = "path/to/save/output" # can be the same as the project path, but make sure to type output_base = project_path
```
*Note*: The first step in Machine Learning is the one that takes the longest time to run.