# FCBI Assignment 1 - Amino Acid Frequency Analysis

## Overview

This project is comprised of 2 parts: 
- **Part 1**: Probability Analysis, Pattern Matching
- **Part 2**: Yeast Proteome Analysis, Comparison of amino-acid frequencies 

## Structure

- `assignment1.py`: Main script for running the analysis
- `plots/`: Output folder for generated figures
- `Assignment 1.docx` or `answers.txt`: Contains the rationale for solving each question, and the calculations used
- `MMG1344H_env.yml`: Conda environment created to perform the computations for this assignment. Python 3.12.9 is used, 
along with basic libraries such as pandas and matplotlib. 

## How to Run

1. Unzip the .zip assignment file 

```bash
pip install matplotlib biopython
```

2. After adding the contents of the .zip file in a folder, the structure should look like this: 

```bash
.
├── Assignment 1.docx
├── FCBI_Assignment_1_Instructions_2025-1.pdf
├── MMG1344H_env.yml
├── README_assignment1.md
├── assignment1.py
├── orf_trans.fasta
├── output.txt
└── plots_a1
    ├── difference_matrix.png
    ├── section4_5_side_by_side.png
    └── section4_dipeptide_freq_matrix.png

2 directories, 10 files
```

3. Create the conda environment 
```bash
conda env create -f MMG1344H_env.yml
```

4. In the `assignment1.py` file, change the following paths at the first section, to match your file structrue 
```python
fasta_path = "/path/to/your/orf_trans.fasta"
assignment_pdf_path = "/path/to/your/FCBI_Assignment_1_Instructions_2025-1.pdf"
```

## Note
The `assignment1.py` script is designed to have its output saved in an output.txt file. You can move the output to the terminal by commenting out lines 14-15 and 421: 
```python
# with open('output.txt', 'w') as file: # line 14
#     sys.stdout = file                 # line 15

# sys.stdout = sys.__stdout__           # line 421
```
Then select the script after the imports to remove the indentation
