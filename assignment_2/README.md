# FCBI Assignment 2 - Sequence Alignment & SNP Population Genetics

## Overview

This project is comprised of 2 parts: 
- **Part 1**: Modifying the Smith-Waterman algorithm for RNA-RNAalignment, at opposite directions with custom nucleotide pair scoring.
- **Part 2**: Interpreting SNPs from the HapMap3 dataset and comparing the FST scores with the findings of Li. et al. 2012

## Structure

- `assignment1.py`: Main script for running the analysis
- `MMG1344H_env.yml`: Conda environment created to perform the computations for this assignment. Python 3.12.9 is used, 
along with basic libraries such as pandas and matplotlib.
- `output.txt`: The program's output to Part 1 for the full test set
- `Assignment_2.docx`: Document with answers to Part 2

## How to Run

1. Unzip the .zip assignment file 

```bash
unzip mariaeleni_fafouti.zip
```

2. After adding the contents of the .zip file in a folder, the structure should look like this: 

```bash
├── Assignment_2.docx
├── MMG1344H_env.yml
├── README.md
├── assignment2.py
└── output.txt

1 directory, 5 files
```

3. Create the conda environment 
```bash
conda env create -f MMG1344H_env.yml
```
4. In the `assignment2.py` file, change the following paths at the first section, to match your file structrue 
```python
fasta_path = "/path/to/your/orf_trans.fasta"
assignment_pdf_path = "/path/to/your/FCBI_Assignment_1_Instructions_2025-1.pdf"
```
5. For Part2, download the HapMap3 data (HapMap3.zip) from: https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3 


6. The program unzips both the HapMap file and the individual hmap files for each population. Once the HapMap3.zip file is decompresed the folder structure should look like this: 
```bash
hapmap3_data
├── HapMap3
│   ├── ASW.hmap.gz
│   ├── CEU.hmap.gz
│   ├── CHB.hmap.gz
│   ├── CHD.hmap.gz
│   ├── GIH.hmap.gz
│   ├── JPT.hmap.gz
│   ├── LWK.hmap.gz
│   ├── MEX.hmap.gz
│   ├── MKK.hmap.gz
│   ├── TSI.hmap.gz
│   └── YRI.hmap.gz
└── __MACOSX
    └── HapMap3

4 directories, 11 files
```

7. The decompressed individual hmap files will be saved in a new folder (`hapmap3_unzipped`) in the same directory. The folder structure should look like this: 
```bash
hapmap3_unzipped
├── ASW.hmap
├── CEU.hmap
├── CHB.hmap
├── CHD.hmap
├── GIH.hmap
├── JPT.hmap
├── LWK.hmap
├── MEX.hmap
├── MKK.hmap
├── TSI.hmap
└── YRI.hmap

1 directory, 11 files
```
This folder will be used for any downstream analyses 

