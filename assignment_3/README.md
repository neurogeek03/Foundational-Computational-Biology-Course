# FCBI Assignment 3 - k-means Clustering and Variant Allele Frequency in Leukemia

## Overview

This project is comprised of 2 parts: 
- **Part 1**: Defining the k-means algorithm as a function and implementing it to cluster 40 different blamstomeres originating from either a two-cell or a 4-cell embryo. 
- **Part 2**: Comparing the Variant Allele Frequencies of cells in leukemia patients before and after treatment, as well as with control data. 

## Structure

- `assignment3.py`: Main script for running the analysis
- `MMG1344H_env.yml`: Conda environment created to perform the computations for this assignment. Python 3.12.9 is used. 
along with basic libraries such as pandas and matplotlib.
- `part1_q2_iterations_heatmap.png`: `heatmap visualization from Part 1 
- `part1_q3_significant_genes_per_cluster.txt`: The program's full output to Part 1, question 3.
- `Part2_q2_variant_rate_fitting.csv`: The program's full output to Part 2, question 2
- `Part2_q4_significant_variants_fdr05.csv`:The program's full output to Part 2, question 4
- `Part2_q6_significant_diagnostic_variants.csv`: The program's full output to Part 2, question 6
- `Assignment 3.docx`: Document with the rationale and answers to all questions of both parts.

## How to Run

1. Unzip the .zip assignment file 

```bash
unzip mariaeleni_fafouti.zip
```

2. After adding the contents of the .zip file in a folder, the structure should look like this: 

```bash
├── Assignment 3.docx
├── MMG1344H_env.yml
├── part1_q2_iterations_heatmap.png
├── part1_q3_significant_genes_per_cluster.txt
├── Part2_q2_variant_rate_fitting.csv
├── Part2_q4_significant_variants_fdr05.csv
├── Part2_q6_significant_diagnostic_variants.csv
└── README.md

1 directory, 8 files
```

3. Create the conda environment 
```bash
conda env create -f MMG1344H_env.yml
```
4. In the `assignment3.py` file, change the following paths at the first section, to match your file structrue. Make sure to have the input data installed and name it `Part1_Biase_2014.csv`.
```python
project_path = "/path/to/your/ass_3"
```

*Note:* You will need to define different paths on certain questions, based on your folder structure and teh way you save outputs. 

5. For Part 2, obtain the input folder, which should look like this:

```bash
├── controls
│   ├── control1.csv
│   ├── control10.csv
│   ├── control11.csv
│   ├── control12.csv
│   ├── control13.csv
│   ├── control14.csv
│   ├── control15.csv
│   ├── control16.csv
│   ├── control2.csv
│   ├── control3.csv
│   ├── control4.csv
│   ├── control5.csv
│   ├── control6.csv
│   ├── control7.csv
│   ├── control8.csv
│   └── control9.csv
├── mutations.csv
└── patients
    ├── Patient_0107.csv
    ├── Patient_0186.csv
    ├── Patient_0247.csv
    ├── Patient_0277.csv
    ├── Patient_0298.csv
    ├── Patient_0488.csv
    ├── Patient_0531.csv
    ├── Patient_0573.csv
    ├── Patient_0574.csv
    ├── Patient_0672.csv
    ├── Patient_0872.csv
    ├── Patient_0967.csv
    ├── Patient_1230.csv
    ├── Patient_1441.csv
    ├── Patient_1589.csv
    ├── Patient_1618.csv
    ├── Patient_1815.csv
    ├── Patient_2056.csv
    ├── Patient_2373.csv
    ├── Patient_2478.csv
    ├── Patient_2570.csv
    ├── Patient_2590.csv
    ├── Patient_2662.csv
    └── Patient_2819.csv

3 directories, 41 files
```

This folder will be used for any downstream analyses 

