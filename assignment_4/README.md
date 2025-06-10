# FCBII Assignment 1 - RNA motif finding program

## Overview

This project is comprised of 1 program that runs 4 different motid-detection models: 
- **HardOOPS**: Assumes one motif per sequence and uses a greedy, deterministic approach to find the most probable motif occurrences.

- **HardZOOPS**: Allows zero or one motif per sequence, combining greedy selection with the possibility of skipping sequences without a clear motif.

- **GibbsZOOPS**: A probabilistic model that allows zero or one motif per sequence and uses Gibbs sampling to iteratively update motif positions and the motif profile.

- **AnnealZOOPS**: Uses simulated annealing to stochastically explore motif placements under a zero-or-one-per-sequence constraint, balancing exploration and convergence.

## Structure

- `assignment_1_pt2.py`: Main script for running the analysis
- `Assignment_1_FCBII.docx`: Document with the rationale about each model and the automation to iterate through the 5 input text files.
- `MMG1344H_env.yml`: Conda environment created to perform the computations for this assignment. Python 3.12.9 is used. 
along with basic libraries such as pandas and matplotlib.
- `test_input_{number}`: There are 5 directories that contain the CLI output for running all models, as well as the plot for log (E,o) per model. The structure of those folders is identical and looks like this: 
```bash
├── AnnealZOOPS
│   └── log_e_history.png
├── GibbsZOOPS
│   └── log_e_history.png
├── HardOOPS
│   └── log_e_history.png
├── HardZOOPS
│   └── log_e_history.png
└── output.txt
```

## How to Run

1. Unzip the .zip assignment file 

```bash
unzip mariaeleni_fafouti.zip
```

2. After adding the contents of the .zip file in a folder, the structure should look like this: 

```bash
├── Assignment_1_FCBII.docx
├── assignment_1_pt2.py
├── MMG1344H_env.yml
├── README.md
├── test_input_1
├── test_input_2
├── test_input_3
├── test_input_4
└── test_input_5

6 directories, 4 files
```

3. Create the conda environment 
```bash
conda env create -f MMG1344H_env.yml
```

4. Make sure that the input files are installed and unzipped. The input folder structure should look like this: 

```bash
├── sample_hardoops.png
├── sample_input.txt
├── sample_output.txt
├── test_input_1.txt
├── test_input_2.txt
├── test_input_3.txt
├── test_input_4.txt
└── test_input_5.txt

1 directory, 8 files
```
4. In the `assignment3.py` file, change the following paths to match your file structrue. 

```python
project_path = "/path/to/your/ass_1_FCB2"
data_path = "path/to/your/ass_1_FCB2_input"
output_base = "path/to/save/output" # can be the same as the project path, but make sure to type output_base = project_path
```

*Note:* You will need to define different paths on certain questions, based on your folder structure and teh way you save outputs. 

5. Most of the code file comprises of functions, while the program is ran at this part: 
```python
# ================ Running the program ================
if __name__ == "__main__":

```
