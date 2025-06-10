import numpy as np
import matplotlib.pyplot as plt
import math
import random
import sys
import os 
import time
import glob
from contextlib import redirect_stdout

# Defining paths 
project_path = '/Users/marlenfaf/Desktop/UofT_PhD/MMG1344H'
data_path = os.path.join(project_path, 'pt2_ass1/FCB2_Asn1_resources')
output_base = os.path.join(project_path, 'Foundational-Computational-Biology-Course/pt2_assignment_1')

input_file = os.path.join(data_path, 'sample_input.txt')



# ================ General Setup ================

NUCLEOTIDES = ['A', 'C', 'G', 'U']

def read_sequences(filepath):
    """Reading input file and gathering all sequences in a list

    Args:
        filepath (string): path to input_data.txt

    Returns:
        list: sequences read from the input file
    """
    sequences = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip().upper()
            if line:
                sequences.append(line)
    return sequences

def initialize_random_offsets(sequences, K):
    """Random initialization of offsets using 1-indexing

    Args:
        sequences (list): list of sequences from input file
        K (int): user-determined length of motif (6 in this assignment)

    Returns:
        list: offsets for the sequences given (certain models give 1 per sequence, 
        others allow certain sequences to not have any offsets)
    """
    offsets = []
    for seq in sequences: 
        max_offset = len(seq) - K + 1
        if max_offset < 1:
            raise ValueError(f"Sequence too short for K={K}: {seq}")
        offset = random.randint(1, max_offset)
        offsets.append(offset)
    return offsets

def extract_kmer(seq, offset, K):
    """Extracting k-mer using 1-indexed offset. This function adapts the 1-indexing 
    used for the offset to python's default 0-indexed string slicing, to extract 
    the k-mer correctly. 

    Args:
        seq (string): 1 sequence out of the list of sequences from input file
        offset (int): 1 offset out of the list of offsets for the sequences given
        K (int): user-determined length of motif (6 in this assignment)

    Returns:
        string: a possible k-mer sequence, composed of different combinations of the 
        4 nucleotides 
    """
    return seq[offset - 1: offset - 1 + K]

def sample_from_distribution(probabilities):
    """Sample an index from a discrete probability distribution, represented as a list
    of probabilities. It picks a random number along the cumulative probability scale. 

    Args:
        probabilities (list): a list of floats, each representing the probability of 
        selecting the corresponding index. 

    Returns:
        list: Relative probability of selecting a specific starting position for a 
            motif in a sequence. Each list corresponds to 1 sequence. 
    """
    r = random.random()
    cumulative = 0.0
    for i,p in enmumerate(probabilities): 
        cumulative += p 
        if r <= cumulative: 
            return i 
    return len(probabilities) - 1 

def build_pfm(sequences, offsets, K, pseudocount=0.25): 
    """Construct a position frequency matrix (PFM) from a set  of sequences and their 
    corresponding motif positions (offsets). The PFM provides information about how 
    often each nucleotide appears at each position of the motif. 

    Args:
        sequences (list): list of sequences from input file
        K (int): user-determined length of motif (6 in this assignment)
        offsets (list): list of integers that serve as the starting points for the k-mer in 
            the sequences.
        pseudocount (float, optional): We initialize the probability of each nucleotide to be 
            present at a given position. At the start, all nucleotides have equal probability (1/4). 

    Returns:
        dictionary: Each key is a nucleotide and the value is a list of counts, one per position
            in the k-mer (motif). Those counts originate from all input sequences. 
    """
    pfm = {nuc: [pseudocount] * K for nuc in NUCLEOTIDES}

    # Iterating over each sequence & extract the motif 
    for seq, offset in zip(sequences, offsets):
        kmer = extract_kmer (seq, offset, K)
        # Keeping count of how many times a base appears and saving it in the PFM
        for i, base in enumerate(kmer): 
            pfm[base][i] += 1
    return pfm

def score_kmer(kmer, pfm, K): 
    """Determine the score of a k-mer using the latest PFM. After obtaining the frequency of a base
    at a specific position in the k-mer, it divides that over the total counts of bases at this same 
    position. Using log conversion turns multiplication into addition which allows us to avoid 
    generating too small probabilities. This facilitates the comparison of probability values 
    across the four models. 

    Args:
        kmer (string): a sequence containing different combinations of nucleotides 
        pfm (dictionary): Each key is a nucleotide and the value is a list of counts, one per position
            in the k-mer (motif). Those counts originate from all input sequences. 
        K (int): user-determined length of motif (6 in this assignment)

    Returns:
        float: log-probability of obsering this k-mer under the model defined by the PFM
    """
    background_prob = 0.25
    score = 0.0
    for i, base in enumerate(kmer):
        col_total = sum(pfm[nuc][i] for nuc in NUCLEOTIDES)
        pfm_prob = pfm[base][i] / col_total
        score += math.log(pfm_prob / background_prob)
    return score

def compute_log_e(sequences, offsets, pfm, K):
    """Computing the log-likelihood of the current motif model M and current set of motif positions o: 
    evaluating how well the current PFM explains the k-mers currently selected from each sequence. It
    extracts and scores all k-mers to see if the model is improving. 

    Args:
        sequences (list): list of sequences from input file
        offsets (list): list of integers that serve as the starting points for the k-mer in 
            the sequences.
        pfm (dictionary): Each key is a nucleotide and the value is a list of counts, one per position
            in the k-mer (motif). Those counts originate from all input sequences. 
        K (int): user-determined length of motif (6 in this assignment)

    Returns:
        float: Represents how well the current PFM fits the selected motifs across all sequences. 
    """
    total_log_prob = 0.0 
    for seq, offset in zip(sequences, offsets): 
        kmer = extract_kmer(seq, offset, K)
        total_log_prob += score_kmer(kmer, pfm, K)
    return total_log_prob

def print_pfm(pfm):
    """Print the position frequency matrix (PFM) as a readable table, showing the probability of each
    nucleotide at each position in the motif. Each row represents the probabilities for 1 nucleotide 
    and each column represents one position in the k-mer.

    Args:
        pfm (dictionary): Each key is a nucleotide and the value is a list of counts, one per position
            in the k-mer (motif). Those counts originate from all input sequences. 
    """
    print("PFM:")
    K = len(next(iter(pfm.values())))
    totals = [sum(pfm[nuc][i] for nuc in NUCLEOTIDES) for i in range(K)]

    # Print header with position numbers
    header = ["Pos"] + [f"{i+1:>6}" for i in range(K)]
    print("".join(header))

    # Print each row with fixed width formatting
    for nuc in NUCLEOTIDES:
        row = [f"{nuc:<4}"] + [f"{pfm[nuc][i] / totals[i]:>6.3f}" for i in range(K)]
        print("".join(row))

def plot_log_e_history(log_e_history, filename="log_e_plot.png"):
    """Plot the log E(M, o) values over iterations and saves the plot to a file.

    Parameters:
        log_e_history (list of float): List of log E(M, o) values for each iteration.
        filename (str): The name of the file to save the plot.
    """
    plt.plot(range(1, len(log_e_history) + 1), log_e_history, marker='o')
    plt.title("Log E(M, o) Over Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Log E(M, o)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# ================ HARDOOPS model ================
def run_hard_oops(sequences, K, max_iter=100):
    """Runs the HARDOOPS model. First, initialize the motif offsets randomly in each sequence. It
    stores the log-likelihood after each iteration in a list to track converegence. For every
    iteration, it buids the PFM from the current offsets. Then, it loops over each input 
    sequence and tries every possible k-mer position, and computes its score. The k-mer score 
    and its offset are saved and the offsets are updated. Then, it calculates the log E (M,o)
    of the current motif model and offsets and saves it in the list initialized at the start.
    After the final iteration, the PFM is rebuilt one last time, using the final offsets.

    Args:
        sequences (list): list of sequences from input file
        K (int): user-determined length of motif (6 in this assignment)
        max_iter (int, optional): Maximum number of iterations performed by the model. Defaults to 100.

    Returns:
        pfm: Final position frequency matrix.
        offsets: List of final offsets (0 if motif not present).
        log_e_history: List of log E values over iterations.
    """
    start_time = time.time()
    offsets = initialize_random_offsets(sequences, K)
    log_e_history = []

    for it in range(max_iter):
        pfm = build_pfm(sequences, offsets, K)
        new_offsets = []

        for seq in sequences:
            best_offset = 1
            best_score = float('-inf')

            for i in range(1, len(seq) - K + 2):  # 1-indexed
                kmer = extract_kmer(seq, i, K)
                score = score_kmer(kmer, pfm, K)
                if score > best_score:
                    best_score = score
                    best_offset = i

            new_offsets.append(best_offset)

        offsets = new_offsets
        log_e = compute_log_e(sequences, offsets, pfm, K)
        log_e_history.append(log_e)

    pfm = build_pfm(sequences, offsets, K)
    runtime_seconds = time.time() - start_time

    # Reporting Result
    print(f"Runtime: {runtime_seconds:.2f} seconds")
    print(f"logE(M,o): {log_e_history[-1]:.3f}")  # log E(M, o)
    print("offset:")
    print(offsets)
    print_pfm(pfm)
    return pfm, offsets, log_e_history

# ================ HARDZOOPS model ================
def run_hard_zoops(sequences, K, max_iter=100):
    """Runs the HardZOOPS model: each sequence may or may not contain a motif. First, it initializes 
    random offsets, while recording whether a sequence has a motif or not (presence). It also 
    iniatializes a list to hold the log E values. Then, it runs the E-M loop in a similar way to HardOOPS, but 
    it ensures to update the presence variable. The PFM is built based on the sequences believed to 
    contain a motif. It loops over all possible substrings of length K and finds the best scoring one
    according to the current PFM. It then records the best k-mer and its offset. Finally, log E 
    is computed using only motif-containing sequences. 

    Args:
        sequences (list): list of sequences from input file
        K (int): user-determined length of motif (6 in this assignment)
        max_iter (int, optional): Maximum number of iterations performed by the model. Defaults to 100.

    Returns:
        pfm: Final position frequency matrix.
        offsets: List of final offsets (0 if motif not present).
        log_e_history: List of log E values over iterations.
    """
    start_time = time.time()
    offsets = initialize_random_offsets(sequences, K)
    presence = [True] * len(sequences)
    log_e_history = []

    for it in range(max_iter):
        pfm = build_pfm([seq for i, seq in enumerate(sequences) if presence[i]],
                        [offset for i, offset in enumerate(offsets) if presence[i]],
                        K)

        new_offsets = []
        new_presence = []
        for seq in sequences:
            best_score = float('-inf')
            best_offset = 1

            for i in range(1, len(seq) - K + 2):
                kmer = extract_kmer(seq, i, K)
                score = score_kmer(kmer, pfm, K)
                if score > best_score:
                    best_score = score
                    best_offset = i

            if best_score > 0:  # heuristic threshold
                new_offsets.append(best_offset)
                new_presence.append(True)
            else:
                new_offsets.append(0)
                new_presence.append(False)

        offsets = new_offsets
        presence = new_presence

        # Compute log E only for sequences with motifs
        pfm = build_pfm([seq for i, seq in enumerate(sequences) if presence[i]],
                        [offset for i, offset in enumerate(offsets) if presence[i] and offset > 0],
                        K)
        log_e = compute_log_e([seq for i, seq in enumerate(sequences) if presence[i] and offsets[i] > 0],
                              [offset for i, offset in enumerate(offsets) if presence[i] and offset > 0],
                              pfm, K)
        log_e_history.append(log_e)

    runtime_seconds = time.time() - start_time

    # Reporting Result
    print(f"Runtime: {runtime_seconds:.2f} seconds")
    print(f"logE(M,o): {log_e_history[-1]:.3f}")
    print("offset:")
    print(offsets)
    print_pfm(pfm)

    return pfm, offsets, log_e_history

# ================ GibbsZOOPS model ================
def run_gibbs_zoops(sequences, K, max_iter=100):
    """
    Runs the GibbsZOOPS model using Gibbs sampling. It initializes the offsets and tracks which
    sequences have a motif (similarly to the previous ZOOPS model). It also initializes a list 
    to hold the log likelihoods. Then, for each iteration it performs the following: it leaves one
    sequence out to remove its current motif offset. It then uses the other sequences to build a PFM
    and then it uses the held-out sequence to score all possible k-mers and then convert those 
    scores into probabilities. Based on those probabilities, it samples a new offset. Before selecting
    a new start position from the distribution, it adds a flat baseline probability to account for 
    the possibility that a sequence does not have a motif. Now the model can either sample a 
    start position or 'no motif' from the distribution. If the latter happens, the offset is assigned
    the value 0 and there is no motif. If the sampled index is > 0, then the offsets are updated 
    accordingly. Finally, log E is computed using only motif-containing sequences. 

    Returns:
        pfm, offsets, log_e_history
    """
    import copy
    start_time = time.time()
    N = len(sequences)
    offsets = [random.randint(1, len(seq) - K + 1) for seq in sequences]
    has_motif = [True for _ in sequences]  # ZOOPS: true means sequence has a motif

    log_e_history = []

    for iteration in range(max_iter):
        i = random.randint(0, N - 1)  # randomly select one sequence to leave out
        # build PFM from all sequences except i, and only those marked as 'has motif'
        pfm = build_pfm([s for j, s in enumerate(sequences) if j != i and has_motif[j]],
                        [o for j, o in enumerate(offsets) if j != i and has_motif[j]],
                        K)

        # Score all k-mers in sequence i and pick one with prob. proportional to score
        probs = []
        for offset in range(1, len(sequences[i]) - K + 2):
            kmer = extract_kmer(sequences[i], offset, K)
            score = score_kmer(kmer, pfm, K)
            probs.append(score)

        probs = [math.exp(p) for p in probs]
        total = sum(probs)

        # Add a "no motif" option for ZOOPS
        probs.append(1.0)  # flat baseline probability
        total += 1.0

        probs = [p / total for p in probs]

        sampled_index = random.choices(range(len(probs)), weights=probs)[0]
        if sampled_index == len(probs) - 1:
            has_motif[i] = False
            offsets[i] = 0
        else:
            has_motif[i] = True
            offsets[i] = sampled_index + 1

        log_e = compute_log_e(sequences, offsets, build_pfm(sequences, offsets, K), K)
        log_e_history.append(log_e)

    pfm = build_pfm(sequences, offsets, K)
    runtime_seconds = time.time() - start_time

    # Reporting Result
    print(f"Runtime: {runtime_seconds:.2f} seconds")
    print(f"logE(M,o): {log_e_history[-1]:.3f}")
    print("offset:")
    print(offsets)
    print_pfm(pfm)

    return pfm, offsets, log_e_history


# ================ AnnealZOOPS model ================
def run_anneal_zoops(sequences, K, max_iter=100, T_init=5.0, T_final=0.1):
    """
    Runs the AnnealZOOPS model using simulated annealing. Similar to GibbsZOOPS, this model 
    assumes that each sequence contains at most one motif (ZOOPS) and tracks motif offsets 
    and which sequences are considered motif-containing. It also initializes a log E history 
    list to track progress over iterations.
    The algorithm starts with a high temperature and gradually cools it down over iterations. 
    At each step, it proposes small random changes to the motif offsets. This may include 
    modifying the start position of a motif in a sequence, or toggling whether a sequence 
    is considered to contain a motif at all. The new configuration is accepted if it improves 
    the overall log likelihood (log E), or probabilistically based on the current temperature 
    if the log likelihood is worse. This allows exploration of the search space and avoids 
    getting trapped in local optima early on.
    As the temperature decreases, the model becomes more selective, eventually settling into 
    a high-scoring motif configuration. The log E score is computed only on motif-containing 
    sequences in each iteration. The final PFM and motif offsets reflect the best configuration 
    encountered during the annealing process.

    Returns:
        pfm, offsets, log_e_history
    """
    import copy
    start_time = time.time()
    offsets = initialize_random_offsets(sequences, K)
    presence = [True] * len(sequences)
    log_e_history = []
    T = T_init

    for it in range(max_iter):
        i = random.randint(0, len(sequences) - 1)

        current_offset = offsets[i]
        current_present = presence[i]

        pfm = build_pfm([seq for j, seq in enumerate(sequences) if j != i and presence[j]],
                        [offsets[j] for j in range(len(sequences)) if j != i and presence[j]],
                        K)

        new_offset = random.randint(1, len(sequences[i]) - K + 1)
        new_present = random.random() < 0.5  # randomly turn motif on/off

        offsets_candidate = offsets[:]
        presence_candidate = presence[:]
        offsets_candidate[i] = new_offset
        presence_candidate[i] = new_present

        pfm_candidate = build_pfm([seq for j, seq in enumerate(sequences) if presence_candidate[j]],
                                  [offsets_candidate[j] for j in range(len(sequences)) if presence_candidate[j]],
                                  K)

        log_e_old = compute_log_e([sequences[j] for j in range(len(sequences)) if presence[j]],
                                  [offsets[j] for j in range(len(sequences)) if presence[j]],
                                  pfm, K)
        log_e_new = compute_log_e([sequences[j] for j in range(len(sequences)) if presence_candidate[j]],
                                  [offsets_candidate[j] for j in range(len(sequences)) if presence_candidate[j]],
                                  pfm_candidate, K)

        delta = log_e_new - log_e_old
        accept = False
        if delta > 0:
            accept = True
        else:
            if random.random() < math.exp(delta / T):
                accept = True

        if accept:
            offsets = offsets_candidate
            presence = presence_candidate

        log_e_history.append(log_e_new if accept else log_e_old)
        T = max(T_final, T * 0.95)  # cooling schedule

    pfm = build_pfm([sequences[j] for j in range(len(sequences)) if presence[j]],
                    [offsets[j] for j in range(len(sequences)) if presence[j]],
                    K)
    runtime_seconds = time.time() - start_time
    
    # Reporting Result
    print(f"Runtime: {runtime_seconds:.2f} seconds")
    print(f"logE(M,o): {log_e_history[-1]:.3f}")
    print("offset:")
    print(offsets)
    print_pfm(pfm)

    return pfm, offsets, log_e_history


# ================ Running the program ================
if __name__ == "__main__":

    random.seed(123)

    # Defining Parameters
    K = 6
    MAX_ITER = 100

    models = [
    ("HardOOPS", run_hard_oops),
    ("HardZOOPS", run_hard_zoops),
    ("GibbsZOOPS", run_gibbs_zoops),
    ("AnnealZOOPS", run_anneal_zoops),
    ]

    input_files = glob.glob(os.path.join(data_path, "test_input_*.txt"))

    for input_file in input_files: 
        sequences = read_sequences(input_file)
        filename = os.path.splitext(os.path.basename(input_file))[0] # to be used for the directory/plot file naming later

        print(f"\n### Processing {filename} ###")

        # output dir creation
        out_dir = os.path.join(output_base, filename)
        os.makedirs(out_dir, exist_ok=True)

        output_path = os.path.join(out_dir, "output.txt")

        with open(output_path, "w") as f:
            with redirect_stdout(f):
                print(f"Input File: {input_file}\nK: {K}\n")

                for model_name, model_func in models:
                    print(f"\n# Running {model_name}...")

                    # Running model
                    pfm, offsets, log_e_history = model_func(sequences, K, max_iter=MAX_ITER)

                    # Report results!
                    print(f"Model: {model_name}")
                    print(f"logE(M,o): {log_e_history[-1]:.3f}")
                    print("Offsets:")
                    print(offsets)
                    print_pfm(pfm)

                    # saving plot separately per each model
                    model_plot_dir = os.path.join(out_dir, model_name)
                    os.makedirs(model_plot_dir, exist_ok=True)
                    plot_path = os.path.join(model_plot_dir, "log_e_history.png")
                    plot_log_e_history(log_e_history, filename=plot_path)


        