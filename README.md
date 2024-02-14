![BuildStatus](https://github.com/AmazonRF/HW5-NW/actions/workflows/pytest.yml/badge.svg?event=push)

# Needleman-Wunsch Sequence Alignment Implementation and Validation

This project includes an implementation of the Needleman-Wunsch algorithm to perform global sequence alignment between two biological sequences. The implementation is provided in Python and utilizes both standard Python libraries and NumPy for efficient matrix operations. Additionally, the project includes validation functions to ensure the accuracy and correctness of the alignment results based on scoring matrices and gap penalties.

## Features

- **Needleman-Wunsch Algorithm Implementation:** Performs global sequence alignment using dynamic programming, producing an optimal alignment between two sequences.
- **Scoring Matrix Support:** Utilizes popular scoring matrices like BLOSUM62 for amino acid alignments, ensuring biologically relevant alignments.
- **Gap Penalty Configuration:** Allows for the specification of gap opening and extension penalties, accommodating different biological scenarios and sequence types.
- **Sequence Input:** Supports sequence inputs directly as strings or from FASTA files, providing flexibility in handling sequence data.

## Getting Started

### Prerequisites

- Python 3.x
- NumPy

### Installation

Ensure you have Python and NumPy installed on your system. If NumPy is not already installed, you can install it using the following command:

```bash
pip install numpy
```

## Usage

### Sequence Initialization

First, prepare your sequences. They can be input directly as string variables or loaded from FASTA files using a helper function:

```python
from align import NeedlemanWunsch, read_fasta

# Direct string input
seqA = "ACDEFGHIKLMNPQRSTVWY"
seqB = "ABCDEFGHIKLMNPQRSTVWY"

# Loading from FASTA files
seq1, _ = read_fasta('path/to/sequence1.fasta')
seq2, _ = read_fasta('path/to/sequence2.fasta')
```
### Performing Sequence Alignment

Instantiate the Needleman-Wunsch class with the desired scoring matrix and gap penalties, then call the `align` method with your sequences:

```python
# Initialize the Needleman-Wunsch aligner
nw_aligner = NeedlemanWunsch(sub_matrix_file="path/to/BLOSUM62", gap_open=-10, gap_extend=-1)

# Perform the alignment
alignment_score, aligned_seq1, aligned_seq2 = nw_aligner.align(seqA, seqB)

# Output the results
print(f"Alignment Score: {alignment_score}")
print(f"Aligned Sequence 1: {aligned_seq1}")
print(f"Aligned Sequence 2: {aligned_seq2}")
```

# Species Order

In this project, different species gene are applied to alignment and analysis base on the sequence similarity with human. The order of species is as follows:

1. `Balaeniceps`
2. `Gallus_gallus`
3. `Mus_musculus`
4. `tursiops`


# Grading
## Code (6 points)
* Pairwise global alignment works properly (6)
    * Correct implementation of Needleman-Wunsch algorithm (4)
    * Produces correct order of species in main.py (1) 
    * Produces correct NW alignment scores in main.py (1)

## Unit tests (3 points)
* `test_nw_alignment` function properly checks that matrices are filled in correctly for alignment of test_seq1.fa and test_seq2.fa (1)
* `test_nw_backtrace` function properly checks that backtrace works correctly (1)
* Ensure functionality with pytest (1)
## Style (1 points)
* Readable code with clear comments and method descriptions (1)
## Extra credit (0.5)
* Github actions/workflow (0.5)
