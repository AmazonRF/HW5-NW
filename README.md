![BuildStatus](https://github.com/AmazonRF/HW5-NW/actions/workflows/pytest.yml/badge.svg?event=push)

# Project 5
Needleman Wunsch Algorithm


# Needleman-Wunsch Sequence Alignment Tests

This document outlines unit tests for validating the Needleman-Wunsch sequence alignment implementation. The tests utilize the BLOSUM62 substitution matrix along with specified gap opening and extension penalties. Test sequences are read from FASTA files, and assertions are made to ensure the alignment scores and the resulting sequences align with expected outcomes.

## Dependencies

The tests require `pytest`, `numpy`, and a `typing` module that implements the Needleman-Wunsch algorithm.


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
