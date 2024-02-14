# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    # Define gap penalties
    gap_open = -10
    gap_extend = -1
    
    # Read sequences from FASTA files
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
   
    # Expected results for the alignment
    score = 4.0
    alignedA = "MYQR"
    alignedB = "M-QR"

    # Perform the alignment using NeedlemanWunsch class
    nw = NeedlemanWunsch(gap_open=gap_open, gap_extend=gap_extend, sub_matrix_file="substitution_matrices/BLOSUM62.mat")
    test_score_12, test_alignedA, test_alignedB = nw.align(seq1, seq2)

    # Assertions to check if the alignment results match the expected outcomes
    assert test_score_12 == score, 'Proposed NW has incorrect expected score'
    assert test_alignedA == alignedA, 'Proposed NW has incorrect aligned output'
    assert test_alignedB == alignedB, 'Proposed NW has incorrect aligned output'
    
    # Assertions to check the integrity of the scoring and traceback matrices
    assert np.all(nw.traceback_matrix == np.array([[0, 2, 2, 2],[1, 3, 2, 2],[1, 1, 3, 2],[1, 1, 3, 3],[1, 1, 3, 3]])), 'Incorrect Traceback Matrix'
    assert np.all(nw.score_matrix == np.array([[  0., -10., -11., -12.],[-10.,   5.,  -6.,  -7.],[-11.,  -6.,   4.,  -7.],[-12.,  -7.,  -1.,   5.],[-13.,  -8.,  -6.,   4.]])), 'Incorrect Traceback Matrix'
    
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """

    # Define gap penalties
    gap_open = -10
    gap_extend = -1

    # Expected results for the alignment
    alignedA = 'MAVHQLIRRP'
    alignedB = 'M---QLIRHP'

    # Read sequences from FASTA files
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    # Perform the alignment using NeedlemanWunsch class
    nw = NeedlemanWunsch(gap_open=gap_open, gap_extend=gap_extend, sub_matrix_file="substitution_matrices/BLOSUM62.mat")
    test_score_34, test_alignedA, test_alignedB = nw.align(seq3, seq4)

    # Assertions to check if the alignment results match the expected outcomes
    assert test_score_34 == 17, 'Proposed NW has incorrect expected score'
    assert test_alignedA == alignedA, 'Proposed NW has incorrect aligned output'
    assert test_alignedB == alignedB, 'Proposed NW has incorrect aligned output'




