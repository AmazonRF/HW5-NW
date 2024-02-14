# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""
        

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        """Perform global sequence alignment using Needleman-Wunsch algorithm."""
        self.seqA = seqA
        self.seqB = seqB
        lenA = len(seqA)
        lenB = len(seqB)

        self.score_matrix = np.zeros((lenA + 1, lenB + 1))
        self.traceback_matrix = np.zeros((lenA + 1, lenB + 1), dtype=np.int8)

        # Initialize gap penalties for the first row and column
        for i in range(1, lenA + 1):
            self.score_matrix[i][0] = self.gap_open + (i - 1) * self.gap_extend
            self.traceback_matrix[i][0] = 1  # Indicates an up direction

        for j in range(1, lenB + 1):
            self.score_matrix[0][j] = self.gap_open + (j - 1) * self.gap_extend
            self.traceback_matrix[0][j] = 2  # Indicates a left direction


        # self.initialize_matrices(len(seqA), len(seqB))
        # self.fill_matrices()
        """Fill score and traceback matrices."""
        for i in range(1, len(self.seqA) + 1):
            for j in range(1, len(self.seqB) + 1):
                """Calculate scores for matching, insertion, and deletion."""
                dict_match_score = self.sub_dict[seqA[i - 1], seqB[j - 1]]
                match_score = self.score_matrix[i - 1, j - 1] + dict_match_score
                if self.traceback_matrix[i - 1, j] == 1:
                    del_score =  self.score_matrix[i - 1, j] + self.gap_extend 
                else: del_score =  self.score_matrix[i - 1, j] + self.gap_open + self.gap_extend
                
                if self.traceback_matrix[i, j - 1] == 2:
                    ins_score = self.score_matrix[i, j - 1] + self.gap_extend
                else: 
                    ins_score = self.score_matrix[i, j - 1] + self.gap_open+ self.gap_extend
                # match, delete, insert = self.calculate_score(i, j)
                best_score = max(match_score, del_score, ins_score)
                self.score_matrix[i, j] = best_score
                
                if best_score == match_score:
                    self.traceback_matrix[i, j] = 3  # Diagonal
                elif best_score == del_score:
                    self.traceback_matrix[i, j] = 1  # Up
                else:
                    self.traceback_matrix[i, j] = 2  # Left
        return self._backtrace()    		
        		    

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        """Reconstruct alignment from the traceback matrix."""
        alignmentA, alignmentB = '', ''
        i, j = len(self.seqA), len(self.seqB)

        while i > 0 or j > 0:
            if self.traceback_matrix[i, j] == 3:  # Diagonal
                alignmentA = self.seqA[i-1] + alignmentA
                alignmentB = self.seqB[j-1] + alignmentB
                i -= 1
                j -= 1
            elif self.traceback_matrix[i, j] == 1:  # Up
                alignmentA = self.seqA[i-1] + alignmentA
                alignmentB = '-' + alignmentB
                i -= 1
            elif self.traceback_matrix[i, j] == 2:  # Left
                alignmentA = '-' + alignmentA
                alignmentB = self.seqB[j-1] + alignmentB
                j -= 1

        return self.score_matrix[-1, -1], alignmentA, alignmentB



def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
