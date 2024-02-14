# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    gap_open = -10
    gap_extend = -1

    nw = NeedlemanWunsch(gap_open=gap_open, gap_extend=gap_extend, sub_matrix_file="substitution_matrices/BLOSUM62.mat")
    score_hg, alignedA, alignedB = nw.align(hs_seq, gg_seq)
    score_hm, alignedA, alignedB = nw.align(hs_seq, mm_seq)
    score_hr, alignedA, alignedB = nw.align(hs_seq, br_seq)
    score_ht, alignedA, alignedB = nw.align(hs_seq, tt_seq)

    

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nameList = ['Gallus_gallus','Mus_musculus','Balaeniceps','tursiops']
    scoreList = np.array([score_hg,score_hm,score_hr,score_ht])

    sorted_indices = np.argsort(scoreList)

    # Use sorted_indices to sort nameList
    sorted_nameList = [nameList[i] for i in sorted_indices]

    print("Original nameList:", nameList)
    print("Sorted nameList:", sorted_nameList)
    

if __name__ == "__main__":
    main()
