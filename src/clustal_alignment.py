"""Multiple sequence alignment (MSA) with clustal method.

Usage:
======
python clustal_alignment.py --input [FASTA file] --gap_score [float value or -8.0 by default] --protein_sequence [boolean value or True by default]
"""

__authors__ = ("Karine Duong")
__contact__ = "karine.duong@etu.u-paris.fr" 
__date__ = "2024-09-05"


import argparse
from Bio import SeqIO
import numpy as np
from pathlib import Path
import pandas as pd

from loguru import logger 


# https://github.com/dmnfarrell/epitopepredict/blob/master/epitopepredict/mhcdata/blosum62.csv
BLOSSUM_MATRIX = pd.read_csv("data/blosum62.csv", index_col=0)

# Source
MATRICE_AN = 0


def extract_sequence_from_fasta(fasta_filename: str, dict_name_seq: dict)-> None:
    """Extract all sequences from the input FASTA file, and store them in a dictionnary.
    
    Parameters
    ----------
        fasta_filename: str
            FASTA filename
        dict_name_seq: dict
            Dictionary to store sequences extracted as value, and sequence name as key.
    """
    records = SeqIO.parse(fasta_filename, "fasta")
    for rec in records:
        dict_name_seq[str(rec.name)] = str(rec.seq)
        
    if len(dict_name_seq) < 3:
        raise Exception("Your FASTA file must have 3 or more sequences.")
    

def needleman_wunsch(seq1: str, seq2: str, gap_score: float, protein_sequence: bool) -> tuple[str, str, float]:
    """Perform pairwise alignment between seq1 and seq2 to find the optimal alignment between theses sequences.
    
    Parameters
    ----------
        seq1: str
            Sequence use in X axis for the alignment.
        seq2: str
            Sequence use in Y axis for the alignment.
        gap_score: float
            Gap score value to calculate and set a value for each case.
        protein_sequence: bool
            Inform the sequence type, to use a appropriate subtitution matrix.
    Return
    ------
        str
            Aligned sequence of seq1, after the optimal alignment was set.
        str
            Aligned sequence of seq2, after the optimal alignment was set.
        float
            Score of the optimal alignment between seq1 and seq2.
    """
    if protein_sequence:
        SUB_MATRIX = BLOSSUM_MATRIX
    else:
        SUB_MATRIX = MATRICE_AN
    
    matrix_NW = np.zeros((len(seq1)+1, len(seq2)+1))
    matrix_backtracking = np.full((len(seq1)+1, len(seq2)+1), '', dtype=str)

    # Initilize the matrix_score
    nb_row, nb_col = matrix_NW.shape
    for i in range(1,nb_row):
        matrix_NW[i,0] = matrix_NW[i-1, 0] + gap_score

    for j in range(1,nb_col):
        matrix_NW[0,j] = matrix_NW[0, j-1] + gap_score
    
    # Fill matrix_score and the backtracking amtrix
    for i in range(1, nb_row):
        for j in range(1, nb_col):
            diag_case = matrix_NW[i-1, j-1] + SUB_MATRIX.loc[seq1[i-1], seq2[j-1]]
            left_case = matrix_NW[i, j-1] + (gap_score)
            top_case = matrix_NW[i-1, j] + (gap_score)
            
            max_value = max(diag_case, left_case, top_case)
            if max_value == diag_case:
                # d as diagonal
                matrix_backtracking[i, j] = "d"
            elif max_value == left_case:
                # l as left
                matrix_backtracking[i, j] = "l"
            else:
                # u as up
                matrix_backtracking[i, j] = "u"
            matrix_NW[i, j] = max_value
        
    # Backtracking to collect the optimal alignment
    seq1_align, seq2_align = "", ""
    i, j = nb_row-1, nb_col-1
    score = matrix_NW[i, j]
    
    while (i>0 or j>0):
        if matrix_backtracking[i, j] == "d":
            seq1_align += seq1[i-1]
            seq2_align += seq2[j-1]
            i -= 1
            j -= 1
        elif matrix_backtracking[i, j] == "u" or j<=0:
            seq1_align += seq1[i-1]
            seq2_align += "_"
            i -= 1
        elif matrix_backtracking[i, j] == "l" or i<=0:
            seq1_align += "_"
            seq2_align += seq2[j-1]
            j -= 1
        else:
            break
        
    return (seq1_align[::-1], seq2_align[::-1], score)


def search_min_ligne_index(matrix: pd.DataFrame, line_name:str=None) -> tuple[float, tuple[str, str]]:
    """Search and return the minimum value in the given matrix. But in the specific line if it's sepecify.
    
    Parameters
    ----------
        matrix: pd.Dataframe
            Distance matrix of differents sequences to do successive branch method.
        line_name: str
            Specific line name to search the minimum value. None by default.

    Return
    ------
        float
            Minimum value finded in the matrix, and in the specific line if specify.
        (str, str)
            Tuple of the line and column name of the minimum value finded.
    """
    if line_name != None:
        val_min = matrix.loc[line_name].min()
        col_name = matrix.loc[line_name][matrix.loc[line_name] == val_min].index[0] if val_min in matrix.loc[line_name].values else None
        return val_min, (line_name, col_name)
    else:
        val_min = matrix.min(skipna=True).min()
        return val_min, matrix.stack().idxmin()


def calc_dist(matrix: pd.DataFrame, row_col: str, dict_dist: dict) -> tuple[str, pd.DataFrame]:
    """Calculate the new distance value between a specific line and column with all the other sequence.
    
    Parameters
    ----------
        matrix: pd.DataFrame
            Distance matrix of differents sequences to do successive branch.
        row_col: str
            Line and column name of the minimum value finded.
        dict_dist: dict
            Dictionary that stock ... as key and ... as value.
    Return
    ------
        str
            Name of the new row and column of the distance we calculate.
        pd.DataFrame
            Matrix with the updated value, row_col ... name removed, row_col line and column name added.
    """
    new_dist_list = []
    name_col_matrix = matrix.columns
    row, col = row_col
    new_row_col = ', '.join(row_col)

    for i in name_col_matrix: 
        if i not in row_col:
            nb_groupe = str(row_col).split(',')
            
            if name_col_matrix.get_loc(i) > name_col_matrix.get_loc(row):
                dist1 = matrix.loc[i, row]
            else:
                dist1 = matrix.loc[row, i]
            
            if name_col_matrix.get_loc(i) > name_col_matrix.get_loc(col):
                dist2 = (matrix.loc[i, col]*(len(nb_groupe)-1))
            else:
                dist2 = (matrix.loc[col, i]*(len(nb_groupe)-1))

            distance = (dist1 + dist2) / len(nb_groupe)
            new_dist_list.append(distance)
            
    if len(name_col_matrix) < 2:
        dict_dist[new_row_col] = search_min_ligne_index(matrix)
        
    matrix = matrix.drop(index=[row, col])
    matrix = matrix.drop(columns=[col, row])
    
    if len(name_col_matrix) > 2:
        matrix.loc[new_row_col] = new_dist_list
        matrix[new_row_col]=[np.nan] * len(matrix)
    
    return new_row_col, matrix


def embranchement_sucessif(matrix: pd.DataFrame) -> dict[str, float]:
    """Do successive branch method from the matrix in the filename.

    Parameters
    ----------
        filename: str
            CSV filename of distance matrix used for successive branch method.
    Return
    ------
        dict[str, float]
            Dictionary that store each cluster and their distance
    """
    dist_min, row_col = search_min_ligne_index(matrix)

    dict_dist = dict()
    dict_dist[row_col] = dist_min

    new_target, matrix = calc_dist(matrix, row_col, dict_dist)

    while len(matrix) > 1:
        dist_min, row_col = search_min_ligne_index(matrix, new_target)
        dict_dist[row_col] = dist_min
        new_target, matrix = calc_dist(matrix, row_col, dict_dist)

    return dict_dist


def needleman_wunsch_after_UPGMA(seq1: list[str], seq2: str, gap_score: float, protein_sequence: bool) -> tuple[str, str, float]:
    """Perform pairwise alignment between seq2 and all the sequence from the cluster seq1 to find the optimal alignment between theses sequences.
    
        Parameters
    ----------
        seq1: list[str]
            Sequences use in X axis for the alignment.
        seq2: str
            Sequence use in Y axis for the alignment.
        gap_score: float
            Gap score value to calculate and set a value for each case.
        protein_sequence: bool
            Inform the sequence type, to use a appropriate subtitution matrix.
    Return
    ------
        str
            Aligned sequence of seq1, after the optimal alignment was set.
        str
            Aligned sequence of seq2, after the optimal alignment was set.
        float
            Score of the optimal alignment between seq1 and seq2.
    """
    if protein_sequence:
        SUB_MATRIX = BLOSSUM_MATRIX
    else:
        SUB_MATRIX = MATRICE_AN
    
    matrix_NW = np.zeros((len(seq1[0])+1, len(seq2)+1))
    matrix_backtracking = np.full((len(seq1[0])+1, len(seq2)+1), '', dtype=str)

    # Initilize the matrix_score
    nb_row, nb_col = matrix_NW.shape
    for i in range(1,nb_row):
        matrix_NW[i,0] = matrix_NW[i-1, 0] + gap_score

    for j in range(1,nb_col):
        matrix_NW[0,j] = matrix_NW[0, j-1] + gap_score
    
    # Fill matrix_score and the backtracking matrix
    for i in range(1, nb_row):
        for j in range(1, nb_col):
            diag_case = matrix_NW[i-1, j-1]
            for x in range(len(seq1)):
                if (seq1[x][i-1] == '_') or  (seq2[j-1] == "_"):
                    diag_case += gap_score
                else:
                    diag_case += BLOSSUM_MATRIX.loc[seq1[x][i-1], seq2[j-1]]
            diag_case = diag_case / len(seq1)
            left_case = matrix_NW[i, j-1] + (gap_score)
            top_case = matrix_NW[i-1, j] + (gap_score)
            
            max_value = max(diag_case, left_case, top_case)
            if max_value == diag_case:
                # d as diagonal
                matrix_backtracking[i, j] = "d"
            elif max_value == left_case:
                # l as left
                matrix_backtracking[i, j] = "l"
            else:
                # u as up
                matrix_backtracking[i, j] = "u"
            matrix_NW[i, j] = max_value
            
    # Backtracking to collect the optimal alignment
    i, j = nb_row-1, nb_col-1
    score = matrix_NW[i, j]
    if len(seq1) > 1:
        seq1_align = [""] * len(seq1)
    else:
        seq1_align = ""
    seq2_align = ""
    
    while (i>0 or j>0):
        if matrix_backtracking[i, j] == "d":
            if len(seq1) > 1:
                for index_list in range(len(seq1)):
                    seq1_align[index_list] += seq1[index_list][i-1]
            else:
                seq1_align += seq1[0][i-1]
            seq2_align += seq2[j-1]
            i -= 1
            j -= 1
        elif matrix_backtracking[i, j] == "u" or j<=0:
            if len(seq1) > 1:
                for index_list in range(len(seq1)):
                    seq1_align[index_list] += seq1[index_list][i-1]
            else:
                seq1_align += seq1[0][i-1]
            seq2_align += "_"
            i -= 1
        elif matrix_backtracking[i, j] == "l" or i<=0:
            if len(seq1) > 1:
                for index_list in range(len(seq1)):
                    seq1_align[index_list] += "_"
            else:
                seq1_align += "_"
            seq2_align += seq2[j-1]
            j -= 1
        else:
            break
    
    if len(seq1) > 1:
        seq1_align_reverse = [s[::-1] for s in seq1_align]
    else:
        seq1_align_reverse = [seq1_align[::-1]]
    
    return (seq1_align_reverse, seq2_align[::-1], score)


def clustal_alignement(filename: str, gap_score: float = -8.0, protein_sequence: bool = True) -> None:
    """Execute the main function to perform the CLUSTAL multiple sequence alignement, with sequences from the input file.

    Parameters
    ----------
        filename: str
            Filepath of a FASTA file
        gap_score: float
            Gap score value. By default at -8.
        protein_sequence: bool
            Inform if sequences from FASTA file are protein or nucleic acids sequence. By default at True.
    """
    logger.info("Debut...")
    dict_pdb_id_seq = dict()
    dict_seq_align = dict()
    
    # Retrieve sequences from FASTA file
    extract_sequence_from_fasta(filename, dict_pdb_id_seq)
    
    # Pairwise alignment
    # And collect their aligned sequence and the score of this optimal alignment
    list_name_seq = list(dict_pdb_id_seq.keys())
    list_seq = list(dict_pdb_id_seq.values())
    for i in range(len(list_seq)):
        for j in range(i + 1, len(list_seq)):            
            seqi_align, seqj_align, score = needleman_wunsch(list_seq[i], list_seq[j], gap_score, protein_sequence)
            
            fusion_name = "--".join([list_name_seq[i], list_name_seq[j]])
            dict_seq_align[fusion_name] = {
                "seq1":list_seq[i],
                "seq2":list_seq[j],
                "seq1_align":seqi_align,
                "seq2_align":seqj_align,
                "score":score,
            }
    
    # Fill the score matrix with the alignment score
    score_matrix = np.full((len(list_name_seq), len(list_name_seq)), np.nan)
    for key, value in dict_seq_align.items():
        name_seq1, name_seq2 = key.split("--")
        index_seq1 = list_name_seq.index(name_seq1)
        index_seq2 = list_name_seq.index(name_seq2)
        score_matrix[index_seq1,index_seq2] = value["score"]
        score_matrix[index_seq2,index_seq1] = value["score"]
        
    # Transform score matrix into distance matrix
    min_value = np.nanmin(score_matrix)
    max_value = np.nanmax(score_matrix)
    
    nb_row, nb_col = score_matrix.shape
    for i in range(0, nb_row):
        for j in range(0, nb_col):
            if i != j:
                score_matrix[i,j] = 1- ((score_matrix[i,j] - min_value) / (max_value - min_value))

    # Perform UPGMA method on this distance matrix
    score_matrix = pd.DataFrame(score_matrix, index=list_name_seq, columns=list_name_seq)
    dict_dist = embranchement_sucessif(score_matrix)
    
    # Performe NW alignment with the order given by UPGMA on all the sequences
    order_UPGMA = list(dict_dist.keys())[-1] 
    list_order_UPGMA = [item.strip() for item in order_UPGMA[0].split(",")]
    list_order_UPGMA.append(order_UPGMA[1])
    
    print(list_order_UPGMA)
    list_seq1, seq2 = [], ""
    for i in range(1, len(list_order_UPGMA)):     
        if len(list_seq1) < 1:
            # If it's my first iteration
            # I begin my MSA with the previous aligned sequence (with NW)
            # depend the order gived by UPGMA
            fusion_name_1 = list_order_UPGMA[i-1]+"--"+list_order_UPGMA[i]
            fusion_name_2 = list_order_UPGMA[i]+"--"+list_order_UPGMA[i-1]
            
            NW_align = dict_seq_align[fusion_name_1] or dict_seq_align[fusion_name_2]
            seq1_align_NW = NW_align["seq1_align"]
            seq2_align_NW = NW_align["seq2_align"]

            list_seq1.append(seq1_align_NW)
            seq2 = seq2_align_NW

        else:
            # Otherwise, my cluster is a list of list_seq1 and seq2 
            # (the result of their alignment, so the result of the previous iteration)
            # and seq2 is the sequence (from the FASTA) depend on the iteration order gived by UPGMA
            list_seq1.append(seq2)
            seq2 = dict_pdb_id_seq[list_order_UPGMA[i]]
        
        # Apply NW alignment between a cluster list_seq1 and the sequence seq2
        list_seq1, seq2, score = needleman_wunsch_after_UPGMA(list_seq1, seq2, gap_score, protein_sequence)
    list_seq1.append(seq2)
    
    for i in range(len(list_order_UPGMA)):
        print(f"{list_order_UPGMA[i]}: {list_seq1[i]}")
    logger.info("Fin...")

def test_exist_type_file(filepath: str) -> str:
    """Check if the given filepath points to an existing structure file (.fasta).

    Parameters
    ----------
        filepath : str
            Path of the file.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not an existing file,
            or if it does not have a '.fasta' extension.

    Returns
    -------
        str
            The validated path.
    """
    filename = Path(filepath)
    if not Path.is_file(filename):
        raise argparse.ArgumentTypeError(f"{filepath} does not exist")

    if filename.suffix != ".fasta":
        raise argparse.ArgumentTypeError(f"{filepath} is not a .fasta file.")
    return filepath


def check_number(score_value: str) -> float:
    """Check if the input gap score is a number

    Parameter
    ---------
        score_value: str
            Gap score value

    Raises
    ------
        argparse.ArgumentTypeError
            If the score given in the input is not a number

    Return
    ------
        float
            Value of input into float
    """
    try:
        score_as_float = float(score_value)
    except argparse.ArgumentTypeError:
        raise argparse.ArgumentTypeError("Argument should be a number")
    return score_as_float


def parse_arg() -> argparse.Namespace:
    """Parse command-line arguments.
    
    Ressources
    ----------
    - https://docs.python.org/3/library/argparse.html

    Return
    ------
        argparse.Namespace
            An object containing the parsed arguments as attributes.
    """
    parser = argparse.ArgumentParser(
        prog="clustal_alignment", 
        description="", 
        usage="clustal_alignment.py [-u] --input fasta_filepath "
    )
    parser.add_argument(
        "--input", 
        type=test_exist_type_file,
        help="Filepath of the fasta file to retrieve the sequence",
        required=True,
    )
    parser.add_argument(
        "--gap_score", 
        type=check_number,
        help="Define a value for the gap score. By default at -8.",
        default=-8.0,
    )
    parser.add_argument(
        "--protein_sequence", 
        help="Inform if the sequence type are a protein sequence or a nucleic acids sequence. Set by default at True",
        default=True,
        action="store_false",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    clustal_alignement(args.input, 
                       args.gap_score,
                       args.protein_sequence)
    