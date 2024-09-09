"""[Description du le but du fichier script]

Usage:
======
[cmt l'utiliser : ligne de commande]
    python nom_de_ce_super_script.py argument1 argument2
    argument1: un entier signifiant un truc
    argument2: une chaîne de caractères décrivant un bidule
"""

__authors__ = ("Karine Duong")
__contact__ = "karine.duong@etu.u-paris.fr" 
__date__ = "2024-09-05"


import argparse
from Bio import SeqIO
import numpy as np
from pathlib import Path
import pandas as pd
import sys

# https://github.com/dmnfarrell/epitopepredict/blob/master/epitopepredict/mhcdata/blosum62.csv
BLOSSUM_MATRIX = pd.read_csv("blosum62.csv", index_col=0)


def extract_sequence_from_fasta(fasta_filename: str, dict_name_seq: dict)-> None:
    records = SeqIO.parse(fasta_filename, "fasta")
    for rec in records:
        dict_name_seq[str(rec.name)] = str(rec.seq)
    

def needleman_wunsch(seq1: str, seq2: str, gap_score: float = -8.0) -> tuple[str, str, float]:
    """_summary_
    
    Parameters
    ----------
        seq1: str
        seq2: str
        gap_score: float
    
    Return
    ------
    
    """    
    matrix_NW = np.zeros((len(seq1)+1, len(seq2)+1))
    matrix_backtracking = np.full((len(seq1)+1, len(seq2)+1), '', dtype=str)

    #initilize the matrix_score
    nb_row, nb_col = matrix_NW.shape
    for i in range(1,nb_row):
        matrix_NW[i,0] = matrix_NW[i-1, 0] + gap_score

    for j in range(1,nb_col):
        matrix_NW[0,j] = matrix_NW[0, j-1] + gap_score
    
    # fill matrix_score
    for i in range(1, nb_row):
        for j in range(1, nb_col):
            diag_case = matrix_NW[i-1, j-1] + BLOSSUM_MATRIX.loc[seq1[i-1], seq2[j-1]]
            left_case = matrix_NW[i, j-1] + (gap_score)
            top_case = matrix_NW[i-1, j] + (gap_score)
            
            max_value = max(diag_case, left_case, top_case)
            if max_value == diag_case:
                matrix_backtracking[i, j] = "diag"
            elif max_value == left_case:
                matrix_backtracking[i, j] = "gauche"
            else:
                matrix_backtracking[i, j] = "haut"
            matrix_NW[i, j] = max_value
        
    # backtracking
    seq1_align, seq2_align = "", ""
    i, j = nb_row-1, nb_col-1
    score = matrix_NW[i, j]
    
    while (i>0 or j>0):
        if matrix_backtracking[i, j] == "d":
            seq1_align += seq1[i-1]
            seq2_align += seq2[j-1]
            i -= 1
            j -= 1
        elif matrix_backtracking[i, j] == "h" or j<=0:
            seq1_align += seq1[i-1]
            seq2_align += "_"
            i -= 1
        elif matrix_backtracking[i, j] == "g" or i<=0:
            seq1_align += "_"
            seq2_align += seq2[j-1]
            j -= 1
        else:
            break
        
    # print(f"{seq1_align[::-1]}\n{seq2_align[::-1]}\n{score}")
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
    """Caculate the new distance value between a specific line and column with all the other sequence.
    
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
        matrix.loc[new_row_col] = new_dist_list # avec loc : cherche dans les lignes
        matrix[new_row_col]=[np.nan] * len(matrix) # sans loc : cherche sur les colonnes
    
    print(matrix)
    return new_row_col, matrix


def embranchement_sucessif(matrix: pd.DataFrame) -> None:
    """Do successive branch method from the matrix in the filename.

    Parameters
    ----------
        filename: str
            CSV filename of distance matrix used for successive branch method.
    Return
    ------
        None
    """
    # matrix = pd.read_csv(filename, index_col=0)
    # matrix = matrix.transpose() # pour avoir le triangle inférieur
    dist_min, row_col = search_min_ligne_index(matrix)

    dict_dist = dict()
    dict_dist[row_col] = dist_min

    new_target, matrix = calc_dist(matrix, row_col, dict_dist)

    while len(matrix) > 1:
        dist_min, row_col = search_min_ligne_index(matrix, new_target)
        dict_dist[row_col] = dist_min
        new_target, matrix = calc_dist(matrix, row_col, dict_dist)

    print(dict_dist)
    return dict_dist


def clustal_alignement(filename: str):
    """_summary_

    Parameters
    ----------
        filename: str
            Filepath given with arg.parse
    """
    dict_pdb_id_seq = dict()
    dict_seq_align = dict()
    
    # Retrieve sequences from FASTA file
    extract_sequence_from_fasta(filename, dict_pdb_id_seq)
    
    # Align each sequence two by two
    # And collect their aligned sequence and the score of the alignment
    list_name_seq = list(dict_pdb_id_seq.keys())
    list_seq = list(dict_pdb_id_seq.values())
    for i in range(len(list_seq)):
        for j in range(i + 1, len(list_seq)):            
            seqi_align, seqj_align, score = needleman_wunsch(list_seq[i], list_seq[j])
            
            fusion_name = "--".join([list_name_seq[i], list_name_seq[j]])
            dict_seq_align[fusion_name] = {
                "seq1":list_seq[i],
                "seq2":list_seq[j],
                "seq1_align":seqi_align,
                "seq2_align":seqj_align,
                "score":score,
            }

    # fill the score matrix with the alignment score
    score_matrix = np.full((len(list_name_seq), len(list_name_seq)), np.nan)
    for key, value in dict_seq_align.items():
        name_seq1, name_seq2 = key.split("--")
        index_seq1 = list_name_seq.index(name_seq1)
        index_seq2 = list_name_seq.index(name_seq2)
        score_matrix[index_seq1,index_seq2] = value["score"]
        score_matrix[index_seq2,index_seq1] = value["score"]
        
    # Transform score matrix into distance matrix
    
        
    # Perform UPGMA method on this distance matrix
    score_matrix = pd.DataFrame(score_matrix, index=list_name_seq, columns=list_name_seq)
    dict_dist = embranchement_sucessif(score_matrix)
    print(dict_dist)


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
        usage="clustal_alignment.py [-h] --input fasta_filepath "
    )
    parser.add_argument(
        "--input", 
        type=test_exist_type_file,
        help="Filepath of the fasta file to retrieve the sequence",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()    
    # extract_sequence_from_fasta(args.input, dict_name_seq)
    # print(dict_name_seq)
    clustal_alignement(args.input)