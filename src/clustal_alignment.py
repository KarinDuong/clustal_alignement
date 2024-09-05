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


import numpy as np
import pandas as pd


def distance_calculate_bis(dist_seq1_seq2, r1, r2, N):
    return (dist_seq1_seq2 - (r1+r2)) / (N-2)


def distance_calculate(dist_matrix):
    dist_dict = dict()
    # d*(S1, S2) = d(S1, S2) - (R1+R2) / (N-2)


def neighbor_joining():
    """_summary_
    
    Parameters
    ----------
    
    Return
    ------
    
    """
    
    # lire le csv / mat de dist
    # with open("matrice_distance_NJ.csv", 'r') as distance_matrix:
    #     print(distance_matrix.read())
    distance_matrix = pd.read_csv("matrice_distance_NJ.csv", index_col=0)
    print(distance_matrix)
    # print(distance_matrix.dtypes)
    
    
    for i in distance_matrix:
        pass
    # boucle 
    #   - calcul dist entre chaque seq : deux à deux 
    #   - fusionner les seq qui ont la dist min : seq-fus
    #   - calcul long des branches entre les seq de dist min avec leur racine
    #   - calcul dist inter-groupe entre seq-fus et les autres seq (les autres val de ma mat reste les memes)
    pass


def upgma():
    # lire mat dist
    distance_matrix = pd.read_csv("matrice_dist_UPGMA.csv", index_col=0)
    print(distance_matrix)
    
    # print(distance_matrix.min(skipna=True).min())
    # print(distance_matrix.stack().idxmin())
    
    dist_dict = dict()
    
    # for row in distance_matrix.values:
    #     for value in row:
    #         print(value)
    for index, row in distance_matrix.iterrows():
        for column in distance_matrix.columns:
            value = row[column]
            if pd.notna(value): print(value)
    
    # fusionner les seq qui ont dist min : seq-fus 
    #   -> dist min/2 = long des branches entre ces deux seq
    # recalculer dist entre seq-fus et les autres seq
    pass


if __name__ == "__main__":
    # neighbor_joining()
    upgma()