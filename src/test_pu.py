import numpy as np
import pandas as pd

def upgma(matrix, labels):
    # Initialisation
    n = len(matrix)
    clusters = {i: [labels[i]] for i in range(n)}
    distances = matrix.copy()

    while len(clusters) > 1:
        # Trouver les deux clusters les plus proches
        min_dist = np.inf
        to_merge = (None, None)
        for i in range(len(distances)):
            for j in range(i):
                if distances[i, j] < min_dist:
                    min_dist = distances[i, j]
                    to_merge = (i, j)

        # Fusionner les deux clusters
        i, j = to_merge
        new_cluster = clusters[i] + clusters[j]
        clusters[i] = new_cluster
        del clusters[j]

        # Mettre à jour la matrice des distances
        new_distances = np.zeros((len(clusters), len(clusters)))
        keys = list(clusters.keys())
        for m in range(len(keys)):
            for n in range(m):
                if keys[m] == i:
                    d1 = len(clusters[keys[n]]) * distances[keys[n], j]
                    d2 = len(clusters[keys[m]]) * distances[keys[m], j]
                    new_distances[m, n] = (d1 + d2) / (len(clusters[keys[m]]) + len(clusters[keys[n]]))
                else:
                    new_distances[m, n] = distances[keys[m], keys[n]]

        distances = new_distances

    return list(clusters.values())

# Exemple de matrice de distances et labels
labels = ['Bsu', 'Bst', 'Lvi', 'Amo', 'Mlu']
data = [
    [0, 0.1715, 0.2147, 0.3091, 0.2326],
    [0, 0, 0.2991, 0.3399, 0.2058],
    [0, 0, 0, 0.2795, 0.3943],
    [0, 0, 0, 0, 0.4289],
    [0, 0, 0, 0, 0]
]

matrix = np.array(data)

# Test de l'algorithme UPGMA
# resultat = upgma(matrix, labels)
# print("Résultat de l'UPGMA:", resultat)


def search_min_ligne(matrix, line_name=None):
    print(type(matrix.index[0]))
    print("matrix.index ", matrix.index)
    print("col_name ", line_name)
    print(matrix)

    if line_name != None:
        val_min = matrix.loc[line_name].min()
        print("line_name rempli ", val_min)
        col_name = matrix.loc[line_name][matrix.loc[line_name] == val_min].index[0] if val_min in matrix.loc[line_name].values else None
        print(line_name, col_name, val_min)
        return val_min, (line_name, col_name)
    else:
        val_min = matrix.min(skipna=True).min()
        print("line_name None ",val_min)
        return val_min, matrix.stack().idxmin()
    
def calc_dist(matrix, row_col, dict_dist):
    name_col_matrix = matrix.columns
    new_dist_list = []

    row, col = row_col

    print(len(matrix))
    print((matrix.shape))
    print("MATRIX debut ------------- \n ",matrix, "\n MATRIX debut------------- \n")
    for i in name_col_matrix: 
        if i not in row_col:
            print("col-i ", col, i)
            print("row-i ", row, i)
            # print("matrix.loc[row, i] [i, row] ",matrix.loc[row, i], matrix.loc[i, row])
            # print("matrix.loc[col, i] [i, col] ",matrix.loc[col, i], matrix.loc[i, col])
            print(name_col_matrix.get_loc(i), name_col_matrix.get_loc(col), name_col_matrix.get_loc(row))
            
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
            print("new_dist_list ", new_dist_list)
            print("\n")

    if len(name_col_matrix) < 2:
        dict_dist[new_row_col] = search_min_ligne(matrix)

    matrix = matrix.drop(index=[row, col])
    matrix = matrix.drop(columns=[col, row])
    
    new_row_col = ', '.join(row_col)
    print("new_row_col ",new_row_col)
    print("new_dist_list ", new_dist_list)
    if len(name_col_matrix) > 2:
        print("len matrix sup 2")
        matrix.loc[new_row_col] = new_dist_list # avec loc : lignes
        matrix[new_row_col] = [np.nan] * len(matrix) # sans loc : colonnes
    # else:
    #     print("len matrix inf 2")
    #     dict_dist[new_row_col] = search_min_ligne(matrix)
        
    # matrix = matrix.drop(index=[row, col])
    # matrix = matrix.drop(columns=[col, row])
        
    print("MATRIX FIN------------- \n ",matrix, "\n MATRIX FIN------------- \n")
    return new_row_col, matrix

matrix = pd.read_csv("matrice_dist_UPGMA.csv", index_col=0)
print(len(matrix))

matrix = matrix.transpose()

dist_min, row_col = search_min_ligne(matrix)

dict_dist = dict()
dict_dist[row_col] = dist_min
print(dict_dist)

new_target, matrix = calc_dist(matrix, row_col, dict_dist)

print("len matrix ", len(matrix))
while len(matrix) > 1:
    print("len matrix debut", len(matrix))
    print("new_target ", new_target)
    dist_min, row_col = search_min_ligne(matrix, new_target)
    dict_dist[row_col] = dist_min
    new_target, matrix = calc_dist(matrix, row_col, dict_dist)
    
    print("len matrix fin", len(matrix))

print(dict_dist)
print(matrix)
