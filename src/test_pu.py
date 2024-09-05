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
resultat = upgma(matrix, labels)
print("Résultat de l'UPGMA:", resultat)
