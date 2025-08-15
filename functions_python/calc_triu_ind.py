import numpy as np

def calc_triu_ind(matrix):
    M, _ = matrix.shape
    return np.triu_indices(M, k=1)
