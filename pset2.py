# encoding: utf-8
# pset2.py

import numpy as np
# don't forget import packages, e.g. scipy
# but make sure you didn't put unnecessary stuff in here

# INPUT : diag_broadcast - list of diagonals value to broadcast,length equal to 3 or 5; n - integer, band matrix shape.
# OUTPUT : L - 2D np.ndarray, L.shape[0] depends on bandwidth, L.shape[1] = n-1, do not store main diagonal, where all ones;                  add zeros to the right side of rows to handle with changing length of diagonals.
#          U - 2D np.ndarray, U.shape[0] = n, U.shape[1] depends on bandwidth;
#              add zeros to the bottom of columns to handle with changing length of diagonals.
def band_lu(diag_broadcast, n): # 5 pts
    # enter your code here
    if len(diag_broadcast) == 3:
        U = np.zeros((2, n))
        L = np.zeros(n - 1)
        U[0][0] = diag_broadcast[1]
        U[1][0:n - 1] = diag_broadcast[2]
        for i in range(1, n):
            L[i - 1] = diag_broadcast[0] / U[0][i - 1]
            U[0][i] = diag_broadcast[1] - L[i - 1] * diag_broadcast[2]
    if len(diag_broadcast) == 5:
        U = np.zeros((3, n))
        L = np.zeros((2, n - 1))
        U[0][0] = diag_broadcast[2]
        U[1][0] = diag_broadcast[3]
        U[2][0:n - 2] = diag_broadcast[4]
        for i in range(1, n):
            L[0][i - 1] = (diag_broadcast[1] - L[1][i - 2] * U[1][i - 2]) / U[0][i - 1]
            L[1][i - 1] = diag_broadcast[0] / U[0][i - 1]
            U[0][i] = diag_broadcast[2] - L[0][i - 1] * U[1][i - 1] - L[1][i - 2] * U[2][i - 2]
            U[1][i] = diag_broadcast[3] - L[0][i - 1] * U[2][i - 1]
            if i >= n - 1:
                U[1][i] = 0
            if i >= n - 1:
                L[1][i - 1] = 0
    U = U.T
    return L, U


# INPUT : rectangular matrix A
# OUTPUT: matrices Q - orthogonal and R - upper triangular such that A = QR
def gram_schmidt_qr(A): # 5 pts
    # your code is here
    return Q, R

# INPUT : rectangular matrix A
# OUTPUT: matrices Q - orthogonal and R - upper triangular such that A = QR
def modified_gram_schmidt_qr(A): # 5 pts
    # your code is here
    return Q, R


# INPUT : rectangular matrix A
# OUTPUT: matrices Q - orthogonal and R - upper triangular such that A=QR
def householder_qr(A): # 7 pts
    # your code is here
    return Q, R


# INPUT:  G - np.ndarray
# OUTPUT: A - np.ndarray (of size G.shape)
def pagerank_matrix(G): # 5 pts
    # enter your code here
    return A


# INPUT:  A - np.ndarray (2D), x0 - np.ndarray (1D), num_iter - integer (positive) 
# OUTPUT: x - np.ndarray (of size x0), l - float, res - np.ndarray (of size num_iter + 1 [include initial guess])
def power_method(A, x0, num_iter): # 5 pts
    # enter your code here
    return x, l, res


# INPUT:  A - np.ndarray (2D), d - float (from 0.0 to 1.0), x - np.ndarray (1D, size of A.shape[0/1])
# OUTPUT: y - np.ndarray (1D, size of x)
def pagerank_matvec(A, d, x): # 2 pts
    # enter your code here
    return y


def return_words():
    # insert the (word, cosine_similarity) tuples
    # for the words 'numerical', 'linear', 'algebra' words from the notebook
    # into the corresponding lists below
    # words_and_cossim = [('word1', 'cossim1'), ...]
    
    numerical_words_and_cossim = []
    linear_words_and_cossim = []
    algebra_words_and_cossim = []
    
    return numerical_words_and_cossim, linear_words_and_cossim, algebra_words_and_cossim