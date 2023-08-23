import numpy as np

def qr_decomposition(matrix):

    m, n = matrix.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    for j in range(n):
        v = matrix[:, j]
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], matrix[:, j])
            v = v - R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]

    return Q, R

def find_eigenvalues(matrix, tolerance=1e-6):

    n = matrix.shape[0]
    eigenvalues = []

    while n > 1:
        Q, R = qr_decomposition(matrix)
        matrix = np.dot(R, Q)
        eigenvalue = matrix[n-1, n-1]
        eigenvalues.append(eigenvalue)
        if np.abs(matrix[:n-1, :n-1] - eigenvalue * np.identity(n-1)).max() < tolerance:
            matrix = matrix[n-1:n, n-1:n]
            n = 1
        else:
            n -= 1

    eigenvalues.append(matrix[0, 0])
    return eigenvalues

matrix = np.array([[1, 7, -1], [-2, 2, -2], [9, -7, 3]])
eigenvalues = find_eigenvalues(matrix)
print(eigenvalues)

