import numpy as np
from numpy.linalg import eigvals

# Residues: 24 units mod 90
R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

# Inversion graph adjacency matrix
n = len(R)
adj = np.zeros((n, n))
for i, u in enumerate(R):
    inv_u = pow(u, -1, 90)
    j = R.index(inv_u)
    adj[i, j] = 1  # Undirected, includes loops if self-inverse

# Group spectrum: eigenvalues of adjacency
spectrum = np.sort(np.real(eigvals(adj)))

print("Group Spectrum (sorted eigenvalues):")
print(spectrum)