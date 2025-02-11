import numpy as np
from scipy.constants import epsilon_0 as ep0

def elec_pot(pos_data, cell_size):
    N = len(pos_data)
    gen_list = np.zeros((N, N))

    for i in range(N):
        for j in range(i, N):
            r = np.linalg.norm(pos_data[i] - pos_data[j])  # distance calculation by norm

            element = sum(1 / ((r + (k + 1) * cell_size) * (4 * np.pi * ep0)) for k in range(10 if r <= 1e-6 else 11))

            gen_list[i, j] = element
            gen_list[j, i] = element  # Symmetry

    return gen_list

