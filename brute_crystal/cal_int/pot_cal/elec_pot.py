import numpy as np
from scipy.constants import epsilon_0 as ep0
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee

def elec_pot(pos_data, cell_size, ortho = False):
    N = len(pos_data)
    gen_list = np.zeros((N, N))

    for i in range(N):
        for j in range(i, N):
            r = np.linalg.norm(pos_data[i] - pos_data[j])  # distance calculation by norm
            
            if ortho:
                element = sum(1 / ((r + (k + 1) * np.average(cell_size)) * (4 * np.pi * ep0)) for k in range(10 if r <= 1e-6 else 11))

            else:
                element = sum(1 / ((r + (k + 1) * cell_size) * (4 * np.pi * ep0)) for k in range(10 if r <= 1e-6 else 11))

            gen_list[i, j] = element
            gen_list[j, i] = element  # Symmetry

    return gen_list



def energy_elec_addup(O, T, Vars, orb_key, matrix, charge, energy, chem):

    # Self interaction term
    for i in range(O):
        for j in range(T):
            energy.add(Vars[j][orb_key[i]] * Vars[j][orb_key[i]] * matrix[i][i] * charge[chem[j]] ** 2 * elec_c ** 2 / eevee)

    # Coulomb interaction term, considering only upper triangle region of matrix
    for i in range(O):
        for j in range(i + 1, O):
            for k1 in range(T):
                for k2 in range(T):
                    interaction = (Vars[k1][orb_key[i]] * Vars[k2][orb_key[j]]* 2 * matrix[i][j] * charge[chem[k1]] * charge[chem[k2]] * elec_c ** 2 / eevee)
                    energy.add(interaction)


