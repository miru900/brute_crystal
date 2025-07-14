from openfermion.ops import FermionOperator
from openfermion.transforms import jordan_wigner
from openfermion.linalg import get_sparse_operator
import numpy as np
from pyscf import fci
from scipy.linalg import eigh

def simulator_init_cons(dist, chem, ion_count, charge, cons):

    grid_size = dist.shape[0]
    norb = len(cons[0]) ; qubo = np.zeros((norb, norb))
    nelec = ion_count[0]

    # variables -> varis and fixed
    varis, fixed, checker = {}, {}, {}
    for i, elem in enumerate(chem):
        for j in range(grid_size):
            if j not in cons[i]:
                continue
            var = f"{elem}_{j}"
            num = grid_size * i + j
            if "Fixed" in cons[i]:
                fixed[var] = num
            else:
                varis[var] = num
            checker[var] = num

    print("=== Variables Debug ===")
    print("varis:", varis)
    print("fixed:", fixed)

    # Objective Function

    for i1 in range(grid_size):
        for j1, elem in enumerate(chem):
            var = f"{elem}_{i1}"
            if var not in checker or var in fixed:
                continue
            idx = varis[var]
            qubo[idx, idx] += dist[i1, i1] * charge[elem]

        for i2 in range(i1 + 1, grid_size):
            for j1, elem in enumerate(chem):
                var1 = f"{elem}_{i1}"
                var2 = f"{elem}_{i2}"
                if var1 not in checker or var2 not in checker:
                    continue
                if var1 in fixed or var2 in fixed:
                    continue
                q_same = 2 * dist[i1, i2] * charge[elem] ** 2
                idx, jdx = varis[var1], varis[var2]
                qubo[idx, jdx] += q_same

    print("=== Objective Function Generated ===")

    return qubo, norb, nelec

def solve_with_pyscf(qubo, norb, nelec):
    h1e_a = np.zeros((norb,norb))
    h1e_b = np.zeros((norb,norb))
    h2e_aa = np.zeros((norb,norb,norb,norb))
    h2e_ab = np.zeros((norb,norb,norb,norb))
    h2e_bb = np.zeros((norb,norb,norb,norb))

    def h1_matrix(h1e, norb, qubo):
        for i in range(norb):
            h1e[i, i] = qubo[i, i]

    def h2_matrix(h2e, norb, qubo):
        for i in range(norb):
            for j in range(i + 1, norb):
                # h2e[j, i, i, j] = -qubo[i, j]
                # h2e[i, j, j, i] = -qubo[i, j]
                # h2e[j, i, j, i] = -qubo[i, j]
                # h2e[i, j, i, j] = -qubo[i, j]
                h2e[i, i, j, j] = qubo[i, j]
                h2e[j, j, i, i] = qubo[i, j]

    # matrix generator
    h1_matrix(h1e_a, norb, qubo)
    h2_matrix(h2e_aa, norb, qubo)

    # tuple generator
    h1 = (h1e_a, h1e_b)
    h2 = (h2e_aa, h2e_ab, h2e_bb)

    # Eigenvalue, Eigenvector calculator
    e, fcivec = fci.direct_uhf.kernel(h1, h2, norb, (nelec, 0), nroots = 10)

    return e, fcivec

        

