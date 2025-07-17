from openfermion.ops import FermionOperator
from openfermion.transforms import jordan_wigner
from openfermion.linalg import get_sparse_operator
import numpy as np
from pyscf import fci
from scipy.linalg import eigh
from pyscf.fci.cistring import addr2str

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
    constant = 0

    for i1 in range(grid_size):
        for j1, elem in enumerate(chem):
            var = f"{elem}_{i1}"
            if var not in checker:
                continue
            if var in fixed:
                constant += dist[i1, i1] * charge[elem] ** 2
                continue
            idx = varis[var]
            qubo[idx, idx] += dist[i1, i1] * charge[elem] ** 2

        for i2 in range(i1 + 1, grid_size):
            for j1, elem in enumerate(chem):
                var1 = f"{elem}_{i1}"
                var2 = f"{elem}_{i2}"
                if var1 not in checker or var2 not in checker:
                    continue
                if var1 in fixed and var2 in fixed:
                    constant += 2 * dist[i1, i2] * charge[elem] ** 2
                    continue
                if var1 in fixed or var2 in fixed:
                    continue
                q_same = 2 * dist[i1, i2] * charge[elem] ** 2
                idx, jdx = varis[var1], varis[var2]
                qubo[idx, jdx] += q_same

                for j2 in range(j1 + 1, len(chem)):
                    var2 = f"{chem[j2]}_{i2}"
                    var3 = f"{chem[j2]}_{i1}"
                    if var2 not in checker or var3 not in checker:
                        continue
                    if var1 in fixed  and var2 in fixed:
                        constant += 4 * dist[i1, i2] * charge[elem] * charge[chem[j2]]
                        continue
                    elif var1 not in fixed and var2 in fixed:
                        q_cross1 = 2 * dist[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        qubo[i1, i1] += q_cross1
                        q_cross2 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        qubo[i2, i2] += q_cross2
                        continue
                    elif var1 in fixed and var2 not in fixed:
                        continue

    print("=== Objective Function Generated ===")
    print(constant)

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
    e, fcivec = fci.direct_uhf.kernel(h1, h2, norb, (nelec, 0), nroots = 924)

    return e, fcivec

def binary_treat(binary, norb):
    binary = binary[2 : ]
    if len(binary) < norb:
        binary = "0" * (norb - len(binary)) + binary

    answer = binary[::-1]

    return answer

def decode_answer(fcivec, norb, nelec):
    bits = []

    for binary in fcivec:
        for i in range(len(binary)):
            if binary[i][0] > 0.99:
                answer = bin(addr2str(norb, nelec, i))
                answer = binary_treat(answer, norb)
                bits.append(answer)

    return bits

def bits_to_answer(bits, cons, chem):
    """
    Decode bitstring back to atom positions
    bits: list of bitlists [[0,1,0,...]]
    cons: constraint list [[0,1,2,...]]
    chem: list of chemical species names
    """
    fix = False
    for pos in cons:
        if isinstance(pos, list) and "Fixed" in pos:
            fix = True
            break

    position = []
    bits = bits[0]
    print("Constraints:", cons)

    for c, bit in zip(cons[0], bits):
        if int(bit):
            position.append(f"{chem[0]}_{c}")

    if fix:
        for idx, pos in enumerate(cons):
            if isinstance(pos, list) and "Fixed" in pos:
                for j in pos[:-1]:  # all except "Fixed"
                    position.append(f"{chem[idx]}_{j}")
            else:
                continue

    return [position]
