from openfermion.ops import FermionOperator
from openfermion.transforms import jordan_wigner
from openfermion.linalg import get_sparse_operator
import numpy as np
from scipy.linalg import eigh


def simulator_init_cons(dist, chem, ion_count, charge, cons):
    grid_size = dist.shape[0]

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
    Hf = FermionOperator('', 0.0)

    for i1 in range(grid_size):
        for j1, elem in enumerate(chem):
            var = f"{elem}_{i1}"
            if var not in checker or var in fixed:
                continue
            idx = varis[var]
            Hf += FermionOperator(f'{idx}^ {idx}', dist[i1, i1] * charge[elem])

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
                Hf += FermionOperator(f'{idx}^ {idx} {jdx}^ {jdx}', q_same)

    print("=== Objective Function Generated ===")

    # Jordan-Wigner transform
    Hq = jordan_wigner(Hf)

    # Convert to dense matrix
    Hs_dense = get_sparse_operator(Hq).toarray()
    dim = Hs_dense.shape[0]
    n_qubits = int(np.log2(dim))
    k = ion_count[0]   # 원하는 particle number

    # Build projector
    basis = []
    for i in range(2**n_qubits):
        bits = format(i, f'0{n_qubits}b')
        if bits.count('1') == k:
            basis.append(i)

    P = np.zeros((len(basis), dim))
    for idx, state in enumerate(basis):
        P[idx, state] = 1

    # Projected Hamiltonian
    H_proj = P @ Hs_dense @ P.T

    # Diagonalize
    e_vals, e_vecs = eigh(H_proj)

    ground_energy = e_vals[0]
    ground_state_proj = e_vecs[:, 0]

    print(f"\n=== Ground energy (k={k}): {ground_energy:.8f} ===")

    # Map back to bitstrings
    energy, bitstrings = [], []
    for coeff, idx in zip(ground_state_proj, basis):
        prob = abs(coeff)**2
        if prob > 1e-4:
            bitstring = format(idx, f'0{n_qubits}b')
            bitlist = [int(b) for b in bitstring]
            energy.append(ground_energy)
            bitstrings.append(bitlist)
            print(f"{bitstring} : {prob:.6f}")

    return energy, bitstrings


def decode_answer_cons(bits, cons, chem):
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
        if bit:
            position.append(f"{chem[0]}_{c}")

    if fix:
        for idx, pos in enumerate(cons):
            if isinstance(pos, list) and "Fixed" in pos:
                for j in pos[:-1]:  # all except "Fixed"
                    position.append(f"{chem[idx]}_{j}")
            else:
                continue

    return [position]

