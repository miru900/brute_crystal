from openfermion.ops import FermionOperator
from openfermion.transforms import jordan_wigner
from openfermion.linalg import get_sparse_operator
from scipy.sparse.linalg import eigsh
import numpy as np

import cal_fermion.cons_matrix as cfcm

def simulator_init_cons(dist, chem, ion_count, charge, cons):
    grid_size = dist.shape[0]
    Q = dist # For Fixed Atoms
    W = dist # For variable atoms, constraints is added to this matrix
    
    W = cfcm.matrix_cons_adder(W, cons, ion_count) # This is only working for only one type of atom, hardcoded!

    # List of Variable atoms and fixed atoms
    varis = {} ; fixed = {} ; checker = {}

    for i, elem in enumerate(chem):
        for j in range(grid_size):
            
            # Only Allowed Variables are added to dictionary.
            if not j in cons[i]:
                pass
            
            else:
                # Fixed and not Fixed variables are distinguished.
                if "Fixed" in cons[i]:
                    var = f"{elem}_{j}"
                    num = grid_size * i + j
                    fixed[var] = num
                    checker[var] = num

                else:
                    var = f"{elem}_{j}"
                    num = grid_size * i + j
                    varis[var] = num
                    checker[var] = num

    # Debug Part
    print("Fermion Variables Debugging")
    print(varis, fixed, checker)

    # Objective Function Generating Part!
    Hf = FermionOperator('', 0.0)

    for i1 in range(grid_size):
        for j1, elem in enumerate(chem):
            var = f"{elem}_{i1}"

            if var not in checker or var in fixed:
                continue
            else:
                idx = varis[var] 
                Hf += FermionOperator(f'{idx}^ {idx}', W[i1, i1] * charge[elem])

        for i2 in range(i1 + 1, grid_size):
            for j1, elem in enumerate(chem):
                var1 = f"{elem}_{i1}" ; var2 = f"{elem}_{i2}"

                if var1 not in checker or var2 not in checker:
                    continue
                elif var1 in fixed or var2 in fixed:
                    pass
                else:
                    q_same = 2 * W[i1, i2] * charge[elem] ** 2
                    idx = varis[var1] ; jdx = varis[var2]
                    Hf += FermionOperator(f'{idx}^ {idx} {jdx}^ {jdx}', q_same)
                    
                for j2 in range(j1 + 1, len(chem)):
                    var3 = f"{chem[j2]}_{i2}"
                    var4 = f"{chem[j2]}_{i1}"
                        
                    if var1 not in checker or var2 not in checker or var3 not in checker or var4 not in checker:
                        continue

                    elif var2 in fixed and var3 in fixed:
                        continue

                    elif var2 not in fixed and var3 in fixed:
                        q_cross1 = 2 * Q[i1, i2] * charge[chem[j1]] * charge[chem[j2]] ; idx = varis[var1]
                        Hf += FermionOperator(f'{idx}^ {idx}', q_cross1)
                        q_cross2 = 2 * Q[i1, i2] * charge[chem[j1]] * charge[chem[j2]] ; idx = varis[var2]
                        Hf += FermionOperator(f'{idx}^ {idx}', q_cross2)
                        continue

                    elif var2 in fixed and var3 not in fixed:
                        q_cross1 = 2 * Q[i1, i2] * charge[chem[j1]] * charge[chem[j2]] ; idx = varis[var3]
                        Hf += FermionOperator(f'{idx}^ {idx}', q_cross1)
                        q_cross2 = 2 * Q[i1, i2] * charge[chem[j1]] * charge[chem[j2]] ; idx = varis[var4] 
                        Hf += FermionOperator(f'{idx}^ {idx}', q_cross2)
                        continue

    print("Objective Function Generated")

    # Jordan-Wigner transform to qubit operator
    Hq = jordan_wigner(Hf)

    # Convert to sparse matrix and diagonalize
    Hs = get_sparse_operator(Hq)
    e_vals, e_vecs = eigsh(Hs, k=1, which='SA')

    energy = []
    bit = []

    # Show ground energy
    print("\nGround energy: {:.8f}".format(e_vals[0]))

    # Compute and show probabilities
    print("\nGround state probabilities:")
    probs = np.abs(e_vecs[:, 0]) ** 2
    n = int(np.log2(len(probs)))
    
    for i, p in enumerate(probs):
        if p > 1e-4:
            bitstring = format(i, f'0{n}b')  # 문자열 예: '1100'
            bitlist = [int(x) for x in bitstring]  # [1,1,0,0]

            # save
            energy.append(e_vals[0])
            bit.append(bitlist)

            # print
            print(f"{bitstring} : {p:.6f}")

    return energy, bits

def decode_answer_cons(bits, cons, chem):
    """
    result: OptimizationResult returned by Qiskit solver
    qp: The QuadraticProgram used (must match variable order)
    """
    fix = False
    # fixed checker:
    for pos in cons:
        if isinstance(pos, list) and "Fixed" in pos:
            fix = True
            break

    position = []
    bits = bits[0]

    for cons, bit in zip(cons[0], bits):
        if bit:
            position.append(f"{chem[0]}_{cons}")

    if fix:
        for idx, pos in enumerate(cons):
            if isinstance(pos, list) and "Fixed" in pos:
                for j in pos[:-1]:  # all except "Fixed"
                    position.append(chem[idx] + "_" + str(j))
            else:
                continue

    return [position]
