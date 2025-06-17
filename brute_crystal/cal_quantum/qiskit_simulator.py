from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit_algorithms import QAOA
from qiskit_algorithms.utils import algorithm_globals
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.optimizers import L_BFGS_B
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Sampler as AerSampler
import os
import numpy as np

def simulator_init(dist, chem, ion_count, charge):

    algorithm_globals.random_seed = 42
    grid_size = dist.shape[0]
    print(f"The length is {grid_size}, this is debug line")
    W = dist
    T = len(chem)

    # Define QUBO problem, variables are now just x{i}, but it will be separated into ion in proceessing part.
    qp = QuadraticProgram()

    for i in chem:
        for j in range(grid_size):
            qp.binary_var(name=f"{i}_{j}")

    # Objective
    linear = {}
    quadratic = {}

    for i1 in range(grid_size):

        for j1 in range(T):
            quadratic[(f"{chem[j1]}_{i1}", f"{chem[j1]}_{i1}")] = W[i1, i1] * charge[chem[j1]] ** 2

        for i2 in range(i1 + 1, grid_size):

            for j1 in range(T):
                quadratic[(f"{chem[j1]}_{i1}", f"{chem[j1]}_{i2}")] = 2 * W[i1, i2] * charge[chem[j1]] ** 2
        
                for j2 in range(j1 + 1, T):
                    quadratic[(f"{chem[j1]}_{i1}", f"{chem[j2]}_{i2}")] = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                    quadratic[(f"{chem[j2]}_{i1}", f"{chem[j1]}_{i2}")] = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]

    qp.minimize(linear=linear, quadratic=quadratic)

    # Constraints
    # First Constraints, there can be only one atom per position -> useless i think, but left as legacy code!!!
    for i in range(grid_size):
        one_atom_one_pos = {}
        for j in range(T):
            var_name = f"{chem[j]}_{i}"
            one_atom_one_pos[var_name] = 1
        # qp.linear_constraint(one_atom_one_pos, sense="==", rhs=1, name=f"uniqueness_{i}")

    # Second Constraint: each atom type i must appear ion_count[i] times

    for i in range(T):
        satisfy_ion_count = {}
        for j in range(grid_size):
            var_name = f"{chem[i]}_{j}"
            satisfy_ion_count[var_name] = 1
        qp.linear_constraint(satisfy_ion_count, sense="==", rhs=ion_count[i], name=f"satisfy_ion_count_{i}")
    
    constant = 0

    return qp, constant

def simulator_init_cons(dist, chem, ion_count, charge, cons):
    algorithm_globals.random_seed = 42
    grid_size = dist.shape[0]
    T = len(chem)
    W = dist

    qp = QuadraticProgram()

    # Helper to check if variable is allowed
    def is_var_allowed(chem_idx, pos):
        return pos in cons[chem_idx]

    def is_fixed(i, cons):
        return "Fixed" in cons[i]

    # Keep track of defined variables
    defined_vars = set()

    # Define binary variables
    for i, elem in enumerate(chem):
        if not is_fixed(i, cons):
            for j in range(grid_size):
                if is_var_allowed(i, j):
                    varname = f"{elem}_{j}"
                    print(varname)
                    qp.binary_var(name=varname)
                    defined_vars.add(varname)

    # Objective
    linear = {}
    quadratic = {}
    constant = 0

    for i1 in range(grid_size):

        for j1 in range(T):

            if not is_var_allowed(j1, i1):
                continue

            if is_fixed(j1, cons):
                continue

            var1 = f"{chem[j1]}_{i1}"
            q = W[i1, i1] * charge[chem[j1]] ** 2
            quadratic[(var1, var1)] = q


        for i2 in range(i1 + 1, grid_size):

            for j1 in range(T):

                if not is_var_allowed(j1, i1):
                    continue
                if not is_var_allowed(j1, i2):
                    continue

                var1 = f"{chem[j1]}_{i1}"
                var2 = f"{chem[j1]}_{i2}"

                if is_fixed(j1, cons):
                    pass
                
                if var1 in defined_vars and var2 in defined_vars:
                    q_same = 2 * W[i1, i2] * charge[chem[j1]] ** 2
                    quadratic[(var1, var2)] = q_same

                for j2 in range(j1 + 1, T):

                    if not is_var_allowed(j2, i2):
                        continue

                    var2 = f"{chem[j2]}_{i2}"
                    var3 = f"{chem[j2]}_{i1}"

                    if is_fixed(j1, cons) and is_fixed(j2, cons):
                        continue

                    elif not is_fixed(j1, cons) and is_fixed(j2, cons):
                        q_cross1 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        linear[(var1)] = q_cross1
                        q_cross2 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        linear[(f"{chem[j1]}_{i2}")] = q_cross2
                        continue

                    elif is_fixed(j1, cons) and not is_fixed(j2, cons):
                        q_cross1 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        linear[(var2)] = q_cross1
                        q_cross2 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        linear[(var3)] = q_cross2
                        continue

                    if var1 in defined_vars and var2 in defined_vars:
                        q_cross = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        quadratic[(var1, var2)] = q_cross
                    if var3 in defined_vars and f"{chem[j1]}_{i2}" in defined_vars:
                        q_cross = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        quadratic[(var3, f"{chem[j1]}_{i2}")] = q_cross

    qp.minimize(linear=linear, quadratic=quadratic)

    # Constraint 2: Each atom type appears ion_count[i] times
    for j in range(1):
        constraint = {}
        for i in range(grid_size):
            if is_var_allowed(j, i):
                var = f"{chem[j]}_{i}"
                if var in defined_vars:
                    constraint[var] = 1
        qp.linear_constraint(constraint, sense="==", rhs=ion_count[j], name=f"satisfy_ion_count_{j}")

    return qp, constant


def simulator_result(qp, constant):

    # Convert to QUBO
    converter = QuadraticProgramToQubo()
    qubo = converter.convert(qp)

    # AerSimulator backend with parallelism
    print("DEBUG:", AerSampler.__module__)
    num_cores = os.cpu_count()
    sampler = AerSampler(run_options={"shots" : 1000000, "max_parallel_shots" : 10000, "max_parallel_threads" : num_cores - 1})

    # QAOA setup
    optimizer = COBYLA()
    qaoa = QAOA(optimizer = optimizer, sampler = sampler, reps = 5)

    # Solve
    solver = MinimumEigenOptimizer(qaoa)
    result = solver.solve(qubo)
    result_sample = result.samples
    obj_value = result.fval + constant

    # Print result
    print("\n✅  QAOA Result (Sampler-based)")
    print(f"Objective value: {result.fval}")
    print(f"Binary solution: {result.x}")
    print(f"Objective Function Constant : {constant}")
    print(f"Status: {result.status.name}")
    print(f"Objective value + Constant : {result.fval + constant}")

    return result, obj_value

def decode_answer(result, chem, dist):
    position = []
    result = result.x
    n = len(dist)

    for idx, i in enumerate(chem):
        for j in range(idx * n, n + idx * n):
            if result[j]:
                position.append(f"{i}_{j - idx * n}")

    return [position]

def decode_answer_cons(result, qp, cons, chem):
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
    binary = result.x  # List of 0/1 results, aligned with qp.variables
    variables = qp.variables

    for var, bit in zip(variables, binary):
        if bit:
            position.append(var.name)

    if fix:
        for idx, pos in enumerate(cons):
            if isinstance(pos, list) and "Fixed" in pos:
                for j in pos[:-1]:  # all except "Fixed"
                    position.append(chem[idx] + "_" + str(j))
            else:
                continue

    return [position]

