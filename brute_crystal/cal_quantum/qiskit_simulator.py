from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit_algorithms import QAOA
from qiskit_algorithms.utils import algorithm_globals
from qiskit.primitives import Sampler  # <-- Use this instead of StatevectorSampler
from qiskit_algorithms.optimizers import COBYLA
import numpy as np

def simulator_init(total_array, type_count, grid_size, ion_count):
    algorithm_globals.random_seed = 42
    n = total_array.shape[0]
    print(f"The length is {n}, this is debug line")
    W = total_array

    # Define QUBO problem, variables are now just x{i}, but it will be separated into ion in proceessing part.
    qp = QuadraticProgram()
    for i in range(n):
        qp.binary_var(name=f"x{i}")

    # Objective
    linear = {}
    quadratic = {}
    for i in range(n):
        if W[i, i] != 0:
            linear[f"x{i}"] = W[i, i]
        for j in range(i + 1, n):
            if W[i, j] != 0:
                quadratic[(f"x{i}", f"x{j}")] = W[i, j]
    
    qp.minimize(linear=linear, quadratic=quadratic)

    # Constraints
    # First Constraints, there can be only one atom per position
    for i in range(grid_size):
        one_atom_one_pos = {}
        for j in range(type_count):
            var_name = f"x{i + j * grid_size}"
            one_atom_one_pos[var_name] = 1
        qp.linear_constraint(one_atom_one_pos, sense="==", rhs=1, name=f"uniqueness_{i}")

    # Second Constraint: each atom type i must appear ion_count[i] times
    for i in range(type_count):
        satisfy_ion_count = {}
        for j in range(grid_size):
            var_name = f"x{i * grid_size + j}"
            satisfy_ion_count[var_name] = 1
        qp.linear_constraint(satisfy_ion_count, sense="==", rhs=ion_count[i], name=f"satisfy_ion_count_{i}")

    return qp

def simulator_result(qp):
    
    # Convert to QUBO
    converter = QuadraticProgramToQubo()
    qubo = converter.convert(qp)

    # QAOA setup
    sampler = Sampler()
    optimizer = COBYLA()
    qaoa = QAOA(optimizer=optimizer, sampler=sampler, reps=1)

    # Solve
    solver = MinimumEigenOptimizer(qaoa)
    result = solver.solve(qubo)
    result_sample = result.samples

    # Print result
    print("\n✅  QAOA Result (Sampler-based)")
    print(f"Objective value: {result.fval}")
    print(f"Binary solution: {result.x}")
    print(f"Status: {result.status.name}")

    return result

def simulator_process(result, grid_size, sort_by="value"):
    result_sample = result.samples
    value_data = {}

    for sample in result_sample:
        x = sample.x
        prob = sample.probability
        fval = sample.fval

        # Convert bitstring to integer
        bitstring = "".join(str(int(b)) for b in x)
        value = int(bitstring, 2)

        if value not in value_data:
            value_data[value] = {
                "probability": 0.0,
                "fvals": []
            }

        value_data[value]["probability"] += prob
        value_data[value]["fvals"].append(fval)

    # 정리된 결과: (value, prob, avg_fval)
    entries = []
    for value, data in value_data.items():
        prob = data["probability"]
        avg_fval = sum(data["fvals"]) / len(data["fvals"])
        entries.append((value, prob, avg_fval))

    # 정렬 방식 선택
    if sort_by == "value":
        sorted_items = sorted(entries, key=lambda x: x[0])
    elif sort_by == "probability":
        sorted_items = sorted(entries, key=lambda x: x[1], reverse=True)
    else:
        raise ValueError("sort_by must be 'value' or 'probability'")

    return sorted_items

def decode_answer(result, type_count, grid_size, type_names=None):
    answer = result.x
    output = []

    if type_names is None:
        type_names = [f"type{i}" for i in range(type_count)]

    for idx, bit in enumerate(answer):
        if bit == 1:
            type_idx = idx // grid_size
            pos_idx = idx % grid_size
            label = f"{type_names[type_idx]}_{pos_idx}"
            output.append(label)

    return output

