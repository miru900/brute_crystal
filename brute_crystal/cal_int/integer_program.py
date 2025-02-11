import gurobipy as gb
import numpy as np
from data.charge import charge
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee

def integer(chem, ion_count, matrix):
    # example " chem = ["Sr", "Ti", "O"], ion_count = [1, 1, 3], grid = grid from open_grid, matrix = matrix from potential calculator
    # orb_key = orb_pos = [0, 1, 2, 3, 4, 5, 6, 7], orb_size = [1, 1, 1, 1, 1, 1, 1, 1], O = N (orbit length = grid square)

    m = gb.Model()
    T = len(chem)
    O = len(matrix)
    orb_key = [i for i in range(O)]
    orb_size = [1] * len(matrix)
    
    # Model & Variables Setting
    m = gb.Model()
    Vars = [[] for i in range(T)]
    
    # Adding Variables and Constraints
    for i in range(O):
        tmp_var = []
        for j in range(T):
            Vars[j] += [m.addVar(vtype = gb.GRB.BINARY, name = chem[j] + "_" + str(orb_key[i]))]
            tmp_var += [Vars[j][-1]]
        m.addConstr(gb.LinExpr([1.0] * T, tmp_var) <=1 )

    for j in range(T):
        tmp = gb.LinExpr()
        for i in range(O):
            tmp.add(Vars[j][i], orb_size[i])
        m.addConstr(tmp == ion_count[j])
    
    print("Variables and Constraints are generated")

    # Making Objective Function, energy, using potential matrix generated from pot_cal
    # i = i1, j1 = j1, j1 = i2, k1 = k1, k2 = j2
    energy = gb.QuadExpr()

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
    
    print("The objetive function was generated")

    # Solving QUBO problem, and return result
    m.setObjective(energy, gb.GRB.MINIMIZE)
    m.optimize()
    min_sym_pos = []
    for v in m.getVars():
        if v.Xn == 1:
            min_sym_pos.append(v.varName)

    return min_sym_pos, m.objVal


