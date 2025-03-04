import gurobipy as gb  ;  import numpy as np  ;  from data.charge import charge
from scipy.constants import elementary_charge as elec_c  ;  from scipy.constants import electron_volt as eevee
from cal_int.pot_cal.elec_pot import energy_elec_addup as elad
from cal_int.pot_cal.ewal_pot import energy_ewal_addup as ewad
from cal_int.pot_cal.buck_pot import energy_buck_addup as buad
from cal_int.pot_cal.buck_pot2 import energy_buck_addup as buad2
import cal_int.constraints as cons

def integer(chem, ion_count, ewal_mat = None, buck_mat = None, info = None):
    # example " chem = ["Sr", "Ti", "O"], ion_count = [1, 1, 3], grid = grid from open_grid, matrix = matrix from potential calculator
    # orb_key = orb_pos = [0, 1, 2, 3, 4, 5, 6, 7], orb_size = [1, 1, 1, 1, 1, 1, 1, 1], O = N (orbit length = grid square)

    T = len(chem)  ;  O = len(ewal_mat)
    orb_key = [i for i in range(O)]  ;  orb_size = [1] * len(ewal_mat)
    
    # Model & Variables Setting
    m = gb.Model() ; Vars = [[] for i in range(T)]
    
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
    
    # add_extra_Cons(m, Vars)

    print("Variables and Constraints are generated")

    # Making Objective Function, energy, using potential matrix generated from pot_cal
    energy1 = gb.QuadExpr()  ;  energy2 = gb.QuadExpr()

    # Potential Addup, remove hashtag u want to calculate
    ewad(O, T, Vars, orb_key, ewal_mat, charge, energy1, chem)
    # buad(O, Vars, orb_key, energy2, chem, info[0], info[1], info[2])
    buad2(O, Vars, orb_key, energy2, chem, buck_mat, info)

    # elad(O, T, Vars, orb_key, matrix, charge, energy1, chem)
    # buad(O, Vars, orb_key, energy2, chem, info[0], info[1], info[2])

    print("The objetive function was generated")

    # Special Constraints for Special Situation
    # cons.NCM_cons(m, Vars)

    # Setting Time limit, Solving QUBO problem, and return result
    m.setParam("TimeLimit", 600)
    m.setObjective(energy1 + energy2, gb.GRB.MINIMIZE)
    m.optimize()
    min_sym_pos = []

    if m.status == gb.GRB.OPTIMAL or m.status == gb.GRB.TIME_LIMIT or gb.GRB.INTERRUPTED:
        for v in m.getVars():
            if v.Xn == 1:
                min_sym_pos.append(v.varName)

    return min_sym_pos, m.objVal


