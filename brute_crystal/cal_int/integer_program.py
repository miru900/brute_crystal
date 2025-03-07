import gurobipy as gb  ;  import numpy as np  ;  from data.charge import charge
from scipy.constants import elementary_charge as elec_c  ;  from scipy.constants import electron_volt as eevee
from cal_int.pot_cal.elec_pot import energy_elec_addup as elad
from cal_int.pot_cal.ewal_pot import energy_ewal_addup as ewad
from cal_int.pot_cal.buck_pot import energy_buck_addup as buad
from cal_int.pot_cal.buck_pot2 import energy_buck_addup as buad2
from cal_int.pot_cal.ase_kim_pot import energy_kim_addup as kiad
import cal_int.constraints as cons

def integer(chem, ion_count, mod, ewal_mat = None, buck_mat = None, kim_mat = None, info = None):
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
    
    print("Variables and Constraints are generated")

    # Potential Addup, remove hashtag u want to calculate
    
    if mod == 4:
        energy1, energy2 = mod_4(m, O, T, Vars, orb_key, ewal_mat, buck_mat, charge, chem, info)
        e_sum = energy1 + energy2
    if mod == 5:
        energy1 = mod_5(m, O, T, Vars, orb_key, ewal_mat, charge, chem, info)
        e_sum = energy1
    if mod == 6:
        energy1, energy2 = mod_6(m, O, T, Vars, orb_key, ewal_mat, buck_mat, charge, chem, info)
        e_sum = energy1 + energy2
    if mod == 7:
        energy1, energy2 = mod_7(m, O, T, Vars, orb_key, ewal_mat, kim_mat, charge, chem, info)
        e_sum = energy1 + energy2

    print("The objetive function was generated")

    # Special Constraints for Special Situation
    # add_cons(m, Vars)

    # Setting Time limit, Solving QUBO problem, and return result
    m.setParam('NoRelHeurTime', 15)
    m.setParam("TimeLimit", 60)
    m.setObjective(e_sum, gb.GRB.MINIMIZE)
    m.optimize()
    min_sym_pos = []

    if m.status == gb.GRB.OPTIMAL or m.status == gb.GRB.TIME_LIMIT or gb.GRB.INTERRUPTED:
        for v in m.getVars():
            if v.Xn == 1:
                min_sym_pos.append(v.varName)

    return min_sym_pos, m.objVal



def mod_4(m, O, T, Vars, orb_key, ewal_mat, buck_mat, charge, chem, info):
    energy1 = gb.QuadExpr()
    ewad(O, T, Vars, orb_key, ewal_mat, charge, energy1, chem)
    energy2 = gb.QuadExpr()
    buad(O, Vars, orb_key, energy2, chem, info[0], info[1], info[2])
    
    return energy1, energy2

def mod_5(m, O, T, Vars, orb_key, ewal_mat, charge, chem, info):
    energy1 = gb.QuadExpr()
    ewad(O, T, Vars, orb_key, ewal_mat, charge, energy1, chem)
    
    return energy1

def mod_6(m, O, T, Vars, orb_key, ewal_mat, buck_mat, charge, chem, info):
    energy1 = gb.QuadExpr()
    ewad(O, T, Vars, orb_key, ewal_mat, charge, energy1, chem)
    energy2 = gb.QuadExpr()
    buad2(O, Vars, orb_key, energy2, chem, buck_mat, info)

    return energy1, energy2

def mod_7(m, O, T, Vars, orb_key, ewal_mat, kim_mat, charge, chem, info):
    energy1 = gb.QuadExpr()
    ewad(O, T, Vars, orb_key, ewal_mat, charge, energy1, chem)
    energy2 = gb.QuadExpr()
    kiad(O, Vars, orb_key, energy2, chem, kim_mat, info)

    return energy1, energy2


def add_cons(m, Vars):
    print("Now adding Constraints")
    cons.NCM_cons(m, Vars)
