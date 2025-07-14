import cal_int.open_matrix as om
import cal_int.pot_cal.ewal_pot2 as pcew2
import cal_pyscf.obj_generator as pyog
from data.charge import charge

from ase import Atoms
import time
import os
from ase.io import write
import numpy as np
import matplotlib.pyplot as plt

def pyscf_simulator(key, settings):

    ortho = isinstance(settings[key]["grid"], list)
    angle = settings[key]["angle"]
    angle_if = isinstance(settings[key]["angle"], list)
    start = time.time()

    # open_matrix.py | open_grid reference | in open_matrix.py, we use chem_compo function as om.og_cc
    chem, ion_count, pos_data1 = om.og_cc(settings[key]["grid"], settings[key]["Symbol"], ortho = ortho)  ;  #print(pos_data1)
    symbol = ""
    for element in chem:
        symbol += element
    grid_size = len(pos_data1)

    # open_matrix | cell matrix generator
    pos_data2 = om.mat_gen(pos_data1, settings[key]["cell_size"], ortho = ortho, angle = angle)  ;  #print(pos_data2)
  
    # Ewald Part
    ewal_mat = pcew2.Ewald(pos_data2, settings[key]["cell_size"], settings[key]["grid"], settings[key]["angle"])
    ewal_mat.vec_gen()  ;  ewal_mat.ewal_mat_gen()  ;  print(ewal_mat.dist)  ;  #print(ewal_mat.position)

    print(f"open_grid ~ generating matrix took {time.time() - start} seconds.")

    # Objective Function Generating Part
    if settings[key]["cons"]:
        qubo, norb, nelec = pyog.simulator_init_cons(ewal_mat.dist, chem, ion_count, charge, settings[key]["cons"])
    else:
        pass

    # Pyscf part
    energy, vector = pyog.solve_with_pyscf(qubo, norb, nelec)
    ge = energy[0]
    print(f"Qubo Energy : {ge:.3f} eevee")
    print(energy)
    print("Number of basis set : ", len(vector[0]))
    print(f"open_grid ~ pyscf matrix diagonalization took {time.time() - start} seconds.")

