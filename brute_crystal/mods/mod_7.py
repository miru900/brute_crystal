import cal_kim.open_grid as og
import cal_int.open_matrix as om
import cal_int.pot_cal.ewal_pot2 as pcew2
import cal_int.pot_cal.ase_kim_pot as pcak
import cal_int.integer_program as ip
from ase import Atoms ; import time ; import os ; from ase.io import write ; import numpy as np

def mod_cell_int_ase_kim(key, settings):

    ortho = isinstance(settings[key]["grid"], list)
    start = time.time()

    # open_matrix.py | open_grid reference
    chem, ion_count, pos_data1 = om.og_cc(settings[key]["grid"], settings[key]["Symbol"], ortho = ortho)  ;  #print(pos_data1)
    symbol = ""

    for element in chem:
        symbol += element

    # open_matrix | cell matrix generator
    pos_data2 = om.mat_gen(pos_data1, settings[key]["cell_size"], ortho = ortho)  ;  #print(pos_data2)
    pos_data3 = og.grid_gen(pos_data1, settings[key]["cell_size"], ortho = ortho)

    # Ewald Part
    ewal_mat = pcew2.Ewald(pos_data2, settings[key]["cell_size"], settings[key]["grid"])
    ewal_mat.vec_gen()  ;  ewal_mat.ewal_mat_gen()  ;  print(ewal_mat.dist)  ;  #print(ewal_mat.position)

    # Ase Part
    kim_mat = pcak.ASE_KIM(pos_data3, settings[key]["cell_size"], settings[key]["grid"], symbol, calc = settings[key]["calc"])  ;  kim_mat.ase_mat_gen()
    ion_pair = kim_mat.info["ion_pair"]  ;   print(kim_mat.dist)

    print(f"open_grid ~ generating matrix took {time.time() - start} seconds.")

    # integer_program.py
    position, op_energy = ip.integer(chem, ion_count, settings[key]["State"] ,ewal_mat.dist, buck_mat = None, kim_mat = kim_mat.dist, info = ion_pair)
    print(f"open_grid ~ integer programming took {time.time() -start} seconds.")

    # Processing answer | splitting answer (Example. Na_1 -> Na, 1)
    print(f"The Calculated objective energy is {op_energy:.3f} eeVee")
    symbol = ""
    pos_num = []
    for pos in position :
        print(f"One of positions {pos}")
        a, b = pos.split("_")
        symbol += a ; b = int(b) ; c = np.array(pos_data1[b])
        pos_num.append(c)

    # processing answer | ASE part
    if ortho:
        cell_size = settings[key]["cell_size"]
        for i, pos in enumerate(pos_num):
            pos = np.multiply(pos, cell_size)
            pos_num[i] = pos
        pos_num = np.array(pos_num)
        atoms = Atoms(symbol, positions = pos_num, cell = settings[key]["cell_size"], pbc=[True, True, True])
        output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)

    else:
        pos_num = np.array(pos_num) ; cell_size = settings[key]["cell_size"] ; pos_num *= cell_size
        atoms = Atoms(symbol, positions = pos_num, cell = [settings[key]["cell_size"]] * 3, pbc=[True, True, True])
        output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")

