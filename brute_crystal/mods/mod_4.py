import cal_int.open_matrix as om
import cal_int.pot_cal.elec_pot as pcep
import cal_int.pot_cal.ewal_pot as pcew
from cal_int.pot_cal.grid_gen import Phase
import cal_int.integer_program as ip
from ase import Atoms ; import time ; import os ; from ase.io import write ; import numpy as np


def mod_cell_int_periodic(key, settings):
    
    ortho = isinstance(settings[key]["grid"], list)
    start = time.time()

    # open_matrix.py | open_grid reference
    chem, ion_count, pos_data1 = om.og_cc(settings[key]["grid"], settings[key]["Symbol"], ortho = ortho)
    symbol = ""
    for element in chem:
        symbol += element
    phase = Phase(symbol)

    # open_matrix | cell matrix generator
    ewal_mat = pcew.get_Ewald(settings[key]["grid"], settings[key]["cell_size"], ortho = ortho)  ;  print(ewal_mat)
    buck_info = [settings[key]["grid"], settings[key]["cell_size"], phase]

    print(f"open_grid ~ generating matrix took {time.time() - start} seconds.")
    
    # integer_program.py
    position, op_energy = ip.integer(chem, ion_count, settings[key]["State"], ewal_mat, buck_mat = None, info = buck_info)
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
        atoms = Atoms(symbol, positions=pos_num, cell= settings[key]["cell_size"], pbc=[True, True, True])
        output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)

    else:
        pos_num = np.array(pos_num) ; cell_size = settings[key]["cell_size"] ; pos_num *= cell_size
        atoms = Atoms(symbol, positions=pos_num, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])
        output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")






    

