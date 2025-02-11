import cal_int.open_matrix as om
import cal_int.pot_cal.elec_pot as pcep
import cal_int.integer_program as ip
from ase import Atoms
import time
import os
from ase.io import write
import numpy as np

def mod_cell_integer(key, settings):
    
    start = time.time()

    # open_matrix.py
    chem, ion_count, pos_data1 = om.og_cc(settings[key]["grid"], settings[key]["Symbol"])
    # print(pos_data1)
    pos_data2 = om.mat_gen(pos_data1, settings[key]["cell_size"])
    # print(pos_data2)
    # print(chem, ion_count)
    print(f" Opening grid information took {time.time() - start} seconds.")


    # pc directory
    elec_mat = pcep.elec_pot(pos_data2, settings[key]["cell_size"] * 1e-10)
    # print(elec_mat)
    print(f"open_grid ~ generating matrix took {time.time() - start} seconds.")
    

    # integer_program.py
    position, op_energy = ip.integer(chem, ion_count, elec_mat)
    print(f"open_grid ~ integer programming took {time.time() -start} seconds.")
    

    # Processing answer
    print(f"The Calculated objective energy is {op_energy:.3f} eeVee")
    symbol = ""
    pos_num = []
    for pos in position :
        print(f"One of positions {pos}")
        a, b = pos.split("_")
        symbol += a ; b = int(b) ; c = np.array(pos_data1[b])
        pos_num.append(c)
    
    pos_num = np.array(pos_num) ; cell_size = settings[key]["cell_size"] ; pos_num *= cell_size

    
    print(pos_num, symbol)
    atoms = Atoms(symbol, positions=pos_num, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])
    output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
    write(file_path, atoms)

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")






    

