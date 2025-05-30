import cal_int.open_matrix as om
import cal_int.pot_cal.elec_pot as pcep
import cal_int.pot_cal.ewal_pot as pcew
import cal_int.pot_cal.ewal_pot2 as pcew2
import cal_int.pot_cal.buck_pot2 as pcbu2
from cal_int.pot_cal.grid_gen import Phase
import cal_int.integer_program as ip
from ase import Atoms ; import time ; import os ; from ase.io import write ; import numpy as np


def mod_cell_int_periodic2(key, settings):
    
    ortho = isinstance(settings[key]["grid"], list)
    angle = settings[key]["angle"]
    angle_if = isinstance(settings[key]["angle"], list)
    start = time.time()

    # open_matrix.py | open_grid reference | in open_matrix.py, we use chem_compo function as om.og_cc
    chem, ion_count, pos_data1 = om.og_cc(settings[key]["grid"], settings[key]["Symbol"], ortho = ortho)  ;  #print(pos_data1)
    symbol = ""
    for element in chem:
        symbol += element

    # open_matrix | cell matrix generator
    pos_data2 = om.mat_gen(pos_data1, settings[key]["cell_size"], ortho = ortho, angle = angle)  ;  #print(pos_data2)

    # Ewald Part
    ewal_mat = pcew2.Ewald(pos_data2, settings[key]["cell_size"], settings[key]["grid"], settings[key]["angle"])
    ewal_mat.vec_gen()  ;  ewal_mat.ewal_mat_gen()  ;  print(ewal_mat.dist)  ;  #print(ewal_mat.position)

    # Buckingham Part
    # buck_mat = pcbu2.Buckingham(pos_data2, settings[key]["cell_size"], settings[key]["grid"], symbol)  ;  buck_mat.operator()
    # ion_pair = buck_mat.info["ion_pair"]  ;   print(buck_mat.dist)

    print(f"open_grid ~ generating matrix took {time.time() - start} seconds.")
    
    # integer_program.py
    position, op_energy, solcount = ip.integer(chem, ion_count, settings[key]["State"], ewal_mat.dist, buck_mat = None, info = None)
    print(f"open_grid ~ integer programming took {time.time() -start} seconds.")
    
    # Processing answer | splitting answer (Example. Na_1 -> Na, 1)
    print(f"The Calculated lowest objective energy is {op_energy[0]:.3f} eeVee")
    
    # Appending positions list of all solutions
    all_positions = []
    all_symbol = []
    all_gridnum = []

    for i in range(solcount):
        symbol = ""
        pos_num = []
        g_num = []
        iter_position = position[i]
        for pos in iter_position:
            print(f"One of positions {pos}")
            a, b = pos.split("_")
            symbol += a ; b = int(b) ; c = np.array(pos_data1[b])
            pos_num.append(c)
            g_num.append(b)

        all_gridnum.append(g_num)
        all_positions.append(pos_num)
        all_symbol.append(symbol)
    
    # print(f"TESTING MOD_5 POSITION LIST : {all_positions}", len(all_positions))

    # processing answer | ASE part
    if angle_if:
        
        for i in range(solcount):
            c_pos = []
            current_grid_num = all_gridnum[i]
            
            for j in current_grid_num:
                c_pos.append(pos_data2[j])

            c_pos = np.array(c_pos)
            c_pos *= 1e10  ;  # print(c_pos)

            atoms = Atoms(all_symbol[i], positions=c_pos, cell = settings[key]["cell_size"] + settings[key]["angle"], pbc=[True, True, True])
            output_dir = f"../results/{settings[key]['Symbol']}" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{i}.cif")
            write(file_path, atoms)
            write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')

        with open(f"../results/{settings[key]['Symbol']}/energies.txt", "w") as enerugi:
            for k in range(solcount):
                enerugi.write(f"{k}, {op_energy[k]} \n")

    elif ortho and not angle_if:

        cell_size = settings[key]["cell_size"]
        for i in range(solcount):
            current_processing_pos = all_positions[i]

            for j, pos in enumerate(current_processing_pos):
                pos = np.multiply(pos, cell_size)
                current_processing_pos[j] = pos

            c_pos = np.array(current_processing_pos)
            atoms = Atoms(all_symbol[i], positions=c_pos, cell= settings[key]["cell_size"], pbc=[True, True, True])
            output_dir = f"../results/{settings[key]['Symbol']}" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{i}.cif")
            write(file_path, atoms)
            write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')

        with open(f"../results/{settings[key]['Symbol']}/energies.txt", "w") as enerugi:
            for k in range(solcount):
                enerugi.write(f"{k}, {op_energy[k]} \n")

    else:
        pos_num = np.array(pos_num) ; cell_size = settings[key]["cell_size"] ; pos_num *= cell_size
        atoms = Atoms(symbol, positions=pos_num, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])
        output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)
        write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")






    

