import cal_int.open_matrix as om
import cal_int.pot_cal.ewal_pot2 as pcew2
import cal_quantum.array_generator as qcag
import cal_quantum.qiskit_simulator as qcqs
from data.charge import charge

from ase import Atoms
import time
import os
from ase.io import write
import numpy as np
import matplotlib.pyplot as plt


def qiskit_simulator(key, settings):

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

    # Qiskit Part
    total_array, type_count, grid_size = qcag.total_array(ewal_mat.dist, chem, charge, grid_size) # Total array generator
    qp = qcqs.simulator_init(total_array, type_count, grid_size, ion_count) # initialization part
    # special constraint part
    result = qcqs.simulator_result(qp) # simulating part
    print(f"Simulating Quantum Computer Finished, time taken from starting is {time.time() - start} seconds.")

    # processing part 1, QAOA result plot
    sorted_items1 = qcqs.simulator_process(result, grid_size)
    value = [sorted_items1[i][0] for i in range(len(sorted_items1))]
    pro = [sorted_items1[i][1] for i in range(len(sorted_items1))]

    sorted_items2 = qcqs.simulator_process(result, grid_size, sort_by = "probability")
    value2 = [sorted_items2[i][0] for i in range(len(sorted_items2))]
    pro2 = [sorted_items2[i][1] for i in range(len(sorted_items2))]
    fval2 = [sorted_items2[i][2] for i in range(len(sorted_items2))]
   
    # processing part 2, bit to answers
    position = qcqs.decode_answer(result, type_count, grid_size, type_names=chem)
    position = [position]
    solcount = len(position)

    # processing part 3, ASE part
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
    
    print(f"The calculated energy is {result.fval:.03f} eevee")

    if angle_if:

        for i in range(solcount):
            c_pos = []
            current_grid_num = all_gridnum[i]

            for j in current_grid_num:
                c_pos.append(pos_data2[j])

            c_pos = np.array(c_pos)
            c_pos *= 1e10  ;  # print(c_pos)

            atoms = Atoms(all_symbol[i], positions=c_pos, cell = settings[key]["cell_size"] + settings[key]["angle"], pbc=[True, True, True])
            output_dir = f"../results/{settings[key]['Symbol']}_quantum" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{i}.cif")
            write(file_path, atoms)
            write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')
            plt.bar(value, pro, color = "hotpink")
            plt.savefig(f"{output_dir}/distribution.png")

    elif ortho and not angle_if:

        cell_size = settings[key]["cell_size"]
        for i in range(solcount):
            current_processing_pos = all_positions[i]

            for j, pos in enumerate(current_processing_pos):
                pos = np.multiply(pos, cell_size)
                current_processing_pos[j] = pos

            c_pos = np.array(current_processing_pos)
            atoms = Atoms(all_symbol[i], positions=c_pos, cell= settings[key]["cell_size"], pbc=[True, True, True])
            output_dir = f"../results/{settings[key]['Symbol']}_quantum" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{i}.cif")
            write(file_path, atoms)
            write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')
            plt.bar(value, pro, color = "hotpink")
            plt.savefig(f"{output_dir}/distribution.png")

    else:
        pos_num = np.array(pos_num) ; cell_size = settings[key]["cell_size"] ; pos_num *= cell_size
        atoms = Atoms(symbol, positions=pos_num, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])
        output_dir = f"../results/{settings[key]['Symbol']}_quantum" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")
        write(file_path, atoms)
        write(file_path.replace(".cif", ".vasp"), atoms, format='vasp')
        plt.bar(value, pro, color = "hotpink")
        plt.savefig(f"{output_dir}/distribution.png")
        result = open(f"{output_dir}/result.txt", "w")
        for i in range(len(value)):
            result.write(f"{value2[i]}, {pro2[i]}, {fval2[i]}\n")
        result.close()

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")

