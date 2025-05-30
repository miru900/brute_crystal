import math
import numpy as np
import matplotlib.pyplot as plt
import time  ;  import os
from ase import Atoms  ;  from ase.calculators.kim import KIM  ;  from ase.visualize import view  ;  from ase.io import write

import cal_kim.calculator as cc
import cal_kim.open_grid as og
import cal_kim.randomizer as rd

angstrom_symbol = "\u212B"

def mod_cell_one(key, settings, sample_size):

    ortho = isinstance(settings[key]["grid"], list)
    start = time.time()
    # open_grid.py
    print("=" * 30 + f" Calculating {key} "+ " | " + "Mod : 1, Calculating only one cel1 size!!! "+ "=" * 30)
    chem, ion, pos_data = og.chem_compo(settings[key]["grid"], settings[key]["Symbol"], ortho = ortho)
    print(f"The Chemical Composition is {chem}")
    print(f"The Ion count is {ion}")
    print(f"Open_grid took time {time.time() - start} seconds")

    # randomizer.py
    all_config = list(rd.gen_config(chem, ion, settings[key]["grid"], sample_size=sample_size, ortho = ortho))
    config_num = len(all_config)
    print(f"Sample Configuration of {key}")
    for config in all_config[:3]:
        print(config)
    print(f"Total Configuration : {config_num}")
    print(f"Open ~ Randomizer took time {time.time() - start} seconds")

    # actual_cases calculation
    i, j = 1, 0
    for indiv_ion in ion:
        i *= math.factorial(indiv_ion)
        j += indiv_ion
    if ortho:
        pass
    else:
        i *= math.factorial((settings[key]["grid"] ** 3) - j)
        actual_cases = math.factorial(settings[key]["grid"] ** 3) // i
        print(f"The actual possible number is {actual_cases}")
        print(f"Open ~ Actual cases took time {time.time() - start} seconds")

    # calculator.py
    best_e, best_c, best_s = cc.calcul(all_config, settings[key]["cell_size"], settings[key]["calc"], pos_data, ortho = ortho)
    print(f"The best energy is {best_e} eeVee")
    
    if settings[key]["view"]:
        if ortho:
            atoms = Atoms(best_s, positions=best_c, cell= settings[key]["cell_size"] , pbc=[True, True, True])  ;  view(atoms)
            output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")  ;  write(file_path, atoms)
        else:
            atoms = Atoms(best_s, positions=best_c, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])  ;  view(atoms)
            output_dir = "../results" ; os.makedirs(output_dir, exist_ok = True) ; file_path = os.path.join(output_dir, f"{settings[key]['Symbol']}.cif")  ;  write(file_path, atoms)

    end = time.time()
    print(f"The Calculation time is {end - start} seconds")

