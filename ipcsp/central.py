import math
import numpy as np
import matplotlib.pyplot as plt
import time
from ase import Atoms
from ase.calculators.kim import KIM
from ase.visualize import view

import calculator as cc
import open_grid as og
import randomizer as rd

angstrom_symbol = "\u212B"
sample_size = 10000

settings = {
    "Pt_1": {
        "State": True, "Symbol": "Pt_4", "grid": 4, "cell_size": 3.92,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "view": True, "cell_vari": False, "vari_range": 0.2
    },
    "Pt3Cu1_1": {
        "State": False, "Symbol": "Pt_3+Cu_1", "grid": 2, "cell_size": 3.92,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "view": True, "cell_vari": True, "vari_range": 0.2
    },
    "SrTiO3_1": {
        "State": False, "Symbol": "Sr_1+Ti_1+O_3", "grid": 2, "cell_size": 5,
        "calc": None, "cell_vari": False
    }
}

print(settings.keys())

if __name__ == "__main__":
    for key in settings.keys():
        if settings[key]["State"] and not settings[key]["cell_vari"]:
            start = time.time()

            # open_grid.py
            print("=" * 20 + f" Calculating {key} " + "=" * 20)
            chem, ion, pos_data = og.chem_compo(settings[key]["grid"], settings[key]["Symbol"])
            print(f"The Chemical Composition is {chem}")
            print(f"The Ion count is {ion}")
            print(f"Open_grid took time {time.time() - start} seconds")

            # randomizer.py
            all_config = list(rd.gen_config(chem, ion, settings[key]["grid"], sample_size=sample_size))
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
            i *= math.factorial((settings[key]["grid"] ** 3) - j)
            actual_cases = math.factorial(settings[key]["grid"] ** 3) // i  # 정수 나눗셈 (//)
            print(f"The actual possible number is {actual_cases}")
            print(f"Open ~ Actual cases took time {time.time() - start} seconds")

            # calculator.py
            best_e, best_c, best_s = cc.calcul(all_config, settings[key]["cell_size"], settings[key]["calc"], pos_data)
            print(f"The best energy is {best_e} eeVee")
            if settings[key]["view"]:
                atoms = Atoms(best_s, positions=best_c, cell=[settings[key]["cell_size"]] * 3, pbc=[True, True, True])
                view(atoms)

            end = time.time()
            print(f"The Calculation time is {end - start} seconds")

        # Cell Variation Mode
        if settings[key]["State"] and settings[key]["cell_vari"]:
            start = time.time()

            print("=" * 20 + f" Calculating {key}, cell size variation mode!!! " + "=" * 20)
            chem, ion, pos_data = og.chem_compo(settings[key]["grid"], settings[key]["Symbol"])
            print(f"The Chemical Composition is {chem}")
            print(f"The Ion count is {ion}")
            print(f"Open_grid took time {time.time() - start} seconds")

            # randomizer.py
            all_config = list(rd.gen_config(chem, ion, settings[key]["grid"], sample_size=sample_size))
            print(f"Sample Configuration of {key}")
            print(f"Total sample : {len(all_config)}")

            # actual_cases calculation
            i, j = 1, 0
            for indiv_ion in ion:
                i *= math.factorial(indiv_ion)
                j += indiv_ion
            i *= math.factorial((settings[key]["grid"] ** 3) - j)
            actual_cases = math.factorial(settings[key]["grid"] ** 3) // i  # 정수 나눗셈 (//)
            print(f"The actual possible number is {actual_cases}")
            print(f"Open ~ actual cases took time {time.time() - start} seconds")

            # calculator.py and cell variation generation
            cell_range = list(np.linspace(
                settings[key]["cell_size"] * (1 - settings[key]["vari_range"]),
                settings[key]["cell_size"] * (1 + settings[key]["vari_range"]),
                20
            ))

            best_best_e = []
            best_best_c = []
            best_best_s = []
            print(cell_range)
            
            # Cell size variation loop
            for size in cell_range:
                best_e, best_c, best_s = cc.calcul(all_config, size, settings[key]["calc"], pos_data)
                best_best_e.append(best_e)
                best_best_c.append(best_c)
                best_best_s.append(best_s)

            # minimum energy
            min_energy = min(best_best_e)
            min_index = best_best_e.index(min_energy)
            best_c_r = best_best_c[min_index]
            best_s_r = best_best_s[min_index]
            best_size = cell_range[min_index]

            print(f"The minimum energy in the cell size range between {cell_range[0]} and {cell_range[-1]} is {min_energy:.3f} eeVee, and the cell size is {best_size:.3f} {angstrom_symbol}")

            # 그래프 출력
            plt.title("The cell size variation and its energy")
            plt.xlabel(f"Cell size, unit: {angstrom_symbol}")
            plt.ylabel("energy, unit : eV")
            plt.plot(cell_range, best_best_e, color="hotpink")
            plt.grid(True)
            plt.show()

            # 3D 시각화
            if settings[key]["view"]:
                atoms = Atoms(best_s_r, positions=best_c_r, cell=[best_size] * 3, pbc=[True, True, True])
                view(atoms)

            end = time.time()
            print(f"The Calculation time is {end - start} seconds")

