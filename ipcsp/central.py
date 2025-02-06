import calculator as cc
import open_grid as og
import randomizer as rd
from ase import Atoms
from ase.visualize import view
import matplotlib.pyplot as plt
import numpy as np

angstrom_symbol = "\u212B"

settings = {
        "Pt_1" : {"State" : False, "Symbol" : "Pt_4", "grid" : 2, "cell_size" : 3.92, "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000", "view" : True, "cell_vari" : True, "vari_range" : 0.2},
        "Pt3Cu1_1" : {"State" : True, "Symbol" : "Pt_3+Cu_1", "grid" : 2, "cell_size" : 3.92, "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000", "view" : True, "cell_vari" : True, "vari_range" : 0.2},
        "SrTiO3_1" : {"State" : False, "Symbol" : "Sr_1+Ti_1+O_3", "grid" : 2, "cell_size" : 5, "calc" : None, "cell_vari" : False}}

print(settings.keys())

if __name__ == "__main__":
    for key in settings.keys():
        if settings[key]["State"] and not settings[key]["cell_vari"] :
            
            # open_grid.py
            print("=" * 20 + f" Calculating {key} " + "=" * 20)
            chem, ion, pos_data = og.chem_compo(settings[key]["grid"], settings[key]["Symbol"])
            print(f"The Chemical Composition is {chem}")
            print(f"The Ion count is {ion}")
            
            # randomizer.py
            all_config = rd.gen_config(chem, ion, settings[key]["grid"])
            config_num = len(all_config)
            print(f"Sample Configuration of {key}")
            for config in all_config[:3]:
                print(config)
            print(f"Total Configuration : {config_num}")

            # calculator.py
            best_e, best_c, best_s = cc.calcul(all_config, settings[key]["cell_size"], settings[key]["calc"], pos_data)
            print(f"The best energy is {best_e} eeVee")
            if settings[key]["view"]:
                atoms = Atoms(best_s, positions = best_c, cell = [settings[key]["cell_size"]] * 3, pbc = [True, True, True])
                view(atoms)
        
        # Cell Variation Mode

        if settings[key]["State"] and settings[key]["cell_vari"] :
            # open_grid.py
            print("=" * 20 + f" Calculating {key}, cell size variation mode!!! " + "=" * 20)
            chem, ion, pos_data = og.chem_compo(settings[key]["grid"], settings[key]["Symbol"])
            print(f"The Chemical Composition is {chem}")
            print(f"The Ion count is {ion}")

            # randomizer.py
            all_config = rd.gen_config(chem, ion, settings[key]["grid"])
            config_num = len(all_config)
            print(f"Sample Configuration of {key}")
            for config in all_config[:3]:
                print(config)
            print(f"Total Configuration : {config_num}")

            # calculator.py and cell variation generation
            cell_range = list(np.linspace (settings[key]["cell_size"] * (1 - settings[key]["vari_range"]), settings[key]["cell_size"] * (1 + settings[key]["vari_range"]), 20))
            best_best_e = []
            print(cell_range)
            for size in cell_range:
                best_e, best_c, best_s = cc.calcul(all_config, size, settings[key]["calc"], pos_data)
                best_best_e.append(best_e)

            min_energy = min(best_best_e)
            min_index = best_best_e.index(min_energy)
            
            print(f"The minimum energy in the cell size range between {cell_range[0]} and {cell_range[-1]} is {min_energy:.3f} eeVee, and the cell size is {cell_range[min_index]:.3f} {angstrom_symbol}")
            plt.title("The cell size variation and its energy")
            plt.xlabel("Cell size, unit: $\AA$")
            plt.ylabel("energy, unit : eV")
            plt.plot(cell_range, best_best_e, color = "hotpink")
            plt.grid(True)
            plt.show()
            if settings[key]["view"]:
                atoms = Atoms(best_s, positions = best_c, cell = [settings[key]["cell_size"]] * 3, pbc = [True, True, True])
                view(atoms)
