import open_grid
import randomizer
from ase import Atoms
from ase.calculators.kim.kim import KIM
from ase.visualize import view

# print(randomizer.all_config)
configs = randomizer.all_config
posdata = open_grid.pos_data

print("The default is Cubic, I will update in later")
cell_size = input("Write Cell size in ANgstrom : ")

calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")

def calcul(configs, cell_size, calc, pbc = [True, True, True]):
    best_energy = 10000000
    symbol = 0
    for config in configs:
        for i in len(config):
            symbol = ""
            pos_tem = []
            if not config[i] == "_":
                symbol += config[i]
                pos_tem.append(posdata[i] * cell_size)
        atoms = Atoms(symbol, positions = pos_tem, cell = [cell_size, cell_size, cell_size], pbc = pbc)
        atoms.calc = calc
        tmp_energy = atoms.get_potential_energy()
        if best_energy > tmp_energy:
            best_energy = tmp_energy
            best_config = pos_tem
            symbol = symbol

    return best_energy, best_config, symbol

best_e, best_p, symbol = calcul(configs, cell_size, calc = calc, pbc = [True, True, True])
print(f"The best energy is {best_e} eeVee")
atoms = Atoms(symbol, positions = best_p, cell = [cell_size, cell_size, cell_size], pbc = [True, True, True])
view(atoms)





