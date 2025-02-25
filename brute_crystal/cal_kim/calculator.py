import cal_kim.open_grid
import cal_kim.randomizer
from ase import Atoms
from ase.calculators.kim import KIM
from ase.visualize import view

def calcul(configs, cell_size, calc, pos_data, pbc=[True, True, True], ortho = False):

    calc_model = KIM(calc)
    best_energy = float("inf")  # Initialize with a large number
    best_config = None
    best_symbols = None

    for config in configs:
        symbols_list = []
        positions = []

        if ortho:
            for i in range(len(config)):
                if config[i] != "_":
                    symbols_list.append(config[i])
                    positions.append([x * cell for x, cell in zip(pos_data[i], cell_size)])

        else:
            for i in range(len(config)):
                if config[i] != "_":
                    symbols_list.append(config[i])  # Collect symbols
                    positions.append([x * cell_size for x in pos_data[i]])  # Scale positions

        # Convert symbols list to string format for ASE
        symbols_string = "".join(symbols_list)  # Convert ['Na', 'Na', 'Cl', 'Cl'] → "NaNaClCl"

        # Create ASE Atoms object
        if ortho:
            atoms = Atoms(symbols_string, positions=positions, cell=cell_size, pbc=pbc)
        else:
            atoms = Atoms(symbols_string, positions=positions, cell=[cell_size]*3, pbc=pbc)
        atoms.calc = calc_model

        # Compute potential energy
        tmp_energy = atoms.get_potential_energy()
        
        # Update best configuration
        if tmp_energy < best_energy:
            best_energy = tmp_energy
            best_config = positions[:]
            best_symbols = symbols_string  # Store best configuration as string

    return best_energy, best_config, best_symbols




