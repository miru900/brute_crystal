import open_grid
import randomizer
from ase import Atoms
from ase.calculators.kim import KIM
from ase.visualize import view

# Get all possible configurations
configs = randomizer.all_config
posdata = open_grid.pos_data

print("The default is Cubic, I will update in later")
cell_size = float(input("Write Cell size in Angstrom : "))  # Ensure float conversion

# Use an appropriate potential
calc = KIM("EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000")

def calcul(configs, cell_size, calc, pbc=[True, True, True]):
    best_energy = float("inf")  # Initialize with a large number
    best_config = None
    best_symbols = None

    for config in configs:
        symbols_list = []
        positions = []

        # Extract symbols and positions
        for i in range(len(config)):
            if config[i] != "_":
                symbols_list.append(config[i])  # Collect symbols
                positions.append([x * cell_size for x in posdata[i]])  # Scale positions

        # Convert symbols list to string format for ASE
        symbols_string = "".join(symbols_list)  # Convert ['Na', 'Na', 'Cl', 'Cl'] → "NaNaClCl"

        # Create ASE Atoms object
        atoms = Atoms(symbols_string, positions=positions, cell=[cell_size]*3, pbc=pbc)
        atoms.calc = calc

        # Compute potential energy
        tmp_energy = atoms.get_potential_energy()
        
        # Update best configuration
        if tmp_energy < best_energy:
            best_energy = tmp_energy
            best_config = positions[:]
            best_symbols = symbols_string  # Store best configuration as string

    return best_energy, best_config, best_symbols

# Run calculation and find the best configuration
best_e, best_p, best_s = calcul(configs, cell_size, calc, pbc=[True, True, True])

# Print results
print(f"The best energy is {best_e} eeVee")
atoms = Atoms(best_s, positions=best_p, cell=[cell_size]*3, pbc=[True, True, True])
view(atoms)

