from ase import Atoms
from ase.calculators.kim import KIM
from data.buck import buck  
from data.radii import radii
from cal_int.pot_cal.ewal_pot2 import Ewald
import numpy as np
import multiprocessing as mp

pbc_t = [True] * 3  

class ASE_KIM(Ewald):

    def __init__(self, position, cell_size, grid, symbol, calc):
        Ewald.__init__(self, position, cell_size, grid)
        self.symbol = symbol
        self.info = buck[self.symbol]
        self.dist = {}
        self.calc_name = calc  # Store the calculator name instead of the object

    def dist_init(self, N):
        for pair in self.info["ion_pair"]:
            self.dist[pair] = np.zeros((N, N))

    def compute_ase_energy(self, pair, pos, cell_size, i_range, calc_name):
        local_dist = np.zeros((len(pos), len(pos)))
        symbol1 = pair.split("+")
        symbol2 = "".join(symbol1)
        symbol3, symbol4 = symbol1[0], symbol1[1]

        calc = KIM(calc_name)  # Initialize KIM inside the worker process
        
        for i in i_range:
            for j in range(i, len(pos)):
                if symbol3 == symbol4:
                    if i == j:
                        atoms = Atoms(symbol3, positions=[pos[i]], cell=np.multiply(cell_size, 1e10), pbc=pbc_t)
                        atoms.calc = calc
                        local_dist[i, j] = atoms.get_potential_energy()
                    elif i != j:
                        atoms = Atoms(symbol2, positions=[pos[i], pos[j]], cell=np.multiply(cell_size, 1e10), pbc=pbc_t)
                        atoms.calc = calc
                        local_dist[i, j] = atoms.get_potential_energy()
                elif symbol3 != symbol4 and i != j:
                    atoms = Atoms(symbol2, positions=[pos[i], pos[j]], cell=np.multiply(cell_size, 1e10), pbc=pbc_t)
                    atoms.calc = calc
                    local_dist[i, j] = atoms.get_potential_energy()

        return local_dist

    def ase_mat_gen(self):
        N = len(self.position)
        self.dist_init(N)
        pos = self.position
        
        pool = mp.get_context("spawn").Pool(mp.cpu_count())  # Use "spawn" to avoid pickle issues
        results = []

        for pair in self.info["ion_pair"]:
            i_ranges = np.array_split(np.arange(N), mp.cpu_count())
            args = [(pair, pos, self.cell_size, i_range, self.calc_name) for i_range in i_ranges]
            results.append(pool.starmap(self.compute_ase_energy, args))

        pool.close()
        pool.join()
        
        for idx, pair in enumerate(self.info["ion_pair"]):
            self.dist[pair] = sum(results[idx])


def energy_kim_addup(N, Vars, o_pos, energy, types, ase, info):
    # addup function is same with buckingham, which different is the matrix of objective function is changed into ase matrix
    for ion_pair in info:
        i_pair, ion_pair = ion_pair, list(ion_pair.split("+")) # The first i_pair is used to extract buckingham matrix...
        if not ((ion_pair[0] in types) and (ion_pair[1] in types)):
            continue

        if ion_pair[0] == ion_pair[1]: # example, ion_pair = ["O", "O"]
            j1 = types.index(ion_pair[0])
            for i1 in range(N):
                for i2 in range(i1, N):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * ase[i_pair][i1, i2])

        else: # example, ion_pair = ["Sr", "Ti"]
            j1 = types.index(ion_pair[0])
            j2 = types.index(ion_pair[1])

            for i1 in range(N):
                for i2 in range(i1 + 1, N):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j2][o_pos[i2]] * ase[i_pair][i1, i2])
                    energy.add(Vars[j2][o_pos[i1]] * Vars[j1][o_pos[i2]] * ase[i_pair][i1, i2])



             

