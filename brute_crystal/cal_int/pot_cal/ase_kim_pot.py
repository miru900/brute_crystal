from ase import Atoms
from ase.calculators.kim import KIM
from data.buck import buck  ;  from data.radii import radii
from cal_int.pot_cal.ewal_pot2 import Ewald
import numpy as np

pbc_t = [True] * 3  ;  pbc_f = [False] * 3

class ASE_KIM(Ewald):

    def __init__(self, position, cell_size, grid, symbol, calc):
        Ewald.__init__(self, position, cell_size, grid) # Caution! : Because we are using ASE, we will not use scaled position which is multiplied by 1e-10, we will use angstrom version ( X : 3.9 * 1e-10, O : 3.9)
        self.symbol = symbol
        self.info = buck[self.symbol] # We are using this because we need ion_pair! If you are using ASE models, write informations, in buck.py file, about molecule like the format written before.
        self.dist = {}
        self.calc = KIM(calc) # Calculation method, same with mod1 and mod2.

    def dist_init(self, N):
        # For each ion pair, we generate N, N empty matrix, this function is used in 3rd line of buck_mat_gen(self)
        for pair in self.info["ion_pair"]:
            self.dist[pair] = np.zeros((N, N))

    def ase_mat_gen(self):
        # programs/lattice_prac/prac5.py -> refer!!!
        # from buck_pot2.py, if the two symbol same, calcul [i, i], if different, treat [i, i] as zero, they are not calculated in ase_pot_addup, because range is i : 0~N and j : i+1 ~ N
        N = len(self.position)  ;  pos = self.position  ;  self.dist_init(N)
        for pair in self.info["ion_pair"]:
            symbol1 = pair.split("+")
            symbol2 = "".join(symbol1)
            symbol3 = symbol1[0]  ;  symbol4 = symbol1[1]
            for i in range(N):
                for j in range(i, N):
                    if symbol3 == symbol4:
                        if i == j:
                            atoms = Atoms(symbol3, positions = [pos[i]], cell = np.multiply( self.cell_size, 1e10 ) , pbc = pbc_t) # Why we multiply 1e10? Because we inherited Ewald, which multiply 1e-10 to settings[key][cell_size]
                            atoms.calc = self.calc  ;  energy = atoms.get_potential_energy()  ;  self.dist[pair][i, j] = energy
                        elif i != j:
                            atoms = Atoms(symbol2, positions = [pos[i], pos[j]], cell = np.multiply( self.cell_size, 1e10 ), pbc = pbc_t)
                            atoms.calc = self.calc  ;  energy = atoms.get_potential_energy()  ;  self.dist[pair][i, j] = energy
                    elif symbol3 != symbol4:
                        if i != j:
                            atoms = Atoms(symbol2, positions = [pos[i], pos[j]], cell = np.multiply( self.cell_size, 1e10 ), pbc = pbc_t)
                            atoms.calc = self.calc  ;  energy = atoms.get_potential_energy()  ;  self.dist[pair][i, j] = energy

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



             

