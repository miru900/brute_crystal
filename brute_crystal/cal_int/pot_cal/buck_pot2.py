import numpy as np  ;  import math
from cal_int.pot_cal.ewal_pot2 import Ewald
from data.buck import buck  ;  from data.radii import radii
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee

class Buckingham(Ewald):

    def __init__(self, position, cell_size, grid, symbol):
        Ewald.__init__(self, position, cell_size, grid)
        self.symbol = symbol
        self.info = buck[self.symbol] # from buck.py buckingham parameter and ion_pair is extracted
        self.dist = {}

    def dist_init(self, N):
        # For each ion pair, we generate N, N empty matrix, this function is used in 3rd line of buck_mat_gen(self)
        for pair in self.info["ion_pair"]:
            self.dist[pair] = np.zeros((N, N))

    def buck_mat_gen(self):
        # Basic setting for matrix generation, for shift vector generation
        N = len(self.position)  ;  realDepth = 4  ;  reciprocalDepth = 4
        pos = self.position  ;  self.dist_init(N) # Empty matrix generation    
        shifts = np.zeros(((2 * realDepth + 1) ** 3 - 1, 3))
        i = 0  ;  tmp = np.array([realDepth, realDepth, realDepth])

        # shift matrix generator
        for shift in np.ndindex(2 * realDepth + 1, 2 * realDepth + 1, 2 * realDepth + 1):
            if shift != (realDepth, realDepth, realDepth):
                shifts[i,] = shift  ;  shifts[i,] = shifts[i,] - tmp  ;  i = i + 1
        shifts = shifts @ self.vecs
    
        # matrix generator for each ion pair made by me!
        for idx, pair in enumerate(self.info["ion_pair"]):
            for i in np.arange(N):
                for j in np.arange(i, N):
                    if i != j:
                        r = np.linalg.norm(pos[i,] - pos[j,])
                        self.dist[pair][i, j] += (self.info["A"][idx]) * math.exp(-r / (self.info["B"][idx] * 1e-10)) - (self.info["C"][idx] * 1e-60) / (r ** 6)

                    for s in np.arange(len(shifts)):
                        r = np.linalg.norm(pos[i,] + shifts[s,] - pos[j,])
                        self.dist[pair][i, j] += (self.info["A"][idx]) * math.exp(-r / (self.info["B"][idx] * 1e-10)) - (self.info["C"][idx] * 1e-60) / (r ** 6)
    
    def close_distance_cal(self, radius_threshold=0.75):
        close_dis = []
        radii_data = radii[self.symbol]
        for pair in self.info["ion_pair"]:
            sp_pair = pair.split("+")  ;  sym1 = sp_pair[0]  ;  sym2 = sp_pair[1]
            c_dis = radius_threshold * (radii_data[sym1] + radii_data[sym2]) * 1e-10
            close_dis.append(c_dis)

        return close_dis

    def buck_mat_gen2(self):
        # matrix generator referred from IPCSP
        N = len(self.position)  ;  self.dist_init(N)  ;  MAX = 300  ;  max_cells = []  ;  closest_distance = self.close_distance_cal()  ;  pos = self.position  ;  print(closest_distance)
        for hi in self.info["E"]:
            max_cells.append(int(np.ceil(hi * 1e-10 / self.cell_size[0])))

        for idx, pair in enumerate(self.info["ion_pair"]):
            for i in np.arange(N):
                for j in np.arange(i, N):
                    for ii in range(-max_cells[idx], max_cells[idx] + 1):
                        for jj in range(-max_cells[idx], max_cells[idx] + 1):
                            for kk in range(-max_cells[idx], max_cells[idx] + 1):
                                if not (ii == 0 and jj == 0 and kk == 0):
                                    r = np.linalg.norm(pos[i,] + self.cell_size * np.array([ii, 0, 0]) + 
                                            self.cell_size * np.array([0, jj, 0]) + self.cell_size * np.array([0, 0, kk]) - pos[j,])

                                    if r <= self.info["E"][idx]:
                                        self.dist[pair][i, j] += (self.info["A"][idx]) * math.exp(-r / (self.info["B"][idx] * 1e-10)) - (self.info["C"][idx] * 1e-60) / (r ** 6)

                                    if r < closest_distance[idx]:
                                        self.dist[pair][i, j] = MAX

                    if np.linalg.norm(pos[i,] - pos[j,]) > 0.001 * 1e-10:
                        r = np.linalg.norm(pos[j,] - pos[i,])

                        if r <= self.info["E"][idx] * 1e-10:
                            self.dist[pair][i, j] += (self.info["A"][idx]) * math.exp(-r / (self.info["B"][idx] * 1e-10)) - (self.info["C"][idx] * 1e-60) / (r ** 6)

                        if r < closest_distance[idx]:
                            self.dist[pair][i, j] = MAX

    def operator(self):
        self.vec_gen()
        self.buck_mat_gen()
        # self.buck_mat_gen2()

def energy_buck_addup(N, Vars, o_pos, energy, types, buck, info):
    # buck is buckingham matrix, which is from self.dist in Buckingham class, info is buckingham ion pair 
    for ion_pair in info:
        i_pair, ion_pair = ion_pair, list(ion_pair.split("+")) # The first i_pair is used to extract buckingham matrix...
        if not ((ion_pair[0] in types) and (ion_pair[1] in types)):
            continue

        if ion_pair[0] == ion_pair[1]: # example, ion_pair = ["O", "O"]
            j1 = types.index(ion_pair[0])
            for i1 in range(N):
                for i2 in range(i1, N):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * buck[i_pair][i1, i2])

        else: # example, ion_pair = ["Sr", "Ti"]
            j1 = types.index(ion_pair[0])
            j2 = types.index(ion_pair[1])

            for i1 in range(N):
                for i2 in range(i1 + 1, N):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j2][o_pos[i2]] * buck[i_pair][i1, i2])
                    energy.add(Vars[j2][o_pos[i1]] * Vars[j1][o_pos[i2]] * buck[i_pair][i1, i2])


