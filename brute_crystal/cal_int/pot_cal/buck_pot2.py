import numpy as np  
import math
import multiprocessing as mp
from cal_int.pot_cal.ewal_pot2 import Ewald
from data.buck import buck  
from data.radii import radii
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee

class Buckingham(Ewald):

    def __init__(self, position, cell_size, grid, symbol):
        Ewald.__init__(self, position, cell_size, grid)
        self.symbol = symbol
        self.info = buck[self.symbol]
        self.dist = {}

    def dist_init(self, N):
        for pair in self.info["ion_pair"]:
            self.dist[pair] = np.zeros((N, N))

    def compute_buckingham(self, pair, pos, shifts, close_dis, i_range, N, MAX, A, B, C):
        local_dist = np.zeros((N, N))
        for i in i_range:
            for j in np.arange(i, N):
                if i != j:
                    r = np.linalg.norm(pos[i,] - pos[j,])
                    if r < close_dis:
                        local_dist[i, j] = MAX
                    else:
                        local_dist[i, j] += A * math.exp(-r / (B * 1e-10)) - (C * 1e-60) / (r ** 6)
                for s in np.arange(len(shifts)):
                    r = np.linalg.norm(pos[i,] + shifts[s,] - pos[j,])
                    if r < close_dis:
                        local_dist[i, j] = MAX
                        break
                    else:
                        local_dist[i, j] += A * math.exp(-r / (B * 1e-10)) - (C * 1e-60) / (r ** 6)
        return local_dist

    def buck_mat_gen(self):
        N = len(self.position)
        realDepth = 4
        reciprocalDepth = 4
        pos = self.position
        self.dist_init(N)
        shifts = np.zeros(((2 * realDepth + 1) ** 3 - 1, 3))
        i = 0  
        tmp = np.array([realDepth, realDepth, realDepth])

        for shift in np.ndindex(2 * realDepth + 1, 2 * realDepth + 1, 2 * realDepth + 1):
            if shift != (realDepth, realDepth, realDepth):
                shifts[i,] = shift - tmp
                i += 1
        shifts = shifts @ self.vecs
        close_dis = self.close_distance_cal()
        MAX = 300
        
        pool = mp.Pool(mp.cpu_count())
        results = []

        for idx, pair in enumerate(self.info["ion_pair"]):
            A, B, C = self.info["A"][idx], self.info["B"][idx], self.info["C"][idx]
            i_ranges = np.array_split(np.arange(N), mp.cpu_count())
            args = [(pair, pos, shifts, close_dis[idx], i_range, N, MAX, A, B, C) for i_range in i_ranges]
            results.append(pool.starmap(self.compute_buckingham, args))

        pool.close()
        pool.join()

        for idx, pair in enumerate(self.info["ion_pair"]):
            self.dist[pair] = sum(results[idx])

    def close_distance_cal(self, radius_threshold=0.75):
        close_dis = []
        radii_data = radii[self.symbol]
        for pair in self.info["ion_pair"]:
            sp_pair = pair.split("+")
            sym1, sym2 = sp_pair[0], sp_pair[1]
            c_dis = radius_threshold * (radii_data[sym1] + radii_data[sym2]) * 1e-10
            close_dis.append(c_dis)
        return close_dis

    def buck_mat_gen2(self):
        N = len(self.position)
        self.dist_init(N)
        MAX = 300
        max_cells = []
        closest_distance = self.close_distance_cal()
        pos = self.position
        
        for hi in self.info["E"]:
            max_cells.append(int(np.ceil(hi * 1e-10 / self.cell_size[0])))

        for idx, pair in enumerate(self.info["ion_pair"]):
            for i in np.arange(N):
                for j in np.arange(i, N):
                    for ii in range(-max_cells[idx], max_cells[idx] + 1):
                        for jj in range(-max_cells[idx], max_cells[idx] + 1):
                            for kk in range(-max_cells[idx], max_cells[idx] + 1):
                                if not (ii == 0 and jj == 0 and kk == 0):
                                    r = np.linalg.norm(
                                        pos[i,] + self.cell_size * np.array([ii, 0, 0]) +
                                        self.cell_size * np.array([0, jj, 0]) +
                                        self.cell_size * np.array([0, 0, kk]) - pos[j,]
                                    )
                                    if r <= self.info["E"][idx]:
                                        self.dist[pair][i, j] += (
                                            self.info["A"][idx] * math.exp(-r / (self.info["B"][idx] * 1e-10))
                                            - (self.info["C"][idx] * 1e-60) / (r ** 6)
                                        )
                                    if r < closest_distance[idx]:
                                        self.dist[pair][i, j] = MAX
                    if np.linalg.norm(pos[i,] - pos[j,]) > 0.001 * 1e-10:
                        r = np.linalg.norm(pos[j,] - pos[i,])
                        if r <= self.info["E"][idx] * 1e-10:
                            self.dist[pair][i, j] += (
                                self.info["A"][idx] * math.exp(-r / (self.info["B"][idx] * 1e-10))
                                - (self.info["C"][idx] * 1e-60) / (r ** 6)
                            )
                        if r < closest_distance[idx]:
                            self.dist[pair][i, j] = MAX

    def operator(self):
        self.vec_gen()
        self.buck_mat_gen()

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


