import numpy as np  
import os  
import math
import multiprocessing as mp
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee
from scipy.constants import epsilon_0 as ep0

scale = (elec_c ** 2)  / (4 * math.pi * ep0 * eevee)

class Ewald:
    def __init__(self, position, cell_size, grid, angle = None):
        self.position = position
        self.angle = angle
        self.angle_T = isinstance(angle, list)
        if isinstance(cell_size, list):
            self.cell_size = np.multiply(cell_size, 1e-10)
            self.grid = grid
        else:
            self.cell_size = np.multiply([cell_size] * 3, 1e-10)
            self.grid = [grid] * 3

    def vec_gen(self):
        if self.angle_T:
            # lattice vector preparation code
            a, b, c = self.cell_size[0], self.cell_size[1], self.cell_size[2]
            sf = np.pi / 180  # convert degrees to radians
            alp, bet, gam = self.angle[0] * sf, self.angle[1] * sf, self.angle[2] * sf

            # Trigonometric terms
            cos_alpha = np.cos(alp)
            cos_beta  = np.cos(bet)
            cos_gamma = np.cos(gam)
            sin_gamma = np.sin(gam)

            # Lattice vectors
            x_vec = [a, 0, 0]
            y_vec = [b * cos_gamma, b * sin_gamma, 0]

            z1 = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
            z2 = np.sqrt(1 - cos_beta**2 - z1**2)
            z_vec = [c * cos_beta, c * z1, c * z2]

            # fractional pos_data -> real cartesian coordinate
            self.vecs = np.array([x_vec, y_vec, z_vec])
            self.cell_volume = abs(np.linalg.det(self.vecs))

        else:
            self.vecs = np.zeros((3, 3))
            for i in range(3):
                self.vecs[i, i] = self.cell_size[i]
            self.cell_volume = abs(np.linalg.det(self.vecs))

        self.reciprocal = np.zeros((3, 3))
        for i in np.arange(3):
            recip_vector = 2 * math.pi * np.cross(self.vecs[(1 + i) % 3,], self.vecs[(2 + i) % 3]) / self.cell_volume
            self.reciprocal[i,] = recip_vector

    def _compute_real_space(self, i, pos, shifts, alpha):
        N = len(pos)
        row = np.zeros(N)
        for j in np.arange(i, N):
            if i != j:
                r = np.linalg.norm(pos[i,] - pos[j,])
                row[j] += math.erfc(alpha * r) / (2 * r)

            for s in range(len(shifts)):
                r = np.linalg.norm(pos[i,] + shifts[s,] - pos[j,])
                row[j] += math.erfc(alpha * r) / (2 * r)
        return i, row

    def _compute_reciprocal_space(self, i, pos, shifts_recip, alpha):
        N = len(pos)
        row = np.zeros(N)
        for j in np.arange(i, N):
            for s in range(len(shifts_recip)):
                k = shifts_recip[s,]
                term = (4 * math.pi ** 2) / np.dot(k, k) * math.exp(-np.dot(k, k) / (4 * alpha ** 2))
                v = pos[j,] - pos[i,]
                term *= math.cos(np.dot(k, v))
                row[j] += term / (2 * math.pi * self.cell_volume)
        return i, row

    def ewal_mat_gen(self, alpha=-1):
        N = len(self.position)
        realDepth, reciprocalDepth = 4, 4
        pos = self.position
        self.dist = np.zeros((N, N))
        if alpha < 0:
            alpha = 2 / (self.cell_volume ** (1.0 / 3))

        tmp = np.array([realDepth, realDepth, realDepth])
        shifts = np.array([shift - tmp for shift in np.ndindex((2 * realDepth + 1,) * 3) if shift != (realDepth, realDepth, realDepth)]) @ self.vecs
        shifts_recip = np.array([shift - tmp for shift in np.ndindex((2 * reciprocalDepth + 1,) * 3) if shift != (reciprocalDepth, reciprocalDepth, reciprocalDepth)]) @ self.reciprocal

        with mp.Pool(mp.cpu_count()) as pool:
            real_results = pool.starmap(self._compute_real_space, [(i, pos, shifts, alpha) for i in range(N)])
            recip_results = pool.starmap(self._compute_reciprocal_space, [(i, pos, shifts_recip, alpha) for i in range(N)])

        for i, row in real_results:
            self.dist[i, i:] += row[i:]
        for i, row in recip_results:
            self.dist[i, i:] += row[i:]

        for i in range(N):
            self.dist[i, i] -= alpha / math.sqrt(math.pi)
        for i in range(N):
            for j in range(i):
                self.dist[i, j] = self.dist[j, i]

        self.dist *= scale

    def operator(self):
        self.vec_gen()
        self.ewal_mat_gen()

def energy_ewal_addup(N, T, Vars, o_pos, dist, charge, energy, chem):
    # np.savetxt('testEwaldNew.out', dist, delimiter=',')
    for i1 in range(N):

        for j1 in range(T):  # self-interaction
            energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i1]] * dist[i1, i1] * charge[chem[j1]] ** 2)

        for i2 in range(i1 + 1, N):
            # print(i2)
            for j1 in range(T):  # pairwise Coulumb
                energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * 2 * dist[i1, i2] * charge[chem[j1]] ** 2)  # i1,i2 have the same type of ion

                for j2 in range(j1 + 1, T):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j2][o_pos[i2]] * 2 * dist[i1, i2] * charge[chem[j1]] * charge[chem[j2]])  # Two different types
                    energy.add(Vars[j2][o_pos[i1]] * Vars[j1][o_pos[i2]] * 2 * dist[i1, i2] * charge[chem[j1]] * charge[chem[j2]])  # Symmetrical case


def energy_ewal_addup_charge2(N, T, Vars, o_pos, dist, charge2, energy, chem):
    # N = grid size, T = ion type count
    # np.savetxt('testEwaldNew.out', dist, delimiter=',')
    for i1 in range(N):

        for j1 in range(T):  # self-interaction
            energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i1]] * dist[i1, i1] * charge2[i1] ** 2)

        for i2 in range(i1 + 1, N):
            # print(i2)
            for j1 in range(T):  # pairwise Coulumb
                energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * 2 * dist[i1, i2] * charge2[i1] * charge2[i2])  # i1,i2 have the same type of ion

                for j2 in range(j1 + 1, T):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j2][o_pos[i2]] * 2 * dist[i1, i2] * charge2[i1] * charge2[i2])  # Two different types
                    energy.add(Vars[j2][o_pos[i1]] * Vars[j1][o_pos[i2]] * 2 * dist[i1, i2] * charge2[i1] * charge2[i2])  # Symmetrical case

def energy_ewal_addup_cons(N, T, Vars, o_pos, dist, charge, energy, chem, cons):

    # variables -> varis and fixed
    grid_size = dist.shape[0]
    varis, fixed, checker = {}, {}, {}
    print(cons)
    for i, elem in enumerate(chem):
        for j in range(grid_size):
            if j not in cons[i]:
                continue
            var = f"{elem}_{j}"
            num = j
            if "Fixed" in cons[i]:
                fixed[var] = num
            else:
                varis[var] = num
            checker[var] = num

    print("=== Variables Debug ===")
    print("varis:", varis)
    print("fixed:", fixed)

    for i1 in range(N):

        for j1 in range(T):
            var = f"{chem[j1]}_{i1}"
            if var not in checker:
                continue
            if var in fixed:
                continue
            energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i1]] * dist[i1, i1] * charge[chem[j1]] ** 2)

        for i2 in range(i1 + 1, N):
            
            for j1 in range(T):
                var1 = f"{chem[j1]}_{i1}"
                var2 = f"{chem[j1]}_{i2}"
                if var1 not in checker or var2 not in checker:
                    continue
                if var1 in fixed and var2 in fixed:
                    continue
                if var1 in fixed or var2 in fixed:
                    continue
                energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * 2 * dist[i1, i2] * charge[chem[j1]] ** 2)  # i1,i2 have the same type of ion

                for j2 in range(j1 + 1, T):
                    var2 = f"{chem[j2]}_{i2}"
                    var3 = f"{chem[j2]}_{i1}"
                    if var2 not in checker or var3 not in checker:
                        continue
                    if var1 in fixed  and var2 in fixed:
                        continue
                    elif var1 not in fixed and var2 in fixed:
                        q_cross1 = 2 * dist[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        energy.add(Vars[j1][o_pos[i1]] * q_cross1)
                        q_cross2 = 2 * W[i1, i2] * charge[chem[j1]] * charge[chem[j2]]
                        energy.add(Vars[j1][o_pos[i2]] * q_cross2)
                        continue
                    elif var1 in fixed and var2 not in fixed:
                        continue


    

