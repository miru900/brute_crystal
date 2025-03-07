import numpy as np  ;  import os  ;  import math
from scipy.constants import elementary_charge as elec_c
from scipy.constants import electron_volt as eevee
from scipy.constants import epsilon_0 as ep0

scale = (elec_c ** 2)  / (4 * math.pi * ep0 * eevee) # scaling factor : e ** 2 / 4 * π * ε_0

class Ewald:

    def __init__(self, position, cell_size, grid):
        # The initialization of position (C#.txt file) and grid size (example, it can be just 2 for cubic, or [3, 2, 2] for orthombic)
        self.position = position
        if isinstance(cell_size, list):
            self.cell_size = np.multiply( cell_size, 1e-10 )  ;  self.grid = grid
        else:
            self.cell_size = np.multiply( [cell_size] * 3, 1e-10 )  ;  self.grid = [grid] * 3

    def vec_gen(self):
        # The lattice vectors and reciprocal vectors are generated in this function, saving it to called instances
        # Lattice vector and cell volume part
        self.vecs = np.zeros( (3, 3) )
        for i in range(3):
            self.vecs[i, i] = self.cell_size[i]
        self.cell_volume = abs(np.linalg.det(self.vecs))

        # Reciprocal vector part
        self.reciprocal = np.zeros( (3, 3) )
        for i in np.arange(3):
            recip_vector = 2 * math.pi * np.cross(self.vecs[(1 + i) % 3,], self.vecs[(2 + i) % 3]) / self.cell_volume
            self.reciprocal[i,] = recip_vector

    def ewal_mat_gen(self, alpha = -1):
        # This part generates matrix which will be used in integer programming (the alpha, objective function = sigma [alpha * pos_1 * pos_2]), self.dist is the return matrix
        # Basic setting for matrix generation
        N = len(self.position)  ;  realDepth = 4  ;  reciprocalDepth = 4
        pos = self.position  ;  self.dist = np.zeros((N, N))  ;  shifts = np.zeros(((2 * realDepth + 1) ** 3 - 1, 3))
        i = 0  ;  tmp = np.array([realDepth, realDepth, realDepth])
        if alpha < 0:
            alpha = 2 / (self.cell_volume ** (1.0 / 3))
        
        # shfit matrix generator
        for shift in np.ndindex(2 * realDepth + 1, 2 * realDepth + 1, 2 * realDepth + 1):
            if shift != (realDepth, realDepth, realDepth):
                shifts[i,] = shift  ;  shifts[i,] = shifts[i,] - tmp  ;  i = i + 1
        shifts = shifts @ self.vecs

        # matrix generator 1 : Ewald real space
        for i in np.arange(N):
            for j in np.arange(i, N):
                if i != j:
                    r = np.linalg.norm(pos[i,] - pos[j,])
                    self.dist[i, j] += math.erfc(alpha * r) / (2 * r)

                for s in np.arange(len(shifts)):
                    r = np.linalg.norm(pos[i,] + shifts[s,] - pos[j,])
                    self.dist[i, j] += math.erfc(alpha * r) / (2 * r)

        # matrix generator 2 : Self interaction term
        for i in np.arange(N):
            self.dist[i, i] = self.dist[i, i] - alpha / math.sqrt(math.pi)

        # matrix generator 3 : Ewald reciprocal space
        shifts_recip = np.zeros(((2 * reciprocalDepth + 1) ** 3 - 1, 3))  ;  i = 0
        tmp = np.array([reciprocalDepth, reciprocalDepth, reciprocalDepth])

        for shift_recip in np.ndindex(2 * reciprocalDepth + 1, 2 * reciprocalDepth + 1, 2 * reciprocalDepth + 1):
            if shift_recip != (reciprocalDepth, reciprocalDepth, reciprocalDepth):
                shifts_recip[i,] = shift_recip  ;  shifts_recip[i,] = shifts_recip[i,] - tmp  ;  i = i + 1
        shifts_recip = shifts_recip @ self.reciprocal

        for i in np.arange(N):
            for j in np.arange(i, N):
                for s in np.arange(len(shifts_recip)):
                    k = shifts_recip[s,]  ;  term = (4 * math.pi ** 2) / np.dot(k, k)  ;  term = term * math.exp(-np.dot(k, k) / (4 * alpha ** 2))
                    v = pos[j,] - pos[i,]  ;  term = term * math.cos(np.dot(k, v))  ;  self.dist[i, j] += term / (2 * math.pi * self.cell_volume)

        # symmetry (upper triangle = lower triangle) and total matrix scaling
        for i in np.arange(N):
            for j in np.arange(i):
                self.dist[i, j] = self.dist[j, i]       
        self.dist = self.dist * scale
        
    def operator(self):
        vec_gen(self)
        ewal_mat_gen(self, alpha = -1)


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





