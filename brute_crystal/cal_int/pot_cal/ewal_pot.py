import numpy as np  ;  import os  ;  import math
from numba import jit, njit  ;  from time import sleep, time  ;  from multiprocessing import Pool
from functools import partial  ;  from pathlib import Path
from cal_int.pot_cal.grid_gen import cubic  ;  from cal_int.pot_cal.grid_gen import generate_orthorhombic as gen_ortho

filedir = Path(__file__).resolve().parent

def QEwald(positions, vecs, reciprocal, cell_volume, alpha=-1):
    """
    positions are relative coordinates
    vecs are the lattice vectors
    reciprocal are the reciprocal lattice vectors. I need them only to use numba, since it doesn't support cross product
    cell_volume is the volume of the cell
    alpha controls the split between the direct and reciprocal sums

    You can tweak realDepth and reciprocalDepth summation constant below in the code
    """
    N = len(positions)
    if alpha < 0:
        alpha = 2 / (cell_volume ** (1.0 / 3))
    realDepth = 4  # 3
    reciprocalDepth = 4  # 3

    pos = positions @ vecs

    d = np.zeros((N, N))

    shifts = np.zeros(((2 * realDepth + 1) ** 3 - 1, 3))
    i = 0
    tmp = np.array([realDepth, realDepth, realDepth])

    for shift in np.ndindex(2 * realDepth + 1, 2 * realDepth + 1, 2 * realDepth + 1):
        if shift != (realDepth, realDepth, realDepth):
            shifts[i,] = shift
            shifts[i,] = shifts[i,] - tmp
            i = i + 1

    shifts = shifts @ vecs

    for i in np.arange(N):
        for j in np.arange(i, N):
            if i != j:
                r = np.linalg.norm(pos[i,] - pos[j,])
                d[i, j] += math.erfc(alpha * r) / (2 * r)

            for s in np.arange(len(shifts)):
                r = np.linalg.norm(pos[i,] + shifts[s,] - pos[j,])
                d[i, j] += math.erfc(alpha * r) / (2 * r)

    # self interaction term
    for i in np.arange(N):
        d[i, i] = d[i, i] - alpha / math.sqrt(math.pi)

    # Ewald reciprocal space

    shifts_recip = np.zeros(((2 * reciprocalDepth + 1) ** 3 - 1, 3))
    i = 0
    tmp = np.array([reciprocalDepth, reciprocalDepth, reciprocalDepth])

    for shift_recip in np.ndindex(2 * reciprocalDepth + 1, 2 * reciprocalDepth + 1, 2 * reciprocalDepth + 1):
        if shift_recip != (reciprocalDepth, reciprocalDepth, reciprocalDepth):
            shifts_recip[i,] = shift_recip
            shifts_recip[i,] = shifts_recip[i,] - tmp
            i = i + 1

    shifts_recip = shifts_recip @ reciprocal

    for i in np.arange(N):
        for j in np.arange(i, N):

            for s in np.arange(len(shifts_recip)):
                k = shifts_recip[s,]
                # k = np.array(shift)@self.reciprocal
                term = (4 * math.pi ** 2) / np.dot(k, k)
                term = term * math.exp(-np.dot(k, k) / (4 * alpha ** 2))
                v = pos[j,] - pos[i,]
                term = term * math.cos(np.dot(k, v))
                d[i, j] += term / (2 * math.pi * cell_volume)

    # Unit conversion
    d = d * 14.399645351950543

    # symmetry completion
    for i in np.arange(N):
        for j in np.arange(i):
            d[i, j] = d[j, i]

    return d

def generate_Ewald(size, ouput_directory, ortho = False):
    """
    Start with the cubic systems
    """
    filename = 'C{size}.npy'.format(size=size)

    # Generate the grid and compute preliminary parameters
    if ortho:
        positions = gen_ortho(size)
        vecs = np.zeros((3, 3))
        for i in range(3):
            vecs[i, i] = size[i]
        cell_volume = 1
        for one_len in size:
            cell_volume = cell_volume * one_len

    else:
        positions = cubic(size)
        vecs = np.zeros((3, 3))
        for i in range(3):
            vecs[i, i] = size
        cell_volume = abs(np.linalg.det(vecs))


    reciprocal = np.zeros((3, 3))

    for i in np.arange(3):
        recip_vector = 2 * math.pi * np.cross(vecs[(1 + i) % 3,], vecs[(2 + i) % 3]) / cell_volume
        reciprocal[i,] = recip_vector

    dist = QEwald(positions, vecs, reciprocal, cell_volume)
    print('Ewald matrix for cubic system of size {size} was generated'.format(size=size))
    print("Its max is ", dist.max())
    print("Its min is ", dist.min())
    # print(dist)

    with open(ouput_directory / filename, 'wb') as outfile:
        np.save(outfile, dist)


def get_Ewald(grid_size, cell_size, ortho = False):
    """
    Generating Ewald matrix. It is done only once when the code is downloaded
    """
    filename = 'C{grid_size}.npy'.format(grid_size=grid_size)

    try:
        with open(filedir / 'Ewald' / filename, 'rb') as f:
            if ortho:
                return np.load(f) * (grid_size / np.average(cell_size))
            else:
                return np.load(f) * (grid_size / cell_size)

    except IOError:
        generate_Ewald(grid_size, filedir / 'Ewald/')
        with open(filedir / 'Ewald' / filename, 'rb') as f:
            return np.load(f) * (grid_size / cell_size)

# ewad(O, T, Vars, orb_key, matrix, charge, energy, chem)

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
