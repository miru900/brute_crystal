import numpy as np
import os
from cal_int.pot_cal.grid_gen import cubic
import math
from numba import jit, njit
from time import sleep, time
from multiprocessing import Pool
from functools import partial

@njit(cache=True)
def BuckinghamTwoIons(pos_i, pos_j, cell_size, A, rho, beta, lo, hi, closest_distance):
    """
    Interaction between two ions under PBC
    pos_i, pos_j belongs to 3D

    TODO:Double check the ranges!
    """
    MAX = 300
    max_cells = int(np.ceil(hi / cell_size))
    energy = 0

    # HACK BEGINS
    # Dirty hack to handle closest distances
    # TODO: remove this
    if rho < 0.00001:
        max_cells = 1

        for i in range(-max_cells, max_cells + 1):
            for j in range(-max_cells, max_cells + 1):
                for k in range(-max_cells, max_cells + 1):
                    if not (i == 0 and j == 0 and k == 0):
                        r = cell_size * np.linalg.norm(pos_j + np.array([i, 0, 0]) +
                                                       np.array([0, j, 0]) + np.array([0, 0, k]) - pos_i)

                        if r < closest_distance:
                            return MAX

        if np.linalg.norm(pos_i - pos_j) > 0.001:  # interaciton within the cell, you can remove it to the row function
            r = cell_size * np.linalg.norm(pos_j - pos_i)

            if r < closest_distance:
                return MAX

        return 0
    # HACK ENDS
    else:
        # interaction with neighbouring cells
        for i in range(-max_cells, max_cells + 1):
            for j in range(-max_cells, max_cells + 1):
                for k in range(-max_cells, max_cells + 1):
                    if not (i == 0 and j == 0 and k == 0):
                        r = cell_size * np.linalg.norm(pos_j + np.array([i, 0, 0]) +
                                                       np.array([0, j, 0]) + np.array([0, 0, k]) - pos_i)
                        if r <= hi:
                            energy += A * math.exp(-1.0 * r / rho) - beta / r ** 6

                        if r < closest_distance:
                            return MAX

        if np.linalg.norm(
                pos_i - pos_j) > 0.001:  # add interaciton within the cell, you can remove it to the row function
            r = cell_size * np.linalg.norm(pos_j - pos_i)
            if r <= hi:
                energy += A * math.exp(-1.0 * r / rho) - beta / r ** 6

            if r < closest_distance:
                return MAX

    return energy

def _Buck_row(pos_i, grid_size, cell_size, A, rho, beta, lo, hi, closest_distance):
    """
    elements is the pair of elements for which we are going to compute the matrix
    """
    result = np.zeros(grid_size ** 3)

    positions = cubic(grid_size)
    for pos_j in range(pos_i, grid_size ** 3):
        # print(BuckinghamTwoIons(pos_i, pos_j, cell_size, A, rho, beta, lo, hi, radius_threshold))
        result[pos_j] = BuckinghamTwoIons(positions[pos_i], positions[pos_j], cell_size, A, rho, beta, lo, hi,
                                          closest_distance)
    return result

def generate_Buck(grid_size, cell_size, phase, output_directory, radius_threshold=0.75, multicpu=False):
    """
    Cubic systems only for the moment
    size is the cube size
    phase
    TODO: compute distances once, e.g. (i1, i2): [cell (0,0,0), cell (1,0,0), cell (0,1,0), cell (0,0,1), ...]
    Then apply the Buckingham function
    """
    for ion_pair, potential_param in phase.buck.items():

        result = np.zeros((grid_size ** 3, grid_size ** 3))
        closest_distance = radius_threshold * (phase.radius[ion_pair[0]] + phase.radius[ion_pair[1]])

        if ion_pair in phase.closest_distance:
            closest_distance = phase.closest_distance[ion_pair]
            print("Closest distance for " + ion_pair[0] + '-' + ion_pair[1] + " was set to be " + str(
                phase.closest_distance[ion_pair]))

        buck_row = partial(_Buck_row, grid_size=grid_size, cell_size=cell_size, A=potential_param['par'][0],
                           rho=potential_param['par'][1], beta=potential_param['par'][2],
                           lo=potential_param['lo'], hi=potential_param['hi'], closest_distance=closest_distance)

        if multicpu:
            with Pool(processes=4) as pool:
                i = 0
                for row in pool.map(buck_row, range(grid_size ** 3)):
                    result[i,] = row
                    i += 1
            # pass
        else:
            i = 0
            for row in map(buck_row, range(grid_size ** 3)):
                result[i,] = row
                i += 1

        # for i in prange(size):
        #     result[i] = test()

        filename = 'C{grid_size}_'.format(grid_size=grid_size) + ion_pair[0] + ion_pair[1] + '_{cell_size}.npy'.format(
            cell_size=cell_size)

        # print(filename)
        # print(ion_pair, '\n', np.around(result, decimals=2))

        with open(output_directory / filename, 'wb') as outfile:
            np.save(outfile, result)

def get_Buck(ion_pair, grid_size, cell_size, phase):
    """
    Loads or generates buckingham matrix
    """
    filename = 'C{grid_size}_'.format(grid_size=grid_size) + ion_pair[0] + ion_pair[1] + '_{cell_size}.npy'.format(
        cell_size=cell_size)

    try:
        with open(phase.location / filename, 'rb') as f:
            return np.load(f)
    except IOError:
        generate_Buck(grid_size, cell_size, phase, phase.location, multicpu=True)
        with open(phase.location / filename, 'rb') as f:
            return np.load(f)

def energy_buck_addup(N, Vars, o_pos, energy, types, grid, cell, phase):
    for ion_pair in phase.buck:
        
        if not ((ion_pair[0] in types) and (ion_pair[1] in types)):
            continue
        
        buck = get_Buck(ion_pair, grid, cell, phase)
        
        if ion_pair[0] == ion_pair[1]:
            j1 = types.index(ion_pair[0])
            
            for i1 in range(N):
                for i2 in range(i1, N):
                    # if i1 == 58 and i2 == 59:
                    #     print(buck[i1, i2])
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j1][o_pos[i2]] * buck[i1, i2])
        
        else:
            j1 = types.index(ion_pair[0])
            j2 = types.index(ion_pair[1])
            
            for i1 in range(N):
                # TODO: Two different ions at the same position is impossible, should we put MAX here?
                # energy.add(Vars[j1][i1] * Vars[j1][i1] * buck[i1, i1])
                
                for i2 in range(i1 + 1, N):
                    energy.add(Vars[j1][o_pos[i1]] * Vars[j2][o_pos[i2]] * buck[i1, i2])
                    energy.add(Vars[j2][o_pos[i1]] * Vars[j1][o_pos[i2]] * buck[i1, i2])

