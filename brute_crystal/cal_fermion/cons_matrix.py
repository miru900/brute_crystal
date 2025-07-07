import numpy as np

def matrix_cons_adder(W, cons, ion_count):

    size = cons[0][-1]

    for i in range(size):
        W[i, i] += 1 - 2 * ion_count[0]
        for j in range(i + 1, size):
            W[i, j] += 2

    return W
    

