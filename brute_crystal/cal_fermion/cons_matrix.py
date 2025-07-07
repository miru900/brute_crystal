import numpy as np
  
def matrix_cons_adder(W, cons, ion_count):

    size = len(cons[0])

    for i in range(size):
        idx = cons[i]
        W[idx, idx] += 1 - 2 * ion_count[0]
        for j in range(i + 1, size):
            jdx = cons[j]
            W[idx, jdx] += 2

    return W

    

