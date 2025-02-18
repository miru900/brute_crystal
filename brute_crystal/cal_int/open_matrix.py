from cal_kim.open_grid import chem_compo as og_cc
# importing file opening and symbol splitting function from cal_kim directory
import numpy as np

def mat_gen(pos_data, cell_size):
    return np.array(pos_data) * cell_size * 1e-10
