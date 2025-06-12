from cal_kim.open_grid import chem_compo as og_cc
# importing file opening and symbol splitting function from cal_kim directory
import numpy as np

def mat_gen(pos_data, cell_size, ortho, angle):

    angle_if = isinstance(angle, list)

    if angle_if:
        # lattice vector preparation code
        a, b, c = cell_size[0], cell_size[1], cell_size[2]
        sf = np.pi / 180  # convert degrees to radians
        alp, bet, gam = angle[0] * sf, angle[1] * sf, angle[2] * sf

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
        l_vector = np.array([x_vec, y_vec, z_vec])
        pos_data_real_cartesian = []
    
        for line in pos_data:
            real_pos = np.dot(line, l_vector)
            real_pos = list(real_pos)
            pos_data_real_cartesian.append(real_pos)

        return np.multiply(pos_data_real_cartesian, 1) * 1e-10

    # pos_data multiplication part, for ortho and simple cubic, just multiplying 1e-10 enough.
    elif ortho and not angle_if:
        return np.multiply(pos_data, cell_size) * 1e-10
    else:
        return np.array(pos_data) * cell_size * 1e-1

def charge_data_to_list(directory):
    direc = f"./data/charge/" + directory
    charge2 = [] # charge will be filled into this list, 1 to last atom

    charge_data = open(direc, "r")
    for i in charge_data:
        i = i.strip()
        try:
            i = float(i)
            i = round(i, 4)
            charge2.append(i)
        except:
            print("Error!")

    print(charge2, len(charge2))

    return charge2


