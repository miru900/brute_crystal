import numpy as np

def total_array(ewal_mat, chem, charge, grid_size):
    
    print(f"{grid_size}, this is grid_size, debugging")
    ion_type_count = len(chem)
    print(f"{ion_type_count}, this is ion_type_count, debugging")
    total_size = grid_size * ion_type_count
    print(f"{total_size}, this is total_size, debugging")

    total_array = np.zeros((total_size, total_size))
    
    # total_array generating part
    for i in range(ion_type_count):
        for j in range(i, ion_type_count):
            single_array = ewal_mat * charge[chem[i]] * charge[chem[j]]

            for k in range(grid_size):
                for l in range(grid_size):
                    total_array[k + i * grid_size, l + j * grid_size] = single_array[k, l]
    
    return total_array, ion_type_count, grid_size


