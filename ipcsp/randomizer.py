from itertools import permutations

def gen_config(chem, ion_count, grid_size):
    # Create a list with all elements
    element_list = []
    for i in range(len(chem)):
        element_list.extend([chem[i]] * ion_count[i])  # Append each element multiple times

    # Calculate number of empty spaces
    total_slots = grid_size ** 3
    empty_slots = total_slots - sum(ion_count)  # Remaining slots are empty

    # Add empty spaces ('_') to the list
    element_list.extend(["_"] * empty_slots)

    # Generate all unique permutations of the element list
    unique_configs = set(permutations(element_list))

    return list(unique_configs)




