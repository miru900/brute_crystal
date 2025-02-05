from itertools import permutations
import open_grid

# Get user input
chem, ion_count = open_grid.chem_compo(open_grid.compo)

def generate_configurations(chem, ion_count, grid_size):
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

# Execute function and print results
grid_size = open_grid.grid_num  # Get grid size
all_config = generate_configurations(chem, ion_count, grid_size)

# Print sample configurations (limit output size)
print("Sample configurations:")
for config in all_config[:3]:  
    print(config)

# Print total number of configurations
print(f"\nTotal Configurations: {len(all_config)}")

