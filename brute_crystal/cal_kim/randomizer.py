import numpy as np
import time

def gen_config(chem, ion_count, grid_size, sample_size=1000000):
    """
    Generate configurations efficiently without storing all permutations.

    Args:
        chem (list): List of element symbols.
        ion_count (list): Count of each element.
        grid_size (int): Total grid size.
        sample_size (int): Number of random configurations to generate.

    Yields:
        tuple: One configuration at a time.
    """

    # Create an array with all elements
    element_list = np.concatenate([np.full(ion_count[i], chem[i], dtype=object) for i in range(len(chem))])

    # Calculate empty slots
    total_slots = grid_size ** 3
    empty_slots = total_slots - sum(ion_count)
    element_list = np.concatenate([element_list, np.full(empty_slots, "_", dtype=object)])

    # Generate random configurations
    for _ in range(sample_size):
        np.random.shuffle(element_list)  # Shuffle in-place for efficiency
        yield tuple(element_list)  # Yield each configuration

if __name__ == "__main__":
    # Define elements and grid size
    start = time.time()

    chem = ["Pt"]
    ion_count = [4]
    grid_size = 4 

    print("Generating Configurations ...")

    # Generate configurations using generator
    test_con = gen_config(chem, ion_count, grid_size, sample_size=1000000)

    # Print the first 5 configurations as a sample
    print("Sample configurations:")
    i = 0
    for _ in test_con:
        i += 1

    end = time.time()
    print(f"Calculation time : {end - start} seconds")
