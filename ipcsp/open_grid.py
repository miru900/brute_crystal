def chem_compo(grid_num, symbol):
    CC = open(f"./data/C{grid_num}.txt")
    pos_data = []
    for line in CC:
        pos_data.append([float(i) for i in line.split()])

    chem = []
    ion_count = []
    indiv_symbol = symbol.split("+")
    compo = len(indiv_symbol)

    for i in range(compo):
        c = indiv_symbol[i]
        b = c.split("_")
        chem.append(b[0])
        ion_count.append(int(b[1]))
    return chem, ion_count, pos_data

# chem, ion_count = chem_compo(compo)
# print(chem, ion_count)
        
if __name__ == "__main__":
    test_num = 2
    test_symbol = "Sr_1+Ti_1+O_3"
    test_chem, test_ion = chem_compo(test_num, test_symbol)
    print(f"The Chemical Composition is {test_chem}")
    print(f"The Ion count is {test_ion}")
    


        
