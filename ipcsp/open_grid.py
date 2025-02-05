grid_num = int(input("Write grid Number : "))

CC = open(f"data/C{grid_num}.txt")
pos_data = []
for line in CC:
    pos_data.append([float(i) for i in line.split()])
# print(pos_data)

print("The example of Composition writing, if SrTiO3 : Sr, Ti, O => Number is 3.")
compo = int(input("Write Number of Composition : "))

def chem_compo(compo):
    print("Write example : SrTiO3 -> Sr_1, Ti_1, O_3")
    chem = []
    ion_count = []
    for i in range(compo):
        a = input("Write the Atom Symbol and ion count : ")
        b = a.split("_")
        chem.append(b[0])
        ion_count.append(int(b[1]))
    return chem, ion_count

# chem, ion_count = chem_compo(compo)
# print(chem, ion_count)
        
        


        
