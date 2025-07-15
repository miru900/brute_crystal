import mods.mod_1 as m1
import mods.mod_2 as m2
import mods.mod_3 as m3
import mods.mod_4 as m4
import mods.mod_5 as m5
import mods.mod_6 as m6
import mods.mod_7 as m7
import mods.mod_8 as m8
import mods.mod_9 as m9
import mods.mod_10 as m10
from settings import settings

angstrom_symbol = "\u212B"
mod_list = [m1.mod_cell_one, m2.mod_cell_vari, m3.mod_cell_integer, m4.mod_cell_int_periodic, m5.mod_cell_int_periodic2, m6.mod_cell_int_periodic3, m7.mod_cell_int_ase_kim, m8.qiskit_simulator, m9.fermion_simulator,
        m10.pyscf_simulator]

iter_num = 1
sample_size = 1

# print(settings.keys())

if __name__ == "__main__":
    for key in settings.keys():
        mod = settings[key]["State"]

        if mod == 0:
            continue

        selected_mod = mod_list[mod - 1]
        
        if mod == 1:
            for i in range(iter_num):
                selected_mod(key, settings, sample_size)

        elif mod == 2:
            selected_mod(key, settings, sample_size)

        else:
            selected_mod(key, settings)
