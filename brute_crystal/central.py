import mods.mod_1 as m1
import mods.mod_2 as m2
import mods.mod_3 as m3

angstrom_symbol = "\u212B"
sample_size = 3812256
mod_list = [m1.mod_cell_one, m2.mod_cell_vari, m3.mod_cell_integer]
iter_num = 1

settings = {
    "Pt_1": {
        "State": 0, "Symbol": "Pt_4", "grid": 4, "cell_size": 3.92, "view" : True,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range": 0.2
    },
    "Pt3Cu1_1": {
        "State": 0, "Symbol": "Pt_3+Cu_1", "grid": 2, "cell_size": 3.92, "view" : True,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range": 0.2
    },
    "C_1": {
        "State": 0, "Symbol": "C_8", "grid": 4, "cell_size": 3.567, "view" : True,
        "calc": "DUNN_WenTadmor_2019v1_C__MO_584345505904_000", 
        "vari_range" : 0.2
    },
    "Pt2Cu1Fe1": {
        "State" : 0, "Symbol": "Pt_2+Cu_1+Fe_1", "grid" : 2, "cell_size" : 3.92, "view" : True,
        "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range" : 0.2},
    "Pt1Cu1Fe1Mo1" : {
        "State" : 0, "Symbol" : "Pt_1+Cu_1+Fe_1+Mo_1", "grid" : 2, "cell_size" :3.92, "view" : True,
        "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range" : 0.3},
    "NaCl_1" : {
        "State" : 3, "Symbol" : "Na_4+Cl_4", "grid" : 2, "cell_size" : 5.64, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",
        "vari_range" : 0.3},
    "SrTiO3_1" : {
        "State" : 0, "Symbol" : "Sr_1+Ti_1+O_3", "grid" : 4, "cell_size" : 3.9, "view" : True,
        "calc" : "Not Determined hehe",
        "vari_range" : 0.3}
}

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

        if mod == 2:
            selected_mod(key, settings, sample_size)

        if mod == 3:
            selected_mod(key, settings)
