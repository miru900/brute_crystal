import math

settings = {
    "NaCl_square" : {
        "State" : 0, "Symbol" : "Na_2+Cl_2", "grid" : [17, 17, 17], "cell_size" : [5.64, 5.64, 5.64], "view" : True,
        "calc" : None, "angle" : None,
        "vari_range" : None},
    "NaCl_cubic" : {
        "State" : 0, "Symbol" : "Na_4+Cl_4", "grid" : [2] * 3, "cell_size" : [5.64] * 3, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None, "cons" : None, "charge" : None,
        "vari_range" : 0.3},
    "SrTiO3_1" : {
        "State" : 0, "Symbol" : "Sr_1+Ti_1+O_3", "grid" : [2, 2, 2], "cell_size" : [3.9] * 3, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None, "charge" : None, "cons" : None,
        "vari_range" : 0.3},
    "NCM_battery" : {
        "State" : 5, "Symbol" : "Li_45+Ni_48+Co_6+Mn_6+O_120", "grid" : [811, 811, 811], "cell_size" : [14.27426, 11.35460, 14.02276], "view" : True, # Li 60 test, 
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" :[89.8172, 90.0286, 120.0661], "charge" : NCM811,
        "vari_range" : 0.2},
    "LiCo_cathode" : {
        "State" : 0, "Symbol" : "Li_12+Co_24+O_48", "grid" : [2025, 6, 5], "cell_size" : [2.81126 * 2, 2.81126 * 2, 13.90946 * 2], "view" : True, # mp-22526 so C22526, for supercell, C22526 * supercell
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : [90, 90, 120], "cons" : [list(range(0, 24)), list(range(24, 48)) + ["Fixed"], list(range(48, 96)) + ["Fixed"]], "charge" : None,
        "vari_range" : 0.2},
    "LiNiAl_Cathode" : {
        "State" : 0, "Symbol" : "Li_3+O_6+Ni_3", "grid" : [2019, 155, 18], "cell_size" : [2.86463, 2.86463, 14.23982], "view" : True, # Al3 doping, for bachelor's degree graduation paper
        "calc" : None, "angle" : [90, 90, 120], "cons" : None, "charge" : None,
        "vari_range" : 0.2}
}

