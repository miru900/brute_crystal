import math

settings = {
    "NaCl_square" : {
        "State" : 0, "Symbol" : "Na_2+Cl_2", "grid" : [17, 17, 17], "cell_size" : [5.64, 5.64, 5.64], "view" : True,
        "calc" : None, "angle" : None,
        "vari_range" : None},
    "NaCl_cubic" : {
        "State" : 8, "Symbol" : "Na_2+Cl_2", "grid" : [2, 2, 2], "cell_size" : [5.64] * 3, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None, "cons" : None,
        "vari_range" : 0.3},
    "SrTiO3_1" : {
        "State" : 0, "Symbol" : "Sr_4+O_4", "grid" : [2, 2, 2], "cell_size" : [3.9] * 3, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None,
        "vari_range" : 0.3},
    "NCM_battery" : {
        "State" : 0, "Symbol" : "Li_3+O_24+Ni_10+Co_1+Mn_1", "grid" : [39, 33, 88], "cell_size" : [10.05628 * 0.5, 5.79685, 14.28257], "view" : True, # original : 24, 48, 20, 2, 2
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None,
        "vari_range" : 0.2},
    "LiCo_cathode" : {
        "State" : 0, "Symbol" : "Li_3+Co_3+O_6", "grid" : [22, 52, 6], "cell_size" : [2.81126, 2.81126, 13.90946], "view" : True, # mp-22526 so C22526, for supercell, C22526 * supercell
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : [90, 90, 120], "cons" : [[0, 1, 2], [3, 4, 5], [6, 7, 8, 9, 10, 11]],
        "vari_range" : 0.2},
    "LiNiAl_Cathode" : {
        "State" : 0, "Symbol" : "Li_3+O_6+Ni_3", "grid" : [2019, 155, 18], "cell_size" : [2.86463, 2.86463, 14.23982], "view" : True, # Al3 doping, for bachelor's degree graduation paper
        "calc" : None, "angle" : [90, 90, 120], 
        "vari_range" : 0.2}
}

