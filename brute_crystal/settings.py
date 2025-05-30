import math

settings = {
    "Pt_1": {
        "State": 0, "Symbol": "Pt_4", "grid": 4, "cell_size": 3.92, "view" : True,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range": 0.2
    },
    "Pt3Cu1_1": {
        "State": 0, "Symbol": "Pt_3+Cu_1", "grid": 6, "cell_size": 3.796, "view" : True,
        "calc": "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range": 0.2
    },
    "C_1": {
        "State": 0, "Symbol": "C_8", "grid": 8, "cell_size": 3.567, "view" : True,
        "calc": "DUNN_WenTadmor_2019v1_C__MO_584345505904_000",
        "vari_range" : 0.2
    },
    "Pt2Cu1Fe1": {
        "State" : 0, "Symbol": "Pt_2+Cu_1+Fe_1", "grid" : 2, "cell_size" : 3.92, "view" : True,
        "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range" : 0.2},
    "Pt1Cu1Fe1Mo1" : {
        "State" : 0, "Symbol" : "Pt_2+Cu_2+Fe_2+Mo_2", "grid" : 4, "cell_size" : 9.455, "view" : True,
        "calc" : "EAM_Dynamo_ZhouJohnsonWadley_2004_CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr__MO_870117231765_000",
        "vari_range" : 0.1},
    "NaCl_1" : {
        "State" : 0, "Symbol" : "Na_4+Cl_4", "grid" : 2, "cell_size" : 5.64, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None,
        "vari_range" : 0.3},
    "SrTiO3_1" : {
        "State" : 8, "Symbol" : "Sr_1+Ti_1+O_3", "grid" : 2, "cell_size" : 3.9, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None,
        "vari_range" : 0.3},
    "NCM_battery" : {
        "State" : 0, "Symbol" : "Li_3+O_24+Ni_10+Co_1+Mn_1", "grid" : [39, 33, 88], "cell_size" : [10.05628 * 0.5, 5.79685, 14.28257], "view" : True, # original : 24, 48, 20, 2, 2
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : None,
        "vari_range" : 0.2},
    "TiZrHfScMo" : {
        "State" : 0, "Symbol" : "Ti_25+Zr_25+Hf_25+Sc_25+Mo_25", "grid" : 10, "cell_size" : 16.625, "view" : True,
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003",
        "vari_range" : 0.2},
    "LiCo_cathode" : {
        "State" : 0, "Symbol" : "Li_3+Co_3+O_6", "grid" : [22, 52, 6], "cell_size" : [2.81126, 2.81126, 13.90946], "view" : True, # mp-22526 so C22526, for supercell, C22526 * supercell
        "calc" : "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003", "angle" : [90, 90, 120],
        "vari_range" : 0.2},
    "LiNiAl_Cathode" : {
        "State" : 0, "Symbol" : "Li_3+O_6+Ni_3", "grid" : [2019, 155, 18], "cell_size" : [2.86463, 2.86463, 14.23982], "view" : True, # Al3 doping, for bachelor's degree graduation paper
        "calc" : None, "angle" : [90, 90, 120], 
        "vari_range" : 0.2}
}

