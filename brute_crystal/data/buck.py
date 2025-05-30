# information about buckingham parameter
# buckingham potential : A * exp(B*r) - C/r**6
# buck = {"symbol" : {"ion_pair" : [ss, ss, ss], "A" : [11, 22, 33], "B" : [44, 55, 66], "C" : [77, 88, 99]}} -> input format.

buck = {
        "SrTiO" : { "ion_pair" : ["Sr+O", "Ti+O", "O+O"], "A" : [1952.39, 4590.7279, 1388.77],
            "B" : [0.33685, 0.261, 0.36262], "C" : [19.22, 0, 175], "D" : [0, 0, 0], "E" : [10.0, 10.0, 10.0]},

        "Pt" : {"ion_pair" : ["Pt+Pt"]},

        "PtCu" : {"ion_pair" : ["Pt+Pt", "Pt+Cu"]},

        "C" : {"ion_pair" : ["C+C"]},

        "NaCl" : {"ion_pair" : ["Na+Na", "Cl+Cl", "Na+Cl"]},

        "PtCuFeMo" : {"ion_pair" : ["Pt+Pt", "Pt+Cu", "Pt+Fe", "Pt+Mo", "Cu+Cu", "Cu+Fe", "Cu+Mo", "Fe+Fe", "Fe+Mo", "Mo+Mo"]},

        "TiZrHfScMo" : {"ion_pair" : ["Ti+Ti", "Ti+Zr", "Ti+Hf", "Ti+Sc", "Ti+Mo", "Zr+Zr", "Zr+Hf", "Zr+Sc", "Zr+Mo", "Hf+Hf",
            "Hf+Sc", "Hf+Mo", "Sc+Sc", "Sc+Mo", "Mo+Mo"]}

            }
