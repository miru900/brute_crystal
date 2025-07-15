def LiCo(m, Vars, cons, ewal_mat):
    if cons:
        grid_size = ewal_mat.shape[0]
        for idx, con in enumerate(cons):
            if "Fixed" in cons:
                con = con[0 : -1]
                for i in con:
                    m.addConstr(Vars[idx][i] == 1)
            else:
                con_set = set(con)
                all_set = [i for i in range(grid_size)] ; all_set = set(all_set)
                fit_set = all_set - con_set ; fit_set = list(fit_set)
                for i in fit_set:
                    m.addConstr(Vars[idx][i] == 0)
                    

    else:
        for i in range(24, 96):
            m.addConstr(Vars[0][i] == 0)
        for i in range(24, 48):
            m.addConstr(Vars[1][i] == 1)
        for i in range(48, 96):
            m.addConstr(Vars[2][i] == 1)
