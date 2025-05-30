def LiCo(m, Vars):
    for i in range(3, 12):
        m.addConstr(Vars[0][i] == 0)
    for i in range(6, 12):
        m.addConstr(Vars[1][i] == 1)
    for i in range(0, 3):
        m.addConstr(Vars[2][i] == 0)
        m.addConstr(Vars[3][i] == 0)
        m.addConstr(Vars[4][i] == 0)

