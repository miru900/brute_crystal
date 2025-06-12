def LiCo(m, Vars):
    for i in range(24, 96):
        m.addConstr(Vars[0][i] == 0)
    for i in range(48, 96):
        m.addConstr(Vars[2][i] == 1)
    for i in range(24, 48):
        m.addConstr(Vars[1][i] == 1)

def NCM811(m, Vars):
    for i in range(60, 240):
        m.addConstr(Vars[0][i] == 0)
    for i in range(60, 108):
        m.addConstr(Vars[1][i] == 1)
    for i in range(108, 114):
        m.addConstr(Vars[2][i] == 1)
    for i in range(114, 120):
        m.addConstr(Vars[3][i] == 1)
    for i in range(120, 240):
        m.addConstr(Vars[4][i] == 1)

