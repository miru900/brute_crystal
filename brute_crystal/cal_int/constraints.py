def add_extra_Cons(m, Vars):
    m.addConstr(Vars[0][0] == 1)
    m.addConstr(Vars[0][21] == 1)
    m.addConstr(Vars[0][33] == 1)
    m.addConstr(Vars[0][63] == 1)

def NCM_cons(m, Vars):
    Li = [i for i in range(24, 96)]
    for i in Li:
        m.addConstr(Vars[0][i] == 0)

    Ni = [i for i in range(24, 48)]  ;  Ni_remove = [25, 26, 30, 41]
    for i in Ni_remove:
        Ni.remove(i)
    for i in Ni:
        m.addConstr(Vars[2][i] == 1)

    Co = [25, 26]  ;  Mn = [30, 41]
    for i in Co:
        m.addConstr(Vars[3][i] == 1)
    for i in Mn:
        m.addConstr(Vars[4][i] == 1)

    O = [i for i in range(48, 96)]
    for i in O:
        m.addConstr(Vars[1][i] == 1)

def NCM_cons8(m, Vars):
    # Li sites
    Li = [i for i in range(24, 96)]
    for i in Li:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[0][i + shift] == 0)
    
    # Ni sites
    Ni = [i for i in range(24, 48)]
    Ni_remove = [25, 26, 30, 41]
    for i in Ni_remove:
        Ni.remove(i)
    for i in Ni:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[2][i + shift] == 1)
    
    # Co and Mn sites
    Co = [25, 26]
    Mn = [30, 41]
    for i in Co:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[3][i + shift] == 1)
    for i in Mn:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[4][i + shift] == 1)
    
    # O sites
    O = [i for i in range(48, 96)]
    for i in O:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[1][i + shift] == 1)

def NCM_cons27(m, Vars):
    # Li sites
    Li = [i for i in range(24, 96)]
    for i in Li:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[0][i + shift] == 0)
    
    # Ni sites
    Ni = [i for i in range(24, 48)]
    Ni_remove = [25, 26, 30, 41]
    for i in Ni_remove:
        Ni.remove(i)
    for i in Ni:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[2][i + shift] == 1)
    
    # Co and Mn sites
    Co = [25, 26]
    Mn = [30, 41]
    for i in Co:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[3][i + shift] == 1)
    for i in Mn:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[4][i + shift] == 1)
    
    # O sites
    O = [i for i in range(48, 96)]
    for i in O:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[1][i + shift] == 1)
