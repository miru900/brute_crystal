def add_extra_Cons(m, Vars):
    m.addConstr(Vars[0][0] == 1)
    m.addConstr(Vars[0][21] == 1)
    m.addConstr(Vars[0][33] == 1)
    m.addConstr(Vars[0][63] == 1)

def NCM_cons_2(m, Vars):
    Li = [i for i in range(12, 48)]
    for i in Li:
        m.addConstr(Vars[0][i] == 0)

    Ni = [i for i in range(12, 24)]  ;  Ni_remove = [13, 17]
    for i in Ni_remove:
        Ni.remove(i)
    for i in Ni:
        m.addConstr(Vars[2][i] == 1)

    Co = [13]  ;  Mn = [17]
    for i in Co:
        m.addConstr(Vars[3][i] == 1)
    for i in Mn:
        m.addConstr(Vars[4][i] == 1)

    O = [i for i in range(24, 48)]
    for i in O:
        m.addConstr(Vars[1][i] == 1)

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

def NCM_fcons8(m, Vars):
    # Li sites
    Li = [i for i in range(24, 96)]
    for i in Li:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[0][i + shift] == 0)
    
    # NCM forbidden site
    NCM1 = [i for i in range(0, 24)]
    NCM2 = [i for i in range(48, 96)]
    NCM = NCM1 + NCM2
    for i in NCM:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[2][i + shift] == 0)
            m.addConstr(Vars[3][i + shift] == 0)
            m.addConstr(Vars[4][i + shift] == 0)

    # O sites
    O = [i for i in range(48, 96)]
    for i in O:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[1][i + shift] == 1)

def NCM_pcons8(m, Vars):
    # Li sites
    # 1/2 : [2, 3, 6, 7, 8, 10, 12, 14, 15, 18, 19, 22]
    # 1/4 : [2, 3, 6, 14, 15, 19]
    # 1/8 : [3, 6, 15]
    Li = [3, 6, 15]
    for i in Li:
        for shift in [0, 96, 192, 288, 384, 480, 576, 672]:
            m.addConstr(Vars[0][i + shift] == 1)

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

def NCM_hcons8(m, Vars):
    # Li sites
    # 1/2 : [2, 3, 6, 7, 8, 10, 12, 14, 15, 18, 19, 22]
    # 1/4 : [2, 3, 6, 14, 15, 19]
    # 1/8 : [3, 6, 15]

    Li_except = [i for i in range(24, 96)]
    for i in Li_except:
        for shift in [96 * j for j in range(8)]:
            m.addConstr(Vars[0][i + shift] == 0)

    Li = [2, 3, 6, 14, 15, 19]
    Li_cons = 0
    for i in Li:
        for shift in [96 * j for j in range(8)]:
            Li_cons += Vars[0][i + shift]
    m.addConstr(Li_cons == 41)

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

def NCM_pcons27(m, Vars):
    # Li sites
    # Li sites
    # 1/2 : [2, 3, 6, 7, 8, 10, 12, 14, 15, 18, 19, 22]
    # 1/4 : [2, 3, 6, 14, 15, 19]
    # 1/8 : [3, 6, 15]
    Li = [3, 6, 15]
    for i in Li:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[0][i + shift] == 1)

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

def NCM_hcons27(m, Vars):
    # Li sites
    # Li sites
    # 1/2 : [2, 3, 6, 7, 8, 10, 12, 14, 15, 18, 19, 22]
    # 1/4 : [2, 3, 6, 14, 15, 19]
    # 1/8 : [3, 6, 15]

    Li_except = [i for i in range(24, 96)]
    for i in Li_except:
        for shift in [96 * j for j in range(27)]:
            m.addConstr(Vars[0][i + shift] == 0)

    Li = [2, 3, 6, 14, 15, 19]
    Li_cons = 0
    for i in Li:
        for shift in [96 * j for j in range(27)]:
            Li_cons += Vars[0][i + shift]
    m.addConstr(Li_cons == 153)


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



"""
def NCM_fcons8(m, Vars):
    # Define site indices
    Li_sites = [i for i in range(24, 96)]
    O_sites = [i for i in range(48, 96)]

    # Ensure Li and O sites are fixed
    for i in Li_sites:
        for shift in range(0, 673, 96):  # Shifts across layers
            m.addConstr(Vars[0][i + shift] == 0)  # Li must not be occupied

    for i in O_sites:
        for shift in range(0, 673, 96):
            m.addConstr(Vars[1][i + shift] == 1)  # O must be present

    # Define all possible sites except Li and O
    all_sites = set(range(24, 96))
    non_Li_O_sites = list(all_sites - set(O_sites))  # Allowed for Ni, Co, Mn

    # Ensure Ni, Co, and Mn are placed in valid sites
    for shift in range(0, 673, 96):
        for i in non_Li_O_sites:
            m.addConstr(Vars[2][i + shift] + Vars[3][i + shift] + Vars[4][i + shift] <= 1)  # Only one per site
"""

