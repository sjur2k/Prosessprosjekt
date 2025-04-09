import constants5
import numpy as np
import MassBalance as mb


# Beregner Ws
def ws(T_inn, p_inn, p_b):
    T_2 = T2_merket(T_inn, p_inn, p_b)
    ws_rev = (constants5.cpg[0]*(T_2-T_inn))
    ws = ws_rev/constants5.eta
    return ws #KJ/KG

# Beregner T2'
def T2_merket(T_inn, p_inn, p_b):
    cp_co2 = 37  # [fra SI, J/K mol^-1]
    k = constants5.gasConst/cp_co2 
    T_b = T_inn*(p_b/p_inn)**k

    return T_b

# beregner den virkelige T2
def T2(ws, T1):
    cp_co2 = 0.868
    T2 = ws/(cp_co2)+T1 # kJ/kg
    return T2

def trestegskompresjon():
    T_inn = constants5.T[8] + 273.15 #[K]
    p_inn = constants5.p[8]
    p_b   = 4 #[bar]

    T_c = 303 #[K]
    p_c = 4 # [bar]
    p_d = 8 # [bar]

    T_e = 303 #[K]
    p_e = 8 # [bar]
    p_f = 20 # [bar]

    T_ut = 303  # [K]
    p_ut = 20  # [bar]

    ws1 = ws(T_inn, p_inn, p_b)
    ws2 = ws(T_c, p_c, p_d)
    ws3 = ws(T_e, p_e, p_f)

    T_b = T2(ws1, T_inn)
    T_d = T2(ws2, T_c)
    T_f = T2(ws3, T_e)

    print(f'T_b: {T_b:.3f} K')
    print(f'T_d: {T_d:.3f} K')
    print(f'T_f: {T_f:.3f} K')
    print(f'ws1: {ws1:.3f} KJ/KG')
    print(f'ws2: {ws2:.3f} KJ/KG')
    print(f'ws3: {ws3:.3f} KJ/KG')

    return ws1+ws2+ws3

m9 = mb.m9

enstegs = ws(constants5.T[8]+273.15, constants5.p[8], 20)
trestegs = trestegskompresjon()

print(f'Enstegskompresjon: {m9*enstegs:.3f} KJ/s')
print(f'Trestegskompresjon: {m9*trestegs:.3f} KJ/s')

