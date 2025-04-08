import constants5
import numpy as np
#from energibalanse import C_p_CO2
from heatcapacity import C_p_CO2
import MassBalance as mb

#Gjerne endre hvis det er en bedre måte å integrere på

def numerisk_integrasjon(f, T_start, T_slutt, n=1000):
    T = np.linspace(T_start, T_slutt, n)
    y = np.array([f(t) for t in T])
    return np.trapezoid(y, T)

def ett_steg_kompresjon():
    T_inn = constants5.T[8] + 273.15  # [K]
    p_inn = constants5.p[8]
    p_b = 20  # [bar]
    T_ut = 303  # [K]
    p_ut = 20  # [bar]

    m9 =  mb.m9

    cp_CO2 = 37  # [SI, J K mol^-1]
    k = constants5.gasConst / cp_CO2  # R/cp_CO2
    T_b = T_inn * (p_b / p_inn) ** k

    C_p_CO2_int = numerisk_integrasjon(C_p_CO2, T_inn, T_b)
    ws_rev = m9 * C_p_CO2_int * (T_b - T_inn)  # [KJ/s]
    ws = ws_rev / constants5.eta

    Q = m9 * C_p_CO2_int * (T_b - T_ut)  # [KJ/s]

    print(f'T_b: {T_b:.3f} K')
    print(f'ws: {ws:.3f} KJ/s')
    print(f'Q: {Q:.3f} KJ/s')
    return ws, Q


def tre_steg_kompresjon():
    T_inn = constants5.T[8] + 273.15  # [K]
    p_inn = constants5.p[8]
    p_b = 4  # [bar]
    T_c = 303  # [K]
    p_c = 4
    p_d = 8
    T_e = 303  # [K]
    p_e = 8  # [bar]
    p_f = 20  # [bar]
    T_ut = 303  # [K]
    p_ut = 20  # [bar]

    m9 =  mb.m9

    cp_CO2 = 37  # [SI, J K mol^-1]
    k = constants5.gasConst / cp_CO2  # R/cp_CO2

    T_b = T_inn * (p_b / p_inn) ** k
    T_d = T_c * (p_d / p_c) ** k
    T_f = T_e * (p_f / p_e) ** k

    C_p_CO2_1 = numerisk_integrasjon(C_p_CO2, T_inn, T_b)
    C_p_CO2_2 = numerisk_integrasjon(C_p_CO2, T_c, T_d)
    C_p_CO2_3 = numerisk_integrasjon(C_p_CO2, T_e, T_f)

    ws1_rev = m9 * C_p_CO2_1 * (T_b - T_inn)
    ws2_rev = m9 * C_p_CO2_2 * (T_d - T_c)
    ws3_rev = m9 * C_p_CO2_3 * (T_f - T_e)
    ws = (ws1_rev + ws2_rev + ws3_rev) / constants5.eta

    Q1 = m9 * C_p_CO2_1 * (T_b - T_c)
    Q2 = m9 * C_p_CO2_2 * (T_d - T_e)
    Q3 = m9 * C_p_CO2_3 * (T_ut - T_f)

    Q = Q1 + Q2 + Q3

    print(f'T_b: {T_b:.3f} K')
    print(f'T_d: {T_d:.3f} K')
    print(f'T_f: {T_f:.3f} K')

    return ws, Q
