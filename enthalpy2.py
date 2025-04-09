import constants5 as const
import scipy as sp
import numpy as np

useApprox = False
if useApprox:
    filename = 'sim_data/massbalance_with_approx.txt'
else:
    filename = 'sim_data/massbalance_no_approx.txt'

m = {}  # Sett feks m3 = m["m3"]
w = {}  # Sett feks wc5 = w["wc5"]

with open(filename, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith("m"):
            key,value = line.split(":")
            m[key] = float(value.split("kg/s")[0].strip())
        elif line.startswith("w"):
            key, value = line.split(":")
            w[key.strip()]= round(float(value.split("%")[0].strip())/100.0,2)

#Ren H20
def cp_h(T):
    T = T + 273.15 # K
    return const.Aw + const.Bw*T + const.Cw*T**2

#Ren MEA
def cp_m(T):
    return const.Aa + const.Ba*T + const.Ca*T**2

#CO2 i MEA-løsning
def cp_c(T):
    return const.Ac + const.Bc*T

#MEA-løsning (aq)
def cp_sol(T, is_rich):
    if is_rich:
        key = "wm4"
    else:
        key = "wm3"
    wm = w[key]
    return (1-wm)*cp_h(T)+wm*cp_m(T)+wm*(1-wm)*(const.As + const.Bs*T + const.Cs*wm*(T-273.15)**(-1.5859))
def choose_cp(T, key):
    key = key[0:2]
    if key == "wc":
        return cp_c(T)
    elif key == "wh":
        return cp_h(T)
    elif key == "wm":
        return cp_sol(T, True)
    elif key == "wn":
        return cp_c(T)
    elif key == "wo":
        return cp_c(T)
    else:
        print("Feil i strømnummer")
        return np.nan
def cp_i(T, i):
    if i<1 or i>9:
        print("Feil i strømnummer")
        return np.nan
    elif i in [1,2,8,9]:
        return 
    else:
        cp = 0
        keys = [f"wc{i}", f"wm{i}"]
        for key in keys:
            cp += w[key]*cp_sol(T, key=="wm4")
        

#Tar T7 inn som Celsiusgrader
def DeltaT_lm_V1(T7):
    DeltaT1 = const.T6 - const.T5
    DeltaT2 = T7 - const.T4
    if DeltaT2 == 0:
        print("DeltaT2 = 0 !!")
        return np.nan
    if 1/1.4 < DeltaT1/DeltaT2 < 1.4:
        return 0.5 * (DeltaT1 + DeltaT2)
    else:
        return (DeltaT1-DeltaT2) / np.log(DeltaT1/DeltaT2)

def V1(vars):
    T7 = vars[0]
    m3 = m["m3"]
    m4 = m["m4"]
    m6 = m["m6"]
    m8 = m["m8"]
    wc3 = w["wc3"]
    wc4 = w["wc4"]
    wc6 = w["wc6"]
    wc8 = w["wc8"]


    deltaT_lm = DeltaT_lm_V1(T7)
return "something"