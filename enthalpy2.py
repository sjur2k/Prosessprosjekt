import constants5 as const
import scipy as sp
import numpy as np
### Constants for the simulation
T1,T2,T3,T4,T5,T6 = const.T[0:6]
T8,T9,T10,T11 = const.T[7:11]

useApprox = True
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

#Ren H20(g)
def cp_h(T):
    T = T + 273.15 # K
    return const.Aw + const.Bw*T + const.Cw*T**2

#Snitt for ren H20(g)
def mean_cp_h(T1,T2):
    T1 = T1 + 273.15 # K
    T2 = T2 + 273.15 # K
    return const.Aw + (const.Bw/2)*(T2**2-T1**2)/(T2-T1) + (const.Cw/3)*((T2**3-T1**3)/(T2-T1))

#Ren MEA
def cp_m(T):
    T = T + 273.15 # K
    return const.Aa + const.Ba*T + const.Ca*T**2

#Snitt for ren MEA
def mean_cp_m(T1,T2):
    T1 = T1 + 273.15 # K
    T2 = T2 + 273.15 # K
    return const.Aa + (const.Ba/2)*(T2**2-T1**2)/(T2-T1) + (const.Ca/3)*((T2**3-T1**3)/(T2-T1))

#CO2 i MEA-løsning
def cp_c(T):
    T = T + 273.15 # K
    return const.Ac + const.Bc*T

#Snitt for CO2 i MEA-løsning
def mean_cp_c(T1,T2):
    T1 = T1 + 273.15 # K
    T2 = T2 + 273.15 # K
    return const.Ac + (const.Bc/2)*(T2**2-T1**2)/(T2-T1)

#MEA-løsning
def cp_sol(T , is_rich): #([Celsius], Bool)
    if is_rich:
        key = "wm4"
    else:
        key = "wm3"
    wm = w[key]
    return (1-wm)*cp_h(T)+wm*cp_m(T)+wm*(1-wm)*(const.As + const.Bs*T + const.Cs*wm*(T)**(-1.5859))

#Snitt for MEA-løsning
def mean_cp_sol(T1,T2, is_rich):
    if is_rich:
        key = "wm4"
    else:
        key = "wm3"
    wm = w[key]
    part1 = (1-wm)*mean_cp_h(T1,T2) # Konverterer til K inni
    part2 = wm*mean_cp_m(T1,T2) # Konverterer til K inni
    T1 = T1 + 273.15 # K
    T2 = T2 + 273.15 # K
    part3 = wm*(1-wm)*(const.As + (const.Bs/2)*(T2**2-T1**2)/(T2-T1) + ((const.Cs*wm)/(-0.5859))*((T2-273.15)**(-0.5859)-(T1-273.15)**(-0.5859))/(T2-T1))
    return part1 + part2 + part3


#Varmekapasitet kJ/(kg K) for strøm i
def cp_i(T, i): #([Celsius], int)
    cp = 0
    if i<1 or i>9:
        print("Feil i strømnummer")
        return np.nan
    elif i in [1,2,8,9]:
        keys=["wc","wh","wn","wo"]
        for j in range(len(keys)):
            key = f"{keys[j]}{i}"
            if key in w:
                cp += w[key]*const.cpg[j] # cpg[j] er omtrent lik i temperaturområdet 45C-107C
    else:
        is_rich = False
        if i in [3,6,7]:
            is_rich = True
        cp += w[f"wc{i}"]*cp_c(T)
        cp += w[f"wm{i}"]*cp_sol(T, is_rich)
    return cp 

#Snitt varmekapasitet mellom T1 og T2 for strøm i
def mean_cp_i(T1,T2,i): #([Celsius], [Celsius], int)
    if i<1 or i>9:
        print("Feil i strømnummer")
        return np.nan
    if i in [1,2,8,9]:
        return cp_i(T1,i) # cpg[j] er omtrent lik i temperaturområdet 45C-107C
    else:
        cp = 0
        is_rich = False
        if i in [3,6,7]:
            is_rich = True
        cp += w[f"wc{i}"]*mean_cp_c(T1,T2)
        cp += w[f"wm{i}"]*mean_cp_sol(T1,T2, is_rich)
        return cp
    
def DeltaT_lm_V1(T7): #([Celsius])
    DeltaT1 = T6 - T5
    DeltaT2 = T7 - T4
    if DeltaT1==0 or DeltaT2 == 0:
        print("DeltaT2 = 0 !!")
        return np.nan
    if DeltaT1/DeltaT2<=0:
        print(f"Invalid: DT2 = T7-T4 = {T7}-{T4} = {DeltaT2} <0")
        return np.nan
    elif 1/1.4 < DeltaT1/DeltaT2 < 1.4:
        return 0.5 * (DeltaT1 + DeltaT2)
    else:
        return (DeltaT1-DeltaT2) / np.log(DeltaT1/DeltaT2)

def V1(vars):
    T7, A = vars

    m4 = m["m4"]
    m6 = m["m6"]
    mean_cp_hot = float(mean_cp_i(T6,T7,6))
    mean_cp_cold = float(mean_cp_i(T4,T5,4))
    
    Q = m4*mean_cp_cold*(T5-T4) #kW
    deltaT_lm = DeltaT_lm_V1(T7) #K
    U = const.U/1000 #kW/(m^2 K)
    
    eq=[0]*2
    eq[0] = T7 - (T6-Q/(m6*mean_cp_hot))
    eq[1] = A  - (Q/(const.U*deltaT_lm))
    return eq

sol = sp.optimize.root(V1, [T6-15, 190]).x
print(sol)
print(mean_cp_i(T5,T6,4))