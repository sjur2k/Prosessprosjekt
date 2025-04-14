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
def mean_cp_h(Ti,Tf):
    if Ti<Tf:
        Ti,Tf = Tf,Ti
    Ti += 273.15 # K
    Tf += 273.15 # K
    return const.Aw + (const.Bw/2)*(Tf**2-Ti**2)/(Tf-Ti) + (const.Cw/3)*((Tf**3-Ti**3)/(Tf-Ti))

#Ren MEA
def cp_m(T):
    T = T + 273.15 # K
    return const.Aa + const.Ba*T + const.Ca*T**2

#Snitt for ren MEA
def mean_cp_m(Ti,Tf):
    if Ti<Tf:
        Ti,Tf = Tf,Ti
    Ti += 273.15 # K
    Tf += 273.15 # K
    return const.Aa + (const.Ba/2)*(Tf**2-Ti**2)/(Tf-Ti) + (const.Ca/3)*((Tf**3-Ti**3)/(Tf-Ti))

#CO2 i MEA-løsning
def cp_c(T):
    T = T + 273.15 # K
    return const.Ac + const.Bc*T

#Snitt for CO2 i MEA-løsning
def mean_cp_c(Ti,Tf):
    if Ti<Tf:
        Ti,Tf = Tf,Ti
    Ti += 273.15 # K
    Tf += 273.15 # K
    return const.Ac + (const.Bc/2)*(Tf**2-Ti**2)/(Tf-Ti)

#MEA-løsning
def cp_sol(T , is_rich): #([Celsius], Bool)
    if is_rich:
        key = "wm4"
    else:
        key = "wm3"
    wm = w[key]
    return (1-wm)*cp_h(T)+wm*cp_m(T)+wm*(1-wm)*(const.As + const.Bs*(T+273.15) + const.Cs*wm*(T)**(-1.5859))

#Snitt for MEA-løsning
def mean_cp_sol(Ti,Tf, is_rich):
    if Ti<Tf:
        Ti,Tf = Tf,Ti
    if is_rich:
        key = "wm4"
    else:
        key = "wm3"
    wm = w[key]
    part1 = (1-wm)*mean_cp_h(Ti,Tf) # Konverterer til K inni
    part2 = wm*mean_cp_m(Ti,Tf) # Konverterer til K inni
    Ti += 273.15 # K
    Tf += 273.15 # K
    part3 = wm*(1-wm)*(const.As + (const.Bs/2)*(Tf**2-Ti**2)/(Tf-Ti) + ((const.Cs*wm)/(-0.5859))*((Tf-273.15)**(-0.5859)-(Ti-273.15)**(-0.5859))/(Tf-Ti))
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
        if i in [4,5]:
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
        if i in [4,5]:
            is_rich = True
        cp += w[f"wc{i}"]*mean_cp_c(T1,T2)
        cp += (1-w[f"wc{i}"])*mean_cp_sol(T1,T2, is_rich)
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
    eq[1] = A  - (Q/(U*deltaT_lm))
    return eq

# Spesifikk entalpi i strøm $i$
def h_i(i):
    h=0
    T_ref = const.Tref_C
    T = const.T[i-1]
    if i < 1 or i > 9:
        print("Feil i strømnummer")
        return np.nan
    elif i in [1,2,8,9]:
        keys=["wc","wh","wn","wo"]
        for j in range(len(keys)):
            key = f"{keys[j]}{i}"
            if key in w:
                h += w[key]*const.hf[j] #Dannelsesentalpi
                h += w[key]*const.cpg[j]*(T-T_ref) #Entalpi temperaturøkning
    else:
        is_rich = False
        if i in [4,5]:
            is_rich = True
        h += (1-w[f"wc{i}"])*(const.hfsol+(T-T_ref)*mean_cp_sol(T_ref,T,is_rich)) #Entalpi for MEA-løsning
        h += w[f"wc{i}"]*(const.hf[0]+(T-T_ref)*mean_cp_c(T_ref,T)) #Entalpi for CO2 i MEA-løsning
        h += w[f"wc{i}"]*const.habs_m*(1000/const.Mw[0]) #Ekstra entalpi for absorpsjon av CO2
    return h

T7, A = sp.optimize.root(V1, [T6-15, 190]).x
const.T[6] = T7
m10 = abs((T3-T7)/(T11-T10)) * (mean_cp_i(T7,T3,7)/mean_cp_h(T10,T11)) * m["m7"]

Q_V3 = m["m8"]*(
    w["wc8"]*mean_cp_c(T8,T9)*(T9-T8)+
    w["wh8"]*(
        mean_cp_h(T8,T9)*(T9-T8) -
        const.dHvap*1000/const.Mw[1]
        )
    )
Q_V4 = (Q_V3 + h_i(5)*m["m5"])-(h_i(6)*m["m6"] + h_i(9)*m["m9"])


print(f"T7 = {round(T7,2)}°C")
print(f"A = {round(A,2)}m^2")
print(f"m10 = m11 = {round(m10,2)} kg/s")
for i in range(9):
    print(f"h_{i+1}: {round(h_i(i+1),2)} kJ/kg")
print(f"Q_V3 = {round(Q_V3,2)} kJ")
print(f"Q_V4 = {round(Q_V4,2)} kJ")
print(mean_cp_c(T8,T9))
print("T_økning strøm 4:",A*DeltaT_lm_V1(T7)*const.U/(1000*m["m4"]*mean_cp_i(T4,T5,4)))
print(mean_cp_i(T6,T7,6))