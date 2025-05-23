import constants5
from scipy.optimize import root
import numpy as np
import MassBalance as mb

useApprox = False
if useApprox:
    filename = 'sim_data/massbalance_with_approx.txt'
else:
    filename = 'sim_data/massbalance_no_approx.txt'

massestrømmer = {}  # Sett feks m3 = massestrømmer["m3"]
massefraksjoner = {}  # Sett feks wc5 = massefraksjoner["wc5"]

with open(filename, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith("m"):
            key,value = line.split(":")
            massestrømmer[key] = float(value.split("kg/s")[0].strip())
        elif line.startswith("w"):
            key, value = line.split(":")
            massefraksjoner[key.strip()]= round(float(value.split("%")[0].strip())/100.0,2)

print("Massestrømmer:\n",massestrømmer,"\n")
print("Massefraksjoner:\n",massefraksjoner)

#numerisk integrasjon
def numerisk_integrasjon(funksjon, T_start, T_slutt, n=1000):
    T_range = np.linspace(T_start, T_slutt, n)
    f_vals = funksjon(T_range)
    return np.trapezoid(f_vals, T_range)

#Varmekapasiteter

#H2O
def C_p_H2O(T):
    T = T + 273.15 # K
    return constants5.Aw + constants5.Bw*T + constants5.Cw*T**2  # kJ/kg

print("C_p_H2O(25) = ", C_p_H2O(25)) # kJ/kg

#Ren MEA
def C_p_MEA(T):
    return constants5.Aa + constants5.Ba*T + constants5.Ca*T**2  # kJ/kg

#CO2 i MEA-løsning
def C_p_CO2(T):
    return constants5.Ac + constants5.Bc*T # kJ/kg

# varmekapasitet for MEA-løsning
def C_p_sol(T):
    cp_H2O = C_p_H2O(T)
    cp_MEA = C_p_MEA(T)
    w_MEA  = constants5.waMEA
    return (1-w_MEA)*cp_H2O + w_MEA*cp_MEA + w_MEA*(1-w_MEA)*(constants5.As + constants5.Bs*T +constants5.Cs*w_MEA*(T-273.15)**-1.5859) # kJ/kg

#DeltaT LM:
def deltaT_LM(TH_inn,TH_ut,TC_inn,TC_ut):
    Delta_T1 = TH_inn - TC_ut
    Delta_T2 = TH_ut - TC_inn
    return (Delta_T2-Delta_T1)/(np.log(Delta_T1/Delta_T2))


def Varme(m_C, Cp_C, TC_ut, TC_inn):
    return m_C * Cp_C * (TC_ut - TC_inn)

#Utregning av A og T7 for varmeveksler
def V_1(unknowns):
    U = constants5.U
    m6 = mb.m6  #mass flow rate 6 (kg/s)
    m4 = mb.m4 #mass flow rate 4 (kg/s)
    wm4 = mb.wm4  #mass fraction of MEA in the solution
    wc4 = 0.10  #mass fraction of CO2 in the solution
    T4 = constants5.T[3] + 273.15 # K
    T5 = constants5.T[4] + 273.15 # K
    T6 = constants5.T[5] + 273.15 # K


    T7 = unknowns[0]
    A = unknowns[1]

    C_p_c = wm4 *numerisk_integrasjon(C_p_sol,T5, T4) + wc4 * numerisk_integrasjon(C_p_CO2, T5, T4)
    #Q = Varme(m4, C_p_c, T5 ,T4 ) # kJ/s mulig denne gir logisk feil siden jeg har gange eq2 med -1
    Q = m4 * C_p_c * (T5 - T4) # kJ/s
    Tdlm = deltaT_LM(T6, T7, T4,T5) 
    

    eq1 = (U * A * Tdlm) - (Q*1000) #J/s
    eq2 = m6 * numerisk_integrasjon(C_p_sol, T4, T7) * 1000 *(T6-T7) - Q * 1000 #Gange med tusen for å få til J/s

    return [eq1,eq2]


A = 200
T7 = 300

guess = [A, T7]
sol = root(V_1, guess)

A = sol[0]
T7 = sol[1]

print(f'A = {A:.3f} m^2')
print(f'T7 = {T7:.3f} K')


def V_2():
    T3 = constants5.T[2] + 273.15 # K
    T10 = constants5.T[9] + 273.15 # K
    T11 = constants5.T[10] + 273.15 # K
    mh = mb.m6 # kg/s

    Cp_H = numerisk_integrasjon(C_p_sol,T7,T3)[0]
    Q = Varme(mh, Cp_H, T3, T7) # kJ/s

    Cp_C = numerisk_integrasjon(C_p_H2O,T10,T11)[0]
    mc = -Q/(Cp_C*(T11-T10)) # kg/s

    return mc,Q

mc = V_2()[0]
Q_V_2 = V_2()[1]

print(f'Massestrøm kald strøm: {mc:.3f} kg/s')
print(f'Totalt energiforbruk fra varmeveksler V_2: {Q_V_2:.3f} kJ/s') #Denne spør de egentlig ikke om, men tar med selv om

def V_3():
    m8 = mb.m8
    wc8 = mb.wc8 
    wh8 = mb.wh8
    M_H2O = constants5.Mw[1] # kg/mol
    T8 = constants5.T[7] # K
    T9 = constants5.T[8] # K
    Cp_CO2 = constants5.cpg[0]
    Cp_H2O = constants5.cpg[1]
    dHfus = constants5.dHfus # kJ/mol
    dHsub = constants5.dHsub # kJ/mol


    q1 = wc8 * Cp_CO2 * (T8 - T9) # kJ/s
    q2 = wh8 * Cp_H2O * (T8 - T9) # kJ/s
    q3 = (wh8 * (dHfus - dHsub)*10e3)/M_H2O # kJ/s

    Q = m8(q1 + q2 + q3) # kJ/s
    return Q


Q_V_3 = V_3()
print(f'Totalt energiforbruk fra varmeveksler V-3: {Q_V_3:.3f} kJ/s')

#Entalpi:

def h_strøm1():
    T1 = constants5.T[0] # Celcius
    Tref = constants5.Tref_C  # Celsius

    w_CO2 = constants5.wCO2 # kg/mol
    w_H2O = constants5.wH2O # kg/mol
    w_N2 = constants5.wN2 # kg/mol
    w_O2 = constants5.wO2 # kg/mol

    Cp_CO2 = constants5.cpg[0] # kJ/kg
    Cp_H2O = constants5.cpg[1] # kJ/kg
    Cp_N2 = constants5.cpg[2] # kJ/kg
    Cp_O2 = constants5.cpg[3] # kJ/kg

    hf_CO2 = w_CO2 * constants5.hf[0]
    hf_H2O = w_H2O * constants5.hf[1]
    #std dannelse entalpi for N2 og O2 = 0

    dH_CO2 = hf_CO2 + w_CO2 * Cp_CO2 * (T1 - Tref) # kJ/kg
    dH_H2O =hf_H2O + w_H2O * Cp_H2O * (T1 - Tref) # kJ/kg
    dH_N2 = w_N2 * Cp_N2 * (T1 - Tref) # kJ/kg
    dH_O2 = w_O2 * Cp_O2 * (T1 - Tref) # kJ/kg

    h1 = dH_CO2 + dH_H2O + dH_N2 + dH_O2 # kJ/kg
    return h1*constants5.m1 #kJ

h1 = h_strøm1()
print(f'Entalpi for strøm 1: {h1:.3f} kJ')

def h_strøm2():
    T2 = constants5.T[1] # Celsius
    Tref = constants5.Tref_C  # Celsius

    w_CO2 =mb.wc2
    w_H2O = mb.wh2
    w_N2 = mb.wn2
    w_O2 = mb.no2

    Cp_CO2 = constants5.cpg[0] # kJ/kg
    Cp_H2O = constants5.cpg[1] # kJ/kg
    Cp_N2 = constants5.cpg[2] # kJ/kg
    Cp_O2 = constants5.cpg[3] # kJ/kg

    hf_CO2 =w_CO2 * constants5.hf[0]
    hf_H2O = w_H2O * constants5.hf[1]
    #std dannelse entalpi for N2 og O2 = 0


    dH_CO2 = hf_CO2 + w_CO2 * Cp_CO2 * (T2 - Tref) # kJ/mol
    dH_H2O =hf_H2O + w_H2O * Cp_H2O * (T2 - Tref) # kJ/mol
    dH_N2 = w_N2 * Cp_N2 * (T2 - Tref) # kJ/mol
    dH_O2 = w_O2 * Cp_O2 * (T2 - Tref) # kJ/mol

    h2 = dH_CO2 + dH_H2O + dH_N2 + dH_O2 # kJ/mol
    return h2*mb.m2 #kJ

h2 = h_strøm2()
print(f'Entalpi for strøm 2: {h2:.3f} kJ')

# h(sol+CO2)(T)
def h_sol_co2_strøm3(T): #Blir den samme for 3,6 og 7

    w_CO2 = mb.wc3
    w_H2O = mb.wh3
    w_MEA = mb.wm3

    d_abs_h_co2 = (constants5.habs_m*1000)/constants5.Mw[0] #Gjør om fra KJ/mol til kJ/kg

    Cp_CO2 = C_p_CO2 # kJ/kg
    Cp_sol = C_p_sol # kJ/kg
    Cp_H2O = C_p_H2O # kJ/kg
    Cp_MEA = constants5.hfsol # kJ/kg

    hf_CO2 =w_CO2 * constants5.hf[0]
    hf_H2O = w_H2O * constants5.hf[1]
    hf_MEA = w_MEA * constants5.hfsol # kJ/kg

    integrert_sol = numerisk_integrasjon(C_p_sol, constants5.Tref_K, T)
    integrert_co2 = numerisk_integrasjon(C_p_CO2, constants5.Tref_K, T)


    dH_CO2 = hf_CO2 + w_CO2 * integrert_co2[0] #KJ/kg
    dH_MEA = hf_MEA + w_MEA * integrert_sol[0] # kJ/kg
    dH_H2O = w_CO2 * d_abs_h_co2

    h3 = dH_CO2 + dH_H2O + dH_MEA# kJ/mol
    return h3*mb.m3

h_sol_co2_strøm3 = h_sol_co2_strøm3(constants5.T[2]+273)
h_sol_co2_strøm6 = h_sol_co2_strøm3(constants5.T[5]+273)
h_sol_co2_strøm7 = h_sol_co2_strøm3(constants5.T[6]+273)
print(f'Entalpi for strøm 3: {h_sol_co2_strøm3:.3f} kJ')
print(f'Entalpi for strøm 6: {h_sol_co2_strøm6:.3f} kJ')
print(f'Entalpi for strøm 7: {h_sol_co2_strøm7:.3f} kJ')

#h(sol+CO2)(T)
def h_sol_co2_strøm4(T): #Blir den samme for 4 og 5
    w_CO2 = mb.wc4
    w_H2O = mb.wh4
    w_MEA = mb.wh4

    d_abs_h_co2 = (constants5.habs_m*1000)/constants5.Mw[0] #Gjør om fra KJ/mol til kJ/kg

    Cp_CO2 = C_p_CO2 # kJ/kg
    Cp_sol = C_p_sol # kJ/kg
    Cp_H2O = C_p_H2O # kJ/kg
    Cp_MEA = constants5.hfsol # kJ/kg

    hf_CO2 =w_CO2 * constants5.hf[0]
    hf_H2O = w_H2O * constants5.hf[1]
    hf_MEA = w_MEA * constants5.hfsol

    integrert_sol = numerisk_integrasjon(C_p_sol, constants5.Tref_K, T)
    integrert_co2 = numerisk_integrasjon(C_p_CO2, constants5.Tref_K, T)

    dH_CO2 = hf_CO2 + w_CO2 * integrert_co2[0] #KJ/kg
    dH_MEA = hf_MEA + w_MEA * integrert_sol[0] # kJ/kg
    dH_H2O = w_CO2 * d_abs_h_co2

    h4 = dH_CO2 + dH_H2O + dH_MEA# kJ/mol
    return h4

h_sol_co2_strøm4 = h_sol_co2_strøm3(constants5.T[3]+273)
h_sol_co2_strøm5 = h_sol_co2_strøm3(constants5.T[4]+273)
print(f'Entalpi for strøm 4: {h_sol_co2_strøm4:.3f} kJ')
print(f'Entalpi for strøm 5: {h_sol_co2_strøm5:.3f} kJ')

def h_strøm8():
    T8 = constants5.T[7] # Celsius
    Tref = constants5.Tref_C  # Celsius

    w_CO2 = mb.wc8
    w_H2O = mb.wh8
    w_N2 = mb.wn8
    w_O2 = mb.wo8

    Cp_CO2 = constants5.cpg[0] # kJ/kg
    Cp_H2O = constants5.cpg[1] # kJ/kg
    Cp_N2 = constants5.cpg[2] # kJ/kg
    Cp_O2 = constants5.cpg[3] # kJ/kg

    hf_CO2 =w_CO2 * constants5.hf[0]
    hf_H2O = w_H2O * constants5.hf[1]
    #std dannelse entalpi for N2 og O2 = 0


    dH_CO2 = hf_CO2 + w_CO2 * Cp_CO2 * (T8 - Tref) # kJ/mol
    dH_H2O =hf_H2O + w_H2O * Cp_H2O * (T8 - Tref) # kJ/mol
    dH_N2 = w_N2 * Cp_N2 * (T8 - Tref) # kJ/mol
    dH_O2 = w_O2 * Cp_O2 * (T8 - Tref) # kJ/mol

    h8 = dH_CO2 + dH_H2O + dH_N2 + dH_O2 # kJ/mol
    return h8*mb.m8

h8 = h_strøm8()
print(f'Entalpi for strøm 8: {h8:.3f} kJ')

def h_strøm9():
    T9 = constants5.T[8] # Celsius
    Tref = constants5.Tref_C  # Celsius

    w_CO2 = constants5.wc9 

    Cp_CO2 = constants5.cpg[0] # kJ/kg

    hf_CO2 = w_CO2 * constants5.hf[0] #kJ/kg


    dH_CO2 = hf_CO2 + w_CO2 * Cp_CO2 * (T9 - Tref) #kJ/kg

    h9 = dH_CO2 
    return h9

h9 = h_strøm9()
print(f'Entalpi for strøm 9: {h9:.3f} kJ')

    

#Koker:
def V_4():
    m5 = mb.m5
    m6 = mb.m6
    m9 =mb.m9

    H5 = m5*h_sol_co2_strøm4(constants5.T[4]+273)
    H6 = m6*h_sol_co2_strøm3(constants5.T[5]+273)
    H9 = m9*h_strøm9()
    Q1 = Q_V_3

    #V3 = H9 + H6 - (H5 + V3)
    Q2 = (H6 + H9 - (H5 + Q1)) *1000 #J/s
    return Q2 
    
Q_V_4 = V_4()
print(f'Totalt energiforbruk Q [kJ/s] i kokeren V-4: {Q_V_4:.3f} J/s')
