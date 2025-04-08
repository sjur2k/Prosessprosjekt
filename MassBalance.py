import scipy.optimize
import numpy as np
import constants5 as const
import WtFrac as wf
import os



# ------------------------------------------------------------
use_mass_approximation = False
#-------------------------------------------------------------



if use_mass_approximation:
    print("Antar m_MEA >> m_CO2\n")
else: print("Antar ikke m_MEA >> m_CO2\n")
print("---------------------------\n")

wc3 = wf.WtFracCO2(const.alpha3, use_mass_approximation) 
wc4 = wf.WtFracCO2(const.alpha4, use_mass_approximation)
wm3 = wf.WtFracMEA(const.alpha3, use_mass_approximation)
wm4 = wf.WtFracMEA(const.alpha4, use_mass_approximation)
wh3 = 1- wm3 - wc3
wh4 = 1- wm4 - wc4

wc6, wm6, wh6 = wc3, wm3, wh3

 

wc8 = (const.xc8*const.Mw[0]) / (const.xc8*const.Mw[0] + (1-const.xc8)*const.Mw[1]) # Molfraksjon -> vektfraksjon
wh8 = 1-wc8

def massBalances(vars):
    # Input
    m2, m4, m6, m8 = vars
    # Output
    eq = [0]*4

    #Vektfraksjoner i strøm 2 
    wn2 = (const.m1*const.wn1)/m2
    wo2 = (const.m1*const.wo1)/m2
    wc2 = ((1-const.wcapture)*const.m1*const.wc1)/m2
    wh2 = 1-(wn2+wo2+wc2)
    
    # Massestrøm fra reflukstank til stripper
    resirkulert_H2O = wh8*m8 

    # CO2-komponentbalanse for absorber
    eq[0] = const.m1*const.wc1 + m6*wc6 - m2*wc2 - m4*wc4 
    
    # Massebalanse for absorber
    eq[1] = const.m1 + m6 - m2 - m4

    # Massebalanse for stripper    
    eq[2] = m4 + resirkulert_H2O - m6 - m8 
    
    # CO2-komponentbalanse for stripper
    eq[3] = m4*wc4 - m6*wc6 - m8*wc8 

    return eq

   
#Initialbetingelser for optimeringsalgoritmen
x0 = [550.0, 1000.0, 1000.0, 50.0]

solution = scipy.optimize.root(massBalances, x0)
m2, m4, m6, m8= solution.x
m9 = const.wcapture*const.m1*const.wc1

wn2 = (const.m1*const.wn1)/m2
wo2 = (const.m1*const.wo1)/m2
wc2 = ((1-const.wcapture)*const.m1*const.wc1)/m2
wh2 = 1-(wn2+wo2+wc2)    

stripper_massebalanse = abs(m4+(1-wc8)*m8-(m6+m8))
absorber_massebalanse = abs(const.m1+m6-(m2+m4))
tot_mass = abs(const.m1-(m2+m9))

# Logging av resultater
print(f"Massebalanse for absorber: {absorber_massebalanse:.2f}")
print(f"Inn i absorber: {const.m1+m6:.2f}")
print(f"Ut av absorber: {m2+m4:.2f}")
print("---------------------------\n")
print(f"Massebalanse for stripper: {stripper_massebalanse:.2f}")
print(f"Inn i stripper: {m4+(1-wc8)*m8:.2f}")
print(f"Ut av stripper: {m6+m8:.2f}")
print("---------------------------\n")
print(f"Total massebalanse: {tot_mass:.2f}")
print(f"Ut av systemet:{m2+m9:.2f}")
print(f"Inn i systemet:{const.m1:.2f}\n")
print("                |            Massefraksjoner       ")
print("     |     m    |        gass         |         væske     ")
print("Strøm|  [kg/s]  | CO2 | H2O | N2 | O2 |  MEA  |  H2O   |  CO2  ")
print(f"  1  |  {const.m1:.2f}  |{const.wc1:.2f} |{const.wh1:.2f} |{const.wn1:.2f}|{const.wo1:.2f}|   0   |   0    |   0   ")
print(f"  2  |  {m2:.2f}  |{wc2:.2f} |{wh2:.2f} |{wn2:.2f}|{wo2:.2f}|   0   |   0    |   0   ")
print(f"  3  |  {m6:.2f}  |  0  |  0  | 0  | 0  | {wm3:.2f}  |  {wh3:.2f}  | {wc3:.2f} ")
if use_mass_approximation:
    print(f"  4  |  {m4:.2f}  |  0  |  0  | 0  | 0  | {wm4:.2f}  |  {wh4:.2f}  | {wc4:.2f} ")
    print(f"  5  |  {m4:.2f}  |  0  |  0  | 0  | 0  | {wm4:.2f}  |  {wh4:.2f}  | {wc4:.2f} ")
else:
    print(f"  4  | {m4:.2f}  |  0  |  0  | 0  | 0  | {wm4:.2f}  |  {wh4:.2f}  | {wc4:.2f} ")
    print(f"  5  | {m4:.2f}  |  0  |  0  | 0  | 0  | {wm4:.2f}  |  {wh4:.2f}  | {wc4:.2f} ")
print(f"  6  |  {m6:.2f}  |  0  |  0  | 0  | 0  | {wm6:.2f}  |  {wh3:.2f}  | {wc6:.2f} ")
print(f"  7  |  {m6:.2f}  |  0  |  0  | 0  | 0  | {wm6:.2f}  |  {wh3:.2f}  | {wc6:.2f} ")
print(f"  8  |   {m8:.2f}  |{wc8:.2f} |{wh8:.2f} | 0  | 0  |   0   |   0    |   0   ")
print(f"  9  |   {m9:.2f}  |{const.wc9:.2f} |  0  | 0  | 0  |   0   |   0    |   0   ")


# Skrive data til fil
# ------------------------------------------------------------
should_write_to_file = False
#------------------------------------------------------------


if should_write_to_file:
    m_vec=[const.m1,m2,m6,m4,m4,m6,m6,m8,m9]
    w_format = ["wc","wh","wn","wo","wm"]
    w_matrix = [[const.wc1,const.wh1,const.wn1,const.wo1,0],
                [wc2, wh2, wn2, wo2, 0],
                [wc3, wh3, 0, 0, wm3],
                [wc4, wh4, 0, 0, wm4],
                [wc4, wh4, 0, 0, wm4],
                [wc6, wh6, 0, 0, wm6],
                [wc6, wh6, 0, 0, wm6],
                [wc8, wh8, 0, 0, 0],
                [1, 0, 0, 0, 0]]

    if use_mass_approximation:
        filename = "massbalance_with_approx.txt"
    else: filename = "massbalance_no_approx.txt"
        
    with open(filename,"w") as f:
        for i in range(len(m_vec)):
            f.write(f"m{i+1}: {m_vec[i]:.2f}\n")
            for j in range(len(w_matrix[i])):
                            if w_matrix[i][j]!=0:
                                f.write(f"{w_format[j]}{i+1}: {w_matrix[i][j]:.2f}\n")
            f.write("\n")
            

