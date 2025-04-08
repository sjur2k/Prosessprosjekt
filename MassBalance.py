import scipy.optimize
import numpy as np
import constants5 as const
import WtFrac as wf



# ------------------------------------------------------------
use_mass_approximation = False 
#-------------------------------------------------------------



if use_mass_approximation:
    print("Antar m_MEA >> m_CO2")
else: print("Antar ikke m_MEA >> m_CO2")

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
    m2, m4, m6, m8 = vars
    eq = [0]*4
    

    #Vektfraksjoner i strøm 2 
    wn2 = (const.m1*const.wn1)/m2
    wo2 = (const.m1*const.wo1)/m2
    wc2 = ((1-const.wcapture)*const.m1*const.wc1)/m2
    wh2 = 1-(wn2+wo2+wc2)
    resirkulert_H2O = wh8*m8 # Fra reflukstank til stripper

    eq[0] = const.m1*const.wc1 + m6*wc6 - m2*wc2 - m4*wc4 # CO2-komponentbalanse for absorberen
    
    eq[1] = const.m1 + m6 - m2 - m4 # Massebalanse for absorber

    eq[2] = m4 + resirkulert_H2O - m6 - m8 # Massebalanse for stripper 
    
    eq[3] = m4*wc4 - m6*wc6 - m8*wc8 # CO2-komponentbalanse for stripperen

    return eq

   

x0 = [550.0, 1000.0, 1000.0, 50.0]  # initial guess for m2, m3, m4, m5, m6, m8, m9

solution = scipy.optimize.root(massBalances, x0)
m2, m4, m6, m8= solution.x
m9 = const.wcapture*const.m1*const.wc1


wn2 = (const.m1*const.wn1)/m2
wo2 = (const.m1*const.wo1)/m2
wc2 = ((1-const.wcapture)*const.m1*const.wc1)/m2
wh2 = 1-(wn2+wo2+wc2)    


stripperen = abs(m4+(1-wc8)*m8-(m6+m8))
absorberen = abs(const.m1+m6-(m2+m4))
tot_mass = abs(const.m1-(m2+m9))
print(f"massebalanse over absorberen: {absorberen:.2f}")
print(f"Inn i absorberen: {const.m1+m6:.2f}")
print(f"Ut av absorberen: {m2+m4:.2f}")
print("---------------------------\n")

print(f"massebalanse over stripperen: {stripperen:.2f}")
print(f"Inn i stripperen: {m4+(1-wc8)*m8:.2f}")
print(f"Ut av stripperen: {m6+m8:.2f}")
print("---------------------------\n")
print(f"Total massebalanse: {tot_mass:.2f}")
print(f"ut av systemet:{m2+m9:.2f}")
print(f"inn i systemet:{const.m1:.2f}\n")


print("                |            massefraksjoner       ")
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
