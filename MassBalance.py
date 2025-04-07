import scipy.optimize
import numpy as np
import constants5 as const
import WtFrac as wf


wc3 = wf.WtFracCO2(const.alpha3) 
wc4 = wf.WtFracCO2(const.alpha4)
wm3 = wm4 = const.waMEA
wh3 = 0.7 - wc3
wh4 = 0.7 - wc4

def massBalances(vars):
    m2, m3, m4, m5, m6, m8, m9 = vars
    
    eq1 = m2 - (const.m1 - m9)
    eq2 = m4 - (const.m1 + m3 - m2)  
    eq3 = m5 - m4  
    eq4 = m6 - (m5 - m9)
    eq5 = m8 - (m9/const.wcapture)
    eq6 = m9 - (const.wcapture * const.wc1 * const.m1)
    eq7 = ((const.m1 * const.wc1 + m3 * wc3 - m4 * wc4)/m2) - const.xc8  
    
    return [eq1, eq2, eq3, eq4, eq5, eq6, eq7]

x0 = np.ones(7)  # initial guess for m2, m3, m4, m5, m6, m8, m9

solution = scipy.optimize.root(massBalances, x0)
m2, m3, m4, m5, m6, m8, m9 = solution.x

wc2 = ((const.m1 * const.wc1) + (m3 * wc3) - (m4 * wc4)) / m2
wh2 = (const.m1 / m2) * const.wh1
wn2 = (const.m1 / m2) * const.wn1
wo2 = (const.m1 / m2) * const.wo1
wt2 = wc2 + wh2 + wn2 + wo2

print("\nMassestrømmer [kg/s]:")
print(f"m2 (gass ut): {m2:.2f}")
print(f"m3 (MEA+CO2 ut): {m3:.2f}")
print(f"m4 (MEA inn): {m4:.2f}")
print(f"m5 (=m4): {m5:.2f}")
print(f"m6: {m6:.2f}")
print(f"m8: {m8:.2f}")
print(f"m9 (fanget CO2): {m9:.2f}")

print("\nMassefraksjoner i strøm 2 (gass ut):")
print(f"CO2: {wc2:.4f}")
print(f"H2O: {wh2:.4f}")
print(f"N2:  {wn2:.4f}")
print(f"O2:  {wo2:.4f}")
print(f"Total: {wt2:.4f}")


