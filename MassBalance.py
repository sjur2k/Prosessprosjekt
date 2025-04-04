import scipy as sp
import numpy as np
import constants5 as const
import WtFrac as wf

def massBalance(vars):
    M_CO2 = const.Mw[0]
    M_H20 = const.Mw[1]
    
    #Strøm2
    wc2, wh2, wn2, wo2, \
    wc3, wh3, \
    wc4, wh4, \
    wc8, wh8, \
    m2, m3, m4, m5, m6, m7, m8, m9 = vars
    
    # mass balance equations
    eq1_MS = m2 - (const.m1-m9)
    eq2_MS = m3 - (m7)
    eq3_MS = m4 - (const.m1+m3-m2)
    eq4_MS = m5 - (m4)
    eq5_MS = m6 - (m5-m9)
    eq6_MS = m7 - (m6)
    eq7_MS = m8 - (m9/wc8)
    eq8_MS = m9 - (const.wcapture*const.wc1*const.m1)
       
    # mass fraction equations
    eq1_MF = wc2 - ((const.m1*const.wc1-m9)/m2)
    eq2_MF = wh2 - (const.m1/m2)*const.wh1
    eq3_MF = wn2 - (const.m1/m2)*const.wn1
    eq4_MF = wo2 - (const.m1/m2)*const.wo1
    eq5_MF = wc3 - (wf.WtFracCO2(const.alpha3))
    eq6_MF = wh3 - (1-const.waMEA-wc3)
    eq7_MF = wc4 - (wf.WtFracCO2(const.alpha4))
    eq8_MF = wh4 - (1-const.waMEA-wc4)
    eq9_MF = wc8 - (1/(1+((1-const.xc8)*M_H20)/(const.xc8*M_CO2)))
    eq10_MF = wh8 - (1-wc8)
    
    return [eq1_MS, eq2_MS, eq3_MS, eq4_MS, eq5_MS, eq6_MS, eq7_MS, eq8_MS,
            eq1_MF, eq2_MF, eq3_MF, eq4_MF, eq5_MF, eq6_MF, eq7_MF, eq8_MF, eq9_MF, eq10_MF]
    
    
Vars=["wc2", "wh2", "wn2", "wo2", \
    "wc3", "wh3", \
    "wc4", "wh4", \
    "wc8", "wh8", \
    "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9"]
sol= sp.optimize.root(massBalance,np.ones(18)).x
for i in range(len(sol)):
    if i<=9:
        print(f"{Vars[i]} = {round(sol[i]*100,1)}%")
    else:
        print(f"{Vars[i]} = {round(sol[i],1)}kg/s")