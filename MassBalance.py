import scipy.optimize
import numpy as np
import constants5 as const
import WtFrac as wf

def massBalance(vars):
    
    
    #Strøm2
    wc2, wh2, wn2, wo2, \
    wc3, wh3, \
    wc4, wh4, \
    wc8, wh8, \
    m2, m3, m4, m5, m6, m7, m8, m9 = vars
    
    # mass balance equations
    eq1_MS = const.m1 -m9-m2
    eq2_MS = m7-m3
    eq3_MS = const.m1+m3-m2-m4
    eq4_MS = m4-m5
    eq5_MS = m5 - m9 - m6
    eq6_MS = m6-m7
    eq7_MS = m9/wc8-m8
    eq7_MC = const.wcapture*const.wc1*const.m1 - m9
    
    
    
    # mass fraction equations
    eq1_MF = (const.m1*const.wc1+m3*wc3-m4*wc4)/m2 - wc2
    eq2_MF = (const.m1/m2)*const.wh1-wh2
    eq3_MF = (const.m1/m2)*const.wn1-wn2
    eq4_MF = (const.m1/m2)*const.wo1-wo2
    eq5_MF = wf.WtFracCO2(const.alpha3)-wc3
    eq6_MF = 1-const.waMEA-wc3-wh3
    eq7_MF = wf.WtFracCO2(const.alpha4)-wc4
    eq8_MF = 1-const.waMEA-wc4-wh4
    eq9_MF = 1/(1+(0.3*const.Mw[1])/(0.7*const.Mw[0]))-wc8
    eq10_MF = 1-wc8-wh8
    
    return [eq1_MS, eq2_MS, eq3_MS, eq4_MS, eq5_MS, eq6_MS, eq7_MS, eq7_MC,
            eq1_MF, eq2_MF, eq3_MF, eq4_MF, eq5_MF, eq6_MF, eq7_MF, eq8_MF, eq9_MF, eq10_MF]
    
    

print(scipy.optimize(massBalance(np.zeros(18))))