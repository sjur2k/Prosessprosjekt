import scipy.optimize
import constants5 as const
import WtFrac as wf

def massBalance(vars):
    
    
    #Strøm2
    wc2, wh2, wn2, wo2, \
    wc3, wh3, wm3,\
    wc4, wh4, wm4, \
    wc5, wh5, wm5, \
    wc6, wh6, wm6, \
    wc7, wh7, wm7, \
    wc8, wh8, \
    wc9, \
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
    
    
    
    
    
    
    