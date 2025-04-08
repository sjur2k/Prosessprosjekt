import scipy.optimize
import numpy as np
import constants5 as const
import WtFrac as wf


wc3 = wf.WtFracCO2(const.alpha3) 
wc4 = wf.WtFracCO2(const.alpha4)
wm3 = wm4 =const.waMEA
wh3 = 0.7 - wc3
wh4 = 0.7 - wc4


# #Vektfraksjoner strøm 3 = strøm 6

# wc3 = wf.WtFracCO2(const.alpha3)
# wm3 = const.waMEA
# wh3 = wc3-wh3
wh6 = wh3
wc6 = wc3


# #Vektfraksjoner strøm 4
# wc4 = wf.WtFracCO2(const.alpha4)


# def liquidMassFracs(alpha):
#     n_MEA = 100.0
#     mass_MEA = n_MEA * const.MwMEA / 1000
    
#     mass_h2o = mass_MEA*(1/0.30 -1)
    
#     n_co2 = alpha*n_MEA    
#     mass_co2 = n_co2 * const.Mw[0] / 1000
    
#     total_mass = mass_MEA + mass_h2o + mass_co2
    
#     wc = mass_co2 / total_mass
#     wh = mass_h2o / total_mass
#     wm = mass_MEA / total_mass
    
#     return wc, wh, wm

# wc3, wh3, wm3 = liquidMassFracs(const.alpha3)
# wc4, wh4, wm4 = liquidMassFracs(const.alpha4)
# wc6, wh6, wm6 = liquidMassFracs(const.alpha3) #same as wc3, wh3, wm3  

wc8 = (const.xc8*const.Mw[0] / 1000) / (const.xc8*const.Mw[0] / 1000 + (1-const.xc8)*const.Mw[1] / 1000)
wh8 = 1-wc8


def massBalances(vars):
    m2, m4, m6, m8 = vars
    eq = [0]*4
    
    #massfractions in stream 2 
    wn2_out = (const.m1*const.wn1)/m2
    wo2_out = (const.m1*const.wo1)/m2
    wc2_out = ((1-const.wcapture)*const.m1*const.wc1)/m2
    wh2_out = 1-(wn2_out+wo2_out+wc2_out)
    
    m9 = const.wcapture*const.m1*const.wc1 #mass flow CO2 in stream 9
    m10 = m8*wh8
    print(f"m10: {m10:.4f}")
    
    
    #massbalances

    
    eq[0] = const.m1*const.wc1+m6*wc6-(m2*wc2_out+m4*wc4) #
    
    eq[1] = const.m1 +m6 -(m2+m4)

    eq[2] = m4 - (m6+m8)
    
    eq[3] = m4*wc4-(m6*wc6+m8*wc8)
    
    
    
    return eq

   

x0 = [550.0, 900.0, 900.0, 60.0]  # initial guess for m2, m3, m4, m5, m6, m8, m9

solution = scipy.optimize.root(massBalances, x0)
m2, m4, m6, m8= solution.x
m9 = const.wcapture*const.m1*const.wc1


wh2 = (const.m1*const.wh1+wh3*m6-wh4*m4)/m2
wn2 = (const.m1 / m2) * const.wn1
wo2 = (const.m1 / m2) * const.wo1
wc2 = 1 - (wh2 + wn2 + wo2)
wt2 = wc2 + wh2 + wn2 + wo2

print("\nMassestrømmer [kg/s]:")
print(f"m2: {m2:.2f}")
print(f"m3: {m6:.2f}")
print(f"m4: {m4:.2f}")
print(f"m5 (=m4): {m4:.2f}")
print(f"m6: {m6:.2f}")
print(f"m7 (=m6): {m6:.2f}")
print(f"m8: {m8:.2f}")
print(f"m9 (fanget CO2): {m9:.2f}")

print("\nMassefraksjoner i strøm 2 (gass ut):")
print(f"CO2: {wc2:.4f}")
print(f"H2O: {wh2:.4f}")
print(f"N2:  {wn2:.4f}")
print(f"O2:  {wo2:.4f}")
print(f"Total: {wt2:.4f}")
print("---------------------------\n")


stripperen = abs(m4- (m6+m8))
absorberen = abs(const.m1+m6-(m2+m4))
tot_mass = abs(const.m1-(m2+m9))
print(f"massebalanse over absorberen: {absorberen:.2f}")
print(f"Inn i absorberen: {const.m1+m6:.2f}")
print(f"Ut av absorberen: {m2+m4:.2f}")
print("---------------------------\n")

print(f"massebalanse over stripperen: {stripperen:.2f}")
print(f"Inn i stripperen: {m4:.2f}")
print(f"Ut av stripperen: {m6+m8:.2f}")
print("---------------------------\n")
print(f"Total massebalanse: {tot_mass:.2f}")
print(f"ut av systemet:{m2+m9:.2f}")
print(f"inn i systemet:{const.m1}")


print(f"wm3: {wm3:.4f}")
print(f"wm4: {wm4:.4f}")

