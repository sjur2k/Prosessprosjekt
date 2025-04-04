import numpy as np
import constants5 as const

c=const.Hc
c=[i*1.0E-3 for i in c]
def K(T,alpha):
    return (c[0]+c[1]*alpha*T**-1)*np.exp(c[2]*alpha**2+c[3]*T**-1+c[4]*alpha*T**-2)
def p_CO2_g(K2,KH,alpha):
    return KH/K2*alpha**2/(1-2*alpha)**2

alpha = np.linspace(0,0.499,200) #Viktig å ikkje la alpha = 0.5 -> divbyzero
T=[273,298,313,333,353,373,393,423]
K2=[3.93E5,3.7E4,1.14E4,2.43E3,5.78E2,2.46E2,4.08E1,6.74]

KH=[]
p_CO2=[]
for i in range(len(T)):
    KH.append([K(T[i],a)for a in alpha])
    p=[]
    for j in range(len(alpha)):
        p.append(p_CO2_g(K2[i],KH[-1][j],alpha[j]))
    p_CO2.append(p)