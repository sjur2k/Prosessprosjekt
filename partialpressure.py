import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",        
    "mathtext.rm": "serif",
    "font.size": 14,
})


c=[4.96563E2,3.41697E5,1.69131,1.47225E3,1.28338E5]
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

plt.figure(figsize=(8,7))
colors = cm.bwr(np.linspace(0.1, 0.9, len(T)))

for i in range(len(p_CO2)):
    plt.plot(alpha,p_CO2[i],color=colors[i],label=f"T = {T[i]}K")
    plt.xlabel(r"'Loading'-faktor, $\alpha$")
    plt.ylabel(r"$p_{\mathrm{CO_2}}$ [Pa]")    
    plt.yscale("log")
    plt.xlim(0,0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.tight_layout()
plt.show()
#plt.savefig("CO2_partialkonsentrasjon.eps",dpi=300)