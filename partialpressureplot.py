import partialpressure as pp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",        
    "mathtext.rm": "serif",
    "font.size": 14,
})

plt.figure(figsize=(8,7))
colors = cm.bwr(np.linspace(0.1, 0.9, len(pp.T)))

for i in range(len(pp.p_CO2)):
    plt.plot(pp.alpha,pp.p_CO2[i],color=colors[i],label=f"T = {pp.T[i]}K")
    plt.xlabel(r"'Loading'-faktor, $\alpha$")
    plt.ylabel(r"$p_{\mathrm{CO_2}}$ [Pa]")    
    plt.yscale("log")
    plt.xlim(0,0.5)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.tight_layout()
plt.show()
#plt.savefig("CO2_partialkonsentrasjon.eps",dpi=300)

p_CO2=[]
alpha=np.linspace(0,0.49,10)
for a in alpha:
    tempvec=[]
    for i in range(len(pp.T)):
        tempvec.append(pp.p_CO2_g(pp.K2[i],pp.K(pp.T[i],a),a))
    p_CO2.append(tempvec)
plt.figure()
for i in p_CO2:
    plt.plot(pp.T,i)
plt.yscale("log")
plt.show()