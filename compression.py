import constants5 as const
import enthalpy as h


# Ettstegs kompresjon
def kompresjon_1():
    Tinn = const.T[8]+273.15 # K
    g = 1.3 # Adiabatisk eksponent for CO2: ca 1.3 ved 20C
    pb = 20 # Bar
    Tb_Kelvin = Tinn*(1+(1/const.eta)*((pb/const.p[8])**((g-1)/g)-1)) # K
    return Tb_Kelvin-273.15 # C

# Trestegs kompresjon
Tinn = const.T[8] + 273.15 # K
Tc_Kelvin = Te_Kelvin = 303
pb = 4 # Bar
pd = 8 # Bar
pf = 20 # Bar
g = 1.3 # Adiabatisk eksponent for CO2: ca 1.3 ved 20C
Tb_Kelvin = Tinn*(1+(1/const.eta)*((pb/const.p[8])**((g-1)/g)-1))
Td_Kelvin = Tc_Kelvin*(1+(1/const.eta)*((pd/pb)**((g-1)/g)-1))
Tf_Kelvin = Te_Kelvin*(1+(1/const.eta)*((pf/pd)**((g-1)/g)-1))
Tb = Tb_Kelvin-273.15
Tc = Tc_Kelvin-273.15
Td = Td_Kelvin-273.15
Te = Te_Kelvin-273.15
Tf = Tf_Kelvin-273.15

P1 = h.m["m9"]*h.mean_cp_c(const.T[8],kompresjon_1())*(kompresjon_1()-const.T[8]) # kJ/kg

P3 = h.m["m9"]*h.mean_cp_c(const.T[8],Tb)*(Tb-const.T[8])
P3 += h.m["m9"]*h.mean_cp_c(Tc,Td)*(Td-Tb)
P3 += h.m["m9"]*h.mean_cp_c(Te,Tf)*(Tf-Te)

print(f"Tb ved 1stegs kompresjon: {round(kompresjon_1(),2)}°C")
print(f"Akselarbeid for 1stegs kompresjon: {round(P1,2)}kW\n")
print(f"Tb ved 3stegs kompresjon: {round(Tb,2)}°C")
print(f"Td ved 3stegs kompresjon: {round(Td,2)}°C")
print(f"Tf ved 3stegs kompresjon: {round(Tf,2)}°C")
print(f"Akselarbeid for 3stegs kompresjon: {round(P3,2)}kW")