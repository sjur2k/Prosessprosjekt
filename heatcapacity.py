import constants5 as const

def Cp_func(T,A,B,C):
    Cp = A + B * (T+273.15) + C * (T + 273.15) ** 2
    return Cp

def Cp_H2O(T):
    return Cp_func(T, const.Aw, const.Bw, const.Cw)
def Cp_MEA(T):
    return Cp_func(T, const.Aa, const.Ba, const.Ca)
def Cp_MEA_H2O(T):
    Cp = (1-const.waMEA)*Cp_H2O(T) + const.waMEA*Cp_MEA(T),\
        +const.waMEA*(1-const.waMEA)*(const.As+const.Bs*(T+273.15)+const.Cs*const.waMEA*(T)**(-1.5859))
    return Cp
def Cp_CO2_ABSORBED(T):
    return Cp_func(T, const.Ac, const.Bc, 0)

def meanCp_H2O(T0,T1):
    return const.Aw+(1/2)*const.Bw*((T1+273.15)**2-(T0+273.15)**2)/(T1-T0)+(1/3)*const.Cw*((T1+273.15)**3-(T0+273.15)**3)/(T1-T0)
def meanCp_MEA(T0,T1):
    return const.Aa+(1/2)*const.Ba*((T1+273.15)**2-(T0+273.15)**2)/(T1-T0)+(1/3)*const.Ca*((T1+273.15)**3-(T0+273.15)**3)/(T1-T0)
def meanCp_MEA_H2O(T0,T1):
    part1 = (1-const.waMEA)*meanCp_H2O(T0,T1)
    part2 = const.waMEA*meanCp_MEA(T0,T1)
    part3 = const.waMEA*(1-const.waMEA)*(const.As+(1/2)*const.Bs*((T1+273.15)**2-(T0+273.15)**2)/(T1-T0)+(const.Cs*const.waMEA/(-0.5859))*(T1**(-0.5859)-T0**(-0.5859))/(T1-T0))
    return part1+part2+part3
def meanCp_CO2_ABSORBED(T0,T1):
    return const.Ac+(1/2)*const.Bc*((T1+273.15)**2-T0+(273.15)**2)/(T1-T0)

print(meanCp_MEA_H2O(40,50))