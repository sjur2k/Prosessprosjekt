import constants5 as const

def Cp(T,A,B,C):
    Cp = A + B * (T+273.15) + C * (T + 273.15) ** 2
    return Cp

def Cp_H2O(T):
    return Cp(T, const.Aw, const.Bw, const.Cw)
def Cp_MEA(T):
    return Cp(T, const.Aa, const.Ba, const.Ca)
def Cp_MEA_H2O(T):
    return Cp(T, const.As, const.Bs, const.Cs)
def Cp_CO2_IN_SOLUTION(T):
    return Cp(T, const.Ac, const.Bc, 0)

def meanCp(T0,T1,A,B,C):
    return A+(1/2)*B*(T1**2-T0**2)/(T1-T0)+(1/3)*C*(T1**3-T0**3)/(T1-T0)

def meanCp_H2O(T0,T1):
    return meanCp(T0,T1,const.Aw,const.Bw,const.Cw)
def meanCp_MEA(T0,T1):
    return meanCp(T0,T1,const.Aa,const.Ba,const.Ca)
def meanCp_MEA_H2O(T0,T1):
    return meanCp(T0,T1,const.As,const.Bs,const.Cs)
def meanCp_CO2_IN_SOLUTION(T0,T1):
    return meanCp(T0,T1,const.Ac,const.Bc,0)

print(meanCp_H2O(5,30))