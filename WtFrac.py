import constants5 as const

def WtFracCO2(loading, approx_bool):
    if approx_bool:
        return 0.3*loading*(const.Mw[0]/const.MwMEA)
    else: return 1/(1+const.MwMEA/(const.Mw[0]*loading*0.3))

def WtFracMEA(loading, approx_bool):
    if approx_bool:
        return const.waMEA
    else: return 1/(1/0.3 + loading*const.Mw[0]/const.MwMEA)