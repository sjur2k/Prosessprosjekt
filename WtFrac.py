import constants5 as const


def WtFracCO2(loading):
    return (loading*const.Mw[0])/((1+ 0.7/0.3)*const.MwMEA)
    
    
print(f"w_co2_3 = {WtFracCO2(0.18):.3f}")
print(f"w_co2_4 = {WtFracCO2(0.47):.3f}")