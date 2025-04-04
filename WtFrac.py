import constants5 as const


def WtFracCO2(loading):
    return (loading*const.M_CO2)/((1+ 0.7/0.3)*const.M_MEA)
    