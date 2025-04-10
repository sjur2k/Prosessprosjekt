import constants5 as const


def WtFracCO2(loading):
    return (loading*(const.Mw[0]))/((1+ 0.7/0.3)*(const.MwMEA))
    