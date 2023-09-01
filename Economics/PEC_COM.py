# Reference: "Technical and economic working domains of industrial heat pumps:
# Part 1 - single stage vapour compression heat pumps"
# by Ommen, Torben Schmidt; Jensen, Jonas Kjær; Markussen, Wiebke Brix;
# Reinholdt, Lars; Elmegaard, Brian.
# In International Journal of Refrigeration, Vol. 55, 2015, p. 168–182.
# - V_dot in [m3/h]
# - PEC_Comp in [€]
def PEC_Com(V_dot, eta_vol_comp):
    PEC_Comp = 19850 * (V_dot / eta_vol_comp / 279.8) ** 0.73 * 1.73 # NZD ,Euro/NZD exchange rate = 1.73
    return PEC_Comp
