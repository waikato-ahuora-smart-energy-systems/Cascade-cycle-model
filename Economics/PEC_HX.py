# Reference: "Technical and economic working domains of industrial heat pumps:
# Part 1 - single stage vapour compression heat pumps"
# by Ommen, Torben Schmidt; Jensen, Jonas Kjær; Markussen, Wiebke Brix;
# Reinholdt, Lars; Elmegaard, Brian.
# In International Journal of Refrigeration, Vol. 55, 2015, p. 168–182.
# - Area in [m2]
# - PEC_HX in [€]
def PEC_HX(area):
    PEC_HX = 15526 * (area / 42) ** 0.8 * 1.73 # NZD ,Euro/NZD exchange rate = 1.73
    return PEC_HX

