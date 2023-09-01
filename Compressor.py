# Determine isentropic efficiency of the compressor
# Reference: J.F. Wang et al.(2018). Heat pump heat recovery options for food industry dryers
def eta_is_comp(PR):
    #n_is = 0.65 + 0.015 * PR - 0.0015 * PR ** 2
    n_is = 0.7
    return n_is