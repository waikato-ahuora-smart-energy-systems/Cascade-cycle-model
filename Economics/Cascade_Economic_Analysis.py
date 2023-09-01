from HP_Model_Input import model
import numpy as np
from Economics.PEC_COM import PEC_Com
from Economics.PEC_HX import PEC_HX

def economic_analysis (Cycle):

    # Input for economic analysis
    #i_eff         = 0.02    # 0.05 [-]Effective interest rate
    i_eff         = 0.05    # 0.05 [-]Effective interest rate based on i = 5.25 % and i_L = 7.25; 12% based on Lana et al.,
    n_lt          = 20     # 20 [-]Lifetime
    OH            = 8000    # [h/year] Operating Hours
    #c_el          = 0.12    # [€/kWh] Specific cost of electricity - Euro price
    c_el          = 0.155   # [NZD/kWh] Specific cost of electricity - New Zealand price based on NZD 40.73 / GJ + carbon price
    c_hsource     = 0       # [€/kWh] Specific cost of heat from heat source
    #c_rev         = 0.04    # [€/kWh] Revenue - sell heat to DH network - - Euro price
    c_rev          = 0.076    # [NZD/kWh] Revenue - sell heat to DH network - New Zealand Price based on 5.157 c/kWh + carbon price
    f_TCI         = 4       # 4 [-] Scaling TCI from PCE (ref. Booster HP case study cost)
    f_OM          = 0.0    # [-] Factor for estimation of O&M cost as one time investment: O&M=TCI*f_OM -  Bejan et al., 1996
    # eta_vol_comp  = 0.9     # [-]Compressor - volumetric efficiency will estimate from the thermodynamic model
    eta_vol_com1 = Cycle.eta_vol_com1
    eta_vol_com2 = Cycle.eta_vol_com2
    eta_vol_com3 = Cycle.eta_vol_com3
    flammable = 'true'      # options = 'true' or 'false'

    # Flammability - check if higher investment on the compressor is required to handle
    # flammability issues (ref. Thermal Design and Optimization)
    if flammable == 'true':
        #f_flammable = 1.2
        f_flammable = 1.0
    else:
        f_flammable = 1.0
    # Cycle specifications and components size
    Qdot_source     = Cycle.Qdot_source              # [W]
    Qdot_sink       = Cycle.Qdot_sink                   # [W]
    W_dot           = Cycle.W_dot                   # [W]
    V_suction_com1  = Cycle.V_suction_com1 * 3600    # [m3/h]
    V_suction_com2  = Cycle.V_suction_com2 * 3600    # [m3/h]
    V_suction_com3  = Cycle.V_suction_com3 * 3600  # [m3/h]
    # Estimate heat transfer area
    htc_eva         = 3000                          # W/m2K @ Zuhlsdorf et al., 2018 Improving the performance of booster heat pumps using zeotropic mixtures
    htc_cond        = 2400                          # W/m2K @ Zuhlsdorf et al., 2018 Improving the performance of booster heat pumps using zeotropic mixtures
    htc_liquid      = 1500                          # W/m2K @ Zuhlsdorf et al., 2018 Improving the performance of booster heat pumps using zeotropic mixtures
    htc_gas         = 250                           # W/m2K Table 1.1 page 8 Fundamental of Heat and Mass Trasfer 6th Edition Incropera, DeWitt, Bergman

    n_states_evap       = model.n_states_evap
    n_states_sh         = model.n_states_sh
    n_states_cond       = model.n_states_cond
    n_states_desup      = model.n_states_desup
    n_states_sub        = model.n_states_sub
    n_states_ihx        = model.n_states_ihx
    n_sates_cascadeHX   = model.n_sates_cascadeHX

    UA_source   = Cycle.UA_source
    UA_cond_top = Cycle.UA_cond_top
    UA_ihx3     = Cycle.UA_ihx3
    UA_ihx2     = Cycle.UA_ihx2
    UA_cax      = Cycle.UA_cax

    U_source    = np.zeros(len(UA_source))
    U_cond_top  = np.zeros(len(UA_cond_top))
    U_ihx3      = np.zeros(len(UA_ihx3))
    U_cax       = np.zeros(len(UA_cax))
    U_ihx2     = np.zeros(len(UA_ihx2))


    # Bottom cycle evaporator area
    for i in range(len(UA_source)):
        if i < n_states_evap:
            U_source[i] = 1/(1/htc_eva+1/htc_liquid) # Evaporation zone
        elif i >= n_states_evap:
            U_source[i] = 1 / (1 / htc_gas + 1 / htc_liquid) # Superheating zone

    A_source = [UA/U for UA,U in zip(UA_source,U_source)]
    A_source_total = sum(A_source)
    print('A_source_total:',A_source_total)

    # Top cycle condenser area
    for i in range(len(UA_cond_top)):
        if i < n_states_sub:
            U_cond_top[i] = 1/(1/htc_liquid+1/htc_gas) # subcooling zone
        elif i >= n_states_sub and i <(n_states_sub+n_states_cond-1):
            U_cond_top[i] = 1/(1/htc_cond+1/htc_gas)   #Evaporation zone
        else:
            U_cond_top = 1/(1/htc_gas + 1/htc_gas)     # Desuperheating zone

    A_cond_top = [UA/U for UA,U in zip(UA_cond_top,U_cond_top)]
    A_cond_top_total = sum(A_cond_top)
    print('A_cond_top_total',A_cond_top_total)

    # Cascade heat exchanger area
    for i in range(len(UA_cax)):
        if i < n_sates_cascadeHX:
            U_cax[i] = 1/(1/htc_eva+1/htc_cond)
        elif i == n_sates_cascadeHX:
            U_cax[i] = 1/(1/htc_gas + 1/htc_eva)
        else:
            U_cax[i] = 1/(1/htc_gas+1/htc_gas)

    A_cax = [UA/U for UA,U in zip(UA_cax,U_cax)]
    A_cax_total = sum(A_cax)
    print('A_cax_total', A_cax_total)

    # IHX3 area
    for i in range(len(UA_ihx3)):
        U_ihx3[i] = 1/(1/htc_liquid + 1/htc_gas)

    A_ihx3 = [UA/U for UA,U in zip(UA_ihx3,U_ihx3)]
    A_ihx3_total = sum(A_ihx3)
    print('A_ihx3_total', A_ihx3_total)

    # IHX2 area
    for i in range(len(UA_ihx2)):
        U_ihx2[i] = 1/(1/htc_gas + 1/htc_liquid)

    A_ihx2 = [UA/U for UA,U in zip(UA_ihx2,U_ihx2)]
    A_ihx2_total = sum(A_ihx2)
    print('A_ihx2_total', A_ihx2_total)
    # =========================COMPRESSOR INVESTMENT====================================================================
    # Compressor 1 - Bottom cycle
    PEC_com1 = f_flammable * PEC_Com(V_suction_com1, eta_vol_com1) # [NZD]
    # Compressor 2 - Bottom cycle
    PEC_com2 = f_flammable * PEC_Com(V_suction_com2, eta_vol_com2) # [NZD]
    #print('PEC_top_comp', PEC_top_comp)
    # Compressor 3 - Top cycle compressor
    PEC_com3 = f_flammable * PEC_Com(V_suction_com3, eta_vol_com3)  # [NZD]
    #print('PEC_bot_comp ', PEC_bot_comp )
    # ========================= HEs INVESTMENT COST ====================================================================
    PEC_Top_Cond = PEC_HX(A_cond_top_total)
    print('PEC_Top_Cond', PEC_Top_Cond)
    PEC_Evap = PEC_HX(A_source_total)
    print('PEC_Evap', PEC_Evap)
    PEC_Cax  = PEC_HX(A_cax_total)
    print('PEC_Cax ', PEC_Cax )
    PEC_IHX3 = PEC_HX(A_ihx3_total)
    print('PEC_IHX3  ', PEC_IHX3 )
    PEC_IHX2 = PEC_HX(A_ihx2_total)
    print('PEC_IHX2  ', PEC_IHX2)
    #===================================================================================================================
    PEC_total = PEC_Top_Cond + PEC_Evap + PEC_Cax + PEC_IHX2 + PEC_IHX3 + PEC_com1 + PEC_com2 + PEC_com3
    print('PEC_total',PEC_total)
    # Total capital investment
    TCI = f_TCI*PEC_total
    print('TCI',TCI)
    TCI_spec = TCI / Qdot_sink * 1000 # [NZD/kW] Specific investment cost per unit of supplied heat
    # Operation and maintenance cost
    OMC_total = f_OM * TCI
    # Capital Recovery Factor
    CRF = i_eff * (1 + i_eff) ** n_lt / ((1 + i_eff) ** n_lt - 1) # Ommen at al., 2015
    print('CRF',CRF)
    # ==================================================================================================================
    # Cash flows
    # Fuel cost HP
    # Annual cash flow for fuel cost electricity of heat pump
    CF_FCel = W_dot * (c_el/3600/1000) * OH * 3600 # NZD/year
    print('CF_FCel',CF_FCel)
    # Annual cash flow for fuel cost for heat source of heat pump
    CF_FChsource = Qdot_source * (c_hsource / 3600 / 1000) * OH * 3600 # [NZD / a]
    print('CF_FChsource',CF_FChsource)
    # Annual cash flow for revenue
    CF_rev  = Qdot_sink * (c_rev/3600/1000) * OH * 3600 #  % [NZD/a]
    print('CF_rev ', CF_rev )
    # Investment Cost cash flow term
    CF_TCI = TCI * CRF # % [NZD / a]
    # Cost for O&M as one time cost at time of investment
    CF_OM = f_OM * TCI * CRF # % [NZD / a]
    # Total annual cash flows
    CF_total = - CF_TCI - CF_OM - CF_FChsource - CF_FCel + CF_rev
    # Net Present Value
    NPV = -TCI - OMC_total + 1 / CRF * (-CF_FCel - CF_FChsource + CF_rev) # NZD
    # Payback time
    PBT = TCI/(CF_rev-CF_FCel-CF_FChsource)
    # Total heat generation
    Energy_Produced_Yearly = OH * Qdot_sink * 1e-6 # MWh
    # Specific cost of heat
    c_H = (CF_FCel + CF_FChsource + CF_TCI) / Energy_Produced_Yearly # [NZD / MWh]
    class Output:
        pass

    output = Output()
    # Outputs of the economic analysis
    output.PEC_com1 = PEC_com1
    output.PEC_com2 = PEC_com2
    output.PEC_com3 = PEC_com3
    output.PEC_Evap = PEC_Evap
    output.PEC_Top_Cond = PEC_Top_Cond
    output.PEC_Cax = PEC_Cax
    output.PEC_IHX2 = PEC_IHX2
    output.PEC_IHX3 = PEC_IHX3

    output.TCI = TCI
    output.TCI_spec = TCI_spec
    output.CF_FCel = CF_FCel
    output.CF_FChsource = CF_FChsource
    output.CF_rev = CF_rev
    output.CF_TCI = CF_TCI
    output.CF_OM = CF_OM
    output.CF_total = CF_total
    output.NPV = NPV
    output.PBT = PBT
    output.Energy_Produced_Yearly = Energy_Produced_Yearly
    output.c_H = c_H
    return output














