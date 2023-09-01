import CoolProp as CP
import numpy
import numpy as np
import math
from timeit import default_timer as timer
from HP_Model_Input import model, fluids, source, ref, sink, comp
from Compressor import eta_is_comp
#import CoolProp.CoolProp


#P.CoolProp .set_config_string(CP.CoolProp .ALTERNATIVE_REFPROP_PATH, r'C:\Program Files (x86)\REFPROP')



#def Cascade_cycle (model, fluids, source, ref, sink, comp,var):
def cascade_cycle (refrigerant,parameters):
    print('RUN CASCADE CYCLE')


    # ___________________________________Heat Pump Parameter__________________________________

    # Discretization
    n_states_evap       = model.n_states_evap
    n_states_sh         = model.n_states_sh
    n_states_cond       = model.n_states_cond
    n_states_desup      = model.n_states_desup
    n_states_sub        = model.n_states_sub
    n_states_ihx        = model.n_states_ihx
    n_sates_cascadeHX   = model.n_sates_cascadeHX

    # ___________________________________Numerical input parameter____________________________

    # Heat pump design models:
    limit_sc = model.limit_sc  # 0.95  Limit above which process is considered as supercritcal
    error_delta_T = model.error_delta_T  # 1e-2 Allowable error on pinch temperature
    iter_max = model.iter_max  # Number of maximum iterations

    tol = model.tol_HP
    step = model.step
    # relaxation      = model.relaxation

    # _________________________________________ Input thermodynamics _______________________________

    T_source_in = source.T_in
    T_source_out = source.T_out_design
    p_source = source.p_in
    delta_T_pinch_source_min = source.delta_T_pinch_min
    delta_T_minSH = ref.delta_T_minSH
    delta_T_SH_IHX1 = parameters[0]
    delta_T_SH_IHX2 = parameters[1]
    Qdot_source = source.Qdot_design
    T_sink_in = sink.T_in
    T_sink_out = sink.T_out_design
    p_sink = sink.p_in
    delta_T_pinch_sink_min = sink.delta_T_pinch_min
    delta_T_pinch_ihx_min = ref.delta_T_pinch_IHX_min
    delta_T_pinch_cascade_HX_min = ref.delta_T_pinch_IHX_min

    #eta_is_comp = comp.eta_is
    eta_motor = comp.eta_motor
    # _________________________________________ Working fluid ____________________________________

    MEDIUM_TOP = refrigerant.top
    MEDIUM_BOT = refrigerant.bot
    SINK = fluids.sink
    SOURCE = fluids.source

    # ===============================Run the heat pump function ==================================================
    #  Calculate constants
    MEDIUM_TOP.update(CP.PT_INPUTS,101325,300)
    T_crit_top = MEDIUM_TOP.T_critical()
    p_crit_top = MEDIUM_TOP.p_critical()

    MEDIUM_BOT.update(CP.PT_INPUTS,101325,300)
    T_crit_bot = MEDIUM_BOT.T_critical()
    p_crit_bot = MEDIUM_BOT.p_critical()

    T_21 	= T_source_in
    T_23  	= T_source_out

    SOURCE.update(CP.PT_INPUTS,p_source,T_source_in)
    h_21 = SOURCE.hmass()
    s_21 = SOURCE.smass()
    cpW  = SOURCE.cpmass()
    muW  = SOURCE.viscosity()
    kW   = SOURCE.conductivity()

    SOURCE.update(CP.PT_INPUTS,p_source,T_source_out)
    h_23 =  SOURCE.hmass()
    s_23 =  SOURCE.smass()

    T_31    = T_sink_in
    T_35    = T_sink_out
    SINK.update(CP.PT_INPUTS,p_sink,T_sink_in)
    h_31    = SINK.hmass()
    s_31    = SINK.smass()
    SINK.update(CP.PT_INPUTS,p_sink,T_sink_out)
    h_35    = SINK.hmass()
    s_35    = SINK.smass()

    m_flow_source = Qdot_source/(h_21-h_23)

    # Allocate Space
    T           = np.zeros(21)
    p           = np.zeros(21)
    h           = np.zeros(21)
    s           = np.zeros(21)
    v           = numpy.zeros(21)

    T_evap      = np.zeros(n_states_evap + n_states_sh - 1)
    T_source    = np.zeros(n_states_evap + n_states_sh - 1)
    delta_T_lm_source = np.zeros(len(T_source)-1)
    UA_source = np.zeros(len(T_source)-1)

    T_sink_top      = np.zeros(n_states_sub+n_states_cond+n_states_desup - 2)
    T_cond_top      = np.zeros(n_states_sub+n_states_cond+n_states_desup - 2)
    delta_T_lm_cond_top = np.zeros(len(T_cond_top)-1)
    UA_cond_top = np.zeros(len(T_cond_top)-1)

    T_ihx1_cold = np.zeros(n_states_ihx)
    T_ihx1_hot  = np.zeros(n_states_ihx)

    T_ihx2_cold = np.zeros(n_states_ihx)
    T_ihx2_hot  = np.zeros(n_states_ihx)
    delta_T_lm_ihx2 = np.zeros(n_states_ihx-1)
    UA_ihx2 = np.zeros(n_states_ihx-1)

    T_sink_ihx3 = np.zeros(n_states_ihx)
    T_cond_ihx3  = np.zeros(n_states_ihx)
    delta_T_lm_ihx3 = np.zeros(n_states_ihx-1)
    UA_ihx3 = np.zeros(n_states_ihx-1)

    #n_sates_cascadeHX = 2
    T_cax_cold = np.zeros(2 * n_sates_cascadeHX)
    T_cax_hot = np.zeros(2 * n_sates_cascadeHX)


    delta_T_lm_cax = np.zeros(len(T_cax_hot)-1)
    UA_cax = np.zeros(len(T_cax_hot)-1)

    # ==================================================================================================================
    # Fix values
    # Bottom cycle Evaporating pressure
    MEDIUM_BOT.update(CP.QT_INPUTS, 0, T_source_out - delta_T_pinch_source_min) # fixed evaporating temperature to 10 deg C
    p_evap_bot = MEDIUM_BOT.p()
    T_intermediate = parameters[2] + 273.15  # T_intermediate # x[2] + 273.15 K
    # Bottom cycle Condensing pressure
    MEDIUM_BOT.update(CP.QT_INPUTS, 0, T_intermediate)
    p_cond_bot = MEDIUM_BOT.p()

    p_mid_bot = p_evap_bot * np.sqrt(p_cond_bot/p_evap_bot) # Medium pressure of the bottom cycle

    # Top cycle evaporating pressure
    MEDIUM_TOP.update(CP.QT_INPUTS, 1, T_intermediate - delta_T_pinch_cascade_HX_min)
    #MEDIUM_TOP.update(CP.QT_INPUTS, 1, T_intermediate - 3)
    p_evap_top = MEDIUM_TOP.p()

    # Disctinction in dew point and bubble point for top cycle evaporator
    MEDIUM_TOP.update(CP.PQ_INPUTS, p_evap_top,1)
    T[12] = MEDIUM_TOP.T()
    h[12] = MEDIUM_TOP.hmass()
    s[12] = MEDIUM_TOP.smass()

    MEDIUM_TOP.update(CP.PQ_INPUTS, p_evap_top,0)
    T21 = MEDIUM_TOP.T()
    h21 = MEDIUM_TOP.hmass()
    s21 = MEDIUM_TOP.smass()

    # Evaporator outlet - Super heater inlet
    MEDIUM_BOT.update(CP.PQ_INPUTS, p_evap_bot,1)
    T[11] = MEDIUM_BOT.T()
    h[11] = MEDIUM_BOT.hmass()
    s[11] = MEDIUM_BOT.smass()

    # Superheater outlet - IHX1 inlet cold side
    T[1] = T_source_in - 5 # fixed T1 = 35 deg C
    MEDIUM_BOT.update(CP.PT_INPUTS, p_evap_bot, T[1])
    h[1] = MEDIUM_BOT.hmass()
    s[1] = MEDIUM_BOT.smass()

    # IHX1 outlet - Com 1 inlet
    T[2] = T[1] + delta_T_SH_IHX1  # delta_T_SH_IHX1  = parameter 1
    MEDIUM_BOT.update(CP.PT_INPUTS, p_evap_bot, T[2])
    h[2] = MEDIUM_BOT.hmass()
    s[2] = MEDIUM_BOT.smass()
    v[2] = 1/MEDIUM_BOT.rhomass()


    # Com 1 outlet - Com 2 inlet
    PR1 = p_mid_bot/p_evap_bot
    eta_is_comp1 = eta_is_comp(PR1)
    MEDIUM_BOT.update(CP.PSmass_INPUTS, p_mid_bot, s[2])
    h_3_is = MEDIUM_BOT.hmass()
    h[3] = (h_3_is - h[2]) / eta_is_comp1 + h[2]
    MEDIUM_BOT.update(CP.HmassP_INPUTS, h[3], p_mid_bot)
    T[3] = MEDIUM_BOT.T()
    s[3] = MEDIUM_BOT.smass()
    v[3] = 1 / MEDIUM_BOT.rhomass()

    # Com 2 outlet - Cascade HX inlet - hot side
    PR2 = p_cond_bot/p_mid_bot
    eta_is_comp2 = eta_is_comp(PR2)
    MEDIUM_BOT.update(CP.PSmass_INPUTS, p_cond_bot, s[3])
    h_4_is = MEDIUM_BOT.hmass()
    h[4] = (h_4_is - h[3]) / eta_is_comp2 + h[3]
    MEDIUM_BOT.update(CP.HmassP_INPUTS, h[4], p_cond_bot)
    T[4] = MEDIUM_BOT.T()
    s[4] = MEDIUM_BOT.smass()
    v[4] = 1 / MEDIUM_BOT.rhomass()

    # Determine minimum condensing pressure of the top cycle to make sure T17 < T16
    T18_1 = T[4]- delta_T_pinch_cascade_HX_min + delta_T_SH_IHX2 + delta_T_pinch_ihx_min # T18 min = T14 + delta_T_pinch_ihx_min
    T18_2 = T_intermediate - delta_T_pinch_cascade_HX_min + delta_T_minSH + delta_T_SH_IHX2 + delta_T_pinch_ihx_min # T[13] min = T12 + minimum SH temperature

    T18   = max(T18_1,T18_2)
    if T18 > T_crit_top:
        print("Minimum top cycle condensing temperature is higher than critical temperature")
    #print('T18', T18 -273.15)
    MEDIUM_TOP.update(CP.QT_INPUTS, 0,T18)
    p_cond_top_min = MEDIUM_TOP.p()
    #print('p_cond_top_min',p_cond_top_min/1e6)
    #print('T18',T18-273.15)
    if p_cond_top_min > p_crit_top:
        print("Minimum top cycle condensing pressure is higher than critical pressure ")

    # Disctinction in Desuperheating, Condensing and Subcooling for Bottom cycle condenser
    MEDIUM_BOT.update(CP.PQ_INPUTS, p_cond_bot, 1)
    h[5] = MEDIUM_BOT.hmass()
    T[5] = MEDIUM_BOT.T()

    MEDIUM_BOT.update(CP.PQ_INPUTS, p_cond_bot, 0)
    h[6] = MEDIUM_BOT.hmass()
    T[6] = MEDIUM_BOT.T()


    # ==================================================================================================================
    # Top cycle condensing pressure guess value, lower, upper bound pressure
    x = np.zeros(1)

    if T_sink_in + delta_T_pinch_sink_min < T_crit_top:
        SubTrans = 'Sub'
        MEDIUM_TOP.update(CP.QT_INPUTS, 0, T_sink_in + delta_T_pinch_sink_min)
        p_min_H = MEDIUM_TOP.p()
        p_max_H = p_crit_top * limit_sc  # limit_sc = 0.95, to make sure satuation line is not too close to the critical point

        try:
            MEDIUM_TOP.update(CP.QT_INPUTS, 1, T_sink_out + delta_T_pinch_sink_min)
            x[0] = MEDIUM_TOP.p()
        except:
            x[0] = p_min_H * 1.1
            print('error 1')

        if x[0] < p_cond_top_min:
            x[0] = p_cond_top_min * 1.05 #To make sure the initial P_cond_top > p_cond_top_min

    else:
        #print('Transcritical guess values')
        SubTrans = 'Trans'
        x[0] = p_crit_top * 1.2
        p_max_H = p_crit_top * 100
        p_min_H = p_crit_top * 1.05


    #===================================================================================================================
    iter = 0
    RelRes      = np.ones(1)
    Res_sink    = np.ones(2)
    y = [0, 1]

    # _________________________________________________________________________________
    #                    Start Newton-Raphson solver
    # __________________________________________________________________________________

    while abs(RelRes) >= tol:
        iter = iter + 1
        for m in range(2):
            p_cond_top = x[0] * (1 + tol * step * y[m]) # step = 1 => deltaP is positive
            #print('Top cycle condensing pressure, Mpa:',p_cond_top/1e6)
            if p_cond_top <= (p_crit_top * limit_sc):
                SubTrans = 'Sub'
            else:
                SubTrans = 'Trans'

            # Cascade HX outlet IHX2 inlet
            T[13] = T[4] - delta_T_pinch_cascade_HX_min
            #print('T[13]',T[13]-273.15)
            MEDIUM_TOP.update(CP.PT_INPUTS, p_evap_top, T[13])
            h[13] = MEDIUM_TOP.hmass()
            s[13] = MEDIUM_TOP.smass()

            #===========================================================================================================
            #Iteration to set pinch at the cascade HX
            error_delta_T_pinch_cascade_HX = 10
            iter_CHX = 1

            while (error_delta_T_pinch_cascade_HX) > error_delta_T:

                # Disctinction in Desuperheating, Condensing and Subcooling for Top cycle condenser
                MEDIUM_TOP.update(CP.PQ_INPUTS, p_cond_top, 1)
                h[16] = MEDIUM_TOP.hmass()
                T[16] = MEDIUM_TOP.T()

                MEDIUM_TOP.update(CP.PQ_INPUTS, p_cond_top, 0)
                h[17] = MEDIUM_TOP.hmass()
                T[17] = MEDIUM_TOP.T()
                #print('h[16] (kJ/kgK):', h[16] / 1000)
                #print('T[16] (deg C):', T[16] -273.15)

                # Compressor 3
                # Iteration to ensure dry compression at outlet of com 3
                error_delta_T_SH = 10
                delta_T_SH = 0
                tic = timer()
                while error_delta_T_SH >= error_delta_T:  # error_delta_T = 1e-2 Allowable error on pinch temperature

                    # IHX2 outlet Compressor 3 inlet
                    T[14] = T[13] + delta_T_SH_IHX2 + delta_T_SH
                    #print('T[14]', T[14] - 273.15)
                    MEDIUM_TOP.update(CP.PT_INPUTS, p_evap_top, T[14])
                    h[14] = MEDIUM_TOP.hmass()
                    s[14] = MEDIUM_TOP.smass()
                    v[14] = 1 / MEDIUM_TOP.rhomass()

                    # Compressor 3 outlet
                    PR3 = p_cond_top/p_evap_top
                    eta_is_comp3 = eta_is_comp(PR3)
                    MEDIUM_TOP.update(CP.PSmass_INPUTS, p_cond_top, s[14])
                    h_15_is = MEDIUM_TOP.hmass()
                    h[15] = (h_15_is - h[14]) / eta_is_comp3 + h[14]
                    MEDIUM_TOP.update(CP.HmassP_INPUTS, h[15], p_cond_top)
                    T[15] = MEDIUM_TOP.T()
                    s[15] = MEDIUM_TOP.smass()
                    v[15] = 1 / MEDIUM_TOP.rhomass()

                    delta_T_SH_summary = [T[14] - T[12], T[15] - T[16]]
                    error_delta_T_SH = delta_T_minSH - min(delta_T_SH_summary) # delta_T_minSH = 5
                    delta_T_SH = delta_T_SH + error_delta_T_SH
                    toc = timer()
                    if toc - tic > 60:
                        print('SH did not converge within 60 s')

                # Subcoler outlet - IHX2 inlet : T17 need to be smaller than T_top_cond

                T[18] = T[14] + delta_T_pinch_ihx_min

                MEDIUM_TOP.update(CP.PT_INPUTS, p_cond_top, T[18])
                h[18] = MEDIUM_TOP.hmass()
                s[18] = MEDIUM_TOP.smass()


                #  IHX2 outlet - TX1 inlet
                h[19] = h[18]-(h[14]-h[13])
                MEDIUM_TOP.update(CP.HmassP_INPUTS, h[19], p_cond_top)
                T[19] = MEDIUM_TOP.T()
                s[19] = MEDIUM_TOP.smass()

                # TX1 outlet - Cascade HX inlet - cold side
                h[20] = h[19]
                MEDIUM_TOP.update(CP.HmassP_INPUTS, h[20], p_evap_top)
                T[20] = MEDIUM_TOP.T()
                s[20] = MEDIUM_TOP.smass()

                x7 = parameters[3]
                # Cascade HX outlet - hot side - IHX3 inlet
                T[7] = T[20] + delta_T_pinch_cascade_HX_min
                if T[5] == T[6] and T21 == T[12]:
                    MEDIUM_BOT.update(CP.PQ_INPUTS, p_cond_bot, x7)
                    h[7] = MEDIUM_BOT.hmass()
                    s[7] = MEDIUM_BOT.smass()
                else:
                    MEDIUM_BOT.update(CP.PT_INPUTS, p_cond_bot, T[7])
                    h[7] = MEDIUM_BOT.hmass()
                    s[7] = MEDIUM_BOT.smass()

                #==========================================================================================================
                # IH3 outlet - IHX1 inlet hot side
                # Determine minimum T[8] to avoid subcool liquid to enter TX2
                MEDIUM_BOT.update(CP.PQ_INPUTS, p_evap_bot,0)
                h_9min = MEDIUM_BOT.hmass()
                h_8min = h_9min + (h[2]-h[1])
                MEDIUM_BOT.update(CP.HmassP_INPUTS, h_8min, p_cond_bot)
                T_8min = MEDIUM_BOT.T()



                if delta_T_SH_IHX1 > 0:
                    T[8] = max(T[2] + delta_T_pinch_ihx_min,T_8min)
                else:
                    T[8] = max(T_sink_in + delta_T_pinch_sink_min,T_8min)

                MEDIUM_BOT.update(CP.PT_INPUTS, p_cond_bot, T[8])
                h[8] = MEDIUM_BOT.hmass()
                s[8] = MEDIUM_BOT.smass()
                #===========================================================================================================
                # Iteration around IHX3 to set pinch
                # To ensure T[7] >= T_32 + delta_T_pinch_IHX_min
                # if T[7] < T_32 + delta_T_pinch_IHX_min: increase T[8]
                # Initial guess value
                iter_IHX3 = 1
                error_delta_T_pinch_ihx3 = 1
                while (error_delta_T_pinch_ihx3) > error_delta_T:
                    # IHX1 outlet hot side - TX2 inlet
                    h[9] = h[8] - (h[2]-h[1])
                    MEDIUM_BOT.update(CP.HmassP_INPUTS, h[9], p_cond_bot)
                    T[9] = MEDIUM_BOT.T()
                    s[9] = MEDIUM_BOT.smass()

                    # TX2 outlet - Evaporator inlet
                    h[10] = h[9]
                    MEDIUM_BOT.update(CP.HmassP_INPUTS, h[10], p_evap_bot)
                    T[10] = MEDIUM_BOT.T()
                    s[10] = MEDIUM_BOT.smass()

                    # Mass flow rate
                    m_flow_bot = Qdot_source / (h[1] - h[10])
                    Q_dot_IHX3 = m_flow_bot * (h[7] - h[8])

                    Q_dot_cascade = m_flow_bot * (h[4] - h[7])

                    m_flow_top = Q_dot_cascade / (h[13] - h[20])

                    Q_dot_des_top = m_flow_top * (h[15] - h[16])
                    Q_dot_cond_top = m_flow_top * (h[16] - h[17])
                    Q_dot_sub_top = m_flow_top * (h[17] - h[18])

                    Q_dot_sink = Q_dot_des_top + Q_dot_cond_top + Q_dot_sub_top + Q_dot_IHX3
                    m_flow_sink = Q_dot_sink / (h_35 - h_31)

                    h_32 = h_31 + Q_dot_IHX3 / m_flow_sink
                    SINK.update(CP.HmassP_INPUTS, h_32, p_sink)
                    T_32 = SINK.T()
                    s_32 = SINK.smass()
                    # Calculation of discretized states on the IHX3
                    h_cond_bot = np.linspace(h[8], h[7], n_states_ihx)
                    h_sink_ihx3 = np.linspace(h_31, h_32, n_states_ihx)

                    # Temperature profile of refrigerant and heat sink fluid in the IHX3
                    #for i in range(1, len(h_cond_bot) - 1):
                    for i in range(n_states_ihx):
                        try:
                            MEDIUM_BOT.update(CP.HmassP_INPUTS, h_cond_bot[i], p_cond_bot)
                            T_cond_ihx3[i] = MEDIUM_BOT.T()
                            SINK.update(CP.HmassP_INPUTS, h_sink_ihx3[i], p_sink)
                            T_sink_ihx3[i] = SINK.T()
                        except:
                            T_cond_ihx3[i] = math.nan
                            T_sink_ihx3[i] = math.nan

                    for item in T_cond_ihx3:
                        if (math.isnan(item)) == True:
                            print('Points did  not converge')

                    # Temperature Differences IHX3
                    delta_T_ihx3        = T_cond_ihx3-T_sink_ihx3
                    error_delta_T_pinch_ihx3 = delta_T_pinch_sink_min - min(delta_T_ihx3) # delta_T_pinch_sink_min = 20
                    if (error_delta_T_pinch_ihx3) <= error_delta_T:
                        break
                    else:
                        T[8] = T[8] + error_delta_T_pinch_ihx3
                        MEDIUM_BOT.update(CP.PT_INPUTS, p_cond_bot, T[8])
                        h[8] = MEDIUM_BOT.hmass()
                        s[8] = MEDIUM_BOT.smass()
                        iter_IHX3 = iter_IHX3 + 1
                        #print('IHX3 iteration:',iter_IHX3)
                    if iter_IHX3 > 100:
                        print('IHX3 did not converge: ', error_delta_T_pinch_ihx3, ' K')

                # Determine the heat transferred at point 12 and point 5
                Q_12 = m_flow_top*(h[12]-h[20])  # W
                Q_5  = m_flow_bot*(h[5]-h[7])    # W
                # Determine the corresponding temperature at Q_5 and Q_12 in the T-H diagram
                h5_top = h[20] + Q_5/m_flow_top
                MEDIUM_TOP.update(CP.HmassP_INPUTS, h5_top, p_evap_top)
                T5_top = MEDIUM_TOP.T()
                h12_bot = h[7] + Q_12/m_flow_bot

                MEDIUM_BOT.update(CP.HmassP_INPUTS, h12_bot, p_cond_bot)
                T12_bot = MEDIUM_BOT.T()
                delta_T_cascade_HX = [T[5]-T5_top,T12_bot-T[12]]
                #print('delta_T_cascade_HX',delta_T_cascade_HX)
                error_delta_T_pinch_cascade_HX = delta_T_pinch_cascade_HX_min -min(delta_T_cascade_HX )
                #print('error_delta_T_pinch_cascade_HX',error_delta_T_pinch_cascade_HX)
                if error_delta_T_pinch_cascade_HX <= error_delta_T:
                    break
                else:
                    T[13] = T[13]-error_delta_T_pinch_cascade_HX
                    MEDIUM_TOP.update(CP.PT_INPUTS, p_evap_top, T[13])
                    h[13] = MEDIUM_TOP.hmass()
                    s[13] = MEDIUM_TOP.smass()
                    iter_CHX = iter_CHX+1
                if iter_CHX > 100:
                    print('Cascade HX did not converge: ', error_delta_T_pinch_cascade_HX, ' K')

            # Heat sink side
            h_33 = h_31 + (Q_dot_IHX3 + Q_dot_sub_top) / m_flow_sink

            h_34 = h_31 + (Q_dot_IHX3 + Q_dot_sub_top + Q_dot_cond_top) / m_flow_sink

            SINK.update(CP.HmassP_INPUTS, h_33, p_sink)
            T_33 = SINK.T()

            SINK.update(CP.HmassP_INPUTS, h_34, p_sink)
            T_34 = SINK.T()


            # Calculation of discretized states on the top cycle condenser
            h_cond_sub       = np.linspace(h[18], h[17], n_states_sub)
            h_sink_sub       = np.linspace(h_32, h_33, n_states_sub)

            h_cond_cond      = np.linspace(h[17], h[16], n_states_cond)
            h_sink_cond      = np.linspace(h_33, h_34, n_states_cond)

            h_cond_desup     = np.linspace(h[16], h[15], n_states_desup)
            h_sink_desup     = np.linspace(h_34, h_35, n_states_desup)

            h_cond_top = np.hstack((h_cond_sub[0:n_states_sub - 1], h_cond_cond[0:n_states_cond - 1], h_cond_desup))
            h_sink_top = np.hstack((h_sink_sub[0:n_states_sub - 1], h_sink_cond[0:n_states_cond - 1], h_sink_desup))

            # Temperature profile of refrigerant and heat sink fluid in the top cycle condenser
            T_cond_top[0] = T[18]
            T_sink_top[0] = T_32
            T_cond_top[-1] = T[15]  # -1 is the index of the last object in a array
            T_sink_top[-1] = T_35

            for i in range(1, len(h_cond_top) - 1):
                try:
                    MEDIUM_TOP.update(CP.HmassP_INPUTS, h_cond_top[i], p_cond_top)
                    T_cond_top[i] = MEDIUM_TOP.T()
                    SINK.update(CP.HmassP_INPUTS, h_sink_top[i], p_sink)
                    T_sink_top[i] = SINK.T()
                except:
                    T_cond_top[i] = math.nan
                    T_sink_top[i] = math.nan

            for item in T_cond_top:
                if (math.isnan(item)) == True:
                    print('Points did  not converge')

            # Temperature Differences in the top cycle condenser
            delta_T_sink_top        = T_cond_top - T_sink_top
            #print('T_cond_top:',T_cond_top-273.15)
            #print('T_sink_top:',T_sink_top-273.15)
            #print('Delta_T_sink_top_condenser:',delta_T_sink_top)

            #Res_sink[m] = delta_T_pinch_sink_min - min(delta_T_sink_top)
            #RelRes = (delta_T_pinch_sink_min - min(delta_T_sink_top)) / delta_T_pinch_sink_min
            Res_sink[m] = delta_T_pinch_sink_min - min(delta_T_sink_top[1:]) # pinch point never at T18
            RelRes = (delta_T_pinch_sink_min - min(delta_T_sink_top[1:])) / delta_T_pinch_sink_min

            #print('RelRes:',RelRes)
            #print('Res_sink:',m, Res_sink[m])

            if abs(RelRes) < tol:
                break

        if iter > iter_max:
            print('Model did not converge with ', iter_max,' iterations')


        if (x[0] - p_cond_top_min)/1e6 <= tol and iter > iter_max: # p_cond_top_min: minimum condensing pressure of the top cycle to make sure T17 < T16 devide buy 1e6 to convert the unit of pressure to MPa
        #if (p_cond_top - p_cond_top_min)/1e6 < tol: # p_cond_top_min: minimum condensing pressure of the top cycle to make sure T17 < T16
            print('Solver reach minimum condensing pressure of the top cycle')

            break # the solution is at p_cond_top = p_cond_top_min*(1+tol) that is to make sure delta_T_sink_min

        if abs(RelRes) >= tol: # update pressure
            #x[0] = x[0]-x[0]*tol/(Res_sink[0]-Res_sink[1])*Res_sink[0] # use when deltaP is negative
            x[0] = x[0] + x[0] * tol / (Res_sink[0] - Res_sink[1]) * Res_sink[0] # use when deltaP is positive
            # Lower bound pressure
            p_min_H = p_cond_top_min * (1 + 0.05 * tol * step) # to make sure P_cond_top always higher than p_cond_top_min to avoid h17<h18

            #p_min_H = p_cond_top_min*(1- 1.05 * tol * step) # step = -1, to make sure P - deltaP still higer than P_min_H
            x[0] = max(x[0], p_min_H)
            #print(' p_min_H', p_min_H/1e6)
            # Upper bound pressure
            x[0] = min(x[0], p_max_H)
        #print('tolerance:',tol)

    # =========================================================================================================
    # Compressor work
    W_dot_comp1 = m_flow_bot * (h[3] - h[2])
    W_dot_comp2 = m_flow_bot * (h[4] - h[3])
    W_dot_comp3 = m_flow_top * (h[15] - h[14])
    W_dot_comp = W_dot_comp1 + W_dot_comp2 + W_dot_comp3
    W_dot = W_dot_comp / eta_motor

    # =========================================================================================================
    # Enthalpy profile of refrigerant on the heat sink side
    Q_dot_sink_ihx3 = m_flow_bot * (h_cond_bot - h[8])
    Q_dot_sink_top_cond = m_flow_top * (h_cond_top - (h[18] - (m_flow_bot / m_flow_top) * (h[7] - h[8])))
    Q_dot_sink_cum = np.hstack((Q_dot_sink_ihx3, Q_dot_sink_top_cond))
    T_cond = np.hstack((T_cond_ihx3, T_cond_top))
    T_sink = np.hstack((T_sink_ihx3, T_sink_top))

    # Calculate UA for top cycle condenser
    Q_dot_cond_top_cum = m_flow_top * (h_cond_top-h[18])
    for i in range(len(h_cond_top)-1):
        if (T_cond_top[i] - T_sink_top[i]) / (T_cond_top[i + 1] - T_sink_top[i + 1]) == 1:
            delta_T_lm_cond_top[i] = T_cond_top[i] - T_sink_top[i]
        else:
            delta_T_lm_cond_top[i] = ((T_cond_top[i] - T_sink_top[i]) - (T_cond_top[i + 1] - T_sink_top[i + 1])) / np.log((T_cond_top[i] - T_sink_top[i]) / (T_cond_top[i + 1] - T_sink_top[i + 1]))
        UA_cond_top[i] = (Q_dot_cond_top_cum[i + 1] - Q_dot_cond_top_cum[i]) / delta_T_lm_cond_top[i]

    # Calculate UA for IHX3
    for i in range(len(h_cond_bot)-1):
        if (T_cond_ihx3[i] - T_sink_ihx3[i]) / (T_cond_ihx3[i + 1] - T_sink_ihx3[i + 1]) == 1:
            delta_T_lm_ihx3[i] = T_cond_ihx3[i] - T_sink_ihx3[i]
        else:
            delta_T_lm_ihx3[i] = ((T_cond_ihx3[i] - T_sink_ihx3[i]) - (T_cond_ihx3[i + 1] - T_sink_ihx3[i + 1])) / np.log((T_cond_ihx3[i] - T_sink_ihx3[i]) / (T_cond_ihx3[i + 1] - T_sink_ihx3[i + 1]))
        UA_ihx3[i] = (Q_dot_sink_ihx3[i + 1] - Q_dot_sink_ihx3[i]) / delta_T_lm_ihx3[i]

    #===================================================================================================================
    # Heat source side
    Q_dot_evap = m_flow_bot * (h[11] - h[10])
    Q_dot_super = m_flow_bot * (h[1] - h[11])

    h_22 = h_23 + Q_dot_evap / m_flow_source
    SOURCE.update(CP.HmassP_INPUTS, h_22, p_source)
    T_22 = SOURCE.T()

    h_evap_evap = np.linspace(h[10], h[11], n_states_evap)
    h_source_evap = np.linspace(h_23, h_22, n_states_evap)

    h_evap_sh = np.linspace(h[11], h[1], n_states_sh)
    h_source_sh = np.linspace(h_22, h_21, n_states_sh)

    h_evap = np.hstack((h_evap_evap[0:n_states_evap - 1], h_evap_sh))
    h_source = np.hstack((h_source_evap[0:n_states_evap - 1], h_source_sh))

    # Temperature Profile -> from 10 to 1

    for i in range(len(h_evap)):
        try:
            MEDIUM_BOT.update(CP.HmassP_INPUTS, h_evap[i], p_evap_bot)
            T_evap[i] = MEDIUM_BOT.T()
            SOURCE.update(CP.HmassP_INPUTS, h_source[i], p_source)
            T_source[i] = SOURCE.T()
        except:
            T_evap[i] = math.nan
            T_source[i] = math.nan

    # delta_T_source = T_source - T_evap
    Q_dot_source_cum = m_flow_bot * (h_evap - h[10])
    # Calculate UA value for the bottom cycle evaporator
    for i in range(len(h_evap)-1):
        if (T_source[i]-T_evap[i]) / (T_source[i+1]-T_evap[i+1]) == 1:
            delta_T_lm_source[i] = (T_source[i]-T_evap[i])
        else:
            delta_T_lm_source[i] = ((T_source[i] - T_evap[i]) - (T_source[i + 1] - T_evap[i + 1])) / np.log((T_source[i] - T_evap[i]) / (T_source[i + 1] - T_evap[i + 1]))

        UA_source[i] = (Q_dot_source_cum[i + 1] - Q_dot_source_cum[i]) / delta_T_lm_source[i]

    # =================================================================================================================
    #h_cax_bot_cond = np.linspace(h[6], h[4], n_sates_cascadeHX)
    #h_cax_top_evap = np.linspace(h[19], h4_top, n_sates_cascadeHX)
    #h_cax_bot_des = np.linspace(h[4], h[3], n_sates_cascadeHX)
    #h_cax_top_sh = np.linspace(h[11], h[12], n_sates_cascadeHX)
    #h_cax_bot = np.hstack((h_cax_bot_cond[0:n_sates_cascadeHX - 1], h_cax_bot_des))
    #h_cax_top = np.hstack((h_cax_top_evap[0:n_sates_cascadeHX - 1], h_cax_top_sh))
    #for i in range(len(h_cax_top)):
    #    try:
    #        MEDIUM_TOP.update(CP.HmassP_INPUTS, h_cax_top[i], p_evap_top)
    #        T_cax_cold[i] = MEDIUM_TOP.T()
    #        MEDIUM_BOT.update(CP.HmassP_INPUTS, h_cax_bot[i], p_cond_bot)
    #        T_cax_hot[i] = MEDIUM_BOT.T()
    #    except:
    #        T_cax_cold[i] = math.nan
    #        T_cax_hot[i] = math.nan
    #Q_dot_cax_bot_cum = m_flow_bot * (h_cax_bot - h[6])
    #Q_dot_cax_top_cum = m_flow_top * (h_cax_top - h[19])




    # Cascade heat exchanger
    if m_flow_bot * (h[5] - h[7]) < m_flow_top * (h[12] - h[20]):
        h_cax_bot_1 = np.linspace(h[7], h[5], n_sates_cascadeHX)
        h_cax_bot_2 = np.linspace(h[5], h12_bot,2)
        h_cax_bot_3 = np.linspace(h12_bot, h[4], n_sates_cascadeHX)
        h_cax_top_1 = np.linspace(h[20], h5_top, n_sates_cascadeHX)
        h_cax_top_2 = np.linspace(h5_top,h[12],2)
        h_cax_top_3 = np.linspace(h[12], h[13], n_sates_cascadeHX)
        #h_cax_top_3 = np.linspace(h[11], h[12], 2)
    elif m_flow_bot * (h[5] - h[7]) == m_flow_top * (h[12] - h[20]):
        h_cax_bot_1 = np.linspace(h[7], h[5], n_sates_cascadeHX)
        h_cax_bot_2 = []
        h_cax_bot_3 = np.linspace(h[5], h[4], n_sates_cascadeHX)
        h_cax_top_1 = np.linspace(h[20], h[12], n_sates_cascadeHX)
        h_cax_top_2 = []
        h_cax_top_3 = np.linspace(h[12], h[13], n_sates_cascadeHX)
    else:
        h_cax_bot_1 = np.linspace(h[7], h12_bot, n_sates_cascadeHX)
        h_cax_bot_2 = np.linspace(h12_bot, h[5],2)
        h_cax_bot_3 = np.linspace(h[5], h[4], n_sates_cascadeHX)
        h_cax_top_1 = np.linspace(h[20], h[12], n_sates_cascadeHX)
        h_cax_top_2 = np.linspace(h[12], h5_top,2)
        h_cax_top_3 = np.linspace(h5_top, h[13], n_sates_cascadeHX)

    if m_flow_bot * (h[5] - h[7]) == m_flow_top * (h[12] - h[20]):
        h_cax_bot = np.hstack((h_cax_bot_1[0:n_sates_cascadeHX - 1], h_cax_bot_3))
        h_cax_top = np.hstack((h_cax_top_1[0:n_sates_cascadeHX - 1], h_cax_top_3))
        T_cax_cold = np.zeros(2 * n_sates_cascadeHX - 1)
        T_cax_hot = np.zeros(2 * n_sates_cascadeHX - 1)
    else:
        h_cax_bot = np.hstack((h_cax_bot_1[0:n_sates_cascadeHX - 1],h_cax_bot_2[0:1], h_cax_bot_3))
        h_cax_top = np.hstack((h_cax_top_1[0:n_sates_cascadeHX - 1], h_cax_top_2[0:1], h_cax_top_3))

    for i in range(len(h_cax_top)):
        try:
            MEDIUM_TOP.update(CP.HmassP_INPUTS, h_cax_top[i], p_evap_top)
            T_cax_cold[i] = MEDIUM_TOP.T()
            MEDIUM_BOT.update(CP.HmassP_INPUTS, h_cax_bot[i], p_cond_bot)
            T_cax_hot[i] = MEDIUM_BOT.T()
        except:
            T_cax_cold[i] = math.nan
            T_cax_hot[i] = math.nan
    Q_dot_cax_bot_cum = m_flow_bot * (h_cax_bot - h[7])

    # Calculate UA for cascade HX

    for i in range(len(h_cax_bot)-1):
        if (T_cax_hot[i]-T_cax_cold[i]) / (T_cax_hot[i+1]-T_cax_cold[i+1]) == 1:
            delta_T_lm_cax[i] = (T_cax_hot[i]-T_cax_cold[i])
        else:
            delta_T_lm_cax[i] = ((T_cax_hot[i] - T_cax_cold[i]) - (T_cax_hot[i + 1] - T_cax_cold[i + 1])) / np.log((T_cax_hot[i] - T_cax_cold[i]) / (T_cax_hot[i + 1] - T_cax_cold[i + 1]))

        UA_cax[i] = (Q_dot_cax_bot_cum[i + 1] - Q_dot_cax_bot_cum[i]) / delta_T_lm_cax[i]

    # IHX2
    h_ihx2_hot = np.linspace(h[19], h[18], n_states_ihx)
    h_ihx2_cold = np.linspace(h[13], h[14], n_states_ihx)
    for i in range(len(h_ihx2_hot)):
        try:
            MEDIUM_TOP.update(CP.HmassP_INPUTS, h_ihx2_hot[i], p_cond_top)
            T_ihx2_hot[i] = MEDIUM_TOP.T()
            MEDIUM_TOP.update(CP.HmassP_INPUTS, h_ihx2_cold[i], p_evap_top)
            T_ihx2_cold[i] = MEDIUM_TOP.T()
        except:
            T_ihx2_hot[i] = math.nan
            T_ihx2_cold[i] = math.nan
    Q_dot_ihx2_cum = m_flow_top * (h_ihx2_hot - h[19])

    # Calculate UA for IHX2
    for i in range(len(h_ihx2_hot)-1):
        if (T_ihx2_hot[i]-T_ihx2_cold[i]) / (T_ihx2_hot[i+1]-T_ihx2_cold[i+1]) == 1:
            delta_T_lm_ihx2[i] = (T_ihx2_hot[i]-T_ihx2_cold[i])
        else:
            delta_T_lm_ihx2[i] = ((T_ihx2_hot[i] - T_ihx2_cold[i]) - (T_ihx2_hot[i + 1] - T_ihx2_cold[i + 1])) / np.log((T_ihx2_hot[i] - T_ihx2_cold[i]) / (T_ihx2_hot[i + 1] - T_ihx2_cold[i + 1]))

        UA_ihx2[i] = (Q_dot_ihx2_cum[i + 1] - Q_dot_ihx2_cum[i]) / delta_T_lm_ihx2[i]

    #===================================================================================================================
    # Calculate out put
    #Dead state
    To = 15 + 273.15  # 15 oC
    Po = 1e5 # 1 bar
    SOURCE.update(CP.PT_INPUTS,Po,To)
    ho = SOURCE.hmass()
    so = SOURCE.smass()
    # Exergy efficiency
    E_dot_sink_in =  m_flow_sink * (h_31-ho - To * (s_31-so)) # J/kg
    E_dot_sink_out = m_flow_sink * (h_35 - ho - To * (s_35 - so))  # J/kg
    E_dot_source_in =  m_flow_source * (h_21-ho - To*(s_21-so)) # J/kg
    E_dot_source_out = m_flow_source * (h_23 - ho - To * (s_23 - so))  # J/k
    exergy_eff = (E_dot_sink_out-E_dot_sink_in)/(E_dot_source_in-E_dot_source_out + W_dot)
    # Lorren efficiency
    T_sink_av = (T_sink_out-T_sink_in)/np.log(T_sink_out/T_sink_in)
    T_source_av = (T_source_in-T_source_out)/np.log(T_source_in/T_source_out)
    COP_Lor = T_sink_av/(T_sink_av-T_source_av)
    #COP
    COP = Q_dot_sink / W_dot
    # Second law efficiency
    n_Lorren = (COP/COP_Lor) * 100
    # Exergy destruction at the top cycle condenser
    E_des_top_cond = m_flow_top * (h[15]-h[18] - To *(s[15] - s[18])) + m_flow_sink * (h_32 - h_35 - To *(s_32 - s_35)) # J/kg
    # Exergy destruction at IHX3
    E_des_IHX3 = m_flow_bot * (h[7] - h[8] - To * (s[7] - s[8])) + m_flow_sink * (h_31 - h_32 - To * (s_31 - s_32))  # J/kg
    # Exergy destruction at cascade HX
    E_des_CHX = m_flow_bot * (h[4] - h[7] - To * (s[4] - s[7])) + m_flow_top * (h[20] - h[13] - To * (s[20] - s[13]))  # J/kg
    # Exergy destruction at bottom cycle Evaporator
    E_des_bot_evap = m_flow_bot * (h[10] - h[1] - To * (s[10] - s[1])) + m_flow_source * (h_21 - h_23 - To * (s_21 - s_23))  # J/kg

    # Volumetric effeciency of compressor
    eta_vol_com1 = 1 - 0.06 * (v[2]/v[3] -1)
    eta_vol_com2 = 1 - 0.06 * (v[3]/v[4] -1)
    eta_vol_com3 = 1 - 0.06 * (v[14]/v[15] -1)



    # _________________________________________________________________________________
    #                    Collect the output
    # _________________________________________________________________________________

    # update state points
    p[0] = p_evap_bot
    p[1] = p_evap_bot
    p[2] = p_evap_bot
    p[3] = p_mid_bot
    p[4] = p_cond_bot
    p[5] = p_cond_bot
    p[6] = p_cond_bot
    p[7] = p_cond_bot
    p[8] = p_cond_bot
    p[9] = p_cond_bot
    p[10] = p_evap_bot
    p[11] = p_evap_bot
    p[12] = p_evap_top
    p[13] = p_evap_top
    p[14] = p_evap_top
    p[15] = p_cond_top
    p[16] = p_cond_top
    p[17] = p_cond_top
    p[18] = p_cond_top
    p[19] = p_cond_top
    p[20] = p_evap_top

    T[0]= T[1]
    h[0] = h[1]
    #for i in range(11):
        #print('Elthalpy at point     ', i, '(kJ/kg)     :', h[i] / 1000)
        #print('Temperature at point  ', i, '(oC)        :', T[i] - 273.15)
    #    MEDIUM_BOT.update(CP.HmassT_INPUTS, h[i], T[i])
    #    s[i] = MEDIUM_BOT.smass() / 1000  # kJ/kgK
        #print('Eltropy at point     ', i, '(kJ/kgK)     :', s[i] / 1000)

    #for i in range(11, 20):
    #    MEDIUM_TOP.update(CP.HmassT_INPUTS, h[i], T[i])
    #    s[i] = MEDIUM_TOP.smass() / 1000  # kJ/kgK

    class Output:
        pass



    output = Output()

    output.Qdot_sink            = Q_dot_sink    # W
    output.Qdot_source          = Qdot_source   # W
    output.mdot_sink            = m_flow_sink
    output.mdot_source          = m_flow_source
    output.mdot_ref_top         = m_flow_top
    output.mdot_ref_bot         = m_flow_bot
    output.p_evap_top           = p_evap_top
    output.p_cond_top           = p_cond_top
    output.p_evap_bot           = p_evap_bot
    output.p_cond_bot           = p_cond_bot
    output.p_mid_bot            = p_mid_bot

    output.UA_source            = UA_source
    output.UA_cond_top          = UA_cond_top
    output.UA_cax               = UA_cax
    output.UA_ihx3              = UA_ihx3
    output.UA_ihx2              = UA_ihx2

    output.T_suction_com1       = T[2] # suction temperature
    output.T_suction_com2       = T[3]  # suction temperature
    output.T_suction_com3       = T[14]  # suction temperature

    output.T_discharge_com2     = T[4]  # discharge temperature
    output.T_discharge_com3     = T[15]  # discharge temperature

    #output.T_ref_out_cond = T[18]  # subcooler outlet
    #output.T_ref_in_eva = T[10]  # evaporator inlet

    output.T_cond_top = T[17]  # condensing temperature of top cycle
    output.eta_vol_com1 = eta_vol_com1
    output.eta_vol_com2 = eta_vol_com2
    output.eta_vol_com3 = eta_vol_com3

    output.eta_is_comp1 = eta_is_comp1
    output.eta_is_comp2 = eta_is_comp2
    output.eta_is_comp3 = eta_is_comp3


    output.specific_volume_2 = v[2]
    output.specific_volume_3 = v[3]
    output.specific_volume_4 = v[4]
    output.specific_volume_14 = v[14]
    output.specific_volume_15 = v[15]

    #MEDIUM_BOT.update(CP.HmassP_INPUTS,h[2], p_evap_bot)
    #output.V_suction_bot     = m_flow_bot/MEDIUM_BOT.rhomass() # [ m3/s]
    #output.VHC_bot = Q_dot_sink / output.V_suction_bot  # [J/m3]
    output.V_suction_com1 = m_flow_bot * v[2]  # [ m3/s]
    output.V_suction_com2 = m_flow_bot * v[3]  # [ m3/s]
    output.V_suction_com3 = m_flow_top * v[14]  # [ m3/s]



    #MEDIUM_TOP.update(CP.HmassP_INPUTS,h[14], p_evap_top)
    #output.V_suction_top     = m_flow_top/MEDIUM_TOP.rhomass() # [ m3/s]
    #output.VHC_top = Q_dot_sink/output.V_suction_top        # [J/m3]



    #output.delta_T_SC_top    = T[17]-T[18] #subcooling_top
    #output.delta_T_SC_bot = T[6] - T[7]  # subcooling_bot
    #output.delta_T_SH    = T[1]-T[11] # superheating

    output.W_dot         = W_dot

    output.exergy_eff = exergy_eff
    output.n_Lorren = n_Lorren
    output.E_des_top_cond = E_des_top_cond/1000
    output.E_des_IHX3 = E_des_IHX3/1000
    output.E_des_CHX  = E_des_CHX/1000
    output.E_des_bot_evap = E_des_bot_evap/1000



    output.Q_dot_sink_cum = Q_dot_sink_cum
    output.Q_dot_source_cum = Q_dot_source_cum
    output.Q_dot_cax_bot_cum = Q_dot_cax_bot_cum
    output.T_source = T_source
    output.T_sink = T_sink
    output.T_evap = T_evap
    output.T_cond = T_cond
    output.T_cax_cold = T_cax_cold
    output.T_cax_hot = T_cax_hot
    output.h_cond = h_cond_top
    output.h =h
    #output.s = s
    output.T =T
    output.p = p
    output.SubTrans = SubTrans
    output.T_31 = T_31
    output.T_32 = T_32
    output.T_33 = T_33
    output.T_34 = T_34
    output.T_35 = T_35

    return [COP,output]
    #return COP
















