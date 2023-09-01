# -------------- EXPLANATION OF THE STATE POINTS -----------------------%
# Simple VAPOUR COMPRESSION cycle, with compressor, condenser, valve, evaporator

# Refrigerant
#     1:   Compressor Inlet
#     2:   Compressor Outlet/Desuperheater Inlet
#     3:   Desuperheater Outlet/Condenser Inlet
#     4:   Condenser Outlet/Subcooler Inlet
#     5:   Subcooler Outlet/Throttling Valve Inlet
#     6:   Throttling Valve Outlet/Evaporator Inlet
#     7:   Evaporator Outlet/Superheater Inlet
#     1:   Superheater Outlet
# Heat Sink
#     31:  Heat Sink Inlet -> Inlet of Subcooler
#     32:  Outlet of Subcooler/Inlet of Condenser
#     33:  Outlet of Condenser/Inlet of Desuperheater
#     34:  Heat Sink Outlet -> Outlet of Desuperheater

# Heat source
#     21:  Heat Source Inlet -> Superheater Inlet
#     22:  Superheater Outlet/Inlet Evaporator
#     23:  Heat Source Outlet -> Outlet Evaporator

# Ref: Benjamin Zühlsdorf , Jonas Kjær Jensen , Brian Elmegaard: Heat pump working fluid selection—economic and thermodynamic comparison of criteria and boundary conditions

import CoolProp as CP
import numpy as np
import math
from timeit import default_timer as timer
from HP_Model_Input import model, fluids, source, ref, sink, comp


def single_state_std_hp(refrigerant):
    print('RUN SINGLE STATE STD CYCLE')
    # ___________________________________Heat Pump Parameter__________________________________

    # Discretization
    n_states_evap = 10
    n_states_sh = 2
    n_states_cond = 10
    n_states_desup = 2
    n_states_sub = 2
    #n_states_ihx = 5

    # ___________________________________Numerical input parameter____________________________

    # Heat pump design models:
    limit_sc = model.limit_sc  # 0.95  Limit above which process is considered as supercritcal
    error_delta_T = model.error_delta_T  # 1e-2 Allowable error on pinch temperature
    iter_max = model.iter_max  # Number of maximum iterations

    tol = model.tol_HP
    step = model.step
    # relaxation      = model.relaxation

    # _________________________________________ Working fluid ____________________________________

    MEDIUM1 = refrigerant.ref_3
    SOURCE = fluids.source
    SINK = fluids.sink
    # _________________________________________ Input thermodynamics _______________________________

    T_source_in = source.T_in
    T_source_out = source.T_out_design
    p_source = source.p_in
    delta_T_pinch_source_min = source.delta_T_pinch_min
    delta_T_minSH = ref.delta_T_minSH
    Qdot_source = source.Qdot_design
    T_sink_in = sink.T_in
    T_sink_out = sink.T_out_design
    p_sink = sink.p_in
    delta_T_pinch_sink_min = sink.delta_T_pinch_min
    #delta_T_pinch_ihx_min = ref.delta_T_pinch_IHX_min

    eta_is_comp = comp.eta_is
    eta_motor = comp.eta_motor

    # _________________________ Run the heat pump function   __________________________________________
    #  Calculate constants
    MEDIUM1.update(CP.PT_INPUTS, 101325, 300)
    T_crit = MEDIUM1.T_critical()
    p_crit = MEDIUM1.p_critical()

    T_21 = T_source_in
    T_23 = T_source_out

    SOURCE.update(CP.PT_INPUTS, p_source, T_source_in)
    h_21 = SOURCE.hmass()
    cpW = SOURCE.cpmass()
    muW = SOURCE.viscosity()
    kW = SOURCE.conductivity()

    SOURCE.update(CP.PT_INPUTS, p_source, T_source_out)
    h_23 = SOURCE.hmass()

    T_31 = T_sink_in
    T_34 = T_sink_out
    SINK.update(CP.PT_INPUTS, p_sink, T_sink_in)
    h_31 = SINK.hmass()
    SINK.update(CP.PT_INPUTS, p_sink, T_sink_out)
    h_34 = SINK.hmass()

    m_flow_source = Qdot_source / (h_21 - h_23)

    # Allocate Space
    T = np.zeros(8)
    p = np.zeros(8)
    h = np.zeros(8)
    Q = np.zeros(8)
    s = np.zeros(8)
    T_evap = np.zeros(n_states_evap+n_states_sh - 1)
    T_source = np.zeros(n_states_evap+n_states_sh - 1)
    T_sink = np.zeros(n_states_sub + n_states_cond + n_states_desup - 2)
    T_cond = np.zeros(n_states_sub + n_states_cond + n_states_desup - 2)

    #===================================================================================================================
    # Guess values - lower, upper bound pressure
    x = np.zeros(2)
    # Evaporating pressure guess value
    #if T_source_in - delta_T_pinch_source_min < T_crit:
    #    MEDIUM1.update(CP.QT_INPUTS, 1, T_source_in-delta_T_pinch_source_min)
    #    x[0]          = MEDIUM1.p()
    #elif T_source_out - delta_T_pinch_source_min < T_crit:
    #    MEDIUM1.update(CP.QT_INPUTS, 0, T_source_out - delta_T_pinch_source_min)
    #    x[0] = MEDIUM1.p()
    #else:
    #    x[0] = p_crit * limit_sc

    if T_source_out - delta_T_pinch_source_min < T_crit:
        MEDIUM1.update(CP.QT_INPUTS, 0, T_source_out - delta_T_pinch_source_min)
        x[0] = MEDIUM1.p()
    elif T_source_in - delta_T_pinch_source_min < T_crit:
        MEDIUM1.update(CP.QT_INPUTS, 1, T_source_in-delta_T_pinch_source_min)
        x[0]          = MEDIUM1.p()
    else:
        x[0] = p_crit * limit_sc

    p_min_L = 0.1 * 1e3  # 100 Pa
    p_max_L = p_crit * limit_sc

    # Condensing pressure guess value

    MEDIUM1.update(CP.PQ_INPUTS, x[0], 1)
    s[0] = MEDIUM1.smass()
    MEDIUM1.update(CP.SmassT_INPUTS, s[0], T_sink_out+delta_T_pinch_sink_min)
    x[1] = MEDIUM1.p()


    if x[1] >= p_crit:
        print('Transcritical guess values')
        SubTrans        = 'Trans'
        p_max_H         = p_crit * 20
        p_min_H = p_crit * 1.05
    else:
        print('Subcritical guess values')
        SubTrans        = 'Sub'
        MEDIUM1.update(CP.QT_INPUTS, 0, T_sink_in + delta_T_pinch_sink_min)
        p_min_H = MEDIUM1.p()
        p_max_H = p_crit * limit_sc  # limit_sc = 0.95, to make sure satuation line is not too close to the critical point
    #===================================================================================================================

    iter = 0
    RelRes = np.ones(2)
    Res_source = np.ones(3)
    Res_sink = np.ones(3)
    y = [0, 1, 0]
    z = [0, 0, 1]
    step = -step  # step =1

    # _________________________________________________________________________________
    #                    Start Newton-Raphson solver
    # __________________________________________________________________________________

    while max(abs(RelRes)) >= tol:
        iter = iter + 1
        print(iter)
        for m in range(3):
            # Definition of pressure levels
            p_evap = x[0] * (1 + tol * step * y[m])
            p_cond = x[1] * (1 + tol * step * z[m])
            # Definition of sub-/transcritical process
            if p_cond <= (p_crit * limit_sc):
                SubTrans = 'Sub'
            else:
                SubTrans = 'Trans'

            # Subcooler
            # Subcooler Outlet
            # Determine the minimum T5 temperature to avoid at expansion valve outlet
            MEDIUM1.update(CP.PQ_INPUTS, p_evap, 0.01)
            h5_min = MEDIUM1.hmass()
            MEDIUM1.update(CP.HmassP_INPUTS, h5_min, p_cond)
            T5_min = MEDIUM1.T()

            T[5] = max(T_sink_in + delta_T_pinch_sink_min, T5_min)
            MEDIUM1.update(CP.PT_INPUTS, p_cond, T[5])
            h[5] = MEDIUM1.hmass()
            Q[5] = MEDIUM1.Q()
            if p_cond < p_crit and Q[5] > 0:
                print('Warning: Two-phase state at condenser outlet')

            # Evaporator Outlet
            MEDIUM1.update(CP.PQ_INPUTS, p_evap, 1)
            h[7] = MEDIUM1.hmass()
            T[7] = MEDIUM1.T()

            # Compressor
            if SubTrans == 'Sub':  # p_cond < p_crit * limit_sc
                # Subcritical Operation
                # Disctinction in Desuperheating, Condensing and Subcooling for Subcritical Process
                MEDIUM1.update(CP.PQ_INPUTS, p_cond, 1)
                h[3] = MEDIUM1.hmass()
                T[3] = MEDIUM1.T()
                MEDIUM1.update(CP.PQ_INPUTS, p_cond, 0)
                h[4] = MEDIUM1.hmass()
                T[4] = MEDIUM1.T()

                # Iteration to ensure dry compression at outlet
                error_delta_T_SH = 10
                delta_T_SH = delta_T_minSH
                tic = timer()
                while error_delta_T_SH >= error_delta_T:  # 1e-2 Allowable error on pinch temperature
                    # Compressor inlet
                    T[1] = T[7] + delta_T_SH
                    MEDIUM1.update(CP.PT_INPUTS, p_evap, T[1])
                    h[1] = MEDIUM1.hmass()
                    s[1] = MEDIUM1.smass()

                    # Compressor outlet
                    MEDIUM1.update(CP.PSmass_INPUTS, p_cond, s[1])
                    h_2_is = MEDIUM1.hmass()
                    h[2] = (h_2_is - h[1]) / eta_is_comp + h[1]
                    MEDIUM1.update(CP.HmassP_INPUTS, h[2], p_cond)
                    T[2] = MEDIUM1.T()
                    delta_T_SH_summary = [T[1] - T[7], T[2] - T[3]]
                    error_delta_T_SH = delta_T_minSH - min(delta_T_SH_summary)
                    delta_T_SH = delta_T_SH + error_delta_T_SH
                    toc = timer()
                    if toc - tic > 60:
                        print('SH did not converge within 60 s')
            else:
                # Transcritical
                T[1] =T[7] + delta_T_minSH
                MEDIUM1.update(CP.PT_INPUTS, p_evap, T[1])
                h[1] = MEDIUM1.hmass()
                s[1] = MEDIUM1.smass()

                MEDIUM1.update(CP.PSmass_INPUTS, p_cond, s[1])
                h_2_is = MEDIUM1.hmass()
                h[2] = (h_2_is - h[1]) / eta_is_comp + h[1]
                MEDIUM1.update(CP.HmassP_INPUTS, h[2], p_cond)
                T[2] = MEDIUM1.T()
                # Equidistant Discretization for Transcritical Process  T
                h[3] = h[2]
                T[3] = T[2]
                h[4] = h[5]
                T[4] = T[5]

            # Expansion valve
            h[6] = h[5]
            MEDIUM1.update(CP.HmassP_INPUTS, h[6], p_evap)
            T[6] = MEDIUM1.T()

            # Mass flow rate
            m_flow = Qdot_source / (h[1] - h[6])
            Q_dot_evap = m_flow * (h[7] - h[6])
            Q_dot_SH = m_flow * (h[1] - h[7])

            h_22 = h_23 + Q_dot_evap / m_flow_source
            SOURCE.update(CP.HmassP_INPUTS, h_22, p_source)
            T_22 = SOURCE.T()

            h_evap = np.append(np.linspace(h[6], h[7], n_states_evap),h[1])
            h_source = np.append(np.linspace(h_23, h_22, n_states_evap),h_21)

            # Temperature Profile -> from 6 to 1
            T_evap[0] = T[6]
            T_source[0] = T_23
            T_evap[-1] = T[1]  # -1 is the index of the last object in a array
            T_source[-1] = T_21
            T_evap[-2] = T[7]
            T_source[-2] = T_22

            for i in range(1,len(h_evap)-2):
                try:
                    MEDIUM1.update(CP.HmassP_INPUTS, h_evap[i], p_evap)
                    T_evap[i] = MEDIUM1.T()
                    SOURCE.update(CP.HmassP_INPUTS, h_source[i], p_source)
                    T_source[i] = SOURCE.T()
                except:
                    T_evap[i] = math.nan
                    T_source[i] = math.nan

            delta_T_source = T_source - T_evap
            # print('Temperature different on the source side:',delta_T_source) ####################
            Q_dot_source_cum = m_flow * (h_evap - h[6])

            # Compressor work
            W_dot_comp = m_flow * (h[2] - h[1])

            # Heat sink side
            Q_dot_sink = m_flow * (h[2] - h[5])
            m_flow_sink = Q_dot_sink / (h_34 - h_31)

            # Heat Flow Rates
            Q_dot_desup = m_flow * (h[2] - h[3])
            Q_dot_cond = max(0, m_flow * (h[3] - h[4]))
            Q_dot_subcool = max(0, m_flow * (h[4] - h[5]))

            if SubTrans == 'Sub':  # % p_cond < p_crit * limit_sc  % < p_max_H
                # Disctinction in Desuperheating, Condensing and Subcooling for Subcritical Process
                h_32 = h_31 + Q_dot_subcool / m_flow_sink
                h_33 = h_31 + (Q_dot_subcool + Q_dot_cond) / m_flow_sink
                SINK.update(CP.HmassP_INPUTS, h_32, p_sink)
                T_32 = SINK.T()
                SINK.update(CP.HmassP_INPUTS, h_33, p_sink)
                T_33 = SINK.T()
            else:
                # Equidistant Discretization for Transcritical Process
                h_32 = h_31
                h_33 = h_34

                T_32 = T_31
                T_33 = T_34

            # Calculation of discretized states on condenser side

            if SubTrans == 'Sub':  # p_cond < p_crit
                h_cond_sub = np.linspace(h[5], h[4], n_states_sub)
                h_sink_sub = np.linspace(h_31, h_32, n_states_sub)

                h_cond_cond = np.linspace(h[4], h[3], n_states_cond)
                h_sink_cond = np.linspace(h_32, h_33, n_states_cond)

                h_cond_desup = np.linspace(h[3], h[2], n_states_desup)
                h_sink_desup = np.linspace(h_33, h_34, n_states_desup)

                h_cond = np.hstack((h_cond_sub[0:n_states_sub - 1], h_cond_cond[0:n_states_cond - 1], h_cond_desup))
                h_sink = np.hstack((h_sink_sub[0:n_states_sub - 1], h_sink_cond[0:n_states_cond - 1], h_sink_desup))

                T_cond[0] = T[5]
                T_sink[0] = T_31
                T_cond[-1] = T[2]  # -1 is the index of the last object in a array
                T_sink[-1] = T_34

                for i in range(1, len(h_cond) - 1):
                    try:
                        MEDIUM1.update(CP.HmassP_INPUTS, h_cond[i], p_cond)
                        T_cond[i] = MEDIUM1.T()
                        SINK.update(CP.HmassP_INPUTS, h_sink[i], p_sink)
                        T_sink[i] = SINK.T()
                    except:
                        T_cond[i] = math.nan
                        T_sink[i] = math.nan

                for item in T_cond:
                    if (math.isnan(item)) == True:
                        print('Points did  not converge')

            else:
                h_cond = np.linspace(h[5], h[2], n_states_sub + n_states_cond + n_states_desup - 2)
                h_sink = np.linspace(h_31, h_34, n_states_sub + n_states_cond + n_states_desup - 2)

                T_cond[0] = T[5]
                T_sink[0] = T_31

                T_cond[-1] = T[2]
                T_sink[-1] = T_34

                for i in range(1, len(h_cond) - 1):
                    try:
                        MEDIUM1.update(CP.HmassP_INPUTS, h_cond[i], p_cond)
                        T_cond[i] = MEDIUM1.T()
                        SINK.update(CP.HmassP_INPUTS, h_sink[i], p_sink)
                        T_sink[i] = SINK.T()
                    except:
                        T_cond[i] = math.nan
                        T_sink[i] = math.nan

            # Temperature Differences
            delta_T_sink = T_cond - T_sink
            # print('Temperature different on the sink side',delta_T_sink) ######
            Q_dot_sink_cum = m_flow * (h_cond - h[5])

            T_sink_out = T_34
            T_source_out = T_23

            Res_source[m] = delta_T_pinch_source_min - min(delta_T_source)
            # print('Temperature redidual on the source side: Res_source',m,':',Res_source[m])
            Res_sink[m] = delta_T_pinch_sink_min - min(delta_T_sink[1:])  # to avoid devide by zero when update condensing pressure
            # print('Temperature redidual on the sink side: Res_sink',m,':',Res_sink[m])

            if m == 0:
                RelRes[0] = (delta_T_pinch_source_min - min(delta_T_source)) / delta_T_pinch_source_min
                RelRes[1] = (delta_T_pinch_sink_min - min(delta_T_sink[1:])) / delta_T_pinch_sink_min

            if max(abs(RelRes)) < tol:
                break

            if iter > iter_max:
                print('Model did not converge with ', iter_max, ' iterations')

        if max(abs(RelRes)) >= tol:  # update pressure
            x[0] = x[0] - x[0] * tol / (Res_source[0] - Res_source[1]) * Res_source[0]
            x[1] = x[1] - x[1] * tol / (Res_sink[0] - Res_sink[2]) * Res_sink[0]
            # Lower bound pressure
            x[0] = max(x[0], p_min_L)
            x[1] = max(x[1], p_min_H)
            # Upper bound pressure
            x[0] = min(x[0], p_max_L)
            x[1] = min(x[1], p_max_H)

    # _________________________________________________________________________________
    #                    Collect the output
    # _________________________________________________________________________________
    # Update the state points
    p[1] = p_evap
    p[2] = p_cond
    p[3] = p_cond
    p[4] = p_cond
    p[5] = p_cond
    p[6] = p_evap
    p[7] = p_evap
    p[0] = p[7]
    h[0] = h[7]
    T[0] = T[7]

    for i in range(8):
        MEDIUM1.update(CP.HmassT_INPUTS, h[i], T[i])
        s[i] = MEDIUM1.smass() / 1000  # kJ/kgK



    class Output:
        pass

    output = Output()

    output.Qdot_cond = Q_dot_sink
    output.mdot_sink = m_flow_sink
    output.mdot_ref_tot = m_flow
    output.mdot_source = m_flow_source
    output.p_evap = p_evap
    output.p_cond = p_cond

    output.T_ref_in_cond = T[2]  # discharge temperature
    output.T_ref_out_cond = T[5]  # subcooler outlet
    output.T_ref_in_eva = T[6]  # evaporator inlet
    output.T_ref_out_eva = T[1]  # suction temperature

    output.delta_T_SC = T[4] - T[5]  # subcooling
    output.delta_T_SH = T[1] - T[7]  # superheating

    output.W_dot = W_dot_comp / eta_motor
    output.COP = Q_dot_sink / output.W_dot

    MEDIUM1.update(CP.HmassP_INPUTS, h[1], p_evap)
    output.V_suction = m_flow / MEDIUM1.rhomass()  # [ m3/s]
    output.VHC = Q_dot_sink / output.V_suction  # [J/m3]
    output.Q_dot_source_cum = Q_dot_source_cum
    output.Q_dot_sink_cum = Q_dot_sink_cum
    output.T_source = T_source
    output.T_sink = T_sink
    output.T_evap = T_evap
    output.T_cond = T_cond
    output.h_cond = h_cond
    output.h = h
    output.s = s
    output.T = T
    output.p = p
    output.SubTrans = SubTrans
    #output.delta_T_pinch_ihx_min = delta_T_pinch_ihx_min
    output.cycle = 'Single state HP Standard'
    output.Qdot_eva = source.Qdot_design
    output.refrigerant = MEDIUM1

    return output






































