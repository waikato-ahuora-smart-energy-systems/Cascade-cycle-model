
import matplotlib.pyplot as plt
import numpy as np
import CoolProp as CP



def TQ_diagram_source(output):

    # Temperature - Heat Diagram for the Source side
    Q_source = output.Q_dot_source_cum * 1e-3
    T_source = output.T_source -273.15
    T_evap = output.T_evap - 273.15
    plt.plot(Q_source, T_source, label = 'T_source', color = 'b')
    plt.plot(Q_source, T_evap, label ='T_evap',color = 'orange')
    plt.xlabel("Transferred Heat, kW")
    plt.ylabel('Temperature, $^\circ$C')
    plt.legend()
    plt.show()

def TQ_diagram_sink(output):

    # Temperature - Heat Diagram for the Sink side
    Q_sink = output.Q_dot_sink_cum * 1e-3
    T_sink = output.T_sink -273.15
    T_cond = output.T_cond - 273.15
    plt.plot(Q_sink, T_sink, label = 'T_sink',color = 'r')
    plt.plot(Q_sink, T_cond, label ='T_cond',color = 'orange')
    plt.xlabel("Transferred Heat, kW")
    plt.ylabel('Temperature, $^\circ$C')
    plt.legend()
    plt.show()

def TQ_diagram_ihx(output):
    # Temperature - Heat Diagram for the IHX
    Q_ihx= output.Q_dot_ihx_cum * 1e-3
    T_ihx_cold = output.T_ihx_cold -273.15
    T_ihx_hot = output.T_ihx_hot- 273.15
    plt.plot(Q_ihx, T_ihx_cold, label = 'T_ihx_cold',color = 'blue')
    plt.plot(Q_ihx, T_ihx_hot, label ='T_ihx_hot',color = 'red')
    plt.xlabel("Transferred Heat, kW")
    plt.ylabel('Temperature, $^\circ$C')
    plt.legend()
    plt.show()

def TQ_diagram_cascade_ihx(output):
    # Temperature - Heat Diagram for the IHX
    Q_cax_bot = output.Q_dot_cax_bot_cum * 1e-3

    T_cax_cold = output.T_cax_cold -273.15
    T_cax_hot = output.T_cax_hot- 273.15
    plt.plot(Q_cax_bot, T_cax_hot, label = 'T_cascade_ihx_bottom_cycle',color = 'red')
    plt.plot(Q_cax_bot, T_cax_cold, label ='T_cascade_ihx_top_cycle',color = 'blue')
    plt.xlabel("Transferred Heat, kW")
    plt.ylabel('Temperature, $^\circ$C')
    plt.legend()
    plt.show()

def logPh_diagram(output):

    MEDIUM1 = output.refrigerant
    p_evap = output.p_evap           # Pa
    p_cond = output.p_cond           # bar
    h = output.h                     # J/kg
    T = output.T                     #K
    p = output.p

    #  Calculate critical point
    MEDIUM1.update(CP.PT_INPUTS,101325,300)
    T_crit = MEDIUM1.T_critical()       # K
    p_crit = MEDIUM1.p_critical()       # Pa
    # Number of calculated points on the saturation line
    n_sat = 50
    # Preallocating space in advance
    T_sat = np.zeros(n_sat)
    p_sat = np.zeros(n_sat)
    h_liq = np.zeros(n_sat)
    h_vap = np.zeros(n_sat)
    T_sat = np.linspace(273.15,T_crit,n_sat)
    for i in range(n_sat):
        MEDIUM1.update(CP.QT_INPUTS,0,T_sat[i])
        h_liq[i] = MEDIUM1.hmass()/1000 # kJ/kg
        p_sat[i] = MEDIUM1.p()/1e5      #bar
        MEDIUM1.update(CP.QT_INPUTS, 1, T_sat[i])
        h_vap[i] = MEDIUM1.hmass()/1000 # kJ/kg

    plt.plot(h_liq, p_sat,color ='black')
    plt.plot(h_vap, p_sat,color ='black')

    # Define pressures of cycle state points

    p = np.append(p,p[0])
    h = np.append(h,h[0])

    plt.plot(h/1000, p*1e-5)
    plt.xlabel("Enthalpy, kJ/kg")
    plt.ylabel('Pressure, bar')
    # convert y-axis to Logarithmic scale
    plt.yscale("log")
    plt.show()

def TS_diagram(output):

    MEDIUM1 = output.refrigerant
    p_evap = output.p_evap      # Pa
    p_cond = output.p_cond      # Pa
    SubTrans = output.SubTrans
    h = output.h                # J/kg
    T = output.T                #K
    s = output.s                # kJ/kgK

    #  Calculate critical point
    MEDIUM1.update(CP.PT_INPUTS,101325,300)
    T_crit = MEDIUM1.T_critical()       # K
    p_crit = MEDIUM1.p_critical()       # Pa
    # Number of calculated points on the saturation line
    n_sat = 50
    # Pre-allocating space in advance
    T_sat = np.zeros(n_sat)
    s_liq = np.zeros(n_sat)
    s_vap = np.zeros(n_sat)
    T_sat = np.linspace(273.15,T_crit,n_sat)
    for i in range(n_sat):
        MEDIUM1.update(CP.QT_INPUTS,0,T_sat[i])
        s_liq[i] = MEDIUM1.smass()/1000 # kJ/kgK
        MEDIUM1.update(CP.QT_INPUTS, 1, T_sat[i])
        s_vap[i] = MEDIUM1.smass()/1000 # kJ/kgK

    plt.plot(s_liq, T_sat-273.15,color ='black')
    plt.plot(s_vap, T_sat-273.15,color ='black')

    if SubTrans == 'Sub':
        s = np.append(s,s[0])
        T = np.append(T,T[0])
        plt.plot(s, T-273.15)
        plt.xlabel("Entropy, kJ/kgK")
        plt.ylabel('Temperature, $^\circ$C')
        plt.show()
    else:
        T_cond = np.flip(output.T_cond)  # K
        h_cond = np.flip(output.h_cond)
        s_cond = np.zeros(len(h_cond))
        for i in range(len(h_cond)):
            MEDIUM1.update(CP.HmassT_INPUTS, h_cond[i], T_cond[i])
            s_cond[i] = MEDIUM1.smass() / 1000  # kJ/kgK

        # IHX
        T_plot = np.hstack((T[0:2], T_cond, T[6:],T[0]))
        s_plot = np.hstack((s[0:2], s_cond, s[6:],s[0]))

        plt.plot(s_plot, T_plot-273.15)
        plt.xlabel("Entropy, kJ/kgK")
        plt.ylabel('Temperature, $^\circ$C')
        plt.show()













