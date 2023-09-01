from Cycles.Single_State_Ihx import single_state_ihx_hp
from Cycles.Single_State_Std import single_state_std_hp
from Cycles.Cascade import cascade_cycle
from Opt_All_Mixture import opt_all_mixture
import Figures.Function_figures as Fig
from xlwt import Workbook
from numpy import asarray
from timeit import default_timer as timer
from numpy import around
from DifferentialEvolution import differential_evolution
import CoolProp as CP
from Economics.Cascade_Economic_Analysis import economic_analysis
import CoolProp.CoolProp
#import os
# Change the RPPREFIX environment variable to where REFPROP is not actually installed
#os.environ['RPPREFIX'] = r'C:\Program Files \REFPROP'

#import CoolProp.CoolProp
#CoolProp.CoolProp .set_config_string(CoolProp.CoolProp .ALTERNATIVE_REFPROP_PATH, r'C:\Program Files (x86)\REFPROP')

class Result:
    pass

class Input:
    pass


HP_DESIGN = Result()
economics = Result()

medium1 = 'Hexane'
medium2 = 'Iso-Butane'
medium3 = 'CarbonDioxide'
medium4 = 'Butane'
medium5 = 'Acetone'
medium6 = 'Neopentane'

MEDIUMPure1 = CP.AbstractState('REFPROP', medium1)
MEDIUMPure2 = CP.AbstractState('REFPROP', medium2)
MEDIUMPure3 = CP.AbstractState('REFPROP', medium3)
MEDIUMPure4 = CP.AbstractState('REFPROP', medium4)
MEDIUMPure5 = CP.AbstractState('REFPROP', medium5)
MEDIUMPure6 = CP.AbstractState('REFPROP', medium6)

refrigerant = Input()
#refrigerant.top          = MEDIUMPure1
refrigerant.top          = MEDIUMPure5  #test
#refrigerant.bot          = MEDIUMPure2
refrigerant.bot          = MEDIUMPure6  #test
refrigerant.ref_3         = MEDIUMPure3 # for single state cycle
refrigerant.ref_4         = MEDIUMPure4

#=======================================================================================================================
#CASCADE CYCLE WITH CONSTANT OPERATING CONDITION
#=======================================================================================================================

# parameters for cascade cycle:
# parameters[0]  delta_T_SH_IHX1 =0
# parameters[1] delta_T_SH_IHX2 = 47 K
# parameters[2] = T_intermediate = 117 deg C
parameters = asarray([0, 0.00426473785591952, 107.364623746215,0.0755442327252178])

def result_cascade(refrigerant,parameters):
    [COP,HP_DESIGN]     = cascade_cycle(refrigerant, parameters)
    economics        = economic_analysis(HP_DESIGN)
    print('Bottom cycle evaporation pressure (bar):             ', HP_DESIGN.p_evap_bot * 1e-5)
    print('Bottom cycle medium pressure (bar):                  ', HP_DESIGN.p_mid_bot * 1e-5)
    print('Bottom cycle condensation pressure (bar):            ', HP_DESIGN.p_cond_bot * 1e-5)
    print('Top cycle evaporation pressure (bar):                ', HP_DESIGN.p_evap_top * 1e-5)
    print('Top cycle condensation pressure (bar):               ', HP_DESIGN.p_cond_top * 1e-5)
    print('COP(design):                                         ', COP)
    print('Total heat sink (kW):                                ', HP_DESIGN.Qdot_sink * 1e-3)
    print('Total heat source (kW):                              ', HP_DESIGN.Qdot_source * 1e-3)
    #print('Bottom cycle suction volume flow rate1 (m3/hr):       ', HP_DESIGN.V_suction_bot * 3600)
    print('Bottom cycle suction volume flow rate (m3/hr):       ', HP_DESIGN.V_suction_com1 * 3600)
    #print('Bottom cycle Volumetric heating capacity (kJ/m3):    ', HP_DESIGN.VHC_bot /1000)
    #print('Top cycle suction volume flow rate1 (m3/hr):          ', HP_DESIGN.V_suction_top * 3600)
    print('Top cycle suction volume flow rate (m3/hr):          ', HP_DESIGN.V_suction_com3 * 3600)
    #print('Top cycle Volumetric heating capacity (kJ/m3):       ', HP_DESIGN.VHC_top /1000)

    print('Compressor 1 volumetric effeciency :          ', HP_DESIGN.eta_vol_com1)
    print('Compressor 2 volumetric effeciency :          ', HP_DESIGN.eta_vol_com2)
    print('Compressor 3 volumetric effeciency :          ', HP_DESIGN.eta_vol_com3)

    h = HP_DESIGN.h
    T = HP_DESIGN.T
    for i in range(len(h)):
        print('Elthalpy at point     ', i, '(kJ/kg)     :', h[i] / 1000)
        print('Temperature at point  ', i, '(oC)        :', T[i] - 273.15)

    print('Delta T at point 15:', T[15] - HP_DESIGN.T_35)
    print('Delta T at point 16:', T[16] - HP_DESIGN.T_34)
    print('Delta T at point 17:', T[17] - HP_DESIGN.T_33)
    print('Delta T at point 7:', T[7] - HP_DESIGN.T_32)
    print('Delta T at point 8:', T[8] - HP_DESIGN.T_31)
    #print('Delta T at point 4:')
    # --------------------------------------------------------------------------------------
    #            Plot Figure
    # --------------------------------------------------------------------------------------
    Fig.TQ_diagram_sink(HP_DESIGN)
    Fig.TQ_diagram_source(HP_DESIGN)
    Fig.TQ_diagram_cascade_ihx(HP_DESIGN)

    # Economic analysis
    print('NPV value: ',economics.NPV)
    print('PBT: ', economics.PBT)
    print('Specific investment cost per unit of supplied heat, NZD/kW',economics.TCI_spec)
    print('Specific cost of heat, NZD / MWh',economics.c_H)
    # Export excel
    wb = Workbook()
    sheet1 = wb.add_sheet('Cascade analysis results')
    sheet1.write(0, 0, 'Top cycle Ref')
    sheet1.write(0, 1, 'Bottom cycle Ref')
    sheet1.write(0, 2, 'deltaT_SH1,oC')
    sheet1.write(0, 3, 'deltaT_SH2, oC')
    sheet1.write(0, 4, 'T_intermediate, oC')
    sheet1.write(0, 5, 'T_cond_top, oC')
    sheet1.write(0, 6, 'COP')
    sheet1.write(0, 7, 'Exergy_eff')
    sheet1.write(0, 8, 'n_Lorren')
    sheet1.write(0, 9, 'E_des_top_cond (kJ/kg)')
    sheet1.write(0, 10, 'E_des_IHX3 (kJ/kg)')
    sheet1.write(0, 11, 'E_des_CHX (kJ/kg)')
    sheet1.write(0, 12, 'E_des_bot_evap (kJ/kg)')
    sheet1.write(0, 13, 'm_bot (Kg/s)')
    sheet1.write(0, 14, 'm_top (Kg/s)')
    sheet1.write(0, 15, 'V_suction_com1 (m3/h)')
    sheet1.write(0, 16, 'V_suction_com2 (m3/h)')
    sheet1.write(0, 17, 'V_suction_com3 (m3/h)')
    sheet1.write(0, 18, 'T_suction_com1 (oC)')
    sheet1.write(0, 19, 'T_suction_com2 (oC)')
    sheet1.write(0, 20, 'T_suction_com3 (oC)')
    sheet1.write(0, 21, 'T_discharge_com2 (oC)')
    sheet1.write(0, 22, 'T_discharge_com3 (oC)')
    sheet1.write(0, 23, 'specific_volume_2 (m3/kg)')
    sheet1.write(0, 24, 'specific_volume_3 (m3/kg)')
    sheet1.write(0, 25, 'specific_volume_4 (m3/kg)')
    sheet1.write(0, 26, 'specific_volume_14 (m3/kg)')
    sheet1.write(0, 27, 'specific_volume_15 (m3/kg)')
    sheet1.write(0, 28, 'eta_vol_com1')
    sheet1.write(0, 29, 'eta_vol_com2')
    sheet1.write(0, 30, 'eta_vol_com3')
    sheet1.write(0, 31, 'p_evap_bot (bar)')
    sheet1.write(0, 32, 'p_mid_bot (bar)')
    sheet1.write(0, 33, 'p_cond_bot (bar)')
    sheet1.write(0, 34, 'p_evap_top (bar)')
    sheet1.write(0, 35, 'p_cond_top (bar)')
    sheet1.write(0, 36, 'Q_dot_sink (kW)')
    sheet1.write(0, 37, 'W_dot (kW)')


    sheet1.write(0, 38, 'PEC_com1 (1000 NZD)')
    sheet1.write(0, 39, 'PEC_com2 (1000 NZD)')
    sheet1.write(0, 40, 'PEC_com3 (1000 NZD)')
    sheet1.write(0, 41, 'PEC_Evap (1000 NZD)')
    sheet1.write(0, 42, 'PEC_Top_Cond (1000 NZD)')
    sheet1.write(0, 43, 'PEC_Cax (1000 NZD)')
    sheet1.write(0, 44, 'PEC_IHX2 (1000 NZD)')
    sheet1.write(0, 45, 'PEC_IHX3 (1000 NZD)')
    sheet1.write(0, 46, 'TCI (1000 NZD)')
    sheet1.write(0, 47, 'TCI_spec (NZD/kW)')
    sheet1.write(0, 48, 'NPV (1000 NZD)')
    sheet1.write(0, 49, 'PBT (years)')
    sheet1.write(0, 50, 'Energy_Produced_Yearly (MWh)')
    sheet1.write(0, 51, 'c_H (NZD / MWh)')
    sheet1.write(0, 52, 'CF_FCel (1000 NZD/year)')
    sheet1.write(0, 53, 'CF_FChsource (1000 NZD/year)')
    sheet1.write(0, 54, 'CF_rev (1000 NZD/year)')
    sheet1.write(0, 55, 'CF_TCI (1000 NZD/year)')
    sheet1.write(0, 56, 'CF_OM (1000 NZD/year)')
    sheet1.write(0, 57, 'CF_total (1000 NZD/year)')
    sheet1.write(0, 58, 'eta_is_com1')
    sheet1.write(0, 59, 'eta_is_com2')
    sheet1.write(0, 60, 'eta_is_com3')



    sheet1.write(1, 0, medium5)
    sheet1.write(1, 1, medium6)
    sheet1.write(1, 2, parameters[0])
    sheet1.write(1, 3, parameters[1])
    sheet1.write(1, 4, parameters[2])

    sheet1.write(1, 5, HP_DESIGN.T_cond_top - 273.15)  # oC
    sheet1.write(1, 6, COP)
    sheet1.write(1, 7, HP_DESIGN.exergy_eff)
    sheet1.write(1, 8, HP_DESIGN.n_Lorren)
    sheet1.write(1, 9, HP_DESIGN.E_des_top_cond)  # (kJ/kg)
    sheet1.write(1, 10, HP_DESIGN.E_des_IHX3)  # (kJ/kg)
    sheet1.write(1, 11, HP_DESIGN.E_des_CHX)  # (kJ/kg)
    sheet1.write(1, 12, HP_DESIGN.E_des_bot_evap)  # (kJ/kg)
    sheet1.write(1, 13, HP_DESIGN.mdot_ref_bot)  # (Kg/s)
    sheet1.write(1, 14, HP_DESIGN.mdot_ref_top)  # (Kg/s)
    sheet1.write(1, 15, HP_DESIGN.V_suction_com1 * 3600)  # (m3/h)
    sheet1.write(1, 16, HP_DESIGN.V_suction_com2 * 3600)  # (m3/h)
    sheet1.write(1, 17, HP_DESIGN.V_suction_com3 * 3600)  # (m3/h)
    sheet1.write(1, 18, HP_DESIGN.T_suction_com1 - 273.15)  # (oC)
    sheet1.write(1, 19, HP_DESIGN.T_suction_com2 - 273.15)  # (oC)
    sheet1.write(1, 20, HP_DESIGN.T_suction_com3 - 273.15)  # (oC)
    sheet1.write(1, 21, HP_DESIGN.T_discharge_com2 - 273.15)  # (oC)
    sheet1.write(1, 22, HP_DESIGN.T_discharge_com3 - 273.15)  # (oC)
    sheet1.write(1, 23, HP_DESIGN.specific_volume_2)  # (m3/kg)
    sheet1.write(1, 24, HP_DESIGN.specific_volume_3)  # (m3/kg)
    sheet1.write(1, 25, HP_DESIGN.specific_volume_4)  # (m3/kg)
    sheet1.write(1, 26, HP_DESIGN.specific_volume_14)  # (m3/kg)
    sheet1.write(1, 27, HP_DESIGN.specific_volume_15)  # (m3/kg)
    sheet1.write(1, 28, HP_DESIGN.eta_vol_com1)
    sheet1.write(1, 29, HP_DESIGN.eta_vol_com2)
    sheet1.write(1, 30, HP_DESIGN.eta_vol_com3)
    sheet1.write(1, 31, HP_DESIGN.p_evap_bot / 1e5)  # (bar)
    sheet1.write(1, 32, HP_DESIGN.p_mid_bot / 1e5)  # (bar)'
    sheet1.write(1, 33, HP_DESIGN.p_cond_bot / 1e5)  # (bar)
    sheet1.write(1, 34, HP_DESIGN.p_evap_top / 1e5)  # (bar)
    sheet1.write(1, 35, HP_DESIGN.p_cond_top / 1e5)  # (bar)
    sheet1.write(1, 36, HP_DESIGN.Qdot_sink/ 1000)  # (kW)
    sheet1.write(1, 37, HP_DESIGN.W_dot / 1000)  # (kW)

    sheet1.write(1, 38, economics.PEC_com1 / 1000)  # (1000 NZD)
    sheet1.write(1, 39, economics.PEC_com2 / 1000)  # (1000 NZD)
    sheet1.write(1, 40, economics.PEC_com3 / 1000)  # (1000 NZD)
    sheet1.write(1, 41, economics.PEC_Evap / 1000)  # (1000 NZD)
    sheet1.write(1, 42, economics.PEC_Top_Cond / 1000)  # (1000 NZD)
    sheet1.write(1, 43, economics.PEC_Cax / 1000)  # (1000 NZD)
    sheet1.write(1, 44, economics.PEC_IHX2 / 1000)  # (1000 NZD)
    sheet1.write(1, 45, economics.PEC_IHX3 / 1000)  # (1000 NZD)
    sheet1.write(1, 46, economics.TCI / 1000)  # (1000 NZD)
    sheet1.write(1, 47, economics.TCI_spec)  # (NZD/kW)
    sheet1.write(1, 48, economics.NPV / 1000)  # (1000 NZD)
    sheet1.write(1, 49, economics.PBT)  # (years)
    sheet1.write(1, 50, economics.Energy_Produced_Yearly)  # (MWh)
    sheet1.write(1, 51, economics.c_H)  # (NZD / MWh)
    sheet1.write(1, 52, economics.CF_FCel / 1000)  # (1000 NZD/year)
    sheet1.write(1, 53, economics.CF_FChsource / 1000)  # (1000 NZD/year)
    sheet1.write(1, 54, economics.CF_rev / 1000)  # (1000 NZD/year)')
    sheet1.write(1, 55, economics.CF_TCI / 1000)  # (1000 NZD/year)')
    sheet1.write(1, 56, economics.CF_OM / 1000)  # (1000 NZD/year)
    sheet1.write(1, 57, economics.CF_total / 1000)  # (1000 NZD/year)
    sheet1.write(1, 58, HP_DESIGN.eta_is_comp1)
    sheet1.write(1, 59, HP_DESIGN.eta_is_comp2)
    sheet1.write(1, 60, HP_DESIGN.eta_is_comp3)

    wb.save('cascade_analysis.xls')

#=======================================================================================================================
# SINGLE STATE STD CYCLE
#=======================================================================================================================


def result_single_state_std(refrigerant):
    HP_DESIGN  = single_state_std_hp(refrigerant)
    print('Evaporation pressure (bar):                  ', HP_DESIGN.p_evap * 1e-5)
    print('Condensation pressure (bar):                 ', HP_DESIGN.p_cond * 1e-5)
    print('COP(design):                                 ', HP_DESIGN.COP)
    print('Total heat flow rate at evaporator (kW):     ', HP_DESIGN.Qdot_eva*1e-3)
    print('Total heat flow rate at condenser (kW):      ', HP_DESIGN.Qdot_cond*1e-3)
    print('Pressure ratio:                              ', HP_DESIGN.p_cond/HP_DESIGN.p_evap)
    print('Evaporator superheating (K):                 ', HP_DESIGN.delta_T_SH)
    print('Subcooling (K):                              ', HP_DESIGN.delta_T_SC)
    print('Suction volume flow rate (m3/hr):            ', HP_DESIGN.V_suction * 3600)
    print('Volumetric heating capacity (kJ/m3):         ', HP_DESIGN.VHC/1000)
    h = HP_DESIGN.h
    T = HP_DESIGN.T
    for i in range(len(h)):
        print('Elthalpy at point     ', i, '(kJ/kg)     :', h[i] / 1000)
        print('Temperature at point  ', i, '(oC)        :', T[i] - 273.15)
    # -------------------------------------------------------------------------------------
    #            Plot Figure
    # --------------------------------------------------------------------------------------
    Fig.TQ_diagram_sink(HP_DESIGN)
    Fig.TQ_diagram_source(HP_DESIGN)
    Fig.logPh_diagram(HP_DESIGN)
    Fig.TS_diagram(HP_DESIGN)

#=======================================================================================================================
# SINGLE STATE IHX CYCLE
#=======================================================================================================================

def result_single_state_ihx(refrigerant):
    HP_DESIGN  = single_state_ihx_hp(refrigerant)
    print('Evaporation pressure (bar):                  ', HP_DESIGN.p_evap * 1e-5)
    print('Condensation pressure (bar):                 ', HP_DESIGN.p_cond * 1e-5)
    print('COP(design):                                 ', HP_DESIGN.COP)
    print('Total heat flow rate at evaporator (kW):     ', HP_DESIGN.Qdot_eva*1e-3)
    print('Total heat flow rate at condenser (kW):      ', HP_DESIGN.Qdot_cond*1e-3)
    print('Pressure ratio:                              ', HP_DESIGN.p_cond/HP_DESIGN.p_evap)
    print('Evaporator superheating (K):                 ', HP_DESIGN.delta_T_SH)
    print('Subcooling (K):                              ', HP_DESIGN.delta_T_SC)
    print('Suction volume flow rate (m3/hr):            ', HP_DESIGN.V_suction * 3600)
    print('Volumetric heating capacity (kJ/m3):         ', HP_DESIGN.VHC/1000)
    print('Minimum temperature different at the IHX, oC:', HP_DESIGN.delta_T_pinch_ihx_min)
    h = HP_DESIGN.h
    T = HP_DESIGN.T
    for i in range(len(h)):
        print('Elthalpy at point     ', i, '(kJ/kg)     :', h[i] / 1000)
        print('Temperature at point  ', i, '(oC)        :', T[i] - 273.15)
    # --------------------------------------------------------------------------------------
    #            Plot Figure
    # --------------------------------------------------------------------------------------
    Fig.TQ_diagram_sink(HP_DESIGN)
    Fig.TQ_diagram_source(HP_DESIGN)
    Fig.TQ_diagram_ihx(HP_DESIGN)
    Fig.logPh_diagram(HP_DESIGN)
    Fig.TS_diagram(HP_DESIGN)


#=======================================================================================================================
# RUN THE CASCADE CYCLE AT DIFFERENT OPERATING CONDITIONS
#=======================================================================================================================

def run_cascade():

    wb = Workbook()
    # Super heat temperature of the IHX1
    sheet1 = wb.add_sheet('deltaT_SH1')
    sheet1.write(0, 0, 'T_SH1, oC')
    sheet1.write(0, 1, 'COP_Cascade')
    k = 0
    parameters[2] = 109.521776387603 # T_intermediate
    parameters[1] = 0 # delta_T_SH_IHX2
    parameters[3] = 0 # vapour quality of point 7
    for i in range(0, 21, 1):
        k = k + 1
        parameters[0] = i  # delta_T_SH_IHX1
        [COP,_] = cascade_cycle(refrigerant,parameters)
        print('COP at T_SH1 = ', i, 'oC:', COP)
        sheet1.write(k, 0, i)
        sheet1.write(k, 1, COP)

    # Superheated temperature of the IHX 2
    sheet2 = wb.add_sheet('deltaT_SH2')
    sheet2.write(0, 0, 'T_SH2, oC')
    sheet2.write(0, 1, 'COP_Cascade')
    k = 0
    parameters[2]  = 109.521776387603 # T_intermediate
    parameters[0] = 0 # delta_T_SH_IHX1 = 0  # K for cascade cycle
    parameters[3] = 0  # vapour quality of point 7
    for i in range(0,51,1):
        k = k+1
        parameters[1] = i   # delta_T_SH_IHX2
        [COP,_] = cascade_cycle(refrigerant,parameters)
        print('COP at T_SH2 = ',i,'oC:', COP)
        sheet2.write(k, 0, i)
        sheet2.write(k, 1, COP)

    # Intermediate temperature
    sheet3 = wb.add_sheet('T_intermediate')
    # sheet1.write(row,col, data, style)
    sheet3.write(0, 0, 'T_intermediate, oC')
    sheet3.write(0, 1, 'COP_Cascade')
    k = 0
    parameters[1] = 0
    parameters[0] = 0  # K for cascade cycle
    parameters[3] = 0  # vapour quality of point 7
    for i in range(90,127,1):
        k = k+1
        parameters[2] = i   # intermediate temperature of the cascade cycle
        [COP,_] = cascade_cycle(refrigerant,parameters)
        print('COP at T_intermediate = ',i,'oC:', COP)
        sheet3.write(k, 0, i)
        sheet3.write(k, 1, COP)

    wb.save('COP_cascade.xls')
#=======================================================================================================================
# PERFORM DIFFERENTIAL EVOLUTION FOR CASCADE CYCLE WITH PURE REFRIGERANT
#=======================================================================================================================
def opt_cascade_pure():
    # define lower and upper bounds for every dimension
    #bounds = asarray([(0, 10),(30, 50),(100,120)])
    bounds = asarray([(0, 10), (0, 50), (90, 120)])
    # parameters[0]  delta_T_SH_IHX1 =0 - 10 K
    # parameters[1] delta_T_SH_IHX2 = 30-50 K
    # parameters[2] = T_intermediate = 110 - 120  deg C
    solution = differential_evolution(refrigerant,bounds)
    print('\nSolution: f([%s]) = %.5f' % (around(solution[0], decimals=5), solution[1]))
    result_cascade(refrigerant,solution[0])

#=======================================================================================================================
# PERFORM DIFFERENTIAL EVOLUTION FOR CASCADE CYCLE WITH MIXTURE REFRIGERANT
#=======================================================================================================================
def opt_cascade_mixed():
    # define lower and upper bounds for every dimension
    bounds = asarray([(0, 10),(0, 50),(100,120)])
    # parameters[0]  delta_T_SH_IHX1 =0 - 10 K
    # parameters[1] delta_T_SH_IHX2 = 30-50 K
    # parameters[2] = T_intermediate = 110 - 120  deg C

    #define mixture for the cascade cycle
    x1 = 0.1
    x2 = 1 - x1
    composition_top = [x1, x2]
    medium_top = medium1 + '&' + medium5
    MEDIUM_TOP = CP.AbstractState('REFPROP', medium_top)
    MEDIUM_TOP.set_mass_fractions(composition_top)
    x3 = 0.1
    x4 = 1 - x3
    composition_bot = [x3, x4]
    medium_bot = medium2 + '&' + medium6
    MEDIUM_BOT = CP.AbstractState('REFPROP', medium_bot)
    MEDIUM_BOT.set_mass_fractions(composition_bot)
    refrigerant.top = MEDIUM_TOP
    refrigerant.bot = MEDIUM_BOT
    solution = differential_evolution(refrigerant,bounds)
    print('\nSolution: f([%s]) = %.5f' % (around(solution[0], decimals=5), solution[1]))
    result_cascade(refrigerant,solution[0])

#=======================================================================================================================
#OPTIMISE THE CASCADE CYCLE FOR DIFFERENT MIXED REFRIGERNAT
#=======================================================================================================================


# RUN THE HEAT PUMP CYCLE
#=======================================================================================================================

#result_cascade(refrigerant,parameters)
#result_single_state_std(refrigerant)
#result_single_state_ihx(refrigerant)
#run_cascade()
#opt_cascade_pure()
#opt_cascade_mixed()
opt_all_mixture()