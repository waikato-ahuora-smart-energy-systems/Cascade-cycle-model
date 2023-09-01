import CoolProp as CP
from numpy.random import rand
from numpy.random import choice
from numpy import asarray
from numpy import clip
from numpy import argmin
from numpy import min
from numpy import around
from Cycles.Cascade import cascade_cycle
from xlwt import Workbook
from DifferentialEvolution import differential_evolution
from Economics.Cascade_Economic_Analysis import economic_analysis



# OPTIMISE THE CASCADE CYCLE FOR DIFFERENT REFRIGERNAT
# =======================================================================================================================

def opt_all_mixture():

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
    sheet1.write(0, 61, 'x7')



    class empty_class:
        pass

    refrigerant = empty_class()
    cascade_results     = empty_class()
    economics   = empty_class()

    #list_medium_bot = ['Butane','Butene','cis-Butene','Isobutane','Isobutene','Neopentane','R1234ze(Z)']
    list_medium_bot = ['Butane']
    list_medium_top = [ 'Acetone','Benzene','Cyclohexane','Cyclopentane','Ethanol','Heptane', 'Hexane','Isohexane', 'Isooctane', 'Methylcyclohexane','Toluene']
    #list_medium_top = ['Acetone','Benzene','Cyclohexane','Cyclopentane','Ethanol','Heptane', 'Hexane','Isohexane', 'Isooctane', 'Methylcyclohexane','Toluene']


    #bounds = asarray([(0, 10), (0, 50), (110, 143)]) #
    #bounds = asarray([(0, 10), (0, 50), (110, 140)])
    k = 0
    for medium_top in list_medium_top:
        MEDIUM_TOP = CP.AbstractState('REFPROP', medium_top)
        for medium_bot in list_medium_bot:
            k = k+1
            MEDIUM_BOT = CP.AbstractState('REFPROP', medium_bot)
            refrigerant.top = MEDIUM_TOP
            refrigerant.bot = MEDIUM_BOT
            MEDIUM_BOT.update(CP.PT_INPUTS, 101325, 300)
            T_crit_bot = MEDIUM_BOT.T_critical() - 273.15 # oC
            #bounds = asarray([(0, 20), (0, 50), (100, T_crit_bot -3),(0,0.5)])
            bounds = asarray([(0, 20), (0, 15), (100, T_crit_bot -3),(0,0.7)])
            # define lower and upper bounds for every dimension

            print(k, ' top cycle refrigerant:', medium_top)
            print(k, ' bottom cycle refrigerant:', medium_bot)

            solution = differential_evolution(refrigerant, bounds)

            print('\nSolution: f([%s]) = %.5f' % (around(solution[0], decimals=5), solution[1]))
            [COP,cascade_results] = cascade_cycle(refrigerant,solution[0])
            economics = economic_analysis(cascade_results)
            sheet1.write(k, 0, medium_top)
            sheet1.write(k, 1, medium_bot)
            sheet1.write(k, 2, solution[0][0])
            sheet1.write(k, 3, solution[0][1])
            sheet1.write(k, 4, solution[0][2])
            sheet1.write(k, 5, cascade_results.T_cond_top - 273.15) # oC
            sheet1.write(k, 6, COP)
            sheet1.write(k, 7, cascade_results.exergy_eff)
            sheet1.write(k, 8, cascade_results.n_Lorren)
            sheet1.write(k, 9, cascade_results.E_des_top_cond) # (kJ/kg)
            sheet1.write(k, 10, cascade_results.E_des_IHX3) # (kJ/kg)
            sheet1.write(k, 11, cascade_results. E_des_CHX) # (kJ/kg)
            sheet1.write(k, 12, cascade_results.E_des_bot_evap) # (kJ/kg)
            sheet1.write(k, 13, cascade_results.mdot_ref_bot)  # (Kg/s)
            sheet1.write(k, 14, cascade_results.mdot_ref_top)  # (Kg/s)
            sheet1.write(k, 15, cascade_results. V_suction_com1 * 3600) # (m3/h)
            sheet1.write(k, 16, cascade_results.V_suction_com2 * 3600) # (m3/h)
            sheet1.write(k, 17, cascade_results.V_suction_com3 * 3600) # (m3/h)
            sheet1.write(k, 18, cascade_results.T_suction_com1 - 273.15) # (oC)
            sheet1.write(k, 19, cascade_results.T_suction_com2 - 273.15) # (oC)
            sheet1.write(k, 20, cascade_results.T_suction_com3 -273.15) # (oC)
            sheet1.write(k, 21, cascade_results.T_discharge_com2 -273.15) # (oC)
            sheet1.write(k, 22, cascade_results.T_discharge_com3 - 273.15) # (oC)
            sheet1.write(k, 23, cascade_results.specific_volume_2) # (m3/kg)
            sheet1.write(k, 24, cascade_results.specific_volume_3) # (m3/kg)
            sheet1.write(k, 25, cascade_results.specific_volume_4) # (m3/kg)
            sheet1.write(k, 26, cascade_results.specific_volume_14) # (m3/kg)
            sheet1.write(k, 27, cascade_results.specific_volume_15) # (m3/kg)
            sheet1.write(k, 28, cascade_results.eta_vol_com1)
            sheet1.write(k, 29, cascade_results.eta_vol_com2)
            sheet1.write(k, 30, cascade_results.eta_vol_com3)
            sheet1.write(k, 31, cascade_results.p_evap_bot/1e5) # (bar)
            sheet1.write(k, 32, cascade_results.p_mid_bot/1e5) # (bar)'
            sheet1.write(k, 33, cascade_results.p_cond_bot/1e5) # (bar)
            sheet1.write(k, 34, cascade_results.p_evap_top/1e5) # (bar)
            sheet1.write(k, 35, cascade_results.p_cond_top/1e5) # (bar)
            sheet1.write(k, 36, cascade_results.Qdot_sink/1000) # (kW)
            sheet1.write(k, 37, cascade_results.W_dot/1000) # (kW)

            sheet1.write(k, 38, economics.PEC_com1/1000) # (1000 NZD)
            sheet1.write(k, 39, economics.PEC_com2/1000) # (1000 NZD)
            sheet1.write(k, 40, economics.PEC_com3/1000) # (1000 NZD)
            sheet1.write(k, 41, economics.PEC_Evap/1000) # (1000 NZD)
            sheet1.write(k, 42, economics.PEC_Top_Cond/1000) # (1000 NZD)
            sheet1.write(k, 43, economics.PEC_Cax/1000) # (1000 NZD)
            sheet1.write(k, 44, economics.PEC_IHX2/1000) # (1000 NZD)
            sheet1.write(k, 45, economics.PEC_IHX3/1000) # (1000 NZD)
            sheet1.write(k, 46, economics.TCI/1000) # (1000 NZD)
            sheet1.write(k, 47, economics.TCI_spec) # (NZD/kW)
            sheet1.write(k, 48, economics.NPV/1000) # (1000 NZD)
            sheet1.write(k, 49, economics.PBT) # (years)
            sheet1.write(k, 50, economics.Energy_Produced_Yearly) # (MWh)
            sheet1.write(k, 51, economics.c_H) # (NZD / MWh)
            sheet1.write(k, 52, economics.CF_FCel/1000) # (1000 NZD/year)
            sheet1.write(k, 53, economics.CF_FChsource/1000) # (1000 NZD/year)
            sheet1.write(k, 54, economics.CF_rev/1000) # (1000 NZD/year)')
            sheet1.write(k, 55, economics.CF_TCI/1000) # (1000 NZD/year)')
            sheet1.write(k, 56, economics.CF_OM/1000) # (1000 NZD/year)
            sheet1.write(k, 57, economics.CF_total/1000) # (1000 NZD/year)
            sheet1.write(k, 58, cascade_results.eta_is_comp1)
            sheet1.write(k, 59, cascade_results.eta_is_comp2)
            sheet1.write(k, 60, cascade_results.eta_is_comp3)
            sheet1.write(k, 61, solution[0][3])

            wb.save('cascade_analysis.xls')

