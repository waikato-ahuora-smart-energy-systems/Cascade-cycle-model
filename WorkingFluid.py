'''
Class of workingfluid and its subclasses
Qun Chen @ Massey University April 2022
Properties from CoolProp as it is free, Low level interface
Will be updated later with RefProp if needed
'''

import CoolProp
import numpy as np
import scipy as sp
import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI


class WorkingFluid:
    '''
    class WorkingFluid - properties and functions related to working fluids.
    Initialised by input string "wfname" for fluid name, variable mf for mass flow in kg/s,
    other parameters in type dict as Pf=, Hf=, Tf=, etc (argument **kwargs).
    This works for all the working fluid initialisation.
    class attribute proplist is a list of strings:
    # C: specific heat, J/kg/K;
    # D: mass density, kg/m3, 1/D is then the specific volume per mass, m3/kg; for humid air, D is for dew point, the mass density will be calculated from specific volume V
    # Ex: specific exergy, J/kg;
    # G: specific Gibbs energy, J/kg; H: specifc enthalpy, J/kg; m
    # K: thermal conductivity, W/m/K;
    # MW: molecular weight, kg/mol
    # P: pressure, Pa; Q: quality; S: specific entropy, J/kg/K; T: temperature, K;
    # V: viscosity, Pa.s; but NOTE: for humid air, V is for mixture volume per kg dry air.

    class attribute allFluids saves all the fluid objects defined in the programme.
    object attribute: fluid - a dictionary containing fluid properties
    '''
    proplist = ['C', 'D', 'G', 'H', 'K', 'MW', 'P', 'Pr', 'Q', 'S', 'T', 'V']
    allFluids = []

    def __init__(self, wfname, mf=1.0, **kwargs):
        self.fluid = {'fluid': wfname, 'mf': mf}  # fluid properties in the form of a dictionary
        for element in WorkingFluid.proplist:
            self.fluid[element + 'f'] = 0.0  # initialised by 0.0
        for key, value in kwargs.items():
            self.fluid[key] = value  # update the item in the fluid dict with inputs
        WorkingFluid.allFluids.append(wfname)  # save the fluid name


class HumidAir(WorkingFluid):
    '''
    class HumidAir, a subclass of WorkingFluid,
    fluid name HumidAir is initialised implicitly, other initialisation inputs are mf, and Tf, Hf, or Wf,etc
    class attribute Patm: the atmospheric pressure, Pa
    '''
    # class attributes
    Patm = 101325.0  # atmospheric pressure, Pa
    HAproplist = [item for item in WorkingFluid.proplist if
                  item != 'Q' or item != 'Pr' or item != 'MW' or item != 'P'] + [
                     'M']  # No Q quality or Pr Prandtl, M for humid air viscosity, Pa.s
    humidprop = ['RH', 'W', 'Tdp',
                 'Twb']  # RH for relative humidity; W for humidity ratio, kg/kg da; Tdp for dew point, K; Twb for wet-bulb temperature, K
    inputlist = ['H', 'T'] + humidprop

    def __init__(self, mf=1.0, **kwargs):
        super().__init__("HumidAir", mf)
        for item in HumidAir.humidprop:  # initialisation of humidity
            self.fluid[item + 'f'] = 0.0
        for key, value in kwargs.items():  # initialise the variables with input arguements in dictionary
            self.fluid[key] = value
        del self.fluid['Qf']  # remove unused variables - quality
        self.fluid['Pf'] = HumidAir.Patm
        self.fluid['MWf'] = 0.02897  # air molar mass, kg/mol
        self.fluid['Prf'] = 0.7  # air Prandtl number
        self.initlist = [prop for prop in HumidAir.inputlist if prop + 'f' in kwargs.keys()]
        for item in HumidAir.HAproplist + HumidAir.inputlist:
            if item != self.initlist[0] and item != self.initlist[1]:
                self.fluid[item + 'f'] = self.fluidprop(item, self.initlist[0], self.initlist[1])
        self.fluid['Exf'] = self.FluidExcergy()

    def humidconvert_T(self, output, sinput, vinput, ta=None):
        '''
        conversion btwn RH, Dp, W and WB, with dry bulb temperature as another input.'''
        if ta is None:
            ta = self.fluid['Tf']
        return HAPropsSI(output, 'T', ta + 273.15, sinput, vinput, 'P', self.fluid['Pf'])

    def humidconvert_H(self, output, sinput, vinput, ha=None):
        '''
        conversion btwn RH, Dp, W and WB, with enthalpy as another input.
        '''
        if ha is None:
            ha = self.fluid['Hf']
        return HAPropsSI(output, 'H', ha, sinput, vinput, 'P', self.fluid['Pf'])

    def fluidprop(self, propout, spropin1, spropin2, vpropin1=None, vpropin2=None):
        '''
        # to look up a property of humid air (defined by "propout" string)
        # with two input conditions spropin1 and spropin2, both strings, with values of vpropin1 and vpropin2,respectively,
        # default values of which are defined in self.fluid
        '''
        if vpropin1 is None:
            vpropin1 = self.fluid[spropin1 + 'f']
        if vpropin2 is None:
            vpropin2 = self.fluid[spropin2 + 'f']
        return HAPropsSI(propout, spropin1, vpropin1, spropin2, vpropin2, 'P', self.fluid['Pf'])

    def FluidExcergy(self, t0=15.0, wa=None, ta=None):
        '''
        to calculate the exergy of humid air
        :param t0: ambient temperature, in degree C
        :param wa: air humidity, kg/kg da
        :param ta: air temperature, in K
        :return: the exergy, J/kg
        '''
        if wa is None:
            wa = self.fluid['Wf']
        if ta is None:
            ta = self.fluid['Tf']
        ha0 = HAPropsSI('H', 'T', t0 + 273.15, 'W', wa, 'P', self.fluid['Pf'])
        sa0 = HAPropsSI('S', 'T', t0 + 273.15, 'W', wa, 'P', self.fluid['Pf'])
        ha = HAPropsSI('H', 'T', ta, 'W', wa, 'P', self.fluid['Pf'])
        sa = HAPropsSI('S', 'T', ta, 'W', wa, 'P', self.fluid['Pf'])
        return (ha - ha0) - t0 * (sa - sa0)


class Refrigerant(WorkingFluid):
    '''
    class Refrigerant, a subclass of WorkingFluid
    Initialised by inputting refrigerant name string, mf- mass flow, and two parameters in the list of
    Hf, Pf, Qf, Sf, Tf
    class attribute propdict contains keywords for CoolProp low level interface
    class attribute inputlist contains the input variables most commonly used.
    '''
    propkeys = ['iCpmass', 'iDmass', 'iGmass', 'iHmass', 'iconductivity', 'imolar_mass', 'iP', 'iPrandtl', 'iQ',
                'iSmass',
                'iT', 'iViscosity']
    propdict = dict(zip(WorkingFluid.proplist, propkeys))
    inputlist = ['H', 'P', 'Q', 'S', 'T']

    def __init__(self, refname, mf=1.0, **kwargs):  # initialise by fluid name, mass flow, kg/s and others
        super().__init__(refname, mf)
        for key, value in kwargs.items():
            self.fluid[key] = value
        self.initlist = [prop for prop in Refrigerant.inputlist if prop + 'f' in kwargs.keys()]
        self.fluid_HEOS = CP.AbstractState('HEOS', self.fluid['fluid'])
        self.criticalPoint()
        self.fluid['Tsat'] = self.Tsaturation()
        self.fluid['Psat'] = self.Psaturation()

    def criticalPoint(self):
        # find the critical point
        self.fluid['Tcrit'] = self.fluid_HEOS.T_critical()
        self.fluid['Pcrit'] = self.fluid_HEOS.p_critical()

    def Tsaturation(self):
        # find the saturation temperature when pressure is known.
        self.fluid_HEOS.update(CoolProp.PQ_INPUTS, self.fluid['Pf'], 0)
        return self.fluid_HEOS.T()

    def Psaturation(self):
        # find the saturation pressure when temperature is known.
        self.fluid_HEOS.update(CoolProp.TQ_INPUTS, self.fluid['Tf'], 0)
        return self.fluid_HEOS.p()

    # The following functions may be integrated into one function with arguments
    def fluidProps_HP(self):
        # find the refrigerant properties based on enthalpy and pressure
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.HmassP_INPUTS, self.fluid['Hf'], self.fluid['Pf'])
        for item in WorkingFluid.proplist:
            if item != 'P' and item != 'H':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def fluidProps_PT(self):
        # find the refrigerant properties based on temperature and pressure
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.PT_INPUTS, self.fluid['Pf'], self.fluid['Tf'])
        for item in WorkingFluid.proplist:
            if item != 'P' and item != 'T':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def fluidProps_PQ(self):
        # find the refrigerant properties based on quality and pressure
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.PQ_INPUTS, self.fluid['Pf'], self.fluid['Qf'])
        for item in WorkingFluid.proplist:
            if item != 'P' and item != 'Q':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def fluidProps_QT(self):
        # find the refrigerant properties based on quality and temperature
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.QT_INPUTS, self.fluid['Qf'], self.fluid['Tf'])
        for item in WorkingFluid.proplist:
            if item != 'T' and item != 'Q':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def fluidProps_PS(self):
        # find the refrigerant properties based on quality and temperature
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.PSmass_INPUTS, self.fluid['Pf'], self.fluid['Sf'])
        for item in WorkingFluid.proplist:
            if item != 'P' and item != 'S':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def fluidProps_HS(self):
        # find the refrigerant properties based on quality and temperature
        self.fluid_HEOS.specify_phase(self.fluid['Pf'], self.fluid['Tf'])
        self.fluid_HEOS.update(CoolProp.HmassSmass_INPUTS, self.fluid['Hf'], self.fluid['Sf'])
        for item in WorkingFluid.proplist:
            if item != 'H' and item != 'S':
                self.fluid[item + 'f'] = self.fluid_HEOS.keyed_output(Refrigerant.propdict[item])

    def phaseSpecify(self, pressure, temperature):
        # to determine the phase by pressure and temperature
        if pressure >= self.fluid['Pcrit']:
            if temperature >= self.fluid['Tcrit']:
                iphase = CoolProp.iphase_supercritical
            else:
                iphase = CoolProp.iphase_supercritical_liquid
        else:
            if temperature >= self.fluid['Tcrit']:
                iphase = CoolProp.iphase_supercritical_gas
            else:
                if pressure > self.fluid['Tsat']:
                    iphase = CoolProp.iphase_liquid
                elif pressure < self.fluid['Psat']:
                    iphase = CoolProp.iphase_gas
                else:
                    iphase = CoolProp.iphase_twophase
        return iphase

    def FluidExcergy(self, t0=15):
        hf0 = CP.PropsSI('H', 'T', t0 + 273.15, 'Q', 0, self.fluid['fluid'])
        sf0 = CP.PropsSI('S', 'T', t0 + 273.15, 'Q', 0, self.fluid['fluid'])
        return (self.fluid['Hf'] - hf0) - t0 * (self.fluid['Sf'] - sf0)


'''
to be further revised
class UWater(Refrigerant):
    # plant utility: cold water
    def __init__(self, mf=1.0,**kwargs):
        super().__init__("water", mf, **kwargs)
        self.fluid['Qf']=0.0
        self.fluid['Tsatf']=self.Tsaturation(self.fluid['Pf'])
        self.fluid['dTsub']=self.fluid['Tsatf']-self.fluid['Tf']
class USteam(Refrigerant):
    # plant utility: steam
    def __init__(self, mf=1.0,**kwargs):
        super().__init__("WATER", mf, **kwargs)
        self.fluid['Qf']=1.0
        self.fluid['Tsatf']=self.Tsaturation(self.fluid['Pf'])
        self.fluid['dTsup']=self.fluid['Tf']-self.fluid['Tsatf']
class RefrigerantMixture(Refrigerant):
    def __init__(self, mixname='Butane[0.5]&Heptane[0.5]',mf=1.0,**kwargs):
        namesplit=mixname.replace('[',' ').replace(']',' ').replace('&', ' ').split()
        self.nRefrig=len(namesplit)/2
        self.refriglist=[]
        self.fractionlist=[]
        self.CAS_reflist=[]
        for n in range(self.nRefrig):
            if n % 2 == 0:
                self.refriglist.append(namesplit[n])
            else:
                self.fractionlist.append(float(namesplit[n]))
        refkeylist=[]
        for n in range(len(self.refriglist)):
            refkeylist.append(['Refrig'+'{}'.format(n)])
            self.CAS_reflist.append(CP.get_fluid_param_string(self.refriglist[n],'CAS'))
        CP.set_config_bool(CP.OVERWRITE_BINARY_INTERACTION, True)
        CP.apply_simple_mixing_rule(self.CAS_reflist[0], self.CAS_reflist[1], 'linear') # Only for two refrigerants!!
        super().__init__(mixname, mf, **kwargs)


class FlueGas(WorkingFluid):
    components=['CO2','H2O','O2','N2']
    def __init__(self, mf=1.0, **kwargs):
        super().__init__('Flue gas', mf)
        pass
    # to be further developed
'''