import CoolProp as CP



class Input:
    pass


model = Input()
# Heat pump design models:
model.limit_sc        = 0.95           # Limit above which process is considered as supercritcal
model.error_delta_T   = 1e-2           # Allowable error on pinch temperature
model.iter_max        = 10             # Number of maximum iterations

#N-R solver of the HP cycle design
model.tol_HP                  = model.error_delta_T
model.step                    = model.error_delta_T*100
#model.relaxation            = 1

# Discretization
model.n_states_evap = 5
model.n_states_sh = 2
model.n_states_cond = 5
model.n_states_desup = 2
model.n_states_sub = 2
model.n_states_ihx = 5
model.n_sates_cascadeHX = 2

# Define the working fluids
# Heat sink/source side
medium_sink = 'air'
SINK = CP.AbstractState('HEOS', medium_sink)     #Using low lever interface to call Refprop/Coolprop for a faster computational time
medium_source = 'water' # water
SOURCE = CP.AbstractState('REFPROP', medium_source)
# Refrigerant side - pure refrigerant
medium1 = 'Hexane'
medium2 = 'Iso-Butane'
medium3 = 'Heptane'
medium4 = 'Pentane'
medium5 = 'Neopentane'
medium6 = 'Acetone'


MEDIUMPure1 = CP.AbstractState('REFPROP', medium1)
MEDIUMPure2 = CP.AbstractState('REFPROP', medium2)
MEDIUMPure3 = CP.AbstractState('REFPROP', medium3)
MEDIUMPure4 = CP.AbstractState('REFPROP', medium4)
MEDIUMPure5 = CP.AbstractState('REFPROP', medium5)
MEDIUMPure6 = CP.AbstractState('REFPROP', medium6)

# Define the fluids struct for the Refprop/Coolprop functions
fluids = Input()
refrigerant = Input()
refrigerant.top          = MEDIUMPure1
refrigerant.bot          = MEDIUMPure2
refrigerant.ref_1        = MEDIUMPure3
refrigerant.ref_2        = MEDIUMPure4
refrigerant.ref_3        = MEDIUMPure5
refrigerant.ref_4        = MEDIUMPure6

fluids.source            = SOURCE
fluids.sink              = SINK

# Define the fixed inputs for cycle design

sink = Input()
source = Input()
ref = Input()
comp = Input()


sink.composition = 1
source.composition = 1


source.T_in = 273.15 + 40  # K 40
source.p_in = 1e5  # Pa
source.T_out_design = 273.15 + 13  # K 25

# Sizing of the HP by fixing the heat source load:
source.Qdot_design = 500 * 1000  # W 500*1000

sink.T_in = 273.15 + 15 # K 50
sink.T_out_design = 273.15 + 200# K 80

if medium_sink == 'water':
    sink.p_in = 30e5  # Pa - no phase change if water used as heat sink fluid, and outlet hot water temperature > 200 oC
elif medium_sink == 'air':
    sink.p_in = 1e5

# Temperature Pinch Differences
source.delta_T_pinch_min = 3        # K
ref.delta_T_pinch_IHX_min = 3       # K
if medium_sink == 'water':
    sink.delta_T_pinch_min = 3      # K
elif medium_sink == 'air':
    sink.delta_T_pinch_min = 20     # K

# Super heat temperature
ref.delta_T_minSH = 5               # K minimum superheat temperature
#ref.delta_T_SH_IHX1 = 0

# compressor
comp.eta_is      = 0.7              # [-] Isentropic Efficiency compressor
comp.eta_motor   = 0.96

