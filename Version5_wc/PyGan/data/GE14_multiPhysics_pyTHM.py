## PyGan script to couple DONJON procedures with pyTHM solution
# Date : 27/08/2026
# Author : Raphael Guasch
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR fuel assembly

from pyTHM.Solver.main import pyTHM_solver 
from pyTHM.WaterProperties.waterProperties import FastIAPWS
from export_results import export_results_to_csv
from setup_case import generate_c2m_procedures, initialise_starterDD_case
from helpers import underRelaxation, convergence_TH, convergence_POW, convergence_pcm, compute_RMS, compute_residuals, quickPlot
from pygan_procedures import initialize_donjon, neutronics_solution
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *

import time


##
########## User input ##########
### Case to be treated : 
assembly_id = "DOM_VAN"

########## Algorithm parameters ###########
start_time = time.time()
solveConduction = True
nIter = 1000 #1000
tol_TH_Tf = 1e-2 # Convergence criterion for TFuel, in K
tol_TH_Tc = 1e-2 # Convergence criterion for TCool, in K
tol_TH_Dc = 1e-2 # Convergence criterion for DCool, in kg/m3
tol_POW = 1e-3 # Convergence criterion for the power axial distribution, in % of nodal power
tol_Keff = 1e-2 # Convergence criterion for the Keff value in pcm

relax_Pow = False # Under relaxation of the Power distribution for the next iteration
Pow_underRelaxationFactor = 0.5 # 0.8, 0.5, 0.2, Under relaxation factor for the power axial distribution to be tested, 0.1 used in Serpent/OpenFoam coupling

relax_TH = False # Under relaxation of the TH fields for the next iteration
TH_underRelaxationFactor = 0.9 # 0.8, 0.5, 0.2, Under relaxation factor for the TH fields to be tested

########## Mesh parameters ###########
zPlotting = [] #If empty, no plotting of the axial distribution of the fields, otherwise, list of the axial positions where the fields are plotted
## Meshing parameters:
If = 8
I1 = 3

# Sensitivity to the meshing parameters
nz = 128 # 32, 64, 128 etc for DOM_VAN and PSZ_DOM_VAN cases 
Ptot = 2.416e6 # 3.6 MW ~ 870 MWth / 240 fuel assemblies for SMR BWR core

########## Choice of Thermalhydraulics correlation ##########
voidFractionCorrel = "GEramp" # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
frfaccorel = "Churchill" # 'base', 'blasius', 'Churchill', 
P2Pcorel = "friedel" # 'base', 'HEM1', 'HEM2', 'MNmodel', 'lockhartMartinelli'
numericalMethod = "FVM" # "FVM": Solves the system using matrix inversion with preconditioning.
                         # "GaussSiedel" : Applies the Gauss-Seidel iterative solver.
                         # "BiCG" : Uses the BiConjugate Gradient method for solving non-symmetric or indefinite matrices.
                         # "BiCGStab" : Applies the BiCGStab (BiConjugate Gradient Stabilized) method to ensure faster and more stable convergence.

########## Thermal hydraulics parameters ##########
## Geometric parameters
if assembly_id == "DOM":
    # Run single channel case with homogenious axial fuel channel geometry
    GE14_fuel_bundle_DOM_geometry = "GEOM_single_bundle_GE14_DOM.yaml"

    fuel_channel_geometry_data, lcm_variables_names, model_initialiser = initialise_starterDD_case(geometry_yaml=GE14_fuel_bundle_DOM_geometry, 
                                                                                              case_name=f"GE14_bundle_{assembly_id}",
                                                                                              number_of_axial_meshes=nz, 
                                                                                              reactor_power=Ptot)
    compo_file_names = ["_CPO_GE14_DOM_DIAG"] # Name of the COMPO object to be used in the neutronics solution
elif assembly_id == "DOM_VAN":
    # Run single channel case with heterogenious axial fuel channel geometry
    GE14_fuel_bundle_DOM_VAN_geometry = "GEOM_single_bundle_GE14_DOM_VAN.yaml"

    fuel_channel_geometry_data, lcm_variables_names, model_initialiser = initialise_starterDD_case(geometry_yaml=GE14_fuel_bundle_DOM_VAN_geometry, 
                                                                                              case_name=f"GE14_bundle_{assembly_id}",
                                                                                              number_of_axial_meshes=nz, 
                                                                                              reactor_power=Ptot)
    compo_file_names = ["_CPO_GE14_DOM_DIAG", "_CPO_GE14_VAN_DIAG"] # Name of the COMPO object to be used in the neutronics solution

else:
    raise ValueError(f"Case {assembly_id} not supported yet")

# recover cle2000 variable names :
fuel_map_name = lcm_variables_names["Fuel_Map"]
matex_name = lcm_variables_names["Material_Indexation"]
track_name = lcm_variables_names["Tracking"]
thdata_name = lcm_variables_names["Thermal_hydraulics_data"]
compo_lcm_names = lcm_variables_names["MultiCompos"]

procedure_names = {"initialisation": "IniDonjon",
                    "cross_section_interpolation": "InterpXS",
                    "neutronics_solve": "NeutroSol"}
initialisation_proc_name = procedure_names["initialisation"]

## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
tInlet = 270 + 273.15 # K, for BWRX-300 SMR core
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1827 kg/s for a whole core with 240 fuel channels
# Assume 88% of the total flow goes to active flow
# Assume that coolant flow rate is evenly split between all fuel channels for now.
massFlowRate = (1827 / 240) * 0.88  # kg/s

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
kClad = 21.5 #21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity


    
########### Print options ###########
FULL_PRINT = True # Print all the fields at each iteration, create figures in tmp/rundir/ ?

############ Nuclear Parameters ###########
# Number of fuel rods and assemblies for a small modular Boiling Water Reactor core

n_rods = 92 # This value is for the ATRIUM-10 fuel assembly (10*10-9), in a GNF2 fuel assembly (BWRX-300) there are 92 fuel rods per assembly
n_assmblies = 240 # This value is for the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"
full_core_power = 870e6 # W, full core thermal power of the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"

volumic_mass_U = 10970 # kg/m3 'Thermophysical Properties of MOX and UO2 Fuels Including the Effects of Irradiation' - ORNL/TM-2000/351, Popov, Carbajo, Ivanov, Yoder. November 2000.

## Fuel assembly parameters :
fuelRadius = fuel_channel_geometry_data["fuel_data"]["fuel_radius"]
height = fuel_channel_geometry_data["fuel_data"]["max_rod_length"]
pitch = fuel_channel_geometry_data["active_flow_data"]["pitch"]
z_mesh = fuel_channel_geometry_data["z_mesh"]
Fuel_volume = n_rods*np.pi*fuelRadius**2*height # m3

Fuel_mass = Fuel_volume*volumic_mass_U*1000 # g
print(f"Fuel mass = {Fuel_mass} kg")
Bundle_volume = Fuel_volume / nz # m3, Bundle <=> 1 axial slice of the fuel channel

### Generate cle2000 procedures

## Call starterDD to produce the relevant c2m procedures.
generate_c2m_procedures(model_initialiser, procedure_names)

########## Fields of the TH problem ##########
TeffFuel = []
Twater = []
rho = []
voidFraction = []
P_ = [] # Pressure field
U_ = [] # Velocity field
H_ = [] # Enthalpy field


Residuals_TeffFuel = []
Residuals_Twater = []
Residuals_rho = []
Residuals_voidFraction = []

####### Power shape to store the axial power shape at each iteration
Keffs = [] # list to store the Keff values at each iteration and check convergence on it
Power_Distribs = []
Relaxed_Power_Distribs = []

Residuals_Power_Distribs = []
Residuals_Relaxed_Power_Distribs = []


######## Creation of results directory ##########
path=os.getcwd()
root_save_dir = f"GE14_multiPhysics_PyGan_{assembly_id}"
a=os.path.exists(root_save_dir)
if a==False:
    os.mkdir(root_save_dir)
print(path)

SAVE_DIR = f"{root_save_dir}/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}"

a=os.path.exists(SAVE_DIR)
if a==False:
    os.makedirs(SAVE_DIR)
########## End creation of results directory ##########

##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

print("$$$ - multiPhysics.py : BEGIN INITIALIZATION- $$$")
# 1.) Initialize the TH parameters inteprolation table :
fast_iapws = FastIAPWS(pOutlet)
# 2.) Initialize guess of sine-shaped power distribution, normalised to Ptot : run pyTHM

THComponentIni =  pyTHM_solver(
        case_name=assembly_id,
        channel_type="square",
        geometric_data = fuel_channel_geometry_data,
        tInlet=tInlet,
        pOutlet=pOutlet,
        qFlow=massFlowRate,
        Powtot=Ptot,
        power_distribution="sine",
        fraction_pow_fuel=1.0,
        k_fuel=kFuel,
        H_gap=Hgap,
        k_clad=kClad,
        I_f=8,  
        I_c=3,
        water_rod=False,
        water_rod_holes=False,
        plot_at_z=[],
        solveConduction=True,
        fast_iapws_table=fast_iapws,
        dt=0,
        t_tot=0,
        frfaccorel=frfaccorel,
        P2Pcorel=P2Pcorel, 
        voidFractionCorrel=voidFractionCorrel, 
        numericalMethod=numericalMethod,
    )

TeffIni, TwaterIni, rhoIni, voidFracIni, Pini, Uini, Hini = THComponentIni.get_TH_parameters() # renamed this function to be clear about what it does 
#                                                           ---> Fuel and coolant temperature + coolant density aren't nuclear parameters per-se
#


TeffFuel.append(TeffIni)
Twater.append(TwaterIni)
rho.append(rhoIni)
voidFraction.append(voidFracIni)
P_.append(Pini)
U_.append(Uini)
H_.append(Hini)

print(f"$$ - FIRST THM resolution iter = 0 : TeffFuel = {TeffIni}, Twater = {TwaterIni}, rho = {rhoIni}")
print(f"After initialization : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

# 2.2) Create LCM object to store the TH data to be used in the neutronics solution
print("Creating initial THData object")
THData = lcm.new('LCM',thdata_name)
THData['TFuelList'] = np.array(np.clip(TeffIni, a_min=550.0, a_max=1999.0), dtype='f')
THData['TCoolList'] = np.array(TwaterIni, dtype='f')
THData['DCoolList'] = np.array(rhoIni/1000, dtype='f')
THData.close() # close without erasing

## 3.) Initializing DONJON solution

Fmap, Matex, Track, Cpos = initialize_donjon(fuel_map_name, matex_name, track_name, compo_lcm_names, compo_file_names, initialisation_proc_name)

# 4.) Iterative scheme for coupled neutronics and thermalhydraulics solution
# Attempting to run multiPhysics scheme.
# Everytime, solve neutron transport, update axial power shape, solve TH, update TH data, and repeat.
powi = Ptot/1e6
# check powi == PFiss
print(f"powi = {powi} MW and PFiss = {Ptot} W")

## Multi-Physics resolution
iter = 0
conv = False
while not conv:
    convKeff = False
    convTH = False
    convPow = False
    iter+=1
    # 4.1) Neutronics solution for TH parameters obtained with initial TH solution
    if iter==1:
        Flux = None 
        Power = None

    Power, Flux, Fmap, Matex, Track = neutronics_solution(procedure_names, iter, Fmap, Matex, THData, Track, Flux, Cpos)
    # 4.2.2) Recover the Power axial distribution obtained from neutronics calculation --> used to update the axial power shape
    Keff = Flux["K-EFFECTIVE"][0] 
    Keffs.append(Keff)
    PowerDistribution = Power["POWER-DISTR"] 
    Flux_value = Power["FLUX"]

    Power_Distribs.append(PowerDistribution) # Store the axial power shape at each iteration
    # 4.3) Update the axial power shape to be used in the TH solution
    print(f"Power distribution : {PowerDistribution}")

    # 4.3.1) Under relaxation of the Power distribution for the next iteration
    if iter > 1:
        Residuals_Power_Distribs.append(compute_residuals(Power_Distribs))
        RMS_power_residual = compute_RMS(Residuals_Power_Distribs[-1])
        print(f"RMS power residual = {RMS_power_residual} at iter {iter}")
        if relax_Pow:
            PowerDistribution_relaxed = underRelaxation(PowerDistribution, Relaxed_Power_Distribs[-1], Pow_underRelaxationFactor)
            Relaxed_Power_Distribs.append(PowerDistribution_relaxed)
            Residuals_Relaxed_Power_Distribs.append(compute_residuals(Relaxed_Power_Distribs))
    # 4.4) Erase the THData object to store the updated TH data
    THData.erase()

    # 5.) TH procedure for updated power shape
    ############# Thermalhydraulic part ##############
    if relax_Pow and iter > 1:
        updated_power_distribution = PowerDistribution_relaxed
    else:
        updated_power_distribution = PowerDistribution

    # 5.1) TH resolution with updated power shape :
    print(f"At iter {iter} : THM resolution with updated power shape")
    THMComponent = pyTHM_solver(
        case_name=assembly_id,
        channel_type="square",
        geometric_data = fuel_channel_geometry_data,
        tInlet=tInlet,
        pOutlet=pOutlet,
        qFlow=massFlowRate,
        Powtot=Ptot,
        power_distribution=updated_power_distribution,
        fraction_pow_fuel=1.0,
        k_fuel=kFuel,
        H_gap=Hgap,
        k_clad=kClad,
        I_f=8,  
        I_c=3,
        water_rod=False,
        water_rod_holes=False,
        plot_at_z=[],
        solveConduction=True,
        fast_iapws_table=fast_iapws,
        dt=0,
        t_tot=0,
        frfaccorel=frfaccorel,
        P2Pcorel=P2Pcorel, 
        voidFractionCorrel=voidFractionCorrel, 
        numericalMethod=numericalMethod,
    )
    # 5.2) recover the new TH data
    TeffTEMP, TwaterTEMP, rhoTEMP, voidFracTEMP, P_TEMP, U_TEMP, H_TEMP = THMComponent.get_TH_parameters()
    
    TeffFuel.append(TeffTEMP)
    Residuals_TeffFuel.append(compute_residuals(TeffFuel))
    Twater.append(TwaterTEMP)
    Residuals_Twater.append(compute_residuals(Twater))
    rho.append(rhoTEMP)
    Residuals_rho.append(compute_residuals(rho))
    voidFraction.append(voidFracTEMP)
    Residuals_voidFraction.append(compute_residuals(voidFraction))
    P_.append(P_TEMP)
    U_.append(U_TEMP)
    H_.append(H_TEMP)


    print(f"THM resolution at iter={iter} : All lists TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

    if FULL_PRINT:
        quickPlot(z_mesh, TeffFuel, "Fuel temperature convergence", "height (m)", "TFuel (K)", "TFuel_convergence.png", path, SAVE_DIR)
        quickPlot(z_mesh, Twater, "Coolant temperature convergence", "height (m)", "TCool (K)", "TCool_convergence.png", path, SAVE_DIR)
        quickPlot(z_mesh, rho, "Coolant density convergence", "height (m)", "DCool (kg/m3)", "DCool_convergence.png", path, SAVE_DIR)

    
    # 5.3) Update the TH data to be used in the next neutronics solution
    
    THData = lcm.new('LCM','THData') # Create new LCM object to store the updated TH data
    THData['TFuelList']    = np.array(np.clip(TeffTEMP, a_min=550.0, a_max=1999.0), dtype='f')
    THData['TCoolList'] = np.array(TwaterTEMP, dtype='f')
    THData['DCoolList'] = np.array(rhoTEMP/1000, dtype='f') # Convert density from kg/m3 to g/cm3
    THData.close() # close without erasing

    # 5.4) Under relaxation of TH fields for the next iteration : check if it is needed ?
    if relax_TH:
        TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], TH_underRelaxationFactor)
        Twater[-1] = underRelaxation(Twater[-1], Twater[-2], TH_underRelaxationFactor)
        rho[-1] = underRelaxation(rho[-1], rho[-2], TH_underRelaxationFactor)

    
    print(f"$$ - END iter = {iter}")
    # 7.) Check for convergence on Keff value : add more converegence criteria for other fields ?
    if iter>1 and convergence_pcm(Keffs[-1], Keffs[-2], tol_Keff):
        print("Convergence on Keff reached after ", iter, " iterations")
        convKeff = True
    if iter>1 and convergence_TH(rho[-1], rho[-2], tol_TH_Dc) and convergence_TH(TeffFuel[-1], TeffFuel[-2], tol_TH_Tf) and convergence_TH(Twater[-1], Twater[-2], tol_TH_Tc):
    #if iter>1 and convergence_TH(rho[-1], rho[-2], tol_TH_Dc) and convergence_TH(Twater[-1], Twater[-2], tol_TH_Tc):
        print("Convergence on TH fields reached after ", iter, " iterations")
        convTH = True
    if iter>1 and convergence_POW(Power_Distribs[-1], Power_Distribs[-2], tol_POW):
        convPOW = True
        print("Convergence on Power reached after ", iter, " iterations")
    if convKeff and convTH and convPOW:
        conv = True
        total_iter = iter
        print("Convergence on Keff, TH fields and Power reached after ", iter, " iterations")
        print("Convergence criteria enforced : tol_Keff = ", tol_Keff, " tol_TH_Tf = ", tol_TH_Tf, " tol_TH_Tc = ", tol_TH_Tc, " tol_TH_Dc = ", tol_TH_Dc, " tol_POW = ", tol_POW)
    if iter == nIter-1:
        total_iter = iter
        print("Convergence not reached after ", iter, " iterations")
        break


print("$$$ - multiPhysics.py : END ITERATIONS - $$$")


# 8.) Check the results
print("$$$ - multiPhysics.py : BEGIN RESULTS - $$$")

print(f"Checking Fields lengths : Keff : {len(Keffs)}, TeffFuel : {len(TeffFuel)}, Twater : {len(Twater)}, rho : {len(rho)}, Power_Distrib : {len(Power_Distribs)}")
print(f"Checking residuals lengths : Residuals_TeffFuel : {len(Residuals_TeffFuel)}, Residuals_Twater : {len(Residuals_Twater)}, Residual_rho : {len(Residuals_rho)}, Residuals_Power_Distrib : {len(Residuals_Power_Distribs)}")

print("$$$ - multiPhysics.py : END RESULTS - $$$")


# 9.) Plot the results

if relax_Pow:
    relaxPOW_id = f"relaxedPOW_{Pow_underRelaxationFactor}"
else:
    relaxPOW_id = "non_relaxedPOW"
if relax_TH:
    relaxTH_id = f"relaxedTH_{TH_underRelaxationFactor}"
else:
    relaxTH_id = "non_relaxedTH"

SAVE_FIG_DIR = f"{SAVE_DIR}/Figures"
SAVE_DATA_DIR = f"{SAVE_DIR}/Data"
print(f"SAVE_FIG_DIR = {SAVE_FIG_DIR}")

a=os.path.exists(SAVE_FIG_DIR)
if a==False:
    os.makedirs(SAVE_FIG_DIR)

print(f"SAVE_DIR = {SAVE_DIR}")
print(f"Keffs = {Keffs}")
# Plot the results
quickPlot(range(len(Keffs)), Keffs, "Keff convergence", "iteration", "Keff", "Keff_convergence.png", path, SAVE_FIG_DIR)

quickPlot(z_mesh, TeffFuel, "Fuel temperature convergence", "height (m)", "TFuel (K)", f"TFuel_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Twater, "Coolant temperature convergence", "height (m)", "TCool (K)", f"TCool_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, rho, "Coolant density convergence", "height (m)", "DCool (kg/m3)", f"DCool_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, voidFraction, "Void fraction convergence", "height (m)", "Void fraction", f"Void_fraction_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Power_Distribs, "Power distribution convergence", "height (m)", "Power (kW)", f"Power_distribution_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, P_, "Pressure convergence", "height (m)", "Pressure (Pa)", f"pressure_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)

quickPlot(z_mesh, Residuals_TeffFuel, "Fuel temperature residuals", "height (m)", "Residuals", f"TFuel_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_Twater, "Coolant temperature residuals", "height (m)", "Residuals", f"TCool_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_rho, "Coolant density residuals", "height (m)", "Residuals", f"DCool_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_voidFraction, "Void fraction residuals", "height (m)", "Residuals", f"Void_fraction_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_Power_Distribs, "Power distribution residuals", "height (m)", "Residuals", f"Power_distribution_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
#print(f"Residuals Power distribution : {Residuals_Power_Distribs}")
#quickPlot(z_mesh, Residuals_Volumic_Powers, "Axial power shape residuals", "height (m)", "Residuals", f"Qfiss_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
print("$$$ - multiPhysics.py : END of PLOTTING - $$$")


# 10.) Save the results for exportation to Serpent/OpenFoam/GeN-Foam
# Save the results in a file

a=os.path.exists(SAVE_DATA_DIR)
if a==False:
    os.makedirs(SAVE_DATA_DIR)
os.chdir(SAVE_DATA_DIR)
print(f"Save_data_dir = {SAVE_DATA_DIR}")
case = assembly_id
TeffFuel = np.array(TeffFuel)
Twater = np.array(Twater)
rho = np.array(rho)

voidFraction = np.array(voidFraction)
Pressures = np.array(P_)
Velocities = np.array(U_)
Enthalpies = np.array(H_)

Power_Distrib = np.array(Power_Distribs)
Residuals_TeffFuel = np.array(Residuals_TeffFuel)
Residuals_Twater = np.array(Residuals_Twater)
Residual_rho = np.array(Residuals_rho)
Residuals_Power_Distrib = np.array(Residuals_Power_Distribs)



np.savetxt(f"TeffFuel_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", TeffFuel[-1])
np.savetxt(f"Keffs_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Keffs)
np.savetxt(f"Twater_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Twater[-1])
np.savetxt(f"rho_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", rho[-1])
np.savetxt(f"Pressure_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Pressures[-1])
np.savetxt(f"Velocity_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Velocities[-1])
np.savetxt(f"Enthalpy_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Enthalpies[-1])
np.savetxt(f"voidFraction_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", voidFraction[-1])
np.savetxt(f"Power_Distrib_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Power_Distrib[-1])
np.savetxt(f"Residuals_TeffFuel_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_TeffFuel[-1])
np.savetxt(f"Residuals_Twater_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_Twater[-1])
np.savetxt(f"Residuals_rho_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residual_rho[-1])
np.savetxt(f"Residuals_voidFraction_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_voidFraction[-1])
np.savetxt(f"Residuals_Power_Distrib_{case}_mesh{nz}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_Power_Distrib[-1])

# 11.) Save power axial distribution and flux axial distribution for the last iteration
os.chdir(path)
shutil.copyfile("Flux01.res", f"{SAVE_DATA_DIR}/Flux01.res")
shutil.copyfile("Flux02.res", f"{SAVE_DATA_DIR}/Flux02.res")
shutil.copyfile("Pdistr.res", f"{SAVE_DATA_DIR}/Pdistr.res")

shutil.copyfile("GE14_multiPhysics_pyTHM.result", f"{SAVE_DATA_DIR}/GE14_multiPhysics_pyTHM_{relaxPOW_id}_{relaxTH_id}.result")

#
print("$$$ - multiPhysics.py : END of EXPORTS - $$$")
print("$$$ - multiPhysics.py : END OF SCRIPT - $$$")

print("$$$ - multiPhysics.py : SANITY CHECK - $$$")

# export to pandas data Frame
results = {"z_mesh": z_mesh, 
            "TFUEL": TeffFuel[-1], 
            "TCOOL": Twater[-1], 
            "DCOOL": rho[-1], 
            "EPS": voidFraction[-1], 
            "PRESSURE": P_[-1], 
            "VELOCITY": U_[-1], 
            "ENTHALPY": H_[-1], 
            "POWER": Power_Distribs[-1], 
            "KEFF": Keffs}
export_results_to_csv(nz, Ptot, "cosine", pitch, height, results, iter, SAVE_DIR)
