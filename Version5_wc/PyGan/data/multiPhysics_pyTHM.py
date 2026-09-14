## PyGan script to couple DONJON procedures with THM solution
# Date : 20/09/2024
# Author : Clément HUET, Raphaël GUASCH
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR pincell
from THM_main import Version5_THM_prototype as THM_prototype
from THM_main import plotting
from export_results import export_results_to_csv, export_to_unconverged
from iapws import IAPWS97
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *

import time



########## Begin helper functions declaration ##########
## Fonction used for the convergence algorithm
def underRelaxation(Field, OldField, underRelaxationFactor):
    print(f"Field = {Field}")
    print(f"OldField = {OldField}")
    underRelaxed_Field = underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField
    print(f"Under relaxed Field = {underRelaxed_Field}")
    return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

def convergence_TH(Field, OldField, tol):
    conv=False
    print("Checking TH convergence")
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    print(f"TH rms error = {rms_error}")
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i]) > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_POW(Field, OldField, tol):
    """
    Tolerance in % of the nodal power, to allow for an easier comparison between different total powers.
    """
    conv=False
    print("Checking POW convergence")
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    print(f"POW rms error = {rms_error}")
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i])/OldField[i] > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_pcm(keffNew, keffOld, tol):
    print(f"keffNew = {keffNew}")
    print(f"keffOld = {keffOld}")
    print(f"keffNew - keffOld = {keffNew - keffOld}")
    conv=False
    if np.abs(keffNew - keffOld)*1e5 > tol:
        conv = False
    else:
        conv = True
    return conv


def compute_RMS(Field):
    return np.sqrt(np.sum(Field**2)/len(Field))

def get_sine_axial_p_form_over_mesh(h, Iz):
    """
    get sine values over a mesh of Iz points along the height h.
    Parameters:
    h (float): Total height of the cylinder [m].
    Iz (int): Number of axial subdivisions.
    """
    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, Iz + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    # Compute sine values at each midpoint
    sine_values = np.sin(np.pi * z_values / h) / np.mean(np.sin(np.pi * z_values / h))  # Normalize to have a mean of 1
    print(sine_values)
    return z_values, sine_values

def get_cosine_axial_p_form_over_mesh(h, Iz):
    """
    get cosine values over a mesh of Iz points along the height h.
    Parameters:
    h (float): Total height of the cylinder [m].
    Iz (int): Number of axial subdivisions.
    """
    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, Iz + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    # Compute cosine values at each midpoint
    cosine_values = np.cos(np.pi * z_values / (2 * h)) / np.mean(np.cos(np.pi * z_values / (2 * h)))  # Normalize to have a mean of 1
    print(cosine_values)
    return z_values, cosine_values

def get_z_mesh_midpoint(h, nz):
    """
    Returns a list of z mesh points based on the height and number of axial mesh points.
    """
    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, nz + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    return z_values

def compute_difference_fields(field):
    """
    Compute the difference between the last two fields of the list field
    field : list : list of the fields to compare
    """
    print(field)
    diff = np.abs(field[-1] - field[-2])
    print(f"Difference between the last two fields = {diff}")
    return diff

def compute_residuals(field):
    """
    Compute the residuals of the field
    field : list : list of the fields to compute the residuals in %
    """
    residuals = (field[-1] - field[-2])*100/field[-2]
    print(f"Residuals of the field = {residuals} %")
    return residuals

######## End helper functions declaration ##########


##
########## User input ##########
########## Algorithm parameters ###########
start_time = time.time()
solveConduction = True
nIter = 1000 #1000
tol_TH_Tf = 1e-3 # Convergence criterion for TFuel, in K
tol_TH_Tc = 1e-3 # Convergence criterion for TCool, in K
tol_TH_Dc = 1e-3 # Convergence criterion for DCool, in kg/m3
tol_POW = 0.01 # Convergence criterion for the power axial distribution, in % of nodal power
tol_Keff =  0.1 # Convergence criterion for the Keff value in pcm

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
Iz1 = 20 # 
#z_nodes_to_test = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200] # list of the axial nodes to test for the meshing sensitivity analysis

#normalization_powers_to_test = [5e3, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3, 40e3] # list of the powers to test for the normalization power sensitivity analysis
Pnorm = 20e3 # W, total fission power in the fuel rod, 40kW, 35kW, 30kW, 25kW, 20kW
# Test powers 40kW, 35kW, 30kW, 25kW, 20kW.
# Suspected limitation of pressure drop correlation at high powers in GeN-Foam comparison (21/12/2024)

power_guess_type = "sine" # "flat", "sine", "cosine"
fraction_pow_fuel = 1.0 # Fraction of the total power that is deposited in the fuel, default value in DONJON is 0.97

########## Choice of Thermalhydraulics correlation ##########
voidFractionCorrel = 'EPRIvoidModel' # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
frfaccorel = "Churchill" # 'base', 'blasius', 'Churchill', 
P2Pcorel = "lockhartMartinelli" # 'base', 'HEM1', 'HEM2', 'MNmodel', 'lockhartMartinelli'
numericalMethod = "BiCG" # "FVM": Solves the system using matrix inversion with preconditioning.
                         # "GaussSiedel" : Applies the Gauss-Seidel iterative solver.
                         # "BiCG" : Uses the BiConjugate Gradient method for solving non-symmetric or indefinite matrices.
                         # "BiCGStab" : Applies the BiCGStab (BiConjugate Gradient Stabilized) method to ensure faster and more stable convergence.

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square" # "square", "cylindrical"
pitch = 1.295e-2 # m : ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 3.8 # m : height : 3.8 m : active core height in BWRX-300 SMR, 1.555 m : for GeNFoam comparison.


## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
tInlet = 270 + 273.15 # K, for BWRX-300 SMR core
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1530 kg/s
massFlowRate = 1530  / (200*91)  # kg/s

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
kClad = 21.5 #21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity
    

############ Nuclear Parameters ###########
# Number of fuel rods and assemblies for a small modular Boiling Water Reactor core

n_rods = 91 # This value is for the ATRIUM-10 fuel assembly (10*10-9), in a GNF2 fuel assembly (BWRX-300) there are 92 fuel rods per assembly
n_assmblies = 240 # This value is for the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"
full_core_power = 870e6 # W, full core thermal power of the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"

volumic_mass_U = 10970 # kg/m3 'Thermophysical Properties of MOX and UO2 Fuels Including the Effects of Irradiation' - ORNL/TM-2000/351, Popov, Carbajo, Ivanov, Yoder. November 2000.

## Fuel rod scale parameters :
Fuel_volume = np.pi*fuelRadius**2*height # m3

#Fuel_rod_volume = np.pi*cladRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_U*1000 # g
print(f"Fuel mass = {Fuel_mass} g")


Bundle_volume = Fuel_volume / Iz1 # m3, Bundle <=> 1 axial slice of the fuel channel

compo_name = "_COMPO_AT10_MPHYS_PT_addeddens" # Name of the COMPO object to be used in the neutronics solution

######## Creation of results directory ##########
path=os.getcwd()
a=os.path.exists(f"multiPhysics_pyTHM/24UOX_cell")
if a==False:
    os.makedirs(f"multiPhysics_pyTHM/24UOX_cell")
print(path)

SAVE_DIR = f"multiPhysics_pyTHM/24UOX_cell/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/"

a=os.path.exists(SAVE_DIR)
if a==False:
    os.makedirs(SAVE_DIR)
########## End creation of results directory ##########

    
print(f"$$ - BEGIN Coupled case initialization : Iz1 = {Iz1}, total normalization power = {Pnorm} W")

########## Fields to store TH quantities for post treatment and analysis ##########
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
axial_p_form_distribs = []
relaxed_axial_p_form_distribs = [] # list to store the relaxed axial power shape at each iteration

Residuals_Power_Distribs = []
Residuals_Relaxed_Power_Distribs = []

##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

# 1.) Guess the axial power shape : used to initialize the TH solution
if power_guess_type == "flat":
    # Flat power density guess
    z_mesh = get_z_mesh_midpoint(height, Iz1)
    axial_p_form = np.ones(Iz1)
elif power_guess_type == "sine":
    # Sine power density guess
    z_mesh, axial_p_form = get_sine_axial_p_form_over_mesh(height, Iz1)  # Sine distribution, normalized to have a mean of 1
    z_mesh, axial_p_form_cos = get_cosine_axial_p_form_over_mesh(height, Iz1)  # Sine distribution, normalized to have a mean of 1
    print(f"z_mesh = {z_mesh}")
    print(f"axial_p_form = {axial_p_form}")
    print(f"axial_p_form_cos = {axial_p_form_cos}")

elif power_guess_type == "cosine":
    # Cosine power density guess
    z_mesh, axial_p_form = get_cosine_axial_p_form_over_mesh(height, Iz1)  # Cosine distribution, normalized to have a mean of 1



# 2.) TH solution for initial guess of power shape : sine shape (if devide by Iz : doesnt work and the initial TH solution isn't in IAPWS domain)
# This can probably be fixed by giving a different initial guess for the total fission power.
## 2.1) Initial thermal hydraulic resolution
THComponentIni = THM_prototype("Initialization of BWR Pincell equivalent channel", canalType, pitch, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, Pnorm, axial_p_form, fraction_pow_fuel, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)

TeffIni, TwaterIni, rhoIni, voidFracIni, Pini, Uini, Hini = THComponentIni.get_TH_parameters()

THM_ini_time = (time.time() - start_time)
print(f"THM_ini_time = {THM_ini_time} s = total elapsed time = {THM_ini_time} s")


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
THData = lcm.new('LCM','THData')
THData['TFuelList']    = np.array(TeffIni, dtype='f')
THData['TCoolList'] = np.array(TwaterIni, dtype='f')
THData['DCoolList'] = np.array(rhoIni/1000, dtype='f')
THData.close() # close without erasing

## 3.) Initializing Neutronics solution
# 3.1) construct the Lifo stack for IniDONJON
ipLifo1=lifo.new()
ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
ipLifo1.pushEmpty("Cpo", "LCM") # Compo
ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
ipLifo1.push(compo_name) # Compo name
ipLifo1.push(int(Iz1)) # Number of axial subdivisions
ipLifo1.push(Fuel_mass) # Mass of the fuel
ipLifo1.push(height) # Height of the fuel rod
ipLifo1.push(pitch) # pitch of the fuel assembly/cell

# 3.2) call IniDONJON Cle-2000 procedure
IniDONJON = cle2000.new('IniDONJON',ipLifo1,1)
IniDONJON.exec()
print("IniDONJON execution completed")

# recover the output LCM objects
Fmap = ipLifo1.node("Fmap")
Matex = ipLifo1.node("Matex")
Cpo = ipLifo1.node("Cpo")
Track = ipLifo1.node("Track")
# State vector and number of parameters can be used to modify the Fmap object directly ?
stateVector = Fmap["STATE-VECTOR"]
mylength = stateVector[0]*stateVector[1]
npar = stateVector[7]

print("Recovered stateVector: ", stateVector)
print("Number of parameters: ", npar)

# empty the Lifo stack for IniDONJON
while ipLifo1.getMax() > 0:
    ipLifo1.pop();

# 4.) Iterative scheme for coupled neutronics and thermalhydraulics solution
# Attempting to run multiPhysics scheme.
# Everytime, solve neutron transport, update axial power shape, solve TH, update TH data, and repeat.
powi = Pnorm/1e6 # Reference at 0.1722 MW from AT10_24UOX test ?

# check powi == Pnorm
print(f"powi = {powi} MW and Pnorm = {Pnorm} W")

# Create Lifo stack for Neutronics solution
ipLifo2 = lifo.new()
Neutronics = cle2000.new('solNeutro',ipLifo2,1)

current_time0 = time.time()
DONJON_ini_time = (current_time0 - THM_ini_time)
total_elapsed_time = (current_time0 - start_time)
print(f"DONJON initialization time = {DONJON_ini_time} s")
print(f"Total elapsed time = {total_elapsed_time} s")

## Multi-Physics solution by Picard iterations
iter = 0
conv = False
while not conv:
    current_time1 = time.time()
    convKeff = False
    convTH = False
    convPow = False
    iter+=1
    # 4.1) Neutronics solution for TH parameters obtained with initial TH solution

    ################## Neutronics part ##################
    # fill the Lifo stack for Neutronics solution
    print(f"$$ - BEGIN iter = {iter}")
    ipLifo2.push(Fmap)
    ipLifo2.push(Matex)
    if iter == 1: # At the first iteration, create empty LCM objects to host the Flux and Power fields
        print("in iter = 1")
        Flux = ipLifo2.pushEmpty("Flux", "LCM")
    else:
        ipLifo2.push(Flux)
        print("Flux and Power at iteration > 1")
    Power = ipLifo2.pushEmpty("Power", "LCM")
    ipLifo2.push(Cpo) # Push COMPO object
    ipLifo2.push(Track) # Push Tracking object for FEM solution : obtained from IniDONJON
    ipLifo2.push(THData)
    ipLifo2.push(iter)
    ipLifo2.push(powi)
    ipLifo2.push(int(Iz1)) 

    # 4.2) Call Neutronics component :
    Neutronics.exec()
    print("solNeutro.c2m execution completed")
    Flux = ipLifo2.node("Flux") # Recover the Flux field
    RecoveredPower = ipLifo2.node("Power") # Recover the Power field (LCM object)
    # 4.2.1) Recover the Keff value
    Keff = Flux["K-EFFECTIVE"][0] 
    Keffs.append(Keff)
    print(f"At iter {iter} : Keff = {Keff}")
    # 4.2.2) Recover the Power axial distribution obtained from neutronics calculation --> used to update the axial power shape
    PowerDistribution = RecoveredPower["POWER-DISTR"] 
    updated_axial_p_form = PowerDistribution / np.mean(PowerDistribution) # This is the axial power shape obtained from the neutronics calculation
    FLUX = RecoveredPower["FLUX"]
    current_time2 = time.time()
    DONJON_iter_time = (current_time2 - current_time1)
    total_elapsed_time = (current_time2 - start_time)
    Power_Distribs.append(PowerDistribution) # Store the axial power shape at each iteration
    # 4.3) Update the axial power shape to be used in the TH solution

    # 4.3.bis) Under relaxation of the Power distribution for the next iteration
    if iter == 1 and relax_Pow:
        axial_p_form_relaxed = underRelaxation(updated_axial_p_form, axial_p_form, Pow_underRelaxationFactor)
        # Under relaxation of the Power distribution?
    elif iter > 1 and relax_Pow:
        axial_p_form_relaxed = underRelaxation(updated_axial_p_form, relaxed_axial_p_form_distribs[-1], Pow_underRelaxationFactor)
    if relax_Pow:
        relaxed_axial_p_form_distribs.append(axial_p_form_relaxed)
        print(f"Relaxed axial power shape : {axial_p_form_relaxed}")

    if iter > 1:
        Residuals_Power_Distribs.append(compute_residuals(Power_Distribs))
        RMS_power_residual = compute_RMS(Residuals_Power_Distribs[-1])

    # 4.4) Erase the THData object to store the updated TH data
    THData.erase()

    # 5.) TH procedure for updated power shape
    ############# Thermalhydraulic part ##############

    current_time3 = time.time()
    time_spent_in_plotting_res_under_relax = (current_time3 - current_time2)
    total_elapsed_time = (current_time3 - start_time)
    print(f"Total elapsed time = {total_elapsed_time} s")
    # 5.1) TH resolution with updated power shape :
    print(f"At iter {iter} : THM resolution with updated power shape")
    THMComponent = THM_prototype("BWR Pincell equivalent channel", canalType, pitch, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, Pnorm, updated_axial_p_form, fraction_pow_fuel, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)    ##### qFiss updated via axial_p_form

    # 5.2) recover the new TH data
    TeffTEMP, TwaterTEMP, rhoTEMP, voidFracTEMP, P_TEMP, U_TEMP, H_TEMP = THMComponent.get_TH_parameters()
    current_time4 = time.time()
    time_spent_in_THM = (current_time4 - current_time3)
    total_elapsed_time = (current_time4 - start_time)
    print(f"Time spent in THM iter {iter} = {time_spent_in_THM} s")
    print(f"TeffTEMP = {TeffTEMP}, TwaterTEMP = {TwaterTEMP}, rhoTEMP = {rhoTEMP}, voidFracTEMP = {voidFracTEMP}, P_TEMP = {P_TEMP}, U_TEMP = {U_TEMP}, H_TEMP = {H_TEMP}")
    print(f"Total elapsed time = {total_elapsed_time} s")
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
    
    # 5.3) Update the TH data to be used in the next neutronics solution
    
    THData = lcm.new('LCM','THData') # Create new LCM object to store the updated TH data
    THData['TFuelList']    = np.array(TeffTEMP, dtype='f') # Store the updated TH data in the LCM object
    THData['TCoolList'] = np.array(TwaterTEMP, dtype='f') # Store the updated TH data in the LCM object
    THData['DCoolList'] = np.array(rhoTEMP/1000, dtype='f') # Convert density from kg/m3 to g/cm3
    THData.close() # close without erasing

    # 5.4) Under relaxation of TH fields for the next iteration : check if it is needed ?
    if relax_TH:
        TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], TH_underRelaxationFactor)
        Twater[-1] = underRelaxation(Twater[-1], Twater[-2], TH_underRelaxationFactor)
        rho[-1] = underRelaxation(rho[-1], rho[-2], TH_underRelaxationFactor)


    # 6.) Empty the ipLifo2 Lifo stack to prepare for the next iteration
    while ipLifo2.getMax() > 0: 
        ipLifo2.pop()
    
    print(f"$$ - END iter = {iter}")
    # 7.) Check for convergence on Keff value : add more converegence criteria for other fields ?
    if iter>1 and convergence_pcm(Keffs[-1], Keffs[-2], tol_Keff):
        print("Convergence on Keff reached after ", iter, " iterations")
        convKeff = True
    if iter>1 and convergence_TH(rho[-1], rho[-2], tol_TH_Dc) and convergence_TH(TeffFuel[-1], TeffFuel[-2], tol_TH_Tf) and convergence_TH(Twater[-1], Twater[-2], tol_TH_Tc):
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
        export_results_to_csv(Iz1, Pnorm, power_guess_type, pitch, height, results, iter, SAVE_DIR)
    if iter == nIter-1:
        total_iter = iter
        print("Convergence not reached after ", iter, " iterations")
        results = {"z_mesh": z_mesh,
                    "TFUEL": TeffFuel[-1], 
                    "TCOOL": Twater[-1], 
                    "DCOOL": rho[-1], 
                    "EPS": voidFraction, # keep a list of void fraction for all iterations, to analyse in case of convergence failure
                    "PRESSURE": P_[-1], 
                    "VELOCITY": U_[-1], 
                    "ENTHALPY": H_[-1], 
                    "POWER": Power_Distribs, # keep a list of power distribution for all iterations, to analyse in case of convergence failure
                    "KEFF": Keffs} # keep a list of Keff for all iterations, to analyse in case of convergence failure
        export_to_unconverged(Iz1, Pnorm, power_guess_type, pitch, height, results, iter, SAVE_DIR)

    current_time5 = time.time()
    time_spent_in_res_and_conv = (current_time5 - current_time4)
    total_elapsed_time = (current_time5 - start_time)
    print(f"Time spent in residuals and convergence check = {time_spent_in_res_and_conv} s")
    print(f"Total elapsed time = {total_elapsed_time} s")
    total_iteration_time = (current_time5 - current_time1)
    print(f"Total iteration time = {total_iteration_time} s for iteration {iter}")

    print("$$$ - multiPhysics_pyTHM : END ITERATIONS - $$$")

# End of the iterative scheme for coupled neutronics and thermalhydraulics solution
print(f"$$ - END Coupled case : Iz1 = {Iz1}, total normalization power = {Pnorm} W")
print("$$$ - multiPhysics_pyTHM : END OF SCRIPT - $$$")

