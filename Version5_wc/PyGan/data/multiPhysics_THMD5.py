## PyGan script to couple DONJON procedures with THM: solution from DJON5 module
# Date : 29/05/2025
# Author : Raphaël GUASCH
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR pincell
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simpson
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *
from D5THM_driver import initializeModel
from D5THM_driver import solNeutronicsDONJON
from D5THM_driver import solTHM
from generate_c2m import generate_c2m_files

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
        if np.abs(Field[i] - OldField[i])*100/OldField[i] > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_POW_RMS(Field, OldField, tol):
    """
    Tolerance in % of the nodal power, to allow for an easier comparison between different total powers.
    """
    rms_error = compute_RMS((Field - OldField)*100/OldField)
    print(f"POW RMS error = {rms_error}")
    if rms_error > tol:
        return False
    else:
        return True

def convegence_TH_RMS(Field, OldField, tol):
    """
    Tolerance in physical units, to allow for an easier physical interpretations.
    """
    rms_error = compute_RMS(Field - OldField)
    print(f"TH RMS error = {rms_error}")
    if rms_error > tol:
        return False
    else:
        return True

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

    return z_values, cosine_values

def get_z_mesh_midpoint(h, nz):
    """
    Returns a list of z mesh points based on the height and number of axial mesh points.
    """
    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, nz + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    return z_values

######## End helper functions declaration ##########


##
########## User input ##########
########## Algorithm parameters ###########
nIter = 5000 #1000
tol_TH_Tf = 1e-3 # Convergence criterion for TFuel, in K
tol_TH_Tc = 1e-3 # Convergence criterion for TCool, in K
tol_TH_Dc = 1e-3 # Convergence criterion for DCool, in kg/m3
tol_POW = 0.01 # Convergence criterion for the power axial distribution, in % of nodal power
tol_Keff =  0.01 # Convergence criterion for the Keff value in pcm


########## Mesh parameters ###########

# Sensitivity to the meshing parameters
# z_nodes_to_test = list of axial nodes to test for spatial convergence 
z_nodes_to_test = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200] # 
z_nodes_to_test = [10,] #
# Testing mesh refinements all the way to 200 nodes is possible in this version of DONJON, in the official version, it is limted to 100.

powers_to_test = [5.0e3, 10.0e3, 15e3, 20.0e3, 25.0e3, 30.0e3, 35.0e3, 40.0e3] # 
powers_to_test = [5e3] # W, total fission power in the fuel rod,
#powers_to_correlations = {5.0e3:["MONO", "BOWR"], 10e3:["MONO", "BOWR", "SAHA"], 20e3:["SAHA", "BOWR"], 25e3:["SAHA"], 30e3:["SAHA"], 35e3:["SAHA"], 40e3:["SAHA"]} # W, total fission power in the fuel rod, used for the correlations
# Issues with sub-cooled boiling correlations, test MONO --> should not be using subcooled boiling model and simply compute PDROP + DFM if properly implementd.
# key : total fission power in the fuel rod, value : list of sub-cooled boiling correlations to be tested for this specific power.
powers_to_correlations = {5.0e3:["SAHA", "BOWR"], 10e3:["SAHA", "BOWR"], 15e3: ["SAHA"], 20e3:["SAHA"], 25e3:["SAHA"], 30e3:["SAHA"], 35e3:["SAHA"], 40e3:["SAHA"]}
THM_under_relaxation_factor = 1.0 # Under relaxation factor for the THM solution, to be used in the underRelaxation function

create_cases = False # True to create the c2m files for the different z_nodes
exec_cases = True # True to execute the cases with the created c2m files

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square" # "square", "cylindrical"
pitch = 1.295e-2 # m : ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 3.8 # m : height : 3.8 m : active core height in BWRX-300 SMR

area = pitch**2 - (np.pi * cladRadius**2) # m2 : cross section area of the fuel assembly
power_guess_type = "sine" # "cosine", "flat" or "sine" for the initial guess of the power shape

## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
tInlet = 270 + 273.15 # K, for BWRX-300 SMR core
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1530 kg/s
massFlowRate = 1530  / (200*91)  # kg/s


############ Nuclear Parameters ###########
# Choice of COMPO object tabulating XS for neutronics simulation

compo_name = "_COMPO_AT10_MPHYS_PT_addeddens" # Name of the COMPO object to be used in the neutronics solution
# Number of fuel rods and assemblies for a small modular Boiling Water Reactor core

volumic_mass_U = 10970 # kg/m3 'Thermophysical Properties of MOX and UO2 Fuels Including the Effects of Irradiation' - ORNL/TM-2000/351, Popov, Carbajo, Ivanov, Yoder. November 2000.

## Fuel rod scale parameters :
Fuel_volume = np.pi*fuelRadius**2*height # m3

#Fuel_rod_volume = np.pi*cladRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_U*1000 # g

######## Creation of results directory ##########

SAVE_DIR = f"multiPhysics_THMD5/multi_physics_cell/test_interp_THMDRV/test_keff_nodal_pow_conv_tighter"

a=os.path.exists(SAVE_DIR)
if a==False:
    os.makedirs(SAVE_DIR)
########## End creation of results directory ##########



### Automating case creation with "generate_c2m.py"
if create_cases:
    for nz in z_nodes_to_test:
        generate_c2m_files(nz, nz_max=200, output_dir="c2m_files", exec_case = False)


if exec_cases:
    for nz in z_nodes_to_test:
        ##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

        # 1.) Guess the axial power shape : used to initialize the TH solution
        if power_guess_type == "flat":
            # Flat power density guess
            z_mesh = get_z_mesh_midpoint(height, nz)
            axial_p_form = np.ones(nz)
        elif power_guess_type == "sine":
            # Sine power density guess
            z_mesh, axial_p_form = get_sine_axial_p_form_over_mesh(height, nz)  # Sine distribution, normalized to have a mean of 1
        elif power_guess_type == "cosine":
            # Cosine power density guess
            z_mesh, axial_p_form = get_cosine_axial_p_form_over_mesh(height, nz)  # Cosine distribution, normalized to have a mean of 1

        for Pnorm in powers_to_test:
            for THM_sub_cooled_correlation in powers_to_correlations[Pnorm]:
                print(f"$$$ - multiPhysics.py : BEGIN EXECUTION FOR nz = {nz}, Pnorm = {Pnorm} W, THM_sub_cooled_correlation = {THM_sub_cooled_correlation} - $$$")
                guessPowForm = lcm.new('LCM','guessPow')
                guessPowForm['pform'] = np.array(axial_p_form, dtype='f')
                guessPowForm.close()
                print(f"Begin test with generated c2m : nz = {nz}, Pnorm = {Pnorm} W")
                # 2.) TH solution for initial guess of power shape : sine shape (if divide by Iz : doesnt work and the initial TH solution isn't in IAPWS domain)
                # This can probably be fixed by giving a different initial guess for the total fission power.
                ## 2.1) Initial thermal hydraulic resolution
                # Call THM: through PyGan
                powi = Pnorm*1e-6 # W to MW
                # Run coupled calculation with generated c2m files
                Fmap, Matex, Cpo, Track, Thm = initializeModel(compo_name, nz, Fuel_mass, height, pitch, powi, guessPowForm, pOutlet, tInlet, massFlowRate, area, fuelRadius, gapRadius, cladRadius, THM_sub_cooled_correlation, THM_under_relaxation_factor)

                stationary_info = Thm["HISTORY-DATA"]["TIMESTEP0000"]["CHANNEL"][0]
                stationary_info.lib();

                # 5.2) recover the usefull TH data
                TFUEL = [np.array(Fmap['PARAM'][0]['P-VALUE'], dtype='f')]
                TCOOL = [np.array(Fmap['PARAM'][1]['P-VALUE'], dtype='f')]
                DCOOL = [np.array(Fmap['PARAM'][2]['P-VALUE'], dtype='f')]
                VCOOL = [np.array(stationary_info["VELOCITIES"], dtype='f')]
                VLIQ  = [np.array(stationary_info["V-LIQ"], dtype='f')]
                VVAP  = [np.array(stationary_info["V-VAP"], dtype='f')]
                PCOOL = [np.array(stationary_info["PRESSURE"], dtype='f')]
                EPS = [np.array(stationary_info["EPSILON"], dtype='f')]
                HCOOL = [np.array(stationary_info["ENTHALPY"], dtype='f')] # Coolant enthalpy
                # 3.) Neutronics solution
                Keffs = []
                Power_distrs = []
                
                iter=1
                conv = False
                while not conv:
                    print(f"$$$ - multiPhysics.py : BEGIN ITERATION {iter} - $$$")
                    if iter==1:    
                        # First iteration, no previous solution
                        # Call solNeutronicsDONJON with empty Flux and Power
                        print("First iteration, no previous solution")
                        Power, Flux, Fmap, Matex, Track = solNeutronicsDONJON(iter, powi, nz, Fmap, Matex, Cpo, Track, Flux=None)
                        ## Fix time stuff later
                    else:
                        Fmap, Thm = solTHM(iter, Fmap, Thm, pOutlet, tInlet, massFlowRate, area, fuelRadius, gapRadius, cladRadius, THM_sub_cooled_correlation, THM_under_relaxation_factor)
                        stationary_info = Thm["HISTORY-DATA"]["TIMESTEP0000"]["CHANNEL"][0]
                        stationary_info.lib();
                        EPS.append(stationary_info["EPSILON"])
                        PCOOL.append(stationary_info["PRESSURE"])
                        VCOOL.append(stationary_info["VELOCITIES"])
                        VLIQ.append(stationary_info["V-LIQ"])
                        VVAP.append(stationary_info["V-VAP"])

                        ## Fix time stuff later
                        
                        Tfuel_iter = np.array(Fmap['PARAM'][0]['P-VALUE'], dtype='f')
                        Tcool_iter = np.array(Fmap['PARAM'][1]['P-VALUE'], dtype='f')
                        Dcool_iter = np.array(Fmap['PARAM'][2]['P-VALUE'], dtype='f')

                        # Store underrelaxed values
                        TFUEL.append(Tfuel_iter)
                        TCOOL.append(Tcool_iter)
                        DCOOL.append(Dcool_iter)

                        # Previous solution available
                        # Call solNeutronicsDONJON with previous Flux and Power and modified Fmap and Matex 
                        Power, Flux, Fmap, Matex, Track = solNeutronicsDONJON(iter, powi, nz, Fmap, Matex, Cpo, Track, Flux)
                        ## Fix time stuff later
                    #print(f"Neutronics_delta_time = {Neutronics_delta_time} s = total elapsed time = {total_time} s")

                    Keff = Flux["K-EFFECTIVE"][0] 
                    Power_distrib = Power["POWER-DISTR"] 
                    Flux_value = Power["FLUX"]
                    Keffs.append(Keff)
                    Power_distrs.append(Power_distrib)
                    
                    if iter>1:
                        #if convergence_pcm(keffNew=Keff, keffOld=Keffs[-2], tol=tol_Keff) and convergence_POW(Field=Power_distrib, OldField=Power_distrs[-2], tol=tol_POW) and convergence_TH(Field=TFUEL[-1], OldField=TFUEL[-2], tol=tol_TH_Tf) and convergence_TH(Field=TCOOL[-1], OldField=TCOOL[-2], tol=tol_TH_Tc) and convergence_TH(Field=DCOOL[-1], OldField=DCOOL[-2], tol=tol_TH_Dc):
                        #if convergence_pcm(keffNew=Keff, keffOld=Keffs[-2], tol=tol_Keff) and convergence_POW_RMS(Field=Power_distrib, OldField=Power_distrs[-2], tol=tol_POW) and convergence_POW(Field=Power_distrib, OldField=Power_distrs[-2], tol=tol_POW):
                        #if convergence_pcm(keffNew=Keff, keffOld=Keffs[-2], tol=tol_Keff) and convergence_POW_RMS(Field=Power_distrib, OldField=Power_distrs[-2], tol=tol_POW):
                        if convergence_pcm(keffNew=Keff, keffOld=Keffs[-2], tol=tol_Keff) and convergence_POW(Field=Power_distrib, OldField=Power_distrs[-2], tol=tol_POW):
                            conv = True
                            print("Convergence reached")
                    iter=iter+1     
                    if iter>nIter:
                        print("Max iterations reached")
                        # Export data from storing unconverged results for analysis
                        vol_node = pitch**2 * height / nz  # m3, volume of each node
                        df_export = pd.DataFrame({
                            'Z Mesh': z_mesh,
                            'TCOMB (K)': np.array(TFUEL[-1]),
                            'TCOOL (K)': np.array(TCOOL[-1]),
                            'DCOOL (Kg/m3)': np.array(DCOOL[-1])*1000, # convert to kg/m3
                            'PCOOL (Pa)': np.array(PCOOL[-1]),
                            'UCOOL (m/s)': np.array(VCOOL[-1]), # add coolant velocity
                            'HCOOL (J/Kg)': np.array(HCOOL[-1]), # add coolant enthalpy
                            'POWER (W)': np.array(Power_distrs[-1]),  # Power distribution along the fuel channel
                            # Store the last 5 power distributions for analysis
                            'Q (W/m3)_last': np.array(Power_distrs[-1]) / vol_node,  # Power density distribution
                            'Q (W/m3)-2': np.array(Power_distrs[-2]) / vol_node,  # Power density distribution
                            'Q (W/m3)-3': np.array(Power_distrs[-3]) / vol_node,  # Power density distribution
                            'Q (W/m3)-4': np.array(Power_distrs[-4]) / vol_node,  # Power density distribution
                            'Q (W/m3)-5': np.array(Power_distrs[-5]) / vol_node,  # Power density distribution
                            # Store the last 5 epsilon values for analysis
                            'VOID_last': np.array(EPS[-1]),  # Void fraction distribution
                            'VOID-2': np.array(EPS[-2]),  # Void fraction distribution
                            'VOID-3': np.array(EPS[-3]),  # Void fraction distribution
                            'VOID-4': np.array(EPS[-4]),  # Void fraction distribution
                            'VOID-5': np.array(EPS[-5]),  # Void fraction distribution
                        })
                        df_export['Keff'] = Keffs[-nz:]  # Add the last nz values of Keff to the DataFrame : look at convergence rate on the keff value
                        df_export['ITERATION'] = iter  # Add the number of iterations
                        # Save the DataFrame to a CSV file
                        df_export.to_csv(f"{SAVE_DIR}/df_mphys_thm_{nz}_h{height*100:.0f}_{power_guess_type}_pow_guess_{int(Pnorm)}W_{THM_sub_cooled_correlation}_THM_relax_{THM_under_relaxation_factor}_unconverged.csv", index=False)
                        Fmap.erase()
                        Matex.erase()
                        break
                    
                    if conv:
                        print(f"Convergence reached after {iter} iterations")
                        print(f"Final Keff = {Keff}")
                        print(f"Final Power distribution = {Power_distrib}")
                        Fmap.erase()
                        Matex.erase()
                    

                        # 5.) Recover converged results and export them in a pandas data frame for cross validation with TPF and Serpent2
                        vol_node = pitch**2 * height / nz  # m3, volume of each node
                        df_export = pd.DataFrame({
                            'Z Mesh': z_mesh,
                            'TCOMB (K)': np.array(TFUEL[-1]),
                            'TCOOL (K)': np.array(TCOOL[-1]),
                            'DCOOL (Kg/m3)': np.array(DCOOL[-1])*1000, # convert to kg/m3
                            'VOID': np.array(EPS[-1]),
                            'PCOOL (Pa)': np.array(PCOOL[-1]),
                            'VCOOL (m/s)': np.array(VCOOL[-1]), # add coolant velocity
                            'VVAP (m/s)': np.array(VVAP[-1]), # add vapor velocity
                            'VLIQ (m/s)': np.array(VLIQ[-1]), # add liquid velocity
                            'HCOOL (J/Kg)': np.array(HCOOL[-1]), # add coolant enthalpy
                            'POWER (W)': np.array(Power_distrs[-1]),  # Power distribution along the fuel channel
                            'Q (W/m3)': np.array(Power_distrs[-1]) / vol_node,  # Power density distribution
                        })
                        df_export['Keff'] = Keffs[-1]  # Add the final Keff value
                        df_export['ITERATION'] = iter  # Add the number of iterations
                        # Save the DataFrame to a CSV file
                        df_export.to_csv(f"{SAVE_DIR}/df_mphys_thm_{nz}_h{height*100:.0f}_{power_guess_type}_pow_guess_{int(Pnorm)}W_{THM_sub_cooled_correlation}_THM_relax_{THM_under_relaxation_factor}.csv", index=False)
                print("$$$ - multiPhysics.py : END OF ITERATIONS - $$$")
    print("$$$ - multiPhysics.py : END OF EXECUTION - $$$")






