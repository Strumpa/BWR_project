## To cleanup main script and make it more readable
# write export results function
# Author : R. Guasch
# Date : 13/06/2025
# Project : BWR-multiphysics / PolyBWR 
# Set of function to format and export results to csv data frames

# Import necessary libraries
import pandas as pd
import numpy as np
import os 
import sys
import shutil

def export_results_to_csv(nz, Pnorm, power_guess_type, pitch, height, results, iterations, output_dir):
    """
    Recover information about the physical case modelled and export converged fields to a csv file
    Parameters:
    nz (int): Number of axial nodes in the mesh
    Pnorm (float): Normalized power level
    power_guess_type (str): Type of power guess used
    pitch (float): Pitch of the fuel channel in meters
    height (float): Height of the fuel channel in meters
    results (dict): Dictionary containing the results of the simulation
    iterations (int): Number of iterations performed to reach convergence
    output_dir (str): Directory where the results will be saved

    """
    # Unpack results
    z_mesh = results['z_mesh']  # Axial mesh points
    TFUEL = results['TFUEL']  # Fuel temperature
    TCOOL = results['TCOOL']  # Coolant temperature
    DCOOL = results['DCOOL']  # Coolant density
    EPS = results['EPS']      # Void fraction
    PCOOL = results['PRESSURE']  # Coolant pressure
    VCOOL = results['VELOCITY']  # Coolant velocity
    HCOOL = results['ENTHALPY']  # Coolant enthalpy
    Power_distrs = results['POWER']  # Power distribution
    Keffs = results['KEFF']  # Keff values
    
    vol_node = pitch**2 * height / nz  # m3, volume of each node
    df_export = pd.DataFrame({
        'Z Mesh': z_mesh,
        'TFUEL (K)': np.array(TFUEL),
        'TCOOL (K)': np.array(TCOOL),
        'DCOOL (Kg/m3)': np.array(DCOOL),
        'VOID': np.array(EPS),
        'PCOOL (Pa)': np.array(PCOOL),
        'UCOOL (m/s)': np.array(VCOOL), # add coolant velocity
        'HCOOL (J/Kg)': np.array(HCOOL), # add coolant enthalpy
        'POWER (W)': np.array(Power_distrs),  # Power distribution along the fuel channel
        'Q (W/m3)': np.array(Power_distrs) / vol_node,  # Power density distribution
    })
    df_export['Keff'] = Keffs[-1]  # Add the final Keff value
    df_export['ITERATION'] = iterations  # Add the number of iterations required to reach convergence
    # Save the DataFrame to a CSV file
    df_export.to_csv(f"{output_dir}/df_mphys_pyTHM_{nz}_h{height*100:.0f}_{power_guess_type}_pow_guess_{int(Pnorm)}W.csv", index=False)
    return

def export_to_unconverged(nz, Pnorm, power_guess_type, pitch, height, results, iterations, output_dir):
    """
    special treatment of exports in the case on an unconverged simulation
    Parameters:
    nz (int): Number of axial nodes in the mesh
    Pnorm (float): Normalized power level
    power_guess_type (str): Type of power guess used
    pitch (float): Pitch of the fuel channel in meters
    height (float): Height of the fuel channel in meters
    results (dict): Dictionary containing the results of the simulation
    iterations (int): Number of iterations performed to reach convergence
    output_dir (str): Directory where the results will be saved
    
    """
    # Unpack results
    z_mesh = results['z_mesh']  # Axial mesh points
    TFUEL = results['TFUEL']  # Fuel temperature
    TCOOL = results['TCOOL']  # Coolant temperature
    DCOOL = results['DCOOL']  # Coolant density
    EPS = results['EPS']      # Void fraction
    PCOOL = results['PRESSURE']  # Coolant pressure
    VCOOL = results['VELOCITY']  # Coolant velocity
    HCOOL = results['ENTHALPY']  # Coolant enthalpy
    Power_distrs = results['POWER']  # Power distribution
    Keffs = results['KEFF']  # Keff values
    
    vol_node = pitch**2 * height / nz  # m3, volume of each node
    df_export = pd.DataFrame({
        'Z Mesh': z_mesh,
        'TCOMB (K)': np.array(TFUEL[-1]),
        'TCOOL (K)': np.array(TCOOL[-1]),
        'DCOOL (Kg/m3)': np.array(DCOOL[-1]),
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
    # Check how to handle keffs in unconverged cases
    # there will be nIter keffs : 1 per iteration
    # df size is nz entries, so we can only add the last nz keffs if nz < nIter else, need to add at max nIter keffs
    if len(Keffs) < nz:
        Keffs = np.concatenate((Keffs, np.full(nz - len(Keffs), np.nan)))  # Pad with NaNs if necessary
    elif len(Keffs) > nz:
        Keffs = Keffs[-nz:]  # Take the last nz values if there are more than nz keffs
    df_export['Keff'] = Keffs  # Add the last nz/nIter values of Keff to the DataFrame : look at convergence rate on the keff value
    df_export['ITERATION'] = iterations  # Add the number of iterations

    # Save the DataFrame to a CSV file
    if not os.path.exists(f"{output_dir}/unconverged"):
        os.makedirs(f"{output_dir}/unconverged")  # Create the unconverged directory if it doesn't exist
    df_export.to_csv(f"{output_dir}/unconverged/df_mphys_pyTHM_{nz}_h{height*100:.0f}_{power_guess_type}_pow_guess_{int(Pnorm)}W_unconverged.csv", index=False)

      
