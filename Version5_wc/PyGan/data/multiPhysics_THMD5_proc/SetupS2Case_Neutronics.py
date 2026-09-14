## Python3 script to export material data to S2 and setup Serpent2 case
# The goal is to be able to setup a Serpent2 case simulating neutron transport
# Using the material distribution obtained from THM: - DONJON5 coupled Picard iterations
# Author : R. Guasch
# Date : 10/10/2024, updated 29/05/2025
# Purpose #1 : generate Serpent2 3D model of pincell with variable fuel and moderator properties
# Purpose #2 : Validation / Cross validation of multi-physics BWR models in the Version5 environment

import numpy as np
import os
import sys
import matplotlib.pyplot as plt


def DensToIsoDens(coolant_density):
    """
    density : density of the material in g/cm3
    isotopic_mass : mass of the isotope in g/mol

    purpose : switch from coolant density to isotopic density used in Serpent2
    """
    # Calculation of moderator data
    # AVOGADRO's number
    A = 6.022094E-1 # Normalizing by 10E-24 to obtain isotopic density in # / b*cm

    M_H2O = 15.9994 + 2.0*1.00794
    N_H2O = coolant_density*A/M_H2O
    N_O = N_H2O
    N_H = 2.0*N_H2O 
    return N_H, N_O


def create_geometry(number_axial_slices, power_scaling_factor, pitch, fuel_radius, gap_radius, clad_radius, height, isGd):
    """Create the geometry for the Serpent2 input file"""
    slice_thinkness = height/number_axial_slices
    x_max = pitch/2
    x_min = -pitch/2
    y_max = pitch/2
    y_min = -pitch/2
    pin_radii = computeSantamarinaradii(fuel_radius, gap_radius, clad_radius, isGd)
    with open(f'AT10_24UOX_MPHYS_h155_mesh{number_axial_slices}_{power_scaling_factor}_mc', 'a') as output_file:
        output_file.write("% --- 1. Pin definition \n")
        for i in range(number_axial_slices):
            if isGd==False:
                output_file.write(f"pin fuel_slice{i+1} \n UOx_A_{i+1} {pin_radii[0]:.4f} \n UOx_B_{i+1} {pin_radii[1]:.4f} \n UOx_C_{i+1}  {pin_radii[2]:.4f} \n UOx_D_{i+1} {pin_radii[3]:.4f} \n gap  {pin_radii[4]:.4f} \n clad {pin_radii[5]:.4f} \n H2O_{i+1} \n")
            else:
                output_file.write(f"pin fuel_slice{i+1} \n Gd_A_{i+1} {pin_radii[0]:.4f} \n Gd_B_{i+1} {pin_radii[1]:.4f} \n Gd_C_{i+1}  {pin_radii[2]:.4f} \n Gd_D_{i+1} {pin_radii[3]:.4f} \n Gd_E_{i+1} {pin_radii[4]:.4f} \n Gd_F_{i+1} {pin_radii[5]:.4f} \n gap {pin_radii[6]:.4f} \n clad {pin_radii[7]:.4f} \n H2O_{i+1} \n")
        output_file.write(f"% --- 2. Surfaces definition \n")
        for i in range(number_axial_slices):
    
            output_file.write(f"surf {i+1} cuboid {x_min} {x_max} {y_min} {y_max} {i*slice_thinkness:.3f} {i*slice_thinkness+slice_thinkness:.3f} \n")
        output_file.write(f"surf {number_axial_slices+1} cuboid {x_min} {x_max} {y_min} {y_max} 0 {height} \n")
        output_file.write(f"% --- 3. Cell definition \n")
        for i in range(number_axial_slices):
            output_file.write(f"cell {i+1} 0 fill fuel_slice{i+1} -{i+1} -{number_axial_slices+1} \n")
        output_file.write(f"cell {number_axial_slices+1} 0 outside {number_axial_slices+1} \n")
    return

def computeSantamarinaradii(outer_fuel_radius, gap_radius, clad_radius, isGd):
    # Helper to defined fuel region radii for fuel pins --> important for evolution calculations
    # UOX : Finer discretization close to the outer radii is important for rim effects : Pu formation
    # Gd2O3 : rim effect + Gd evolution = challenge for deterministic methods, need 6 sub regions.
    """
    r_out = float, fuel radius
    isGd = bool, adapting discretization to Gd pin --> 6 radial sub regions instead.
    A. Santamarina recommandations :
    volumes for UOX pins : 50%, 80%, 95% and 100%
    volumes for Gd2O3 pins : 20%, 40%, 60%, 80%, 95% and 100%
    """
    if isGd==False:
        pin_radii=[(0.5**0.5)*outer_fuel_radius, (0.8**0.5)*outer_fuel_radius, (0.95**0.5)*outer_fuel_radius, outer_fuel_radius, gap_radius, clad_radius]
        
    else :
        pin_radii=[(0.2**0.5)*outer_fuel_radius, (0.4**0.5)*outer_fuel_radius, (0.6**0.5)*outer_fuel_radius, 
                    (0.8**0.5)*outer_fuel_radius, (0.95**0.5)*outer_fuel_radius, outer_fuel_radius, gap_radius, clad_radius]
    return pin_radii

def create_material_volumes(number_axial_slices, power_scaling_factor, Teff_fuel, TWater, iso_dens_O, iso_dens_H, isGd, XS_lib):
    """
    Generate Sepent2 material volumes
    Teff_fuel : list of effective temperatures in fuel
    TWater : list of temperatures in water
    iso_dens_O : list of oxygen isotopic densities (in water) obtained from DensToIsoDens(rho)
    iso_dens_H : list of hydrogen isotopic densities (in water) obtained from DensToIsoDens(rho)
    isGd : bool, True if Gd is present in the fuel pin --> then radially subdivide the fuel pin in 6 regions, if False 4 regions.
    """
    # check if XS_lib is valid
    if XS_lib not in ["PyNjoy2016", "oldlib"]:
        print("Error: invalid XS library selected")
        sys.exit(1)

    if isGd == False:
            fuel_names = ["UOx_A", "UOx_B", "UOx_C", "UOx_D"]
    else:
        fuel_names = ["Gd_A", "Gd_B", "Gd_C", "Gd_D", "Gd_E", "Gd_F"]

    if XS_lib == "PyNjoy2016":
        clad_lib_id = "05c"
        gap_lib_id = "05c"
        temp_tabulation_points = [293.0, 550.0, 900.0, 1200.0, 1999.0]
    elif XS_lib == "oldlib":
        clad_lib_id = "03c"
        gap_lib_id = "06c"
        temp_tabulation_points = [300.0, 600.0, 900.0, 1200.0, 1500, 1800.0]

    with open(f'AT10_24UOX_MPHYS_h155_mesh{number_axial_slices}_{power_scaling_factor}_mc', 'a') as output_file:
        
        output_file.write(f"% --- 4. Materials definition \n")
        for i in range(len(Teff_fuel)):
            for fuel_name in fuel_names:
                output_file.write(f"mat {fuel_name}_{i+1} sum burn 1 tmp {Teff_fuel[i]:.2f} \n")
                if XS_lib == "PyNjoy2016":
                    if 293.0<Teff_fuel[i]<550.0:
                        temp_suffix = "02c"
                    elif 550.0<=Teff_fuel[i]<900.0:
                        temp_suffix = "05c"
                    elif 900.0<=Teff_fuel[i]<1200.0:
                        temp_suffix = "09c"
                    elif 1200.0<=Teff_fuel[i]<1999.0:
                        temp_suffix = "12c"
                    elif 1999.0<=Teff_fuel[i]:
                        temp_suffix = "20c"
                elif XS_lib == "oldlib":
                    if 300.0<Teff_fuel[i]<600.0:
                        temp_suffix = "03c"
                    elif 600.0<=Teff_fuel[i]<900.0:
                        temp_suffix = "06c"
                    elif 900.0<=Teff_fuel[i]<1200.0:
                        temp_suffix = "09c"
                    elif 1200.0<=Teff_fuel[i]<1500.0:
                        temp_suffix = "12c"
                    elif 1500.0<=Teff_fuel[i]<1800.0:
                        temp_suffix = "15c"
                    elif 1800.0<=Teff_fuel[i]:
                        temp_suffix = "18c"
                output_file.write(f" 92234.{temp_suffix}   5.15910E-06 \n 92235.{temp_suffix}   5.67035E-04 \n 92238.{temp_suffix}   2.27631E-02 \n 8016.{temp_suffix}   4.66705E-02 \n")
        for i in range(len(TWater)):
            if 293.0<TWater[i]<550.0:
                    temp_suffix = "02c"
            elif 550.0<=TWater[i]<900.0:
                temp_suffix = "05c"
            elif 900.0<=TWater[i]<1200.0:
                temp_suffix = "09c"
            elif 1200.0<=TWater[i]<1999.0:
                temp_suffix = "12c"
            elif 1999.0<=TWater[i]:
                temp_suffix = "20c"
            output_file.write(f"mat H2O_{i+1} sum tmp {TWater[i]:.02f} moder lwtr 1001 \n")
            output_file.write(f" 1001.{temp_suffix}   {iso_dens_H[i]:E} \n 8016.{temp_suffix}   {iso_dens_O[i]:E} \n")
        output_file.write("\n")
        output_file.write("% helium \n mat gap sum tmp 750.0 \n")
        output_file.write(f"2004.{gap_lib_id} 1.50456e-04 \n")
        output_file.write("% clad 559K \n mat clad -6.55  tmp 559 \n")
        output_file.write(f"28060.{clad_lib_id} -1.472672e-04\n")
        output_file.write(f"40094.{clad_lib_id} -1.755744e-01\n")
        output_file.write(f"26058.{clad_lib_id} -3.949349e-06\n")
        output_file.write(f"50112.{clad_lib_id} -1.325870e-04\n")
        output_file.write(f"50116.{clad_lib_id} -2.058423e-03\n")
        output_file.write(f"40096.{clad_lib_id} -2.888887e-02\n")
        output_file.write(f"24050.{clad_lib_id} -4.173687e-05\n")
        output_file.write(f"40092.{clad_lib_id} -1.695587e-01\n")
        output_file.write(f"50114.{clad_lib_id} -9.182450e-05\n")
        output_file.write(f"50117.{clad_lib_id} -1.096647e-03\n")
        output_file.write(f"40091.{clad_lib_id} -1.097236e-01\n")
        output_file.write(f"50122.{clad_lib_id} -6.894095e-04\n")
        output_file.write(f"8017.{clad_lib_id} -5.058471e-07\n")
        output_file.write(f"50124.{clad_lib_id} -8.762916e-04\n")
        output_file.write(f"24054.{clad_lib_id} -2.453364e-05\n")
        output_file.write(f"28062.{clad_lib_id} -2.109210e-05\n")
        output_file.write(f"24053.{clad_lib_id} -9.673587e-05\n")
        output_file.write(f"50115.{clad_lib_id} -4.771906e-05\n")
        output_file.write(f"28061.{clad_lib_id} -6.508475e-06\n")
        output_file.write(f"28064.{clad_lib_id} -5.544228e-06\n")
        output_file.write(f"24052.{clad_lib_id} -8.369936e-04\n")
        output_file.write(f"8016.{clad_lib_id} -1.249494e-03\n")
        output_file.write(f"40090.{clad_lib_id} -4.976045e-01\n")
        output_file.write(f"50120.{clad_lib_id} -4.771539e-03\n")
        output_file.write(f"28058.{clad_lib_id} -3.695879e-04\n")
        output_file.write(f"26054.{clad_lib_id} -7.621504e-05\n")
        output_file.write(f"26057.{clad_lib_id} -2.916498e-05\n")
        output_file.write(f"26056.{clad_lib_id} -1.240671e-03\n")
        output_file.write(f"50118.{clad_lib_id} -3.487981e-03\n")
        output_file.write(f"50119.{clad_lib_id} -1.247577e-03\n")
        output_file.write("\n")
        output_file.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n % \n % Isotope temperature at all heights \n % \n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

        for i in range(len(Teff_fuel)):
            if XS_lib == "PyNjoy2016":
                if 293.0<Teff_fuel[i]<550.0:
                    temp_suffix = "02c"
                elif 550.0<=Teff_fuel[i]<900.0:
                    temp_suffix = "05c"
                elif 900.0<=Teff_fuel[i]<1200.0:
                    temp_suffix = "09c"
                elif 1200.0<=Teff_fuel[i]<1999.0:
                    temp_suffix = "12c"
                elif 1999.0<=Teff_fuel[i]:
                    temp_suffix = "20c"
            elif XS_lib == "oldlib":
                if 300.0<Teff_fuel[i]<600.0:
                    temp_suffix = "03c"
                elif 600.0<=Teff_fuel[i]<900.0:
                    temp_suffix = "06c"
                elif 900.0<=Teff_fuel[i]<1200.0:
                    temp_suffix = "09c"
                elif 1200.0<=Teff_fuel[i]<1500.0:
                    temp_suffix = "12c"
                elif 1500.0<=Teff_fuel[i]<1800.0:
                    temp_suffix = "15c"
                elif 1800.0<=Teff_fuel[i]:
                    temp_suffix = "18c"
            output_file.write(f"mat U235_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92235.{temp_suffix} 1.0\n")
            output_file.write(f"mat U238_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92238.{temp_suffix} 1.0\n")
            output_file.write(f"mat U234_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92234.{temp_suffix} 1.0\n")
        output_file.write("\n")

def create_detectors(number_axial_slices, power_scaling_factor, isGd):
    if isGd == False:
        fuel_names = ["UOx_A", "UOx_B", "UOx_C", "UOx_D"]
    else:
        fuel_names = ["Gd_A", "Gd_B", "Gd_C", "Gd_D", "Gd_E", "Gd_F"]
    with open(f'AT10_24UOX_MPHYS_h155_mesh{number_axial_slices}_{power_scaling_factor}_mc', 'a') as output_file:
        output_file.write("\n")
        output_file.write(f"% --- 5. Detectors definition \n")
        output_file.write(f"% --- energy grid for the detectors (1g) \n")
        output_file.write("ene 3 1 1.1E-11 1.9640E+1 \n")
        output_file.write(f"% --- energy grid for the detectors (2g) \n")
        output_file.write("ene 2 1 1.1E-11 6.25E-7 1.9640E+1 \n")
        output_file.write("\n")


        for i in range(number_axial_slices):
            output_file.write(f"det _FLUX_2G_{i+1} \n") 
            output_file.write(f"de 2 \n")
            output_file.write(f"dc {i+1}  \n")
            output_file.write("\n")

            output_file.write(f"det _FUEL_1G_{i+1} \n")
            output_file.write(f"de 3 \n")
            output_file.write(f"dr 102 U235_{i+1} \n")
            output_file.write(f"dr 102 U238_{i+1} \n")
            output_file.write(f"dr 102 U234_{i+1} \n")
            output_file.write(f"dr -6 U235_{i+1} \n")
            output_file.write(f"dr -6 U238_{i+1} \n")
            output_file.write(f"dr -6 U234_{i+1} \n")
            output_file.write(f"dr 301 U235_{i+1} \n")
            output_file.write(f"dr 301 U238_{i+1} \n")
            output_file.write(f"dr 301 U234_{i+1} \n")
            output_file.write(f"dc {i+1} \n")
            output_file.write("\n")

    return

def create_header(number_axial_slices, power_scaling_factor, XS_lib):
    # check if XS_lib is valid
    if XS_lib not in ["PyNjoy2016", "oldlib"]:
        print("Error: invalid XS library")
        sys.exit(1)
    with open(f'AT10_24UOX_MPHYS_h155_mesh{number_axial_slices}_{power_scaling_factor}_mc', 'w') as output_file:
        output_file.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n")
        output_file.write("% Author  : R.Guasch                                                                                              % \n")
        output_file.write("% Date    : 2024/10/23                                                                                            % \n")
        output_file.write("% Title   : AT10_24UOx                                                                                            % \n")
        output_file.write(f"% Purpose : Serpent2 code for square BWR pin cell, 3D model with {number_axial_slices} axial slices with varying TF,TC,DC          % \n")
        output_file.write("% Code    : Serpent2                                                                                              % \n")
        output_file.write("% Version : 2.1.31                                                                                                % \n")
        output_file.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n")
        output_file.write("\n")
        output_file.write("% --- Problem title: \n")
        output_file.write(f"set title \"AT10_24UOx Pin cell benchmark for MPHYS validation\" \n")
        output_file.write("\n")
        output_file.write("% --- Cross section library file path: \n")
        if XS_lib == "oldlib":
            output_file.write('set acelib "/home/p117902/Serpent2/xs/jeff311/sss_jeff311u.data" \n')
        elif XS_lib == "PyNjoy2016":
            output_file.write('set acelib "/home/p117902/Serpent2/xs/jeff311/JEFF-311_pynjoy2016.xsdata" \n')
        output_file.write("\n")
        output_file.write("% --- Decay data library: \n")
        output_file.write('set declib "/home/p117902/Serpent2/xs/jeff311/sss_jeff311.dec" \n')
        output_file.write("\n")
        output_file.write("% --- Neutron-induced fission yield library: \n")
        output_file.write('set nfylib "/home/p117902/Serpent2/xs/jeff311/sss_jeff311.nfy"')
        output_file.write("\n")
    return

def create_calculation_options(number_axial_slices, power_scaling_factor, neutron_population, criticality_cycles, discard_cycles, ures, XS_lib):
    with open(f'AT10_24UOX_MPHYS_h155_mesh{number_axial_slices}_{power_scaling_factor}_mc', 'a') as output_file:
        output_file.write("\n")
        output_file.write("% --- 5. Serpent parameters \n")
        output_file.write(f"set ures {ures} \n")
        output_file.write(f"% --- 5.1 Thermal scattering data for light water: \n")
        if XS_lib == "PyNjoy2016":
            output_file.write("therm lwtr lwtr.05t \n")
        elif XS_lib == "oldlib":
            output_file.write("%therm lwtr lwj3.11t \n")
        output_file.write("\n")
        output_file.write("% --- 5.2 Reflective boundary condition on +/- X/Y and void in +/-Z: \n")
        output_file.write("set bc 2 2 1 \n")
        output_file.write("\n")
        output_file.write("% --- 5.3 Neutron population and criticality cycles: \n")
        output_file.write(f"set gcu 0 \n")
        output_file.write(f"set nfg 2 0.625E-6 \n")
        output_file.write("%set nbuf 10.0 \n")
        output_file.write(f"set pop {neutron_population} {criticality_cycles} {discard_cycles} \n")
        output_file.write("\n")
        output_file.write("% --- 5.4 Output options: \n")
        output_file.write("plot 3 500 500 0.0 \n")
        output_file.write("plot 3 500 500 1 \n")
        output_file.write("plot 3 500 500 100 \n")
        output_file.write("plot 3 500 500 190 \n")
        output_file.write("plot 3 500 500 200 \n")
        output_file.write("plot 3 500 500 300 \n")
        output_file.write("\n")
        output_file.write("plot 1 500 500 \n")
        output_file.write("mesh 3 500 500 \n")
        output_file.write("\n")

        return

def add_moderator_data():
    """
    Specify the moderator data for the Serpent2 case :
    find the right lwtr thermal scattering data file and interpolate the data accordingly following card structure :
    % --- 5.1 Thermal scattering data for light water:
    lwtr <temp> lwtr.<suffix1>t lwtr.<suffix2>t
    or 
    lwtr lwtr.<suffix>t if interpolation is not needed / possible.
    """
    return    
### Clean this up and make it functional for easy interfacing with Version5 - cross validation workflow
TH_underRelaxationFactor = 0.1
Pow_underRelaxationFactor = 0.5 # Under relaxation factor for the power axial distribution to be tested, 0.1 used in Serpent/OpenFoam coupling
relax_TH = False # Under relaxation of the TH fields for the next iteration
relax_Pow = True # Under relaxation of the Power distribution for the next iteration

voidFractionCorrel = 'EPRIvoidModel' # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
frfaccorel = "blasius" # 'base', 'blasius', 'Churchill', 
P2Pcorel = "lockhartMartinelli" # 'base', 'HEM1', 'HEM2', 'MNmodel', 'lockhartMartinelli'

number_axial_slices_all_cases = [160]
#number_axial_slices_all_cases = [70]
#number_axial_slices = 10

pow_scaling_factor_all_cases = [1]
#pow_scaling_factor = 1
neutrons_per_cycle = 50000
num_cycles = 10000
cycles_to_discard = 500
ures_activation_option = 1
cross_section_library = "PyNjoy2016"

height = 155.0 # 380.0
if height == 155.0:
    read_id = "h1555"
elif height == 380.0:
    read_id = "h380"

if relax_Pow:
    relaxPOW_id = f"relaxedPOW_{Pow_underRelaxationFactor}"
else:
    relaxPOW_id = "non_relaxedPOW"
if relax_TH:
    relaxTH_id = f"relaxedTH_{TH_underRelaxationFactor}"
else:
    relaxTH_id = "non_relaxedTH"

for number_axial_slices in number_axial_slices_all_cases:
    for pow_scaling_factor in pow_scaling_factor_all_cases:
        Teff_fuel = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{number_axial_slices}_{pow_scaling_factor}/Data/TeffFuel_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_{relaxPOW_id}_{relaxTH_id}.txt")
        #print(T_data)
        TWater = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{number_axial_slices}_{pow_scaling_factor}/Data/Twater_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_{relaxPOW_id}_{relaxTH_id}.txt")
        DWater = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{number_axial_slices}_{pow_scaling_factor}/Data/rho_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_{relaxPOW_id}_{relaxTH_id}.txt")

        DWater = DWater*1e-3 # kg/m^3 to g/cm^3
        iso_dens_H = []
        iso_dens_O = []
        for rho in DWater:
            iso_dens_H.append(DensToIsoDens.DensToIsoDens(rho)[0])
            iso_dens_O.append(DensToIsoDens.DensToIsoDens(rho)[1])


        # check number of elements in each array
        if len(Teff_fuel) != len(TWater) or len(TWater) != len(DWater) or len(DWater) != len(iso_dens_H) or len(iso_dens_H) != len(iso_dens_O):
            print("Error: unable to generate Serpent2 case")
            sys.exit(1)
        else:
            
            # Create Serpent2 case header
            create_header(number_axial_slices=len(Teff_fuel), 
                        power_scaling_factor=pow_scaling_factor,
                        XS_lib = cross_section_library)

            # Create Serpent2 geometry input
            create_geometry(number_axial_slices=len(Teff_fuel), 
                            power_scaling_factor=pow_scaling_factor, 
                            pitch=1.295, 
                            fuel_radius=0.4435, 
                            gap_radius=0.4520, 
                            clad_radius= 0.5140, 
                            height=height, isGd=False)

            # Create Serpent2 material volumes
            create_material_volumes(len(Teff_fuel), 
                                    pow_scaling_factor, 
                                    Teff_fuel, TWater, 
                                    iso_dens_O, 
                                    iso_dens_H, 
                                    isGd=False,
                                    XS_lib = cross_section_library)

            # Create Serpent2 detectors
            create_detectors(number_axial_slices=len(Teff_fuel), 
                            power_scaling_factor=pow_scaling_factor, 
                            isGd=False)
            
            # Create Serpent2 calculation options
            create_calculation_options(number_axial_slices=len(Teff_fuel), 
                                    power_scaling_factor=pow_scaling_factor, 
                                    neutron_population=neutrons_per_cycle, 
                                    criticality_cycles=num_cycles, 
                                    discard_cycles=cycles_to_discard, ures=ures_activation_option,
                                    XS_lib = cross_section_library)

    






