# Script used to create instances of THM classes and solve convection + conduction in fuel rod.
# author : R. Guasch
# Purpose : prototyping for further developments in THM module of Donjon5

from THM_MONO import FDM_HeatConductionInFuelPin as FDM_Fuel
from THM_MONO import FVM_ConvectionInCanal_MONO as FVM_Canal
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

#Cas 1 : base parameters
# Parameters used to create object from FDM_HeatConductioninFuelpin class
Qfiss = 0.3e9 # W/m^3
fuel_radius = 5.6e-3 # m
gap_width = 0.54e-3 # m
clad_width = 0.38e-3 # m
k_fuel = 5 # W/m/K
H_gap = 10000 # W/m^2/K
k_clad = 10 # W/m/K
I_f = 8
I_c = 3

# Paramters used to create object from FVM_ConvectioninCanal class
canal_type = "cylindrical"
canal_width = 0.5e-3 # m
Lf = 2 # m
T_in = 500 # K
Q_flow = 7000 # kg/m^2/s
P_cool = 10.8 #MPa
I_z = 10

initial_water_z0 = IAPWS97(P=P_cool,T=T_in)
h_ini = initial_water_z0.h 
z=0


def run_Conduction_In_Fuel(fuel_radius, I_f, gap_width, clad_width, I_c, z, Qfiss, k_fuel, k_clad, H_gap, T_surf):

    heat_conduction = FDM_Fuel(fuel_radius, I_f, gap_width, clad_width, I_c, Qfiss, k_fuel, k_clad, H_gap, z, T_surf)

    for i in range(1,heat_conduction.N_node-1):

        if i<heat_conduction.I_f-1: # setting Aij and Di values for nodes inside the fuel 
            heat_conduction.set_ADi_cond(i, 
                                -heat_conduction.get_Di_half(i-1), 
                                heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                -heat_conduction.get_Di_half(i), 
                                heat_conduction.deltaA_f*heat_conduction.Qfiss)
        elif i==heat_conduction.I_f-1: # setting Aij and Di values for last fuel element
            heat_conduction.set_ADi_cond(i,
                            -heat_conduction.get_Di_half(i-1),
                            (heat_conduction.get_Di_half(i-1)+heat_conduction.get_Ei_gap()),
                            -heat_conduction.get_Ei_gap(),
                            heat_conduction.deltaA_f*heat_conduction.Qfiss)
        elif i==heat_conduction.I_f: # setting Aij and Di values first fuel / gap interface
            heat_conduction.set_ADi_cond(i, 
                                -heat_conduction.get_Ei_gap(), 
                                heat_conduction.get_Ei_gap()+heat_conduction.get_Gi(), 
                                -heat_conduction.get_Gi(), 
                                0)
        elif i==heat_conduction.I_f+1: # setting Aij and Di values second gap / clad interface
            heat_conduction.set_ADi_cond(i, 
                                -heat_conduction.get_Gi(), 
                                heat_conduction.get_Fi_gap()+heat_conduction.get_Gi(), 
                                -heat_conduction.get_Fi_gap(), 
                                0)
        elif i>heat_conduction.I_f+1 : # setting Aij and Di for all elements in the clad, apart from the last one
            heat_conduction.set_ADi_cond(i, 
                                -heat_conduction.get_Di_half(i-1), 
                                heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                -heat_conduction.get_Di_half(i), 
                                0)
    A0,Am1 = np.zeros(heat_conduction.N_node), np.zeros(heat_conduction.N_node) 
    A0[:2] = [heat_conduction.get_Di_half(0), -heat_conduction.get_Di_half(0)]
    Am1[-2:] = [-heat_conduction.get_Di_half(heat_conduction.N_node-2), heat_conduction.get_Di_half(heat_conduction.N_node-2)+heat_conduction.get_Ei_clad()]
    D0 = heat_conduction.deltaA_f*heat_conduction.Qfiss
    Dm1 = heat_conduction.get_Ei_clad()*heat_conduction.T_surf
    print(f"Ei_clad = {heat_conduction.get_Ei_clad()}")
    print(f"T_surf = {heat_conduction.T_surf}")
    heat_conduction.set_CL_cond(A0, Am1, D0, Dm1)
   
    for row in heat_conduction.A:
        line = "[  "
        for elem in row:
            line+=f"{elem:.3f}   "
        line += "  ]\n"
        print(line)
    heat_conduction.solve_T_in_pin()
    return heat_conduction




convection_test = FVM_Canal(Lf, T_in, Q_flow, P_cool, I_z, canal_type, rf=fuel_radius, rc=fuel_radius+gap_width+clad_width, rw=fuel_radius+gap_width+clad_width+canal_width)
convection_test.set_Fission_Power(0.3e9, variation_type="constant")
for i in range(1,convection_test.N_vol-1):
    convection_test.set_ADi_conv(i,
                            ci=-1,
                            ai=1,
                            bi=0,
                            di = convection_test.q_fluid[i]*convection_test.dz/(convection_test.Q_flow*convection_test.A_canal))
A0,Am1 = np.zeros(convection_test.N_vol), np.zeros(convection_test.N_vol)
A0[0] = 1
D0 = convection_test.h_z0 + convection_test.q_fluid[i]*convection_test.dz/(2*convection_test.Q_flow*convection_test.A_canal)
Am1[-2:]=[-1, 1]
Dm1 = convection_test.q_fluid[i]*convection_test.dz/(convection_test.Q_flow*convection_test.A_canal)
convection_test.set_CL_conv(A0,Am1,D0,Dm1)
for row in convection_test.A:
        line = "[  "
        for elem in row:
            line+=f"{elem:.3f}   "
        line += "  ]\n"
        print(line)
convection_test.h_z = convection_test.solve_h_in_canal()
print(convection_test.h_z)


Tsurf = convection_test.compute_T_surf()
print(f"T_surf = {Tsurf}")
print(f"T_water = {convection_test.T_water}")
temp_distrib = []
for axial_plane_nb in range(convection_test.N_vol):
    z = convection_test.z_mesh[axial_plane_nb]
    T_surf = convection_test.T_surf[axial_plane_nb]
    Qfiss = convection_test.get_Fission_Power()[axial_plane_nb]
    print(f"Qfiss={Qfiss}")
    temp_distrib.append(run_Conduction_In_Fuel(fuel_radius, I_f, fuel_radius+gap_width, fuel_radius+gap_width+clad_width, I_c, z, Qfiss, k_fuel, k_clad, H_gap, T_surf))
    print(f"temperature disrtib at first z element is {temp_distrib[0].T_distrib}")
    print(f"computed temp distrib at z={z}")
    if z == 0.7:
        fig,ax = plt.subplots(dpi=200)
        print(f"Calculation radial mesh is = {np.sqrt(2*temp_distrib[-1].A_calculation_mesh)}")
        print(f"plot mesh is {temp_distrib[-1].plot_mesh}")
        ax.scatter(temp_distrib[-1].plot_mesh, temp_distrib[-1].T_distrib, marker = "x", s=5, label="Radial temperature distribution in Fuel rod.")
        #ax.scatter(np.sqrt(2*temp_distrib[-1].A_calculation_mesh), temp_distrib[-1].T_distrib, marker = "x", s=5, label="Radial temperature distribution in Fuel rod.")
        ax.legend(loc = "best")
        ax.grid()
        ax.set_xlabel(f"Radial position in {temp_distrib[-1].plotting_units}")
        ax.set_ylabel(f"Temperature in K")
        ax.set_title(f"Temperature distribution in fuel rod at z = {z}, case 1")
        fig.savefig(f"Case1_Figure_plane{axial_plane_nb}")
        colors = ["red", "yellow", "green", "blue"]
        temp_distrib[-1].extend_to_canal_visu(rw = convection_test.wall_dist, Tw = convection_test.T_water[axial_plane_nb])
        print(f"T_surf = {convection_test.T_surf[axial_plane_nb]} K and T_water = {convection_test.T_water[axial_plane_nb]} K")
        print(f"radii at bounds {temp_distrib[-1].radii_at_bounds}")
        fig_filled, axs = plt.subplots(dpi=200)
        for i in range(len(temp_distrib[-1].physical_regions_bounds)-1):
            axs.fill_between(x=temp_distrib[-1].radii_at_bounds, y1=temp_distrib[-1].T_distrib[1:], y2=400*np.ones(len(temp_distrib[-1].radii_at_bounds)),where=(temp_distrib[-1].radii_at_bounds>=temp_distrib[-1].physical_regions_bounds[i])&(temp_distrib[-1].radii_at_bounds<=temp_distrib[-1].physical_regions_bounds[i+1]), color = colors[i])
        
        #axs.fill_betweenx(y=temp_distrib[-1].T_distrib[1:], x1=temp_distrib[-1].radii_at_bounds, where=(temp_distrib[-1].radii_at_bounds<temp_distrib[-1].r_f), facecolor='red')
        #axs.fill_betweenx(convection_test.T_surf[axial_plane_nb])
        ax.legend(loc = "best")
        ax.grid()
        ax.set_xlabel(f"Radial position in {temp_distrib[-1].plotting_units}")
        ax.set_ylabel(f"Temperature in K")
        ax.set_title(f"Temperature distribution in fuel rod at z = {z}, case 1")
        fig_filled.savefig(f"Case1_Figure_plane{axial_plane_nb}_colors")

