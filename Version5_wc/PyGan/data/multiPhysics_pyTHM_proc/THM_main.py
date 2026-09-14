#Used to run the THM prototype class and compare the results with a reference THM_DONJON case.
#Authors : Clement Huet, Raphael Guasch


from THM_conduction import FDM_HeatConductionInFuelPin as FDM_Fuel
from THM_convection import DFMclass
from THM_DONJON_parser import THM_DONJON_parser
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt
import pandas as pd
from openpyxl import Workbook
import os
import re

class Version5_THM_prototype:
    def __init__(self, case_name, canal_type,
                 canal_radius, fuel_radius, gap_radius, clad_radius, fuel_rod_length, tInlet, pOutlet, qFlow, Powtot, axial_p_form, fraction_pow_fuel,
                 k_fuel, H_gap, k_clad, I_z, I_f, I_c, plot_at_z, solveConduction,
                 dt, t_tot, frfaccorel = 'base', P2Pcorel = 'base', voidFractionCorrel = 'GEramp', numericalMethod= 'FVM'):
        """
        Main constructor for THM case, first set of parameters correspond to canal properties, second set to fuel/gap/clad properties
        The structure followed is : 
        In FVM_ConvectionInCanal class : use a finite volume method to solve heat convection in the canal, then use the Dittus-Boelter correlation to obtain the convective heat transfer coef 
        between the water and the fuel rod's outer surface. This allows to solve for the temperature at this outer surface. 
        Then in the FDM_HeatConductionInFuelPin class, solve for the heat conduction using MCFD method. Compute temperature at the center of the fuel rod.
        Options to plot results can be activated giving an array of z values at which the results should be plotted.
        """
        self.name = case_name
        # time atributes to prepare for transient simulations
        self.t0 = 0
        self.dt = dt
        self.t_end = t_tot

        # canal attributes

        self.r_w = canal_radius # outer canal radius (m) if type is cylindrical, if type = square rw is the diameter of inscribed circle in the square canal, ie half the square's side.
        self.canal_type = canal_type # cylindrical or square, used to determine the cross sectional flow area in the canal and the hydraulic diameter
        self.Lf = fuel_rod_length # fuel rod length in m
        

        self.tInlet = tInlet
        self.qFlow = qFlow #  mass flux in kg/s, assumed to be constant along the axial profile.
        self.I_z = I_z # number of mesh elements on axial mesh
        self.rhoInlet = 1000
        self.pOutlet =  pOutlet #Pa
        if canal_type == "cylindrical":
            self.flowArea = np.pi*self.r_w**2 - np.pi*fuel_radius**2
        else:
            self.flowArea = self.r_w**2 - fuel_radius**2
            
        self.uInlet = self.qFlow / (self.rhoInlet*self.flowArea) #m/s

        self.Powtot = Powtot # Total reactor power in W
        self.axial_pow_form = axial_p_form # axial power form factors, representing the power distribution along the axial dimension of the fuel rod, used to compute the fission power in the fuel rod.
        self.Fpow = fraction_pow_fuel # fraction of the total power that is deposited in the fuel, used to compute the fission power in the fuel rod.

        self.r_f = fuel_radius # fuel pin radius in meters
        self.gap_r = gap_radius # gap radius in meters, used to determine mesh elements for constant surface discretization
        self.clad_r = clad_radius # clad radius in meters, used to determine mesh elements for constant surface discretization
        self.k_fuel = k_fuel # thermal conductivity coefficient in fuel W/m/K
        self.H_gap = H_gap # Heat transfer coefficient through gap W/m^2/K
        self.k_clad = k_clad # thermal conductivity coefficient in clad W/m/K
        self.I_f = I_f # number of mesh elements in the fuel
        self.I_c = I_c # number of mesh elements in clad

        self.frfaccorel = frfaccorel # friction factor correlation
        self.P2Pcorel = P2Pcorel # pressure drop correlation
        self.voidFractionCorrel = voidFractionCorrel # void fraction correlation
        self.numericalMethod = numericalMethod # numerical method used to solve the convection problem in the canal
        #Poro = 0.5655077285

        self.plot_results = plot_at_z
        self.solveConduction = solveConduction

        print(f"$$$---------- THM: prototype, case treated : {self.name}.")
        if self.dt == 0 :
            self.transient = False
            print("$$$---------- THM: prototype, steady state case.")

        # Prepare and solve 1D heat convection along the z direction in the canal.
        print("$$---------- Calling DFM class.")
        print(f"Setting up heat convection solution along the axial dimension. zmax = {self.Lf} m with {self.I_z} axial elements.")
        # Create an object of the class DFMclass
        print(f'self.I_z: {self.I_z}')
        print(f'self.qFlow: {self.qFlow}')
        print(f'self.pOutlet: {self.pOutlet}')
        print(f'self.Lf: {self.Lf}')
        print(f'self.r_f: {self.r_f}')
        print(f'self.clad_r: {self.clad_r}')
        print(f'self.r_w: {self.r_w}')
        print(f'self.Dz: {self.Lf/self.I_z}')
        print(f'self.dt: {self.dt}')
        print(f'Courant number: {self.uInlet*self.dt/(self.Lf/self.I_z)}')
        print(f"Numerical Method {numericalMethod}")
        self.convection_sol = DFMclass(self.canal_type, self.I_z, self.tInlet, self.qFlow, self.pOutlet, self.Lf, self.r_f, self.clad_r, self.r_w, self.numericalMethod, self.frfaccorel, self.P2Pcorel, self.voidFractionCorrel, dt = self.dt, t_tot = self.t_end)
        print(f'Hydraulic diameter: {self.convection_sol.D_h}')
        # Set the fission power in the fuel rod
        self.convection_sol.set_Fission_Power(self.Powtot, self.axial_pow_form, self.Fpow) # set the fission power in the fuel rod, given the total power, the axial power form factors and the fraction of power deposited in the fuel
        # Resolve the DFM
        self.convection_sol.resolveDFM()
        #print(f'Pressure: {self.convection_sol.P[-1]} Pa')
        #print(f'Enthalpy: {self.convection_sol.H[-1]} J/kg')
        #print(f'Void fraction: {self.convection_sol.voidFraction[-1]}')
        #print(f'Density: {self.convection_sol.rho[-1]} kg/m^3')
        #print(f'DV: {self.convection_sol.DV}')
        if self.solveConduction:
            self.Tsurf = self.convection_sol.compute_T_surf()
            #print(f'Temperature at the surface: {self.Tsurf} K')
            #print(f'Temperature of water: {self.convection_sol.T_water} K')

        if self.solveConduction:
            # Prepare and solve 1D radial heat conduction in the fuel rod, given a Clad surface temperature as a bondary condition 
            self.SetupAndSolve_Conduction_at_all_z() # creates a list of Temperature distributions in the fuel rod given a surface temperature computed by solving the conection problem
            self.get_TFuel_rowlands() # compute and store in the T_eff_fuel attribute the effective fuel temperature given by the Rowlands formula
            self.get_Tfuel_surface() # store in the T_fuel_surface attribute the fuel surface temperature computed

            # extend to Twater : adding a mesh point corresponding to the middle of the canal in the plotting array, add rw to the bounds array and add Twater to the results array
            for index_z in range(len(self.convection_sol.z_mesh)):
                self.T_distributions_axial[index_z].extend_to_canal_visu(rw = self.convection_sol.wall_dist, Tw = self.convection_sol.T_water[index_z])
                
            if self.plot_results:
                for z_val in self.plot_results:
                    self.plot_Temperature_at_z(z_val)
    
    def set_transitoire(self, t_tot, Tini, dt):
        self.t_tot, self.dt = t_tot, dt           
        self.N_temps = round(self.t_tot / self.dt) # pas de temps (timesteps), il faut etre un nombre entier
        self.T = np.zeros((self.N_temps+1, self.N_vol)) # tableau 2D de temperature. 
        self.T[0] = Tini # Tini est une liste
        return


    def SetupAndSolve_Conduction_at_all_z(self, transient = False):
        self.T_distributions_axial = []
        for axial_plane_nb in range(self.convection_sol.nCells):
            z = self.convection_sol.z_mesh[axial_plane_nb]
            T_surf = self.convection_sol.T_surf[axial_plane_nb]
            Qfiss = self.convection_sol.get_Fission_Power()[axial_plane_nb]
            self.T_distributions_axial.append(self.run_Conduction_In_Fuel_at_z(z,Qfiss,T_surf, transient))

        return
    
    def run_Conduction_In_Fuel_at_z(self,z,Qfiss_z,T_surf_z, transient = False):
        print(f"$$---------- Setting up FDM_HeatConductionInFuelPin class for z = {z} m, Qfiss(z) = {Qfiss_z} W/m^3 and T_surf(z) = {T_surf_z} K")
        heat_conduction = FDM_Fuel(self.r_f, self.I_f, self.gap_r, self.clad_r, self.I_c, Qfiss_z, self.k_fuel, self.k_clad, self.H_gap, z, T_surf_z)
        if transient:
            print(f"Error: transient case not implemented yet")
        else:
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
                                    (heat_conduction.get_Di_half(i-1)+heat_conduction.get_Ei_fuel()),
                                    -heat_conduction.get_Ei_fuel(),
                                    heat_conduction.deltaA_f*heat_conduction.Qfiss)
                elif i==heat_conduction.I_f: # setting Aij and Di values first fuel / gap interface
                    heat_conduction.set_ADi_cond(i, 
                                        -heat_conduction.get_Ei_fuel(), 
                                        heat_conduction.get_Ei_fuel()+heat_conduction.get_Gi(), 
                                        -heat_conduction.get_Gi(), 
                                        0)
                elif i == heat_conduction.I_f+1: # setting Aij and Di values second gap / clad interface
                    heat_conduction.set_ADi_cond(i, 
                                        -heat_conduction.get_Gi(), 
                                        heat_conduction.get_Fi_gap()+heat_conduction.get_Gi(), 
                                        -heat_conduction.get_Fi_gap(), 
                                        0)
                elif i == heat_conduction.I_f+2: # Treating the first clad element interface with the gap.
                    heat_conduction.set_ADi_cond(i,
                                                -heat_conduction.get_Ei_gap(),
                                                (heat_conduction.get_Di_half(i)+heat_conduction.get_Ei_gap()),
                                                -heat_conduction.get_Di_half(i),
                                                0)
                elif i>heat_conduction.I_f+2 : # setting Aij and Di for all elements in the clad, apart from the last one
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
            heat_conduction.set_CL_cond(A0, Am1, D0, Dm1)
            print(f"$---------- Solving for T(r) using the Finite Difference Method, at z = {z}.")
            heat_conduction.solve_T_in_pin()
        return heat_conduction
    
    def get_TFuel_rowlands(self):
        self.T_eff_fuel = np.zeros(self.convection_sol.nCells)
        for i in range(len(self.T_distributions_axial)):
            self.T_distributions_axial[i].compute_T_eff()
            T_eff_z = self.T_distributions_axial[i].T_eff
            self.T_eff_fuel[i] = T_eff_z
        return
    
    def get_Tfuel_surface(self):
        self.T_fuel_surface = np.zeros(self.convection_sol.nCells)
        for i in range(len(self.T_distributions_axial)):
            if len(self.T_distributions_axial[i].T_distrib) == self.T_distributions_axial[i].N_node:
                T_surf_fuel_z = self.T_distributions_axial[i].T_distrib[self.I_f]
            else:
                T_surf_fuel_z = self.T_distributions_axial[i].T_distrib[self.I_f+1]
            self.T_fuel_surface[i] = T_surf_fuel_z
        return

    def plot_Temperature_at_z(self, z_val):
        print(f"$$---------- Plotting Temperature distribution in rod + canal z = {z_val} m")

        print(f"z_val is {z_val}")
        print(f"z_mesh is {self.convection_sol.z_mesh}")
        if z_val in self.convection_sol.z_mesh:
            plane_index = int(np.where(self.convection_sol.z_mesh==z_val)[0][0])
            Temperature_distrib_to_plot = self.T_distributions_axial[plane_index].T_distrib
            plotting_mesh = self.T_distributions_axial[plane_index].plot_mesh
            radii_at_bounds = self.T_distributions_axial[plane_index].radii_at_bounds
            physical_regions_bounds = self.T_distributions_axial[plane_index].physical_regions_bounds
            plotting_units = self.T_distributions_axial[plane_index].plotting_units
            Tsurf = self.convection_sol.T_surf[plane_index]
            Twater = self.convection_sol.T_water[plane_index]
            Tcenter = self.T_distributions_axial[plane_index].T_center
        else: # Interpolate between nearest z values to obtain Temperature distribution at a given z.
            second_plane_index = np.where(self.convection_sol.z_mesh>z_val)[0][0]
            first_plane_index = second_plane_index-1
            plane_index = (first_plane_index+second_plane_index)/2
            print(f"plane index used is {plane_index}")
            plotting_mesh = self.T_distributions_axial[first_plane_index].plot_mesh
            radii_at_bounds = self.T_distributions_axial[first_plane_index].radii_at_bounds
            physical_regions_bounds = self.T_distributions_axial[first_plane_index].physical_regions_bounds
            Temperature_distrib_to_plot = self.T_distributions_axial[first_plane_index].T_distrib+(z_val-self.convection_sol.z_mesh[first_plane_index])*(self.T_distributions_axial[second_plane_index].T_distrib-self.T_distributions_axial[first_plane_index].T_distrib)/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            plotting_units = self.T_distributions_axial[first_plane_index].plotting_units
            Tsurf = self.convection_sol.T_surf[first_plane_index] + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.convection_sol.T_surf[second_plane_index]-self.convection_sol.T_surf[first_plane_index])/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            Twater = self.convection_sol.T_water[first_plane_index] + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.convection_sol.T_water[second_plane_index]-self.convection_sol.T_water[first_plane_index])/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            Tcenter = self.T_distributions_axial[first_plane_index].T_center + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.T_distributions_axial[second_plane_index].T_center-self.T_distributions_axial[first_plane_index].T_center)/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])

        
        
        if (isinstance(plane_index, int)):
            plane_index_print = plane_index
        else:
            plane_index_print = str(plane_index).split(".")[0]+str(plane_index).split(".")[1]
        print(f"at z = {z_val}, temp distrib is = {Temperature_distrib_to_plot}")
        print(f'z_axis is {plotting_mesh}')
        colors = ["lime", "bisque", "chocolate", "royalblue"]
        labels = ["Fuel", "Gap", "Clad", "Water"]
        fig_filled, axs = plt.subplots()
        for i in range(len(physical_regions_bounds)-1):
            axs.fill_between(x=radii_at_bounds, y1=(Tcenter+50)*np.ones(len(radii_at_bounds)), y2=(Twater-50)*np.ones(len(radii_at_bounds)),where=(radii_at_bounds>=physical_regions_bounds[i])&(radii_at_bounds<=physical_regions_bounds[i+1]), color = colors[i], label = labels[i])
        axs.scatter(plotting_mesh, Temperature_distrib_to_plot, marker = "D", color="black",s=10, label="Radial temperature distribution in Fuel rod.")
        axs.legend(loc = "best")
        axs.grid()
        axs.set_xlabel(f"Radial position in {plotting_units}")
        axs.set_ylabel(f"Temperature in K")
        axs.set_title(f"Temperature distribution in fuel rod at z = {z_val}, {self.name}")
        fig_filled.savefig(f"{self.name}_Figure_plane{plane_index_print}_colors")
        plt.show()

    def plotColorMap(self):
        print("$$---------- Plotting colormap of temperature distribution in fuel rod.")
        fig, ax = plt.subplots()
        T = []
        for i in range(len(self.T_distributions_axial)):
            T.append(self.T_distributions_axial[i].T_distrib)
        R = self.T_distributions_axial[0].plot_mesh
        Z = self.convection_sol.z_mesh
        print(f'T: {T}')
        print(f'R: {R}')
        print(f'Z: {Z}')
        plt.xlabel('Rayon (mm)')
        plt.ylabel('Hauteur (m)')
        plt.title('Temperature (K) en fonction du rayon et de la hauteur')
        plt.pcolormesh(R, Z, T, cmap = 'plasma') 
        plt.colorbar()
        plt.show()


    def compare_with_THM_DONJON(self, THM_DONJON_path, visu_params):
        """
        Function used to compare the results obtained in the current instance with a reference DONJON/THM.
        """
        self.reference_case = THM_DONJON_parser(THM_DONJON_path, self.t0, self.dt, self.t_end, self.convection_sol.nCells, self.Lf)
        self.visu_TFuel = visu_params[0]
        self.visu_TFuelSurface = visu_params[1]
        self.visu_TWater = visu_params[2]
        self.visu_deltaTFuel = visu_params[3]
        self.visu_deltaTFuelSurface = visu_params[4]
        self.visu_deltaTWater = visu_params[5]
        print(f"$$---------- Comparing results with DONJON/THM: results from {THM_DONJON_path}.")
        #print(f"The fuel efective teperature from reference results is : {self.reference_case.TCOMB}")
        # Analyzing relaive errors on the fuel effective temperature, fuel surface temperature and water temperature
        if self.dt == 0:
            TEFF_FUEL = np.flip(self.reference_case.TCOMB[0])
            TSURF_FUEL = np.flip(self.reference_case.TSURF[0])
            TCOOL = np.flip(self.reference_case.TCOOL[0])
            # plotting the relative errors on properties in steady state case
            if self.visu_TFuel:
                print("$ --- Visualizing TFUEL")
                fig, ax = plt.subplots(dpi=200)
                ax.plot(self.convection_sol.z_mesh, self.T_eff_fuel, label="THM prototype")
                ax.plot(self.convection_sol.z_mesh, TEFF_FUEL, label="DONJON/THM")
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Efective fuel temperature in K")
                ax.set_title("Effective fuel temperature comparison")        
                ax.legend(loc="best")
                fig.savefig(f"{self.name}_Teff_fuel_comparison_THM_DONJON")

            
            if self.visu_deltaTFuel:
                print("$ --- Visualizing error on TFUEL")
                fig, ax = plt.subplots(dpi=200)
                error_Tfuel = []
                for i in range(len(self.convection_sol.z_mesh)):
                    error_Tfuel.append((self.T_eff_fuel[i]-TEFF_FUEL[i])*100/TEFF_FUEL[i])
                ax.plot(self.convection_sol.z_mesh, error_Tfuel, label="THM prototype", marker="x", linestyle="--", markersize=5)
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Relative error on effective fuel temperature (%)")
                ax.set_title("Error on effective fuel temperature (Prototype vs DONJON)")        
                ax.legend(loc="best")
                ax.grid()
                fig.savefig(f"{self.name}_error_Teff_fuel_THM_DONJON")
            
            
            if self.visu_TFuelSurface:
                print("$ --- Visualizing TSURF")
                fig, ax = plt.subplots(dpi=200)
                ax.plot(self.convection_sol.z_mesh, self.T_fuel_surface, label="THM prototype")
                ax.plot(self.convection_sol.z_mesh, TSURF_FUEL, label="DONJON/THM")
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Fuel surface temperature in K")
                ax.set_title("Fuel surface temperature comparison")        
                ax.legend(loc="best")
                fig.savefig(f"{self.name}_Tsurf_comparison_THM_DONJON")
            
            if self.visu_deltaTFuelSurface:
                print("$ --- Visualizing error on TSURF")
                fig, ax = plt.subplots(dpi=200)
                error_Tsurf = []
                for i in range(len(self.convection_sol.z_mesh)):
                    error_Tsurf.append((self.T_fuel_surface[i]-TSURF_FUEL[i])*100/TSURF_FUEL[i])
                ax.plot(self.convection_sol.z_mesh, error_Tsurf, label="THM prototype", marker="x", linestyle="--", markersize=5)
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Relative error on fuel surface temperature (%)")
                ax.set_title("Error on fuel surface temperature (Prototype vs DONJON)")        
                ax.legend(loc="best")
                fig.savefig(f"{self.name}_error_Tsurf_THM_DONJON")

            if self.visu_TWater:
                print("$ --- Visualizing TCOOL")
                fig, ax = plt.subplots(dpi=200)
                ax.plot(self.convection_sol.z_mesh, self.convection_sol.T_water, label="THM prototype")
                ax.plot(self.convection_sol.z_mesh, TCOOL, label="DONJON/THM")
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Coolant temperature in K")
                ax.set_title("Coolant temperature comparison")        
                ax.legend(loc="best")
                fig.savefig(f"{self.name}_Twater_comparison_THM_DONJON")
            
            if self.visu_deltaTWater:
                print("$ --- Visualizing error on TCOOL")
                fig, ax = plt.subplots(dpi=200)
                error_Twater = []
                for i in range(len(self.convection_sol.z_mesh)):
                    error_Twater.append((self.convection_sol.T_water[i]-TCOOL[i])*100/TCOOL[i])
                ax.plot(self.convection_sol.z_mesh, error_Twater, label="THM prototype", marker="x", linestyle="--", markersize=5)
                ax.set_xlabel("Axial position in m")
                ax.set_ylabel("Relative error on coolant temperature (%)")
                ax.set_title("Error on coolant temperature (Prototype vs DONJON)")        
                ax.legend(loc="best")
                fig.savefig(f"{self.name}_error_Twater_THM_DONJON")

            plt.show()
        #print('TCOMB',TCOMB)
        return


    def get_TH_parameters(self):
        if self.solveConduction:
            print(f'Tfuel : {self.T_eff_fuel} K')
            print(f'Twater : {self.convection_sol.T_water} K')
            print(f'Water void fraction: {self.convection_sol.voidFraction[-1]} K')
            print(f'Water density: {self.convection_sol.rho[-1]} kg/m^3')
            
            return np.array(self.T_eff_fuel), np.array(self.convection_sol.T_water), np.array(self.convection_sol.rho[-1]), np.array(self.convection_sol.voidFraction[-1]), np.array(self.convection_sol.P[-1]), np.array(self.convection_sol.U[-1]), np.array(self.convection_sol.H[-1])
        
        else: 
            print(f'Twater : {self.convection_sol.T_water} K')
            print(f'Water void fraction: {self.convection_sol.voidFraction[-1]} K')
            print(f'Water density: {self.convection_sol.rho[-1]} kg/m^3')
            
            return [0], self.convection_sol.T_water, self.convection_sol.voidFraction[-1], np.array(self.convection_sol.rho[-1]), self.convection_sol.P[-1], self.convection_sol.U[-1], self.convection_sol.H[-1]


    def plotThermohydraulicParameters(self, visuParam):
        
        if visuParam[0]:
            fig1, ax1 = plt.subplots()
            if self.solveConduction:
                ax1.plot(self.convection_sol.z_mesh, self.T_eff_fuel, label="Fuel temperature")
            ax1.plot(self.convection_sol.z_mesh, self.convection_sol.T_water, label="Coolant temperature")
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Temperature in K")
            ax1.set_title("Temperature distribution pincell")
            ax1.legend(loc="best")

        if visuParam[1]:
            fig2, ax2 = plt.subplots()
            ax2.plot(self.convection_sol.z_mesh, self.convection_sol.voidFraction[-1], label="Void fraction")
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Void fraction")
            ax2.set_title("Void fraction distribution in coolant canal")
            ax2.legend(loc="best")

        if visuParam[2]:
            fig3, ax3 = plt.subplots()
            ax3.plot(self.convection_sol.z_mesh, self.convection_sol.rho[-1], label="Density")
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Density in kg/m^3")
            ax3.set_title("Density distribution in coolant canal")
            ax3.legend(loc="best")

        if visuParam[3]:
            fig4, ax4 = plt.subplots() 
            ax4.plot(self.convection_sol.z_mesh, self.convection_sol.P[-1], label="Pressure")
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Pressure in Pa")
            ax4.set_title("Pressure distribution in coolant canal")

        if visuParam[4]:
            fig5, ax5 = plt.subplots()
            ax5.plot(self.convection_sol.z_mesh, self.convection_sol.U[-1], label="Enthalpy")
            ax5.set_xlabel("Axial position in m")
            ax5.set_ylabel("Velocity in m/s")
            ax5.set_title("Velocity distribution in coolant canal")

        plt.show()
        return
    

class plotting:
    def __init__(self, caseList):
        self.caseList = caseList

    def writeResults(self, filename):
    # Create a new Excel workbook
        temperature_cases, void_fraction_cases, density_cases, pressure, velocity, enthalpy_cases = [self.caseList[0].convection_sol.z_mesh],[self.caseList[0].convection_sol.z_mesh],[self.caseList[0].convection_sol.z_mesh],[self.caseList[0].convection_sol.z_mesh],[self.caseList[0].convection_sol.z_mesh],[self.caseList[0].convection_sol.z_mesh]
        for case in self.caseList:
            temperature_cases.append(list(case.convection_sol.T_water))
            void_fraction_cases.append(list(case.convection_sol.voidFraction[-1]))
            density_cases.append(list(case.convection_sol.rho[-1]))
            pressure.append(list(case.convection_sol.P[-1]))
            velocity.append(list(case.convection_sol.U[-1]))
            enthalpy_cases.append(list(case.convection_sol.H[-1]))
            parameters = [['waterRadius', 'fuelRadius','gapRadius','cladRadius','height','pOutlet', 'tInlet', 'u_inlet', 'qFlow', 'Iz1', 'qFiss', 'Dz', 'frfaccorel', 'P2Pcorel', 'voidFractionCorrel', 'numericalMethod']]
            for case in self.caseList:
                parameters.append([case.r_w, case.r_f, case.gap_r, case.clad_r, case.Lf, case.pOutlet, case.tInlet, case.uInlet, case.qFlow, case.I_z, case.Qfiss, case.convection_sol.DV, case.frfaccorel, case.P2Pcorel, case.voidFractionCorrel, case.numericalMethod])

        print(len(temperature_cases))

        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            #Write parameters for each case to a sheet
            parameters_df = pd.DataFrame(parameters).T
            parameters_df.columns = ['Parameter'] + [f"Case {i}" for i in range(len(self.caseList))]
            parameters_df.to_excel(writer, sheet_name='Parameters', index=False)

            #Write temperature data to a sheet
            temperature_df = pd.DataFrame(temperature_cases).T
            temperature_df.columns = ['Axial position'] + [f"Temperature (K) Case {i}" for i in range(1,len(temperature_cases))]
            temperature_df.to_excel(writer, sheet_name='Temperature', index=False)

            #Write void fraction data to a sheet
            void_fraction_df = pd.DataFrame(void_fraction_cases).T
            void_fraction_df.columns =  ['Axial position']+[f"Void Fraction Case {i}" for i in range(1,len(void_fraction_cases))]
            void_fraction_df.to_excel(writer, sheet_name='Void Fraction', index=False)

            #Write density data to a sheet
            density_df = pd.DataFrame(density_cases).T
            density_df.columns =  ['Axial position']+[f"Density (kg/m3) Case {i}" for i in range(1,len(density_cases))]
            density_df.to_excel(writer, sheet_name='Density', index=False)

            #Write pressure data to a sheet
            pressure_df = pd.DataFrame(pressure).T
            pressure_df.columns = ['Axial position']+[f"Pressure (Pa) Case {i}" for i in range(1,len(pressure))]
            pressure_df.to_excel(writer, sheet_name='Pressure', index=False)

            #Write velocity data to a sheet
            velocity_df = pd.DataFrame(velocity).T
            velocity_df.columns =  ['Axial position']+[f"Velocity (m/s) Case {i}" for i in range(1,len(velocity))]
            velocity_df.to_excel(writer, sheet_name='Velocity', index=False)

            #Write enthalpy data to a sheet
            enthalpy_df = pd.DataFrame(enthalpy_cases).T
            enthalpy_df.columns = ['Axial position']+[f"Enthalpy (kg/m3) Case {i+1}" for i in range(1,len(enthalpy_cases))]
            enthalpy_df.to_excel(writer, sheet_name='Enthalpy', index=False)
            
    def compute_error(self, GenFoamPathCase, compParam, genFoamVolumeFraction):

        # Read the Excel file
        df = pd.read_excel(GenFoamPathCase)

        # Create empty lists for each column
        columns = df.columns.tolist()
        data = [[] for _ in columns]

        # Iterate over each row and append values to the corresponding list
        for index, row in df.iterrows():
            for i, col in enumerate(columns):
                data[i].append(row[col])

        for i in range(len(data[7])):
            data[7][i] = (1/(1-genFoamVolumeFraction)) * data[7][i]

        genfoamCASE = [data[0], data[3], data[7], data[3], data[1], data[5]]
        print(genfoamCASE[2])

        Tw_error, voidFraction_error, pressure_error, velocity_error, = [], [], [], []
        
        for i in range(len(self.caseList)):
            Tw_error.append([])
            voidFraction_error.append([])
            pressure_error.append([])
            velocity_error.append([])
            for j in range(len(self.caseList[i].convection_sol.z_mesh)):
                jGF = int(j*len(genfoamCASE[0])/len(self.caseList[i].convection_sol.z_mesh))
                Tw_error[i].append(100*abs(self.caseList[i].convection_sol.T_water[j] - genfoamCASE[1][jGF])/genfoamCASE[1][jGF])
                voidFraction_error[i].append(100*abs(self.caseList[i].convection_sol.voidFraction[-1][j] - genfoamCASE[2][jGF]))
                pressure_error[i].append(100*abs(self.caseList[i].convection_sol.P[-1][j] - genfoamCASE[4][jGF])/genfoamCASE[4][jGF])
                velocity_error[i].append(100*abs(self.caseList[i].convection_sol.U[-1][j] - genfoamCASE[5][jGF])/genfoamCASE[5][jGF])

        Tw_error_RMS, voidFraction_error_RMS, pressure_error_RMS, velocity_error_RMS, = [], [], [], []

        for i in range(len(self.caseList)):
            Tw_error_RMS.append([])
            voidFraction_error_RMS.append([])
            pressure_error_RMS.append([])
            velocity_error_RMS.append([])
            for j in range(len(self.caseList[i].convection_sol.z_mesh)):
                jGF = int(j*len(genfoamCASE[0])/len(self.caseList[i].convection_sol.z_mesh))
                Tw_error_RMS[i].append(self.caseList[i].convection_sol.T_water[j] - genfoamCASE[1][jGF])
                voidFraction_error_RMS[i].append(self.caseList[i].convection_sol.voidFraction[-1][j] - genfoamCASE[2][jGF])
                pressure_error_RMS[i].append(self.caseList[i].convection_sol.P[-1][j] - genfoamCASE[4][jGF])
                velocity_error_RMS[i].append(self.caseList[i].convection_sol.U[-1][j] - genfoamCASE[5][jGF])

        if compParam == 'numericalMethod':
            models = [self.caseList[i].numericalMethod for i in range(len(self.caseList))]
            fig1, ax1 = plt.subplots()
            for i in range(len(self.caseList)):
                ax1.plot(self.caseList[0].convection_sol.z_mesh, Tw_error[i], label=self.caseList[i].numericalMethod)
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Error (%)")
            ax1.set_title("Error in temperature distribution")
            ax1.legend(loc="best")

            fig2, ax2 = plt.subplots()
            for i in range(len(self.caseList)):
                ax2.plot(self.caseList[0].convection_sol.z_mesh, voidFraction_error[i], label=self.caseList[i].numericalMethod)
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Error (%)")
            ax2.set_title("Error in void fraction distribution")
            ax2.legend(loc="best")

            fig3, ax3 = plt.subplots()
            for i in range(len(self.caseList)):
                ax3.plot(self.caseList[0].convection_sol.z_mesh, pressure_error[i], label=self.caseList[i].numericalMethod)
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Error (%)")
            ax3.set_title("Error in pressure distribution")
            ax3.legend(loc="best")

            fig4, ax4 = plt.subplots()
            for i in range(len(self.caseList)):
                ax4.plot(self.caseList[0].convection_sol.z_mesh, velocity_error[i], label=self.caseList[i].numericalMethod)
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Error (%)")
            ax4.set_title("Error in velocity distribution")
            ax4.legend(loc="best")

            plt.show()

        elif compParam == 'voidFractionCorrel':
            models = [self.caseList[i].voidFractionCorrel for i in range(len(self.caseList))]
            fig1, ax1 = plt.subplots()
            for i in range(len(self.caseList)):
                ax1.plot(self.caseList[0].convection_sol.z_mesh, Tw_error[i], label=self.caseList[i].voidFractionCorrel)
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Error (%)")
            ax1.set_title("Error in temperature distribution")
            ax1.legend(loc="best")

            fig2, ax2 = plt.subplots()
            for i in range(len(self.caseList)):
                ax2.plot(self.caseList[0].convection_sol.z_mesh, voidFraction_error[i], label=self.caseList[i].voidFractionCorrel)
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Error (%)")
            ax2.set_title("Error in void fraction distribution")
            ax2.legend(loc="best")

            fig3, ax3 = plt.subplots()
            for i in range(len(self.caseList)):
                ax3.plot(self.caseList[0].convection_sol.z_mesh, pressure_error[i], label=self.caseList[i].voidFractionCorrel)
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Error (%)")
            ax3.set_title("Error in pressure distribution")
            ax3.legend(loc="best")

            fig4, ax4 = plt.subplots()
            for i in range(len(self.caseList)):
                ax4.plot(self.caseList[0].convection_sol.z_mesh, velocity_error[i], label=self.caseList[i].voidFractionCorrel)
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Error (%)")
            ax4.set_title("Error in velocity distribution")
            ax4.legend(loc="best")

            plt.show()

        elif compParam == 'frfaccorel':
            models = [self.caseList[i].frfaccorel for i in range(len(self.caseList))]
            fig1, ax1 = plt.subplots()
            for i in range(len(self.caseList)):
                ax1.plot(self.caseList[0].convection_sol.z_mesh, Tw_error[i], label=self.caseList[i].frfaccorel)
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Error (%)")
            ax1.set_title("Error in temperature distribution")
            ax1.legend(loc="best")

            fig2, ax2 = plt.subplots()
            for i in range(len(self.caseList)):
                ax2.plot(self.caseList[0].convection_sol.z_mesh, voidFraction_error[i], label=self.caseList[i].frfaccorel)
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Error (%)")
            ax2.set_title("Error in void fraction distribution")
            ax2.legend(loc="best")

            fig3, ax3 = plt.subplots()
            for i in range(len(self.caseList)):
                ax3.plot(self.caseList[0].convection_sol.z_mesh, pressure_error[i], label=self.caseList[i].frfaccorel)
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Error (%)")
            ax3.set_title("Error in pressure distribution")
            ax3.legend(loc="best")

            fig4, ax4 = plt.subplots()
            for i in range(len(self.caseList)):
                ax4.plot(self.caseList[0].convection_sol.z_mesh, velocity_error[i], label=self.caseList[i].frfaccorel)
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Error (%)")
            ax4.set_title("Error in velocity distribution")
            ax4.legend(loc="best")

            plt.show()
        
        elif compParam == 'nCells':
            models = [self.caseList[i].convection_sol.nCells for i in range(len(self.caseList))]
            fig1, ax1 = plt.subplots()
            for i in range(len(self.caseList)):
                ax1.plot(self.caseList[0].convection_sol.z_mesh, Tw_error[i], label=self.caseList[i].convection_sol.nCells)
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Error (%)")
            ax1.set_title("Error in temperature distribution")
            ax1.legend(loc="best")

            fig2, ax2 = plt.subplots()
            for i in range(len(self.caseList)):
                ax2.plot(self.caseList[0].convection_sol.z_mesh, voidFraction_error[i], label=self.caseList[i].convection_sol.nCells)
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Error (%)")
            ax2.set_title("Error in void fraction distribution")
            ax2.legend(loc="best")

            fig3, ax3 = plt.subplots()
            for i in range(len(self.caseList)):
                ax3.plot(self.caseList[0].convection_sol.z_mesh, pressure_error[i], label=self.caseList[i].convection_sol.nCells)
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Error (%)")
            ax3.set_title("Error in pressure distribution")
            ax3.legend(loc="best")

            fig4, ax4 = plt.subplots()
            for i in range(len(self.caseList)):
                ax4.plot(self.caseList[0].convection_sol.z_mesh, velocity_error[i], label=self.caseList[i].convection_sol.nCells)
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Error (%)")
            ax4.set_title("Error in velocity distribution")
            ax4.legend(loc="best")

            plt.show()

        elif compParam == 'P2Pcorrel':
            models = [self.caseList[i].P2Pcorel for i in range(len(self.caseList))]
            fig1, ax1 = plt.subplots()
            for i in range(len(self.caseList)):
                ax1.plot(self.caseList[0].convection_sol.z_mesh, Tw_error[i], label=self.caseList[i].P2Pcorel)
            ax1.set_xlabel("Axial position in m")
            ax1.set_ylabel("Error (%)")
            ax1.set_title("Error in temperature distribution")
            ax1.legend(loc="best")

            fig2, ax2 = plt.subplots()
            for i in range(len(self.caseList)):
                ax2.plot(self.caseList[0].convection_sol.z_mesh, voidFraction_error[i], label=self.caseList[i].P2Pcorel)
            ax2.set_xlabel("Axial position in m")
            ax2.set_ylabel("Error (%)")
            ax2.set_title("Error in void fraction distribution")
            ax2.legend(loc="best")

            fig3, ax3 = plt.subplots()
            for i in range(len(self.caseList)):
                ax3.plot(self.caseList[0].convection_sol.z_mesh, pressure_error[i], label=self.caseList[i].P2Pcorel)
            ax3.set_xlabel("Axial position in m")
            ax3.set_ylabel("Error (%)")
            ax3.set_title("Error in pressure distribution")
            ax3.legend(loc="best")

            fig4, ax4 = plt.subplots()
            for i in range(len(self.caseList)):
                ax4.plot(self.caseList[0].convection_sol.z_mesh, velocity_error[i], label=self.caseList[i].P2Pcorel)
            ax4.set_xlabel("Axial position in m")
            ax4.set_ylabel("Error (%)")
            ax4.set_title("Error in velocity distribution")
            ax4.legend(loc="best")

            plt.show()

        
        temperature_errors = {
                'mean': [],
                'min': [],
                'max': []
            }
        pressure_errors = {
                'mean': [],
                'min': [],
                'max': []
            }
        velocity_errors = {
                'mean': [],
                'min': [],
                'max': []
            }
        voidFraction_errors = {
                'mean': [],
                'min': [],
                'max': []
            }
        

        for i in range(len(self.caseList)):

            print(f'Case {i}')
            voidFraction_error[i] = self.cleanList(voidFraction_error[i])
            Tw_error[i] = self.cleanList(Tw_error[i])
            pressure_error[i] = self.cleanList(pressure_error[i])
            velocity_error[i] = self.cleanList(velocity_error[i])
            voidFraction_error_RMS[i] = self.cleanList(voidFraction_error_RMS[i])
            Tw_error_RMS[i] = self.cleanList(Tw_error_RMS[i])
            pressure_error_RMS[i] = self.cleanList(pressure_error_RMS[i])
            velocity_error_RMS[i] = self.cleanList(velocity_error_RMS[i])

            meanVoid = 0
            meanTemp = 0
            meanPressure = 0
            meanVelocity = 0

            for j in range(len(voidFraction_error[i])):
                meanVoid += (voidFraction_error[i][j]**2)
            for j in range(len(Tw_error[i])):
                meanTemp += (Tw_error[i][j]**2)
            for j in range(len(pressure_error[i])):
                meanPressure += (pressure_error[i][j]**2)
            for j in range(len(velocity_error[i])):
                meanVelocity += (velocity_error[i][j]**2)

            meanVoid = np.sqrt(meanVoid/len(voidFraction_error))
            meanTemp = np.sqrt(meanTemp/len(Tw_error))
            meanPressure = np.sqrt(meanPressure/len(pressure_error))
            meanVelocity = np.sqrt(meanVelocity/len(velocity_error))

            print(f'Case {i}')
            print(f'mean voidFraction: {voidFraction_error[i]}')
            print(f'Void Fraction error moyenne: {meanVoid}, erreur max: {np.max(voidFraction_error[i])}, erreur min: {np.min(voidFraction_error[i])}')
            print(f'Temperature error moyenne: {meanTemp}, erreur max: {np.max(Tw_error[i])}, erreur min: {np.min(Tw_error[i])}')
            print(f'Pressure error moyenne: {meanPressure}, erreur max: {np.max(pressure_error[i])}, erreur min: {np.min(pressure_error[i])}')
            print(f'Velocity error moyenne: {meanVelocity}, erreur max: {np.max(velocity_error[i])}, erreur min: {np.min(velocity_error[i])}')

            voidFraction_errors['mean'].append(np.mean(voidFraction_error[i]))
            voidFraction_errors['max'].append(np.max(voidFraction_error[i]))
            voidFraction_errors['min'].append(np.min(voidFraction_error[i]))
            temperature_errors['mean'].append(np.mean(Tw_error[i]))
            temperature_errors['max'].append(np.max(Tw_error[i]))
            temperature_errors['min'].append(np.min(Tw_error[i]))
            pressure_errors['mean'].append(np.mean(pressure_error[i]))
            pressure_errors['max'].append(np.max(pressure_error[i]))
            pressure_errors['min'].append(np.min(pressure_error[i]))
            velocity_errors['mean'].append(np.mean(velocity_error[i]))
            velocity_errors['max'].append(np.max(velocity_error[i]))
            velocity_errors['min'].append(np.min(velocity_error[i]))

        # Création des graphiques pour pression, température et vitesse
        self.plot_error_graph(models, voidFraction_errors, temperature_errors, pressure_errors, velocity_errors)
        plt.show()

    def cleanList(self, data):
        # Convert to numpy array if it's not already
        data = np.array(data)
        # Filtrer les NaN et les valeurs infinies
        cleaned_data = data[np.isfinite(data)]
        return cleaned_data

    def plotSimple(self):
        fig1, ax1 = plt.subplots()
        for i in range(len(self.caseList)):
            ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].voidFractionCorrel)
        ax1.set_xlabel("Axial position in m")
        ax1.set_ylabel("Temperature in K")
        ax1.set_title("Temperature distribution in pincell")
        ax1.legend(loc="best")


        fig2, ax2 = plt.subplots()
        for i in range(len(self.caseList)):
            ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].voidFractionCorrel)
        ax2.set_xlabel("Axial position in m")
        ax2.set_ylabel("Void fraction")
        ax2.set_title("Void fraction distribution in coolant canal")
        ax2.legend(loc="best")


        fig3, ax3 = plt.subplots()
        for i in range(len(self.caseList)):
            ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].voidFractionCorrel)
        ax3.set_xlabel("Axial position in m")
        ax3.set_ylabel("Density in kg/m^3")
        ax3.set_title("Density distribution in coolant canal")
        ax3.legend(loc="best")


        fig4, ax4 = plt.subplots()
        for i in range(len(self.caseList)):
            ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].voidFractionCorrel)
        ax4.set_xlabel("Axial position in m")
        ax4.set_ylabel("Pressure in Pa")
        ax4.set_title("Pressure distribution in coolant canal")
        ax4.legend(loc="best")
        plt.show()

        fig5, ax5 = plt.subplots()
        for i in range(len(self.caseList)):
            ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].voidFractionCorrel)
        ax5.set_xlabel("Axial position in m")
        ax5.set_ylabel("Velocity in m/s")
        ax5.set_title("Velocity distribution in coolant canal")
        ax5.legend(loc="best")
        plt.show()


    # Fonction pour créer des graphiques pour les différentes variables
    def plot_error_graph(self, models, void_fraction_errors, temperature_errors, pressure_errors, velocity_errors):
        fig, ax = plt.subplots()

        # Positions pour chaque groupe de barres
        width = 0.2  # Largeur des barres
        x = np.arange(len(models))  # Positions des modèles

        # Barres pour les erreurs de fraction de vide
        bars_void_fraction = ax.bar(x - 1.5*width, void_fraction_errors['mean'], width, 
                                    yerr=[void_fraction_errors['mean'], [void_fraction_errors['max'][i] - void_fraction_errors['mean'][i] for i in range(len(models))]],
                                    capsize=5, label='Fraction de vide', color='skyblue')

        # Barres pour les erreurs de température
        bars_temperature = ax.bar(x - 0.5*width, temperature_errors['mean'], width, 
                                yerr=[temperature_errors['mean'], [temperature_errors['max'][i] - temperature_errors['mean'][i] for i in range(len(models))]],
                                capsize=5, label='Température', color='lightcoral')

        # Barres pour les erreurs de pression
        bars_pressure = ax.bar(x + 0.5*width, pressure_errors['mean'], width, 
                            yerr=[pressure_errors['mean'], [pressure_errors['max'][i] - pressure_errors['mean'][i] for i in range(len(models))]],
                            capsize=5, label='Pression', color='lightgreen')

        """ # Barres pour les erreurs de vitesse
        bars_velocity = ax.bar(x + 1.5*width, velocity_errors['mean'], width, 
                            yerr=[velocity_errors['mean'], [velocity_errors['max'][i] - velocity_errors['mean'][i] for i in range(len(models))]],
                            capsize=5, label='Vitesse', color='orange') """

        # Ajout des labels et titre
        ax.set_xlabel('Modèles')
        ax.set_ylabel('Erreurs mean/max (%)')
        ax.set_title('Comparaison des erreurs de fraction de vide, température, pression et vitesse')
        ax.set_xticks(x)
        ax.set_xticklabels(models)
        ax.legend(loc="upper left")

        # Affichage du graphique
        plt.tight_layout()

    def plotComparison(self, compParam, visuParam):
        if compParam == 'voidFractionCorrel':
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].voidFractionCorrel)
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title("Temperature distribution in pincell")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].voidFractionCorrel)
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title("Void fraction distribution in coolant canal")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].voidFractionCorrel)
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title("Density distribution in coolant canal")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].voidFractionCorrel)
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title("Pressure distribution in coolant canal")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].voidFractionCorrel)
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title("Velocity distribution in coolant canal")
                ax5.legend(loc="best")

            fig6, ax6 = plt.subplots()
            for i in range(len(self.caseList)):
                ax6.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.get_Fission_Power(), label="Fission power")
            ax6.set_xlabel("Axial position in m")
            ax6.set_ylabel("Fission power in W/m^3")
            ax6.set_title("Fission power distribution in fuel rod")

            
            fig7, ax7 = plt.subplots()
            for i in range(len(self.caseList)):
                ax7.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.H[-1], label=self.caseList[i].voidFractionCorrel)
            ax7.set_xlabel("Axial position in m")
            ax7.set_ylabel("Enthalpy in K")
            ax7.set_title("Enthalpy distribution in pincell")
            ax7.legend(loc="best")

            plt.show()
    
        elif compParam == 'frfaccorel':
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].frfaccorel)
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title("Temperature distribution in pincell")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].frfaccorel)
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title("Void fraction distribution in coolant canal")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].frfaccorel)
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title("Density distribution in coolant canal")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].frfaccorel)
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title("Pressure distribution in coolant canal")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].frfaccorel)
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title("Velocity distribution in coolant canal")
                ax5.legend(loc="best")

            fig6, ax6 = plt.subplots()
            for i in range(len(self.caseList)):
                ax6.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.get_Fission_Power(), label="Fission power")
            ax6.set_xlabel("Axial position in m")
            ax6.set_ylabel("Fission power in W/m^3")
            ax6.set_title("Fission power distribution in fuel rod")

            
            fig7, ax7 = plt.subplots()
            for i in range(len(self.caseList)):
                ax7.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.H[-1], label=self.caseList[i].frfaccorel)
            ax7.set_xlabel("Axial position in m")
            ax7.set_ylabel("Enthalpy in kJ/kg")
            ax7.set_title("Enthalpy distribution in pincell")
            ax7.legend(loc="best")

            plt.show()
    
        elif compParam == 'P2Pcorrel':
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].P2Pcorel)
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title("Temperature distribution in pincell")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].P2Pcorel)
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title("Void fraction distribution in coolant canal")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].P2Pcorel)
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title("Density distribution in coolant canal")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].P2Pcorel)
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title("Pressure distribution in coolant canal")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].P2Pcorel)
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title("Velocity distribution in coolant canal")
                ax5.legend(loc="best")

            fig6, ax6 = plt.subplots()
            for i in range(len(self.caseList)):
                ax6.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.get_Fission_Power(), label="Fission power")
            ax6.set_xlabel("Axial position in m")
            ax6.set_ylabel("Fission power in W/m^3")
            ax6.set_title("Fission power distribution in fuel rod")

            
            fig7, ax7 = plt.subplots()
            for i in range(len(self.caseList)):
                ax7.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.H[-1], label=self.caseList[i].P2Pcorel)
            ax7.set_xlabel("Axial position in m")
            ax7.set_ylabel("Enthalpy in K")
            ax7.set_title("Enthalpy distribution in pincell")
            ax7.legend(loc="best")

            plt.show()

        elif compParam == 'numericalMethod':
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    print(f'Twater: { self.caseList[i].convection_sol.T_water}')
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].numericalMethod)
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title("Temperature distribution in pincell")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].numericalMethod)
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title("Void fraction distribution in coolant canal")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].numericalMethod)
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title("Density distribution in coolant canal")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].numericalMethod)
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title("Pressure distribution in coolant canal")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].numericalMethod)
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title("Velocity distribution in coolant canal")
                ax5.legend(loc="best")

            if visuParam[5]:
                fig6, ax6 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax6.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.get_Fission_Power(), label="Fission power")
                ax6.set_xlabel("Axial position in m")
                ax6.set_ylabel("Fission power in W/m^3")
                ax6.set_title("Fission power distribution in fuel rod")

            if visuParam[6]:
                fig7, ax7 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax7.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.H[-1], label=self.caseList[i].numericalMethod)
                ax7.set_xlabel("Axial position in m")
                ax7.set_ylabel("Enthalpy in K")
                ax7.set_title("Enthalpy distribution in pincell")
                ax7.legend(loc="best")

            plt.show()

    
    def GenFoamComp(self, GenFoamPathCase, compParam, visuParam, genFoamVolumeFraction):

        # Read the Excel file
        df = pd.read_excel(GenFoamPathCase)

        # Create empty lists for each column
        columns = df.columns.tolist()
        data = [[] for _ in columns]

        # Iterate over each row and append values to the corresponding list
        for index, row in df.iterrows():
            for i, col in enumerate(columns):
                data[i].append(row[col])

        for i in range(len(data[7])):
            data[7][i] = (1/(1-genFoamVolumeFraction)) * data[7][i]

        genfoamCASE = [data[0], data[3], data[7], data[4], data[1], data[5]]

        if compParam == 'voidFractionCorrel':
            title = f"Methode numérique: {self.caseList[0].numericalMethod}, \n Correlation multiplicateur biphasique: {self.caseList[0].convection_sol.P2Pcorel}, \n Correlation facteur de friction: {self.caseList[0].convection_sol.frfaccorel}"
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].voidFractionCorrel)
                ax1.step(genfoamCASE[0], genfoamCASE[1], label="GenFoam")
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.grid()
                ax1.set_title(title)
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].voidFractionCorrel)
                ax2.step(genfoamCASE[0], genfoamCASE[2], label="GenFoam")
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.grid()
                ax2.set_title(title)
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].voidFractionCorrel)
                ax3.plot(genfoamCASE[0], genfoamCASE[3], label="GenFoam")
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title(title)
                ax3.grid()
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].voidFractionCorrel)
                ax4.step(genfoamCASE[0], genfoamCASE[4], label="GenFoam")
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title(title)
                ax4.grid()
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].voidFractionCorrel)
                ax5.step(genfoamCASE[0], genfoamCASE[5], label="GenFoam")
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title(title)
                ax5.grid()
                ax5.legend(loc="best")

            plt.show()

        if compParam == 'frfaccorel':
            title = f"Methode numérique: {self.caseList[0].numericalMethod}, \n Correlation multiplicateur biphasique: {self.caseList[0].convection_sol.P2Pcorel}, \n Correlation void fraction: {self.caseList[0].convection_sol.voidFractionCorrel}"
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].frfaccorel)
                ax1.plot(genfoamCASE[0], genfoamCASE[1], label="GenFoam")
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title(f"{title}")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].frfaccorel)
                ax2.plot(genfoamCASE[0], genfoamCASE[2], label="GenFoam")
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title(f"{title}")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].frfaccorel)
                #ax3.plot(genfoamCASE[0], genfoamCASE[3], label="GenFoam")
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title(f"{title}")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].frfaccorel)
                ax4.plot(genfoamCASE[0], genfoamCASE[4], label="GenFoam")
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title(f"{title}")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].frfaccorel)
                ax5.plot(genfoamCASE[0], genfoamCASE[5], label="GenFoam")
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title(f"{title}")
                ax5.legend(loc="best")

            fig7, ax7 = plt.subplots()
            for i in range(len(self.caseList)):
                ax7.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.H[-1], label=self.caseList[i].frfaccorel)
            ax7.set_xlabel("Axial position in m")
            ax7.set_ylabel("Enthalpy in kJ/kg")
            ax7.set_title("Enthalpy distribution in pincell")
            ax7.legend(loc="best")
            
            

            plt.show()

        if compParam == 'numericalMethod':
            title = f"Correlation fracteur de friction: {self.caseList[0].convection_sol.frfaccorel}, \n Correlation multiplicateur biphasique: {self.caseList[0].convection_sol.P2Pcorel}, \n Correlation void fraction: {self.caseList[0].convection_sol.voidFractionCorrel}"
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].numericalMethod)
                ax1.plot(genfoamCASE[0], genfoamCASE[1], label="GenFoam")
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title(f"{title}")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].numericalMethod)
                ax2.plot(genfoamCASE[0], genfoamCASE[2], label="GenFoam")
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title(f"{title}")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].numericalMethod)
                #ax3.plot(genfoamCASE[0], genfoamCASE[3], label="GenFoam")
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title(f"{title}")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].numericalMethod)
                ax4.plot(genfoamCASE[0], genfoamCASE[4], label="GenFoam")
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title(f"{title}")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].numericalMethod)
                ax5.plot(genfoamCASE[0], genfoamCASE[5], label="GenFoam")
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title(f"{title}")
                ax5.legend(loc="best")

            plt.show()

        if compParam == 'P2Pcorrel':
            title = f"Correlation fracteur de friction: {self.caseList[0].convection_sol.frfaccorel}, \n Méthode numérique: {self.caseList[0].convection_sol.numericalMethod}, \n Correlation void fraction: {self.caseList[0].convection_sol.voidFractionCorrel}"
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].P2Pcorel)
                ax1.step(genfoamCASE[0], genfoamCASE[1], label="GenFoam")
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title(f"{title}")
                ax1.grid()
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].P2Pcorel)
                ax2.step(genfoamCASE[0], genfoamCASE[2], label="GenFoam")
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title(f"{title}")
                ax2.grid()
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].P2Pcorel)
                #ax3.plot(genfoamCASE[0], genfoamCASE[3], label="GenFoam")
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title(f"{title}")
                ax3.grid()
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].P2Pcorel)
                ax4.step(genfoamCASE[0], genfoamCASE[4], label="GenFoam")
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.grid()
                ax4.set_title(f"{title}")

                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.step(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].P2Pcorel)
                ax5.step(genfoamCASE[0], genfoamCASE[5], label="GenFoam")
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title(f"{title}")
                ax5.grid()
                ax5.legend(loc="best")

            plt.show()

        if compParam == 'nCells':
            title = f"Correlation fracteur de friction: {self.caseList[0].convection_sol.frfaccorel}, \n Méthode numérique: {self.caseList[0].convection_sol.numericalMethod}, \n Correlation void fraction: {self.caseList[0].convection_sol.voidFractionCorrel}"
            if visuParam[0]:
                fig1, ax1 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax1.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.T_water, label=self.caseList[i].convection_sol.nCells)
                ax1.plot(genfoamCASE[0], genfoamCASE[1], label="GenFoam")
                ax1.set_xlabel("Axial position in m")
                ax1.set_ylabel("Temperature in K")
                ax1.set_title(f"{title}")
                ax1.legend(loc="best")

            if visuParam[1]:
                fig2, ax2 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax2.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.voidFraction[-1], label=self.caseList[i].convection_sol.nCells)
                ax2.plot(genfoamCASE[0], genfoamCASE[2], label="GenFoam")
                ax2.set_xlabel("Axial position in m")
                ax2.set_ylabel("Void fraction")
                ax2.set_title(f"{title}")
                ax2.legend(loc="best")

            if visuParam[2]:
                fig3, ax3 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax3.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.rho[-1], label=self.caseList[i].convection_sol.nCells)
                #ax3.plot(genfoamCASE[0], genfoamCASE[3], label="GenFoam")
                ax3.set_xlabel("Axial position in m")
                ax3.set_ylabel("Density in kg/m^3")
                ax3.set_title(f"{title}")
                ax3.legend(loc="best")

            if visuParam[3]:
                fig4, ax4 = plt.subplots() 
                for i in range(len(self.caseList)):
                    ax4.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.P[-1], label=self.caseList[i].convection_sol.nCells)
                ax4.plot(genfoamCASE[0], genfoamCASE[4], label="GenFoam")
                ax4.set_xlabel("Axial position in m")
                ax4.set_ylabel("Pressure in Pa")
                ax4.set_title(f"{title}")
                ax4.legend(loc="best")

            if visuParam[4]:
                fig5, ax5 = plt.subplots()
                for i in range(len(self.caseList)):
                    ax5.plot(self.caseList[i].convection_sol.z_mesh, self.caseList[i].convection_sol.U[-1], label=self.caseList[i].convection_sol.nCells)
                ax5.plot(genfoamCASE[0], genfoamCASE[5], label="GenFoam")
                ax5.set_xlabel("Axial position in m")
                ax5.set_ylabel("Velocity in m/s")
                ax5.set_title(f"{title}")
                ax5.legend(loc="best")

            plt.show()

        