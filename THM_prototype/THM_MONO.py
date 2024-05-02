# Python3 class part of THM_prototype
# uses : - setting up mesh centered finite difference solution to heat conduction in fuel rod, FDM_ConductionInFuelPin
#        - setting up finite volume discretization for solving heat convection in coolant, FVM_ConvectionInCanal
# Author : R. Guasch
# technical documentation : "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018
# document available at http://merlin.polymtl.ca/downloads/thm.pdf

import numpy as np
from iapws import IAPWS97


class FDM_HeatConductionInFuelPin:
    def __init__(self, r_fuel, I_f, gap_r, clad_r, I_c, Qfiss, kf, kc, Hgap, z, T_surf):
        # Physical prameters
        self.r_f = r_fuel # fuel pin radius in meters
        self.I_f = I_f # number of mesh elements in the fuel
        self.gap_r = gap_r # gap radius in meters, used to determine mesh elements for constant surface discretization
        self.clad_r = clad_r # clad radius in meters, used to determine mesh elements for constant surface discretization

        self.I_c = I_c # number of mesh elements in clad
        self.Qfiss = Qfiss # Fission power density in W/m^3
        self.Hg = Hgap # Heat transfer coefficient through gap W/m^2/K
        self.z = z # corresponding height in m corresponding to axial discretization used in FVM_ConvectionInCanal class
        self.T_surf = T_surf # Boundary condition outer clad surface temperature computed from FVM_ConvectionInCanal class 
        self.kf = kf
        self.kc = kc
        # compute relevant quantities to initialise object
        self.N_node = I_f + I_c +2
        self.A = np.eye(self.N_node)
        self.D = np.zeros(self.N_node)
        self.compute_Area_meshes()
        self.compute_radii()
        self.initialize_ks()

        self.initialise_plotting_mesh("m")
        self.physical_regions_bounds = [0, self.r_f, self.gap_r, self.clad_r]
        
        
    def initialize_ks(self):
        # this array is probaby not needed here as one might assume that kf and kc
        # are constant in fuel/clad. I wanted to keep and option to let them vary according to temperature at node as conductive properties might differ when high temperature gradients are present.
        self.k = np.zeros(self.N_node-1) # associate a k to the center of each mesh element
        for i in range(len(self.k)):
            if i <self.I_f:
                self.k[i]=self.kf
            elif i == self.I_f:
                self.k[i]=0 # in gap!
            elif i<self.N_node:
                self.k[i]=self.kc
        return


    def compute_Area_meshes(self):
        """
        building necessary meshes for solving the heat conduction equation on the constant area discretization 
        """
        self.A_mesh_bounds = []
        self.A_mesh_centers = []
        A_f = self.r_f**2/2
        A_gf =self.gap_r**2/2
        A_cgf = self.clad_r**2/2
        self.deltaA_f = A_f / self.I_f # base assumption is that delta_A is constant in each region --> delta A fuel = constant in fuel, delta A clad = constant in clad and 1 delta A gap.
        for i in range(self.I_f+1):
            self.A_mesh_bounds.append(i*self.deltaA_f)
        for i in range(self.I_f):
            self.A_mesh_centers.append(i*self.deltaA_f+self.deltaA_f/2)
    
        self.deltaA_g = A_gf-A_f
        self.A_mesh_bounds.append(self.A_mesh_bounds[-1]+self.deltaA_g)
        self.A_mesh_centers.append(self.A_mesh_centers[-1]+self.deltaA_f/2+self.deltaA_g/2) # last center in fuel + half of the fuel area step to get to the last fuel bound + half of the gap area step to get to the center of the gap
        self.deltaA_c = (A_cgf-A_gf)/self.I_c
        for i in range(self.I_c):
            self.A_mesh_bounds.append(self.A_mesh_bounds[-1]+self.deltaA_c)
        for i in range(self.I_c):
            if i==0:
                self.A_mesh_centers.append(self.A_mesh_centers[-1]+self.deltaA_c/2+self.deltaA_g/2)
            else:
                self.A_mesh_centers.append(self.A_mesh_centers[-1]+self.deltaA_c)
        self.A_mesh_centers = np.array(self.A_mesh_centers)
        self.A_mesh_bounds = np.array(self.A_mesh_bounds)
        self.A_calculation_mesh = np.zeros(self.N_node)
        for i in range(self.N_node):
            if i < self.I_f:
                self.A_calculation_mesh[i] = i*self.deltaA_f + self.deltaA_f/2
            elif i == self.I_f:
                self.A_calculation_mesh[i] = A_f
            elif i == self.I_f+1:
                self.A_calculation_mesh[i] = A_f + self.deltaA_g
            elif i > self.I_f+1:
                self.A_calculation_mesh[i] = A_gf + self.deltaA_c/2 + (i-(self.I_f+2))*self.deltaA_c

        return
    
    def compute_radii(self):
        self.radii_at_centers = np.sqrt(self.A_mesh_centers*2)
        self.radii_at_bounds = np.sqrt(self.A_mesh_bounds*2)
        return
    
    def get_Di_half(self,i):
        if i > self.I_f+1:
            i=i-1
            Di_half = 4*self.A_mesh_bounds[i+1]/((self.deltaA_c/self.k[i])+(self.deltaA_c/self.k[i+1]))
        else:
            Di_half = 4*self.A_mesh_bounds[i+1]/((self.deltaA_f/self.k[i])+(self.deltaA_f/self.k[i+1]))
        return Di_half
    
    def get_Ei_fuel(self):
        Ei_half = 4*self.A_mesh_bounds[self.I_f]*self.k[self.I_f-1]/self.deltaA_f
        return Ei_half
    
    def get_Ei_gap(self):
        Ei_half = 4*self.A_mesh_bounds[self.I_f+1]*self.k[self.I_f+1]/self.deltaA_c
        return Ei_half
    
    def get_Ei_clad(self):
        Ei_half = 4*self.A_mesh_bounds[-1]*self.k[-1]/self.deltaA_c
        return Ei_half
    
    def get_Fi_gap(self):
        Fi_half = 4*self.A_mesh_bounds[self.I_f+1]*self.k[self.I_f+1]/self.deltaA_c
        return Fi_half
    
    def get_Gi(self):
        return self.Hg*self.radii_at_centers[self.I_f]
    
    def set_ADi_cond(self, i, ci, ai, bi, di):
        # create lines for the tri-diagonal entries in
        self.A[i, i-1:i+2] = [ci, ai, bi]
        self.D[i] = di
        return
    
    def set_CL_cond(self, A0, Am1, D0, Dm1):
        # I_and_half is the index for the I+1/2 element of th mesh which corresponds to the last point in the fuel.
        # conditions aux limites
        # A0 = A[0], Am1 = A[-1], A moins 1, 
        # D0 = D[0], Dm1 = D[-1], D moins 1.
        self.A[0], self.A[-1] = A0, Am1
        self.D[0], self.D[-1] = D0, Dm1
        return
    
    def solve_T_in_pin(self):
        for row in self.A:
            line = "[  "
            for elem in row:
                line+=f"{elem:.3f}   "
            line += "  ]\n"
            print(line)
        self.T_distrib = np.linalg.solve(self.A, self.D)
        self.compute_T_center()
        T_distrib_with_center = np.zeros(self.N_node+1)
        T_distrib_with_center[0] = self.T_center
        for i in range(1,self.N_node+1):
            T_distrib_with_center[i] = self.T_distrib[i-1]
        self.T_distrib = T_distrib_with_center

    def compute_T_center(self):
        # using equation (13) of "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018 to compute T_(3/2).
        T_3_2 = (self.deltaA_f*self.k[0]*self.T_distrib[0]+self.deltaA_f*self.k[1]*self.T_distrib[1])/(self.deltaA_f*self.k[0]+self.deltaA_f*self.k[1])
        # using equation (46) of "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018 to compute T_center.
        self.T_center = 2*self.T_distrib[0] - T_3_2
        return
    def compute_T_eff(self):
        # using equation (45) of "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018 to compute T_eff
        # This corresponds to the simplified correlation, the so-called Rowlands formula :
        if len(self.T_distrib) == self.N_node:
            self.T_eff = (5/9)*self.T_distrib[self.I_f] + (4/9)*self.T_center
        elif len(self.T_distrib) == self.N_node+1:
            self.T_eff = (5/9)*self.T_distrib[self.I_f+1] + (4/9)*self.T_center
        return
    
    def initialise_plotting_mesh(self, unit):
        """
        unit = "m" or "mm"
        this only affects the visualization/plotting units. However all input quantities have to be in MKS units to be consistent with the solver. 
        """
        # building plot mesh in radial units taking the central temperature into account.
        self.plot_mesh = np.zeros(self.N_node+1)
        self.plot_mesh[0] = 0
        self.plotting_units = unit
        for i in range(1,self.N_node+1):
            if unit == "m":
                self.plot_mesh[i] = np.sqrt(2*self.A_calculation_mesh[i-1])
            elif unit == "mm":
                self.plot_mesh[i] = np.sqrt(2*self.A_calculation_mesh[i-1])*1e3
    
    def extend_to_canal_visu(self, rw, Tw):
        A_w = rw**2/2
        deltA_w = A_w - self.A_mesh_bounds[-1] 
        w_center = np.sqrt(2*(self.A_mesh_bounds[-1] + deltA_w/2)) 
        self.plot_mesh = np.append(self.plot_mesh,[w_center])
        self.physical_regions_bounds.append(rw)        
        self.T_distrib = np.append(self.T_distrib, [Tw])
        self.radii_at_bounds = np.append(self.radii_at_bounds, [rw])

        
class FVM_ConvectionInCanal_MONO:
    def __init__(self, Lf, T_in, Q_flow, P_cool, I_z, canal_type, rf, rc, rw):
        """
        Lf = fuel rod length in m
        T_in = inlet water temperature K
        Q_flow = mass flux in kg/m^2/s, assumed to be constant along the axial profile.
        P_cool = coolant pressure in MPa, assumed to be constant along the axial profile.
        I_z = number of mesh elements on axial mesh
        canal_type = cylindrical or square, used to determine the cross sectional flow area in the canal and the hydraulic diameter
        rf = fuel rod radius
        rc, rw = outer clad radius (m), outer canal radius (m) if type is cylindrical, if type = square rw is the radius of inscribed circle in the square canal, ie half the square's side.

        Important note : cross sectional flow area is assumed to constant along z axis. Would be interesting to expand this to treat variations in flow area as in :
        THE MODELING OF ADVANCED BWR FUEL DESIGNS WITH THE NRC FUEL DEPLETION CODES PARCS/PATHS - A. WYSOCKI et al. August 2014
        """
        # Physical parameters
        self.Lf = Lf
        self.T_in = T_in
        self.Q_flow = Q_flow
        self.P_cool = P_cool
        self.fuel_radius = rf
        self.clad_radius = rc
        self.wall_dist = rw
        print(f"wall dist = {rw}")
        initial_water_state_z0 = IAPWS97(P=P_cool,T=T_in)
        self.h_z0 = initial_water_state_z0.h*10**3 # returs enthalpy at z0 for given (P,T), this assumes 1 phase liquid water at z=0, converted to J/kg
        print(f"first enthalpy is {self.h_z0}")
        self.canal_type = canal_type
        # Calculating mesh parameters.
        self.N_vol = I_z
        self.dz = self.Lf/self.N_vol
        self.A, self.D = np.eye(self.N_vol), np.zeros(self.N_vol)
        self.z_mesh = np.linspace(0.5*self.dz, self.Lf-0.5*self.dz, self.N_vol) # creating z mesh from volume center to volume center. 
        
        # calculation A_canal and DH depending on the canal type : idea would be to generalize to square section "CARCEL" to use results in DONJON5.
        if self.canal_type == "cylindrical":
            self.A_canal = np.pi*self.wall_dist**2 - np.pi*self.clad_radius**2
            self.P_wetted = 2*np.pi*(self.wall_dist + self.clad_radius)
            print(f"Perimeter of clad = {2*np.pi*self.clad_radius}, perimeter of canal = {2*np.pi*self.wall_dist}")
        elif self.canal_type == "square":
            self.A_canal = (2*self.wall_dist)**2-np.pi*self.clad_radius**2
            self.P_wetted = 2*(4*self.wall_dist+np.pi*self.clad_radius)
            print(f"Perimeter of clad = {2*np.pi*self.clad_radius}, perimeter of canal = {4*self.wall_dist}")
        self.DH = 4*self.A_canal / self.P_wetted # DH = 4*A/P where A = wetted cross sectional area and P is wetted perimeter = perimeter of fuel rod in contact with coolant + perimeter of canal in contact with coolant.
        print(f"The calculated A_canal is : {self.A_canal} m^2, DH is {self.DH} m")
        print(f"The wetted perimeter is {self.P_wetted} m, calculated in a {self.canal_type} type geometry.")
        


        return
    
    def set_ADi_conv(self, i, ci, ai, bi, di):
        self.A[i, i-1:i+2] = [ci, ai, bi]
        self.D[i] = di
        return
    
    def set_CL_conv(self, A0, Am1, D0, Dm1):
        self.A[0], self.A[-1] = A0, Am1
        self.D[0], self.D[-1] = D0, Dm1
        return
    def set_Fission_Power(self, amplitude, variation_type):
        """
        option to set fission source axial profile : 
        amblitude : fission power from fuel per unit volume W/m^3
        variation_type : string used to allow for constant, sinusoial or cosinusoidal axial profiles ("constant", "sine", "cosine" keywords allowed)
        """
        self.Power_profile = np.ones(self.N_vol)
        #self.q_fluid = np.ones(self.N_vol)
        if variation_type == "constant":
            self.Power_profile = amplitude*self.Power_profile
            self.q_fluid = np.pi*self.fuel_radius**2*self.Power_profile
        elif variation_type == "sine":
            for i in range(self.N_vol):
                self.Power_profile[i] = amplitude*np.sin(np.pi*self.z_mesh[i]/self.Lf)
            self.q_fluid = np.pi*self.fuel_radius**2*self.Power_profile
            print("POWER PROFILE")
            print(self.Power_profile)
        elif variation_type == "cosine":
            print("Keyword for cosine axial variation of fuel power not implemented yet")
        
        return
    
    def get_Fission_Power(self):
        """
        function to retrieve a given source term from the axial profile used to model fission power distribution in the fuel rod
        """
        return self.Power_profile
    
    def solve_h_in_canal(self):
        for row in self.A:
            line = "[  "
            for elem in row:
                line+=f"{elem:.3f}   "
            line += "  ]\n"
            print(line)
        self.h_z = np.linalg.solve(self.A, self.D)
        return self.h_z
    

        
    def compute_T_surf(self):
        self.T_surf = np.zeros(self.N_vol)
        self.Hc = np.zeros(self.N_vol)
        self.T_water = np.zeros(self.N_vol)
        self.T_water[0] = self.T_in
        for i in range(self.N_vol):
            Pr_number = IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.Prandt
            Re_number = self.Q_flow*self.DH/(IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.mu)
            k_fluid = IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.k
            print(f"At axial slice = {i}, computed Reynold # = {Re_number}, computed Prandt # = {Pr_number}, k_fluid = {k_fluid}")
            self.Hc[i] = (0.023)*(Pr_number)**0.4*(Re_number)**0.8*k_fluid/self.DH
            self.T_water[i] = IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).T
            self.T_surf[i] = (self.q_fluid[i]/(2*np.pi*self.clad_radius)/self.Hc[i]+self.T_water[i])
    
        return self.T_surf
    

    