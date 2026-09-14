"""
This program simulates radial heat conduction in a nuclear fuel pin using the finite difference method (FDM).
The code models the temperature distribution across the fuel pin, considering the fuel, gap, and cladding regions.
technical documentation : "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018
document available at http://merlin.polymtl.ca/downloads/thm.pdf
Author : R. Guasch, C. Huet

Key features:
- Discretization of the radial domain into finite difference mesh elements.
- Calculation of thermal conductivity for each region (fuel, gap, and cladding).
- Solver for the temperature distribution using a tridiagonal matrix approach.
- Boundary condition application and computation of center and effective temperatures.
- Optional extension to visualize the entire system, including the coolant channel.

This program follows standard finite difference approaches to solve the heat conduction equation in cylindrical coordinates, with assumptions of constant properties per region (fuel, gap, cladding), and fixed boundary conditions for outer surface temperature.
"""

import numpy as np

class FDM_HeatConductionInFuelPin:

    """
    Parameters:
    - `r_fuel` : Radius of the fuel pin (m).
    - `I_f`    : Number of mesh elements in the fuel.
    - `gap_r`  : Radius of the gap between fuel and cladding (m).
    - `clad_r` : Radius of the cladding (m).
    - `I_c`    : Number of mesh elements in the cladding.
    - `Qfiss`  : Fission power density in the fuel (W/m^3).
    - `kf`     : Thermal conductivity of the fuel (W/m/K).
    - `kc`     : Thermal conductivity of the cladding (W/m/K).
    - `Hgap`   : Heat transfer coefficient through the gap (W/m^2/K).
    - `z`      : Height corresponding to axial discretization (for potential coupling with axial convection models).
    - `T_surf` : Outer cladding surface temperature, obtained from convection models.

    Methods:
    - `compute_Area_meshes`: Calculates mesh boundaries and areas for each region.
    - `compute_radii`: Converts areas into radial positions for mesh centers and boundaries.
    - `solve_T_in_pin`: Solves the system of equations for the temperature distribution.
    - `compute_T_center`: Computes the fuel center temperature using established formulas.
    - `compute_T_eff`: Computes the effective fuel temperature using simplified correlations.
    - `initialize_ks`: Initializes the thermal conductivity values for each mesh region.
    - `set_ADi_cond`: Sets the tri-diagonal matrix entries for internal nodes.
    - `set_CL_cond`: Applies boundary conditions at the first and last nodes.
    - `extend_to_canal_visu`: Extends the mesh and temperature distribution for visualization, including coolant.
    - `initialise_plotting_mesh`: Prepares the radial mesh for plotting the temperature profile.
    """

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
            #print(line)
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

