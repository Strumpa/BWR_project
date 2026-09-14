# Python3 class for coupling neutron transport and thermal-hydraulics
# Author: R. Guasch
# Date: 24/09/2024
# Purpose : automate the coupling of neutron transport and thermal-hydraulics to perform a multi-physics simulation
# Usage: import in multiPhysics.py

import numpy as np
import matplotlib.pyplot as plt
import os
from THM_main import Version5_THM_prototype as THM_prototype
import lifo
import lcm
import cle2000


class coupled_scheme:
    def __init__(self, case_name, compo_name, neutron_transport_data, thermal_hydraulics_data, coupling_parameters, save_data, save_figures):
        """
        Contructor of the class coupling_scheme :

        :param case_name : string : name of the case
        :param compo_name : string : name of the component
        :param neutron_transport_data: dictionary with the neutron transport data
        :param thermal_hydraulics_data: dictionary with the thermal-hydraulics data
        :param coupling_parameters: dictionary with the coupling parameters

        :param save_data : list of fields to save for export, last list is a list of the iteration numbers to save data for
        :param save_figures : list of fields to save figures of, last list is a list of the iteration numbers to save figures for
        

        """
        self.case_name = case_name
        self.compo_name = compo_name

        self.neutron_transport_data = neutron_transport_data
        self.TH_data = thermal_hydraulics_data
        self.coupling_parameters = coupling_parameters

        self.save_data = save_data
        self.save_figures = save_figures

        # Initialize the class attributes
        self.Keffs = []
        self.Power_Distribs = []
        self.Volumic_Power_Distribs = []

        self.Residuals_Power_Distrib = []
        
        self.Teff_list = []
        self.TCool_list = []
        self.rhoCool_list = []
        self.voidFraction_list = []

        self.Residual_TFuel_list = []
        self.Residual_TCool_list = []
        self.Residual_rhoCool_list = []


        # Recover coupling scheme parameters
        self.convergence = False
        self.maxIter = self.coupling_parameters["maxIter"]
        self.tolPOW = self.coupling_parameters["tolPOW"] 
        self.tolTH = self.coupling_parameters["tolTH"]
        self.POWunderRelaxationFactor = self.coupling_parameters["POWunderRelaxationFactor"]
        self.THunderRelaxationFactor = self.coupling_parameters["THunderRelaxationFactor"]
        if self.coupling_parameters["convergenceCriterion"] == "POW":
            self.convergenceCriterion = "POWField"
        elif self.coupling_parameters["convergenceCriterion"] == "TH":
            self.convergenceCriterion = "THField"
        
        if self.POWunderRelaxationFactor == 0:
            self.relaxPOW = False
        else:
            self.relaxPOW = True
        if self.THunderRelaxationFactor == 0:
            self.relaxTH = False
        else:
            self.relaxTH = True

        #### Begin initialization of the coupled scheme
        print("$$$ --- Initialization of the coupled scheme")
        # Guess the axial power shape
        self.z_mesh, self.qFiss_init = self.guessAxialPowerShape(self.TH_data["qFiss_init"], self.TH_data["Iz1"], self.TH_data["height"], self.TH_data["fuelRadius"])

        # Initialization of TH solver
        self.InitializeTH_solver()

        # Prepare LCM data for DONJON
        self.InitializeLCM_THData()

        # Initialization of Neutron Transport solver
        self.InitializeDONJON_solver()
        #### End initialization of the coupled scheme
        print("$$$ --- Initialization of the coupled scheme completed")

        #### Begin the coupling iterations
        print("$$$ --- Start coupling iterations")
        self.multiphysics_iterations()
        

    def InitializeTH_solver(self):
        """
        Initialize the thermal-hydraulics solver
        """
        self.canalType = self.TH_data["canalType"]
        self.waterRadius = self.TH_data["waterRadius"]
        self.fuelRadius = self.TH_data["fuelRadius"]
        self.gapRadius = self.TH_data["gapRadius"]
        self.cladRadius = self.TH_data["cladRadius"]
        self.height = self.TH_data["height"]
        self.tInlet = self.TH_data["tInlet"]
        self.pOutlet = self.TH_data["pOutlet"]
        self.massFlowRate = self.TH_data["massFlowRate"]
        #qFiss_init = self.TH_data["qFiss_init"]
        self.kFuel = self.TH_data["kFuel"]
        self.Hgap = self.TH_data["Hgap"]
        self.kClad = self.TH_data["kClad"]
        self.Iz1 = self.TH_data["Iz1"]
        self.If = self.TH_data["If"]
        self.I1 = self.TH_data["I1"]
        self.zPlotting = self.TH_data["zPlotting"]
        self.solveConduction = self.TH_data["solveConduction"]
        self.frfaccorel = self.TH_data["frfaccorel"]
        self.P2Pcorel = self.TH_data["P2Pcorel"]
        self.voidFractionCorrel = self.TH_data["voidFractionCorrel"]
        self.numericalMethod = self.TH_data["numericalMethod"]


        self.THMComponent = THM_prototype("Initialization of BWR Pincell equivalent canal", self.canalType, self.waterRadius, self.fuelRadius, self.gapRadius, self.cladRadius, 
                            self.height, self.tInlet, self.pOutlet, self.massFlowRate, self.qFiss_init, self.kFuel, self.Hgap, self.kClad, self.Iz1, self.If, self.I1, self.zPlotting, 
                            self.solveConduction, dt = 0, t_tot = 0, frfaccorel = self.frfaccorel, P2Pcorel = self.P2Pcorel, voidFractionCorrel = self.voidFractionCorrel, 
                            numericalMethod = self.numericalMethod)
        
        self.TFuel_init, self.TCool_init, self.rhoCool_init = self.THMComponent.get_TH_parameters()
        self.Teff_list.append(self.TFuel_init)
        self.TCool_list.append(self.TCool_init)
        self.rhoCool_list.append(self.rhoCool_init)

        return 
    
    def InitializeLCM_THData(self):
        """
        Initialize the data for the Lumped Capacitance Model for the thermal-hydraulics solver
        """
        THData = lcm.new('LCM','THData')
        THData['TFuelList']    = np.array(self.TFuel_init, dtype='f')
        THData['TCoolList'] = np.array(self.TCool_init, dtype='f')
        THData['DCoolList'] = np.array(self.rhoCool_init/1000, dtype='f') # kg/m3 to g/cm3
        THData.close() # close without erasing the data
        self.THData = THData

        return

    def InitializeDONJON_solver(self):
        """
        Initialize the neutron transport solver DONJON calling IniDONJON.c2m
        """
        ## 3.) Initializing Neutronics solution
        # 3.1) construct the Lifo stack for IniDONJON
        ipLifo1=lifo.new()
        ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
        ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
        ipLifo1.pushEmpty("Cpo", "LCM") # Compo
        ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
        ipLifo1.push(self.THData) # Thermal Hydraulic data for initialization
        ipLifo1.push(self.compo_name) # Compo name

        # 3.2) call IniDONJON Cle-2000 procedure
        IniDONJON = cle2000.new('IniDONJON',ipLifo1,1)
        IniDONJON.exec()
        print("IniDONJON execution completed")
        # recover the output LCM objects
        self.Fmap = ipLifo1.node("Fmap")
        self.Matex = ipLifo1.node("Matex")
        self.Cpo = ipLifo1.node("Cpo")
        self.Track = ipLifo1.node("Track")
        self.InitTHData = ipLifo1.node("THData") # Recover the TH data after the initialization is this needed ? 
        #       --> THData is just used to initialize the Fuel Map object but not modified so can use initial THData object
        print(f"Initial THData object TFuel = : {self.InitTHData['TFuelList']}, TCool = {self.InitTHData['TCoolList']}, DCool = {self.InitTHData['DCoolList']}")

        # empty the Lifo stack for IniDONJON
        while ipLifo1.getMax() > 0:
            ipLifo1.pop();
        
        # Create Lifo stack for Neutronics solution
        self.ipLifo2 = lifo.new()
        self.Neutronics = cle2000.new('Neutronics',self.ipLifo2,1)
    
        return

    def multiphysics_iterations(self):
        """
        Iterate between the neutron transport and thermal-hydraulics solvers until convergence
        """
        iter = 0
        while not self.convergence:
            iter+=1
            print(f"$$$ --- Begin Iteration {iter}")
            # 1.) Compute the neutron flux/power with DONJON
            self.ComputeNeutronFlux()
            self.UpdateFluxAndPowerProfile()
            # 6.) Empty the ipLifo2 Lifo stack to prepare for the next iteration
            while self.ipLifo2.getMax() > 0:
                print("Clearing ipLifo2 stack at iter = ", iter) 
                self.ipLifo2.pop();

            # 2.) Update the thermal-hydraulics with the axial power shape
            self.ComputeTHFields()
            self.UpdateTHFields()

            # 3.) Check convergence
            self.convergence = self.CheckConvergence()
            if self.convergence:
                print(f"Convergence reached in {iter} iterations")
                break

            if iter == self.maxIter:
                print(f"Maximum number of iterations ({self.maxIter}) reached")
                break
            print(f"$$$ --- End Iteration {iter}")

    def ComputeNeutronFlux(self, iter):
        # 4.1) Neutronics solution for TH parameters obtained with initial TH solution

        ################## Neutronics part ##################
        # fill the Lifo stack for Neutronics solution
        print(f"$$ - BEGIN iter = {iter}")
        self.ipLifo2.push(self.Fmap);
        self.ipLifo2.push(self.Matex);
        if iter == 1: # At the first iteration, create empty LCM objects to host the Flux and Power fields
            print("in iter = 1")
            self.Flux = self.ipLifo2.pushEmpty("Flux", "LCM")
        else:
            self.ipLifo2.push(self.Flux)
            print("Flux and Power at iteration > 1")
        self.Power = self.ipLifo2.pushEmpty("Power", "LCM")
        self.ipLifo2.push(self.Cpo) # Push COMPO object
        self.ipLifo2.push(self.Track) # Push Tracking object for FEM solution : obtained from IniDONJON
        
        if iter == 1:
            self.ipLifo2.push(self.InitTHData)
            print("Pushing THData object at iter = 1")
        else:
            self.ipLifo2.push(self.UpdatedTHData)
            print("Pushing THData object")

        #ipLifo2.push(UpdatedTHData) # Push the TH data obtained from the TH solution
        self.ipLifo2.push(iter)
        self.ipLifo2.push(self.powi) 

        # 4.2) Call Neutronics component :
        print("call Neutronics procedure at iter=", iter)
        self.Neutronics.exec()
        print("Neutronics.c2m execution completed")

    def UpdateFluxAndPowerProfile(self):
        """
        Update the power profile with the power obtained from the neutron transport solver
        """
        Flux = self.ipLifo2.node("Flux") # Recover the Flux field
        RecoveredPower = self.ipLifo2.node("Power") # Recover the Power field (LCM object)
        # 4.2.1) Recover the Keff value
        Keff = Flux["K-EFFECTIVE"][0] 
        self.Keffs.append(Keff)
        print(f"At iter {iter} : Keff = {Keff}")
        # 4.2.2) Recover the Power axial distribution obtained from neutronics calculation --> used to update the axial power shape
        PowerDistribution = RecoveredPower["POWER-DISTR"]
        print(f"Power distribution : {PowerDistribution} kW")
        if iter > 1:
            if self.relaxPOW:
                # Under-relaxation of the power shape
                PowerDistribution = self.underRelaxation(PowerDistribution, self.Power_Distrib[-1], self.POWunderRelaxationFactor)
                self.Power_Distrib.append(PowerDistribution) # Store the axial power shape at each iteration
                print(f"Under-relaxed power distribution : {PowerDistribution} kW")
            else:
                self.Power_Distrib.append(PowerDistribution)
            self.Residuals_Power_Distrib.append(self.compute_residuals(self.Power_Distrib))
        # 4.2.3) Update the power profile
        self.qFiss = PowerDistribution*1000 # Updating the axial power shape, converting to W from kW, and dividing by the bundle volume to get W/m3
        print(f"Updated axial power/vol shape : {self.qFiss}")
        

        return
    
    def ComputeTHFields(self):
        """
        Compute the thermal-hydraulics fields with the updated axial power shape
        """
        # 5.) TH procedure for updated power shape
        ############# Thermalhydraulic part ##############
        # 5.1) TH resolution with updated power shape :
        self.THMComponent = THM_prototype(f"{self.case_name} equivalent canal", self.canalType, self.waterRadius, self.fuelRadius, self.gapRadius, self.cladRadius, 
                                self.height, self.tInlet, self.pOutlet, self.massFlowRate, self.qFiss, self.kFuel, self.Hgap, self.kClad, self.Iz1, self.If, self.I1, self.zPlotting, 
                                self.solveConduction, dt = 0, t_tot = 0, frfaccorel = self.frfaccorel, P2Pcorel = self.P2Pcorel, voidFractionCorrel = self.voidFractionCorrel, 
                                numericalMethod = self.numericalMethod)    ##### qFiss updated!
        
        self.TFuel, self.TCool, self.rhoCool = self.THMComponent.get_TH_parameters()
        if self.relaxTH:
            # Under-relaxation of the TH fields
            self.TFuel = self.underRelaxation(self.TFuel, self.Teff_list[-1], self.THunderRelaxationFactor)
            self.TCool = self.underRelaxation(self.TCool, self.TCool_list[-1], self.THunderRelaxationFactor)
            self.rhoCool = self.underRelaxation(self.rhoCool, self.rhoCool_list[-1], self.THunderRelaxationFactor)
        self.Teff_list.append(self.TFuel)
        self.TCool_list.append(self.TCool)
        self.rhoCool_list.append(self.rhoCool)

        self.Residual_TFuel_list.append(self.compute_residuals(self.Teff_list))
        self.Residual_TCool_list.append(self.compute_residuals(self.TCool_list))
        self.Residual_rhoCool_list.append(self.compute_residuals(self.rhoCool_list))

        return

    def UpdateTHFields(self):
        """
        Update the thermal-hydraulics fields with the updated power shape
        """
        # 5.2) Update the TH fields with the updated power shape
        self.UpdatedTHData = lcm.new('LCM','THData')
        self.UpdatedTHData['TFuelList']    = np.array(self.TFuel, dtype='f')
        self.UpdatedTHData['TCoolList'] = np.array(self.TCool, dtype='f')
        self.UpdatedTHData['DCoolList'] = np.array(self.rhoCool/1000, dtype='f')

        return

    def CheckConvergence(self):
        """
        Check the convergence of the coupled scheme
        """
        # 6.) Check convergence
        if self.convergenceCriterion == "POWField":
            conv = self.convergence(self.qFiss, self.qFiss_init, self.tolPOW)
        elif self.convergenceCriterion == "THField":
            conv = self.convergence(self.TFuel, self.TFuel_init, self.tolTH)
        return conv
    ########## Begin helper functions declaration ##########
    ## Fonction used for the convergence algorithm
    def underRelaxation(self, Field, OldField, underRelaxationFactor):
        return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

    def convergence(self, Field, OldField, tol, norm):
        conv=False
        for i in range(len(Field)):
            if np.abs(Field[i] - OldField[i]) > tol:
                conv = False
                break
            else:
                conv = True
        return conv

    def guessAxialPowerShape(self, Ptot, h, r, num_control_volumes):
        """
        Compute the power density in a cylinder with a sine-shaped profile along the axial direction.

        Parameters:
        Ptot (float): Total power released in the cylinder [W].
        h (float): Total height of the cylinder [m].
        r (float): Radius of the cylinder [m].
        num_control_volumes (int): Number of control volumes along the z-axis.

        Returns:
        z_values (numpy array): Array of z positions (center of each control volume) along the height of the cylinder.
        power_density (numpy array): Array of power densities corresponding to each control volume.
        q0 (float): The peak power density constant.
        """
        
        # Cross-sectional area of the cylinder
        A_cross_section = np.pi * r**2

        # Function for the power density profile along the axial direction
        def q(z, q0):
            return q0 * np.sin(np.pi * z / h)

        # Calculate q0 to ensure total power matches Ptot
        q0 = Ptot / (A_cross_section * (2 * h / np.pi))  # derived from the integral

        # Compute the boundaries and midpoints of each control volume
        z_boundaries = np.linspace(0, h, num_control_volumes + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes

        # Compute power density at each midpoint
        power_density = q(z_values, q0)

        return z_values, power_density, q0


    def compute_difference_fields(self,field):
        """
        Compute the difference between the last two fields of the list field
        field : list : list of the fields to compare
        """
        print(field)
        diff = np.abs(field[-1] - field[-2])
        print(f"Difference between the last two fields = {diff}")
        return diff

    def compute_residuals(self,field):
        """
        Compute the residuals of the field
        field : list : list of the fields to compute the residuals
        """
        residuals = (field[-1] - field[-2])*100/field[-2]
        print(f"Residuals of the field = {residuals} %")
        return residuals

    def quickPlot(self, x, y, title, xlabel, ylabel, saveName, path, SAVE_DIR):
        fig,ax = plt.subplots()
        if len(y) == len(x) and len(y[0]) == 1:
            print("Plotting a scatter plot")
            print(f"y = {y}")
            print(f"x = {x}")
            print(f"len(y) = {len(y)}")
            print(f"len(x) = {len(x)}")
            ax.scatter(x, y)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            os.chdir(SAVE_DIR)
            fig.savefig(saveName)
            os.chdir(path)
        else:
            for i in range(len(y)):
                if i%5==0:
                    data = y[i]
                    ax.plot(x, data, '2-',linewidth=1, label=f"iteration {i}")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.set_title(title)
            os.chdir(SAVE_DIR)
            fig.savefig(saveName)
            os.chdir(path)
        return
    ######## End helper functions declaration ##########
