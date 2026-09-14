# This file contains the implementation of the drift flux model for the THM prototype
# This class models the dynamic and steady-state behavior of a two-phase flow system in different geometries (square and cylindrical channels). 
# It discretizes the channel geometry and sets up the necessary fields for fluid flow, pressure, enthalpy, and void fraction in a thermal-hydraulic model. 
# It also includes methods for transient and steady-state simulations using various numerical techniques. The class supports setting up fission power, initializing 
# the flow fields, creating systems of equations for both steady and transient analysis, and solving them using finite volume method (FVM). Visualization tools 
# are provided to track residuals and convergence during iterative solving.

# Authors : Clement Huet
# Python3 class part of THM_prototype
# uses : - Drift flux model for two-phase flow
#        - Finite volume method for discretization of the conservation equations
#        - IAPWS97 for water properties
#        - THM_linalg for numerical resolution, it include a newton simple iteration method, a Gauss Siedel method, a BiCGStab method and a BiCG method
#        - THM_waterProp for water properties, the calculation of the void fraction, the calculation of the friction factor and the two-phase mltp depend on correlations
#        - THM_plotting for plotting

import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt
from THM_linalg import FVM
from THM_linalg import numericalResolution
from THM_waterProp import statesVariables
import cProfile

class DFMclass():
    def __init__(self, canal_type, nCells, tInlet, qFlow, pOutlet, height, fuelRadius, cladRadius, pitch,  numericalMethod, frfaccorel, P2P2corel, voidFractionCorrel, dt = 0, t_tot = 0, D_h = 0, volumetricArea = 0):
        
        """
        Attributes:
        - nFaces: Number of discretized cells.
        - pOutlet (Pa), tInlet (K): Inlet velocity, outle t pressure, and inlet enthalpy.
        - qFlow: Mass flow rate of the fluid (kg/s).
        - height (m), fuelRadius (m), cladRadius (m): Geometry of the channel (length, fuel, and clad radii).
        - pitch: Channel width or distance, depending on the geometry.
        - canalType: Geometry type of the channel, either 'square' or 'cylindrical'.
        - numericalMethod: Chosen method for numerical resolution (e.g., Gauss-Seidel, FVM, BiGStab).
        - voidFractionCorrel, frfaccorel, P2Pcorel: Correlations used for void fraction and other flow properties.
        - dt, t_tot: Time-step and total simulation time for transient analysis.
        """

        """
        Methods:
        - __init__(...): Initializes the class with geometric, inlet, outlet, and other user-specified parameters. It also sets up physical constants and mesh properties.
        - set_Fission_Power(Ptot, axial_p_form, Fpow): Sets the fission power source distribution term for the system.
        - get_Fission_Power(): Returns the source term (fission power) distribution along the channel.
        - setInitialFields(): Initializes the primary variables (velocity, pressure, enthalpy, void fraction) for steady-state simulation and updates flow properties.
        - createSystem(): Constructs the system of equations for solving the steady-state two-phase flow problem using the finite volume method.
        - createSystemTransient(): Sets up the system of equations for transient simulation.
        - calculateResiduals(): Calculates the residuals for velocity, pressure, and void fraction, and monitors the convergence.
        - testConvergence(k): Checks if the solution has converged based on residuals after iteration k.
        - residualsVisu(): Updates and visualizes residuals during the iterative solving process.
        - resolveDFM(): Main function that orchestrates the simulation by calling initializations, solving the system, and managing iterations and convergence criteria.
        - plotResults(): Plots the results of the simulation, including velocity, pressure, enthalpy, and void fraction profiles.
        - setInitialFieldsTransient(): Initializes the fields for transient simulation.
        - compute_T_surf(): Computes the surface temperature on the clad on the enthalpy profile.
        - sousRelaxation(): Implements under-relaxation for the solution update.
        - mergeVAR(): Merges the variables for the system of equations.
        - splitVAR(): Splits the variables for the system of equations.
        - createBoundaryEnthalpy(): Sets the boundary saturation lines for enthalpy.
        """ 

        if __name__ == "__main__":
            cProfile.run('main()', filename='profiling_result.prof')


        #user choice
        self.frfaccorel = frfaccorel
        self.P2Pcorel = P2P2corel
        self.numericalMethod = numericalMethod
        self.voidFractionCorrel = voidFractionCorrel
        self.voidFractionEquation = 'base'

        self.nCells = nCells
        self.nFaces = self.nCells + 1 #Number of faces
        self.pOutlet = pOutlet
        self.tInlet = tInlet

        #calculate temporary hInlet
        pressureDrop = 186737 #Pa/m
        falsePInlet = pOutlet + height * pressureDrop
        self.hInlet = IAPWS97(T = self.tInlet, P = falsePInlet * 10**(-6)).h*1000 #J/kg
        #print(f'hInlet: {self.hInlet}')

        #Geometry parameters
        self.height = height #m
        self.fuelRadius = fuelRadius #External radius of the fuel m
        self.cladRadius = cladRadius #External radius of the clad m
        self.pitch = pitch
        self.wall_dist = pitch
        self.canalType = canal_type

        if self.canalType == 'square':
            self.flowArea = self.pitch ** 2 - np.pi * self.cladRadius ** 2
        elif self.canalType == 'cylindrical':
            self.waterGap = self.pitch -  self.cladRadius #Gap between the clad and the water m
            self.waterRadius =  self.pitch #External radius of the water m
            self.flowArea = np.pi * self.waterRadius ** 2 - np.pi * self.cladRadius ** 2


        #calculate temporary uInlet
        self.qFlow = qFlow #kg/s
        self.rhoInlet = IAPWS97(T = self.tInlet, P = falsePInlet*10**(-6)).rho #kg/m3
        self.uInlet = self.qFlow / (self.flowArea * self.rhoInlet) #m/s
        #print('Velocity at the inlet: ', self.uInlet)
        #print(f'uInlet: {self.uInlet}')

        self.DV = (self.height/self.nCells) * self.flowArea #Volume of the control volume m3
        
        if self.canalType == 'square':
            
            self.Dh =  4 * self.flowArea / ( 2*np.pi * self.cladRadius) #2*self.pitch +
            #print(f'Dh: {self.Dh}')
        elif self.canalType == 'cylindrical':
            self.Dh = 4 * self.flowArea / (np.pi * self.waterRadius*2 + np.pi * self.cladRadius*2)


        self.Dz = self.height/self.nCells #Height of the control volume m
        self.z_mesh = np.linspace(0, self.height, self.nFaces)
        self.epsilonTarget = 0.1
        self.K_loss = 0.0

        #Porous media parameters
        self.porosity = 1
        self.poro = []
        for i in range(self.nFaces):
            if i > self.nFaces/4:
                self.poro.append(self.porosity)
            else:
                self.poro.append(1)

        self.D_h = []
        self.areaMatrix = []
        for i in range(self.nFaces):
            self.areaMatrix.append(self.poro[i] * self.flowArea)
            self.D_h.append(self.Dh * self.poro[i]**2)


        self.epsInnerIteration = 1e-4
        self.maxInnerIteration = 1000
        if self.numericalMethod == 'BiCGStab' or self.numericalMethod == 'BiCG':
            self.sousRelaxFactor = 0.8
        else:
            self.sousRelaxFactor = 1
        self.epsOuterIteration = 1e-4
        self.maxOuterIteration = 1000

        #Universal constant
        self.g = 9.81 #m/s^2
        self.R = 8.314 #J/(mol*K)

        #residuals
        self.EPSresiduals = []
        self.rhoResiduals = []
        self.rhoGResiduals = []
        self.rhoLResiduals = []
        self.xThResiduals = []
        self.UResiduals = []
        self.Iteration = []
        self.I = []

        self.hlSat = []
        self.hgSat = []

        #Transient parameters
        self.dt = dt
        if dt != 0:
            self.t_tot = t_tot
            self.timeList = np.arange(0, self.t_tot, self.dt)
            self.timeCount = 0
    
    #Fission power
    def set_Fission_Power(self, Ptot, axial_p_forms, Fpow):
        print("Setting fission power distribution...")
        print(f'Ptot: {Ptot} W, axial_p_forms: {axial_p_forms}, Fpow: {Fpow}')
        dv = self.pitch**2*(self.height / self.nCells) #m
        print(f"dv: {dv} m3, pitch: {self.pitch} m, height: {self.height} m, nCells: {self.nCells}, flowArea: {self.flowArea} m2")
        Power_dist = Ptot * axial_p_forms / self.nCells #W Axial power distribution
        linear_powers = Power_dist
        assembly_section = self.pitch**2
        print(f'assembly_section: {assembly_section} m2, fuelRadius: {self.fuelRadius} m, cladRadius: {self.cladRadius} m')
        proportion_fuel = assembly_section / (np.pi * self.fuelRadius**2) #m2
        fraction_in_fuel = Fpow*proportion_fuel #Fraction of the total power released in fuel
        fraction_in_coolant = (1.0-Fpow)*proportion_fuel # Fraction of the total power released in coolant
        if self.dt == 0:
            self.QFUEL = Power_dist*fraction_in_fuel / dv #W/m3
            phi2 = 0.5*self.QFUEL*self.fuelRadius**2 / self.cladRadius
            self.q__ = 2 * np.pi * self.cladRadius / (self.flowArea) * phi2 # W/m3 in the coolant
            print(f'q__: {self.q__}')
        if self.dt != 0: # This option is not suppoted yet, dummy initialization for now
            t_final_q = 0
            self.q__ = np.zeros((len(self.timeList), len(axial_p_forms)))
            for t, time in enumerate(self.timeList):
                if time < t_final_q:
                    self.q__[t] = (self.timeList[t]/t_final_q)*Power_dist*fraction_in_fuel / dv #W/m3
                else:
                    self.q__[t] = Power_dist*fraction_in_fuel / dv #W/m3
            
            print(f'q__: {self.q__}')
            plt.plot(self.timeList, self.q__[:, -1])
            plt.show()

    #Recovers the fission power distribution
    def get_QFUEL(self):
        """
        function to retrieve the fission power density distribution in the fuel rod
        """
        return self.QFUEL
    
    def get_Fission_Power(self):
        """
        function to retrieve a given source term from the axial profile used to model fission power distribution in the fuel rod
        """
        return self.q__
        
    #Define the initial fields for the first iteration
    def setInitialFields(self):
        
        #Steady state
        if self.dt == 0:
        
            self.U = [np.ones(self.nFaces)*self.uInlet] #
            self.P = [np.ones(self.nFaces)*self.pOutlet] #
            self.H = [np.ones(self.nFaces)*self.hInlet] #
            self.voidFraction = [np.array([i*self.epsilonTarget/self.nFaces for i in range(self.nFaces)])]

            updateVariables = statesVariables(self.U[-1], self.P[-1], self.H[-1], self.voidFraction[-1], self.D_h, self.areaMatrix, self.DV, self.voidFractionCorrel, self.frfaccorel, self.P2Pcorel, self.Dz, self.q__, self.qFlow, self.fuelRadius,  self.pitch/2)
            updateVariables.createFields()
                
            self.xTh = [np.ones(self.nFaces)]
            self.rhoL= [updateVariables.rholTEMP]
            self.rhoG = [updateVariables.rhogTEMP]
            self.rho = [updateVariables.rhoTEMP]
            self.Dhfg = [updateVariables.DhfgTEMP]
            self.f = [updateVariables.fTEMP]
            self.areaMatrix_1 = [updateVariables.areaMatrix_1TEMP]
            self.areaMatrix_2 = [updateVariables.areaMatrix_2TEMP]
            self.areaMatrix = updateVariables.areaMatrixTEMP
            self.Vgj = [updateVariables.VgjTEMP]
            self.C0 =[updateVariables.C0TEMP]
            self.VgjPrime = [updateVariables.VgjPrimeTEMP]

        #Transient
        else:
            if self.timeCount == 0:
                self.U = [self.velocityList[self.timeCount]]
                self.P = [self.pressureList[self.timeCount]]
                self.H = [self.enthalpyList[self.timeCount]]
                self.voidFraction = [self.voidFractionList[self.timeCount]]

                updateVariables = statesVariables(self.U[-1], self.P[-1], self.H[-1], self.voidFraction[-1], self.D_h, self.areaMatrix, self.DV, self.voidFractionCorrel, self.frfaccorel, self.P2Pcorel, self.Dz, self.q__, self.qFlow, self.fuelRadius, self.pitch/2)
                updateVariables.createFields()

                self.xTh = [np.ones(self.nFaces)]
                self.rhoL= [updateVariables.rholTEMP]
                self.rhoG = [updateVariables.rhogTEMP]
                self.rho = [updateVariables.rhoTEMP]
                self.Dhfg = [updateVariables.DhfgTEMP]
                self.f = [updateVariables.fTEMP]
                self.areaMatrix_1 = [updateVariables.areaMatrix_1TEMP]
                self.areaMatrix_2 = [updateVariables.areaMatrix_2TEMP]
                self.areaMatrix = updateVariables.areaMatrixTEMP
                self.Vgj = [updateVariables.VgjTEMP]
                self.C0 =[updateVariables.C0TEMP]
                self.VgjPrime = [updateVariables.VgjPrimeTEMP]

                self.xThList[self.timeCount] = self.xTh[-1]
                self.rhoList[self.timeCount] = self.rho[-1]
                self.rhoGList[self.timeCount] = self.rhoG[-1]
                self.rhoLList[self.timeCount] = self.rhoL[-1]
                self.DhfgList[self.timeCount] = self.Dhfg[-1]
                self.fList[self.timeCount] = self.f[-1]
                self.areaMatrix_1List[self.timeCount] = self.areaMatrix_1[-1]
                self.areaMatrix_2List[self.timeCount] = self.areaMatrix_2[-1]
                self.areaMatrixList[self.timeCount] = self.areaMatrix
                self.VgjList[self.timeCount] = self.Vgj[-1]
                self.C0List[self.timeCount] = self.C0[-1]
                self.VgjPrimeList[self.timeCount] = self.VgjPrime[-1]
            
            if self.timeCount != 0:
                self.U = [self.velocityList[self.timeCount-1]]
                self.P = [self.pressureList[self.timeCount-1]]
                self.H = [self.enthalpyList[self.timeCount-1]]
                self.voidFraction = [self.voidFractionList[self.timeCount-1]]
                self.xTh = [self.xThList[self.timeCount-1]]
                self.rhoL = [self.rhoLList[self.timeCount-1]]
                self.rhoG = [self.rhoGList[self.timeCount-1]]
                self.rho = [self.rhoList[self.timeCount-1]]
                self.Dhfg = [self.DhfgList[self.timeCount-1]]
                self.f = [self.fList[self.timeCount-1]]
                self.areaMatrix_1 = [self.areaMatrix_1List[self.timeCount-1]]
                self.areaMatrix_2 = [self.areaMatrix_2List[self.timeCount-1]]
                self.areaMatrix = self.areaMatrixList[self.timeCount-1]
                self.Vgj = [self.VgjList[self.timeCount-1]]
                self.C0 = [self.C0List[self.timeCount-1]]
                self.VgjPrime = [self.VgjPrimeList[self.timeCount-1]]

    #Create the matrix for the velocity and pressure coupling resolution equation system
    def createSystemVelocityPressure(self):
        
        U_old = self.U[-1]
        P_old = self.P[-1]
        H_old = self.H[-1]
        epsilon_old = self.voidFraction[-1]
        rho_old = self.rho[-1]
        rho_g_old = self.rhoG[-1]
        rho_l_old = self.rhoL[-1]
        areaMatrix = self.areaMatrix
        areaMatrix_old_1 = self.areaMatrix_1[-1]
        areaMatrix_old_2 = self.areaMatrix_2[-1]
        Dhfg = self.Dhfg[-1]
        x_th_old = self.xTh[-1]
        f = self.f[-1]
        V_gj_old = self.Vgj[-1]
        Vgj_prime = self.VgjPrime[-1]
        C0 = self.C0[-1]

        VAR_old = self.mergeVar(U_old,P_old)
        rho_old = self.mergeVar(rho_old, rho_old)
        rho_g_old = self.mergeVar(rho_g_old, rho_g_old)
        rho_l_old = self.mergeVar(rho_l_old, rho_l_old)
        epsilon_old = self.mergeVar(epsilon_old, epsilon_old)
        areaMatrix = self.mergeVar(areaMatrix, areaMatrix)
        areaMatrix_old_1 = self.mergeVar(areaMatrix_old_1, areaMatrix_old_1)
        areaMatrix_old_2 = self.mergeVar(areaMatrix_old_2, areaMatrix_old_2)
        V_gj_old = self.mergeVar(V_gj_old, V_gj_old)
        Vgj_prime = self.mergeVar(Vgj_prime, Vgj_prime)
        Dhfg = self.mergeVar(Dhfg, Dhfg)
        C0 = self.mergeVar(C0, C0)
        x_th_old = self.mergeVar(x_th_old, x_th_old)
        
        VAR_VFM_Class = FVM(A00 = 1, A01 = 0, Am0 = 0, Am1 = 1, D0 = self.uInlet, Dm1 = self.pOutlet, N_vol = 2*self.nFaces, H = self.height)
        VAR_VFM_Class.boundaryFilling()

        for i in range(1, 2*self.nFaces-1):
            #Inside the velocity submatrix
            if i < self.nFaces-1:
                VAR_VFM_Class.set_ADi(i, ci = - rho_old[i-1]*areaMatrix[i-1],
                ai = rho_old[i]*areaMatrix[i],
                bi = 0,
                di =  0)
            elif i == self.nFaces-1:
                VAR_VFM_Class.set_ADi(i, 
                ci = - rho_old[i-1]*areaMatrix[i-1],
                ai = rho_old[i]*areaMatrix[i],
                bi = 0,
                di =  0)

            #Inside the pressure submatrix
            elif i == self.nFaces:
                DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 * areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 * areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
                VAR_VFM_Class.set_ADi(self.nFaces, 
                ci = 0,
                ai = - areaMatrix[i],
                bi = areaMatrix[i+1],
                di = - (((rho_old[i+1]+ rho_old[i])* self.g/2) * self.DV * ((self.poro[i%self.nFaces]+ self.poro[(i+1)%self.nFaces])/2) / 2) + DI)
            
                VAR_VFM_Class.fillingOutsideBoundary(i, i-self.nFaces,
                ai = - rho_old[i]*VAR_old[i-self.nFaces]*areaMatrix_old_2[i],
                bi = rho_old[i+1]*VAR_old[i-self.nFaces+1]*areaMatrix_old_1[i+1])

            elif i > self.nFaces and i < 2*self.nFaces-1:
                DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 * areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 * areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
                VAR_VFM_Class.set_ADi(i, ci = 0,
                ai = - areaMatrix[i],
                bi = areaMatrix[i+1],
                di = - (((rho_old[i+1]+ rho_old[i])* self.g/2) * self.DV * ((self.poro[i%self.nFaces]+ self.poro[(i+1)%self.nFaces])/2)/ 2) + DI)
            
                VAR_VFM_Class.fillingOutsideBoundary(i, i-self.nFaces,
                ai = - rho_old[i]*VAR_old[i-self.nFaces]*areaMatrix_old_2[i],
                bi = rho_old[i+1]*VAR_old[i+1-self.nFaces]*areaMatrix_old_1[i+1])

        
        self.FVM = VAR_VFM_Class

    #Create the enthalpy matrix resolution equation system
    def createSystemEnthalpy(self):

        U_old = self.U[-1]
        P_old = self.P[-1]
        H_old = self.H[-1]
        epsilon_old = self.voidFraction[-1]
        rho_old = self.rho[-1]
        rho_g_old = self.rhoG[-1]
        rho_l_old = self.rhoL[-1]
        areaMatrix = self.areaMatrix
        Dhfg = self.Dhfg[-1]
        x_th_old = self.xTh[-1]
        f = self.f[-1]
        V_gj_old = self.Vgj[-1]
        Vgj_prime = self.VgjPrime[-1]
        C0 = self.C0[-1]


        i = -1
        DI = (1/2) * (P_old[i-1]*areaMatrix[i-1] - P_old[i]*areaMatrix[i]) * ((U_old[i]+ ((epsilon_old[i] * (rho_l_old[i] - rho_g_old[i]) * V_gj_old[i])/ rho_old[i]))+ (U_old[i-1]+ ((epsilon_old[i-1] * (rho_l_old[i-1] - rho_g_old[i-1]) * V_gj_old[i-1])/ rho_old[i-1]) ) )
        DI2 = - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]*areaMatrix[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]*areaMatrix[i-1]/rho_old[i-1])
        DM1 = self.q__[i-1] * self.DV * (self.poro[i]) + DI + DI2
        VAR_VFM_Class = FVM(A00 = 1, A01 = 0, Am0 = - rho_old[-2] * U_old[-2] * areaMatrix[-2], Am1 = rho_old[-1] * U_old[-1] * areaMatrix[-1], D0 = self.hInlet, Dm1 = DM1, N_vol = self.nFaces, H = self.height)
        VAR_VFM_Class.boundaryFilling()
        for i in range(1,self.nFaces -1):
            #Inside the enthalpy submatrix
            DI = (1/2) * (P_old[i-1]*areaMatrix[i-1] - P_old[i]*areaMatrix[i]) * ((U_old[i]+ ((epsilon_old[i] * (rho_l_old[i] - rho_g_old[i]) * V_gj_old[i])/ rho_old[i]))+ (U_old[i-1]+ ((epsilon_old[i-1] * (rho_l_old[i-1] - rho_g_old[i-1]) * V_gj_old[i-1])/ rho_old[i-1]) ) )
            DI2 = - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]*areaMatrix[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]*areaMatrix[i-1]/rho_old[i-1])
            VAR_VFM_Class.set_ADi(i, ci =  - rho_old[i-1] * U_old[i-1] * areaMatrix[i-1],
                ai = rho_old[i] * U_old[i] * areaMatrix[i],
                bi = 0,
                di =  self.q__[i-1] * self.DV * (self.poro[i]) + DI2 + DI)
        
        self.FVM = VAR_VFM_Class

    #Create the enthalpy matrix for the transient resolution equation system
    def createSystemEnthalpyTransient(self):

        U_old = self.U[-1]
        P_old = self.P[-1]
        H_old = self.H[-1]
        epsilon_old = self.voidFraction[-1]
        rho_old = self.rho[-1]
        rho_g_old = self.rhoG[-1]
        rho_l_old = self.rhoL[-1]
        areaMatrix = self.areaMatrix
        Dhfg = self.Dhfg[-1]
        x_th_old = self.xTh[-1]
        f = self.f[-1]
        V_gj_old = self.Vgj[-1]
        Vgj_prime = self.VgjPrime[-1]
        C0 = self.C0[-1]


        i = -1
        DI = (1/2) * (P_old[i]*areaMatrix[i] - P_old[i-1]*areaMatrix[i-1]) * ((U_old[i]+ ((epsilon_old[i] * (rho_l_old[i] - rho_g_old[i]) * V_gj_old[i])/ rho_old[i]))+ (U_old[i-1]+ ((epsilon_old[i-1] * (rho_l_old[i-1] - rho_g_old[i-1]) * V_gj_old[i-1])/ rho_old[i-1]) ) )
        DI2 = - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]*areaMatrix[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]*areaMatrix[i-1]/rho_old[i-1])
        DT1 = - (self.pressureList[self.timeCount][i%self.nFaces] * self.areaMatrix[i] - P_old[i] * areaMatrix[i])*(self.Dz/self.dt) + (self.rhoList[self.timeCount][i%self.nFaces] * self.enthalpyList[self.timeCount][i%self.nFaces] * areaMatrix[i] * (self.Dz / self.dt))
        DM1 = self.q__[self.timeCount][i-1] * self.DV * (self.poro[i]) + DI + DI2 + DT1
        VAR_VFM_Class = FVM(A00 = 1, A01 = 0, Am0 = - rho_old[-2] * U_old[-2] * areaMatrix[-2] + rho_old[-2] * areaMatrix[-2] * (self.Dz / self.dt), Am1 = rho_old[-1] * U_old[-1] * areaMatrix[-1], D0 = self.hInlet, Dm1 = DM1, N_vol = self.nFaces, H = self.height)
        VAR_VFM_Class.boundaryFilling()
        for i in range(1,self.nFaces -1):
            #Inside the enthalpy submatrix
            DI = (1/2) * (P_old[i]*areaMatrix[i] - P_old[i-1]*areaMatrix[i-1]) * ((U_old[i]+ ((epsilon_old[i] * (rho_l_old[i] - rho_g_old[i]) * V_gj_old[i])/ rho_old[i]))+ (U_old[i-1]+ ((epsilon_old[i-1] * (rho_l_old[i-1] - rho_g_old[i-1]) * V_gj_old[i-1])/ rho_old[i-1]) ) )
            DI2 = - (epsilon_old[i]*rho_l_old[i]*rho_g_old[i]*Dhfg[i]*V_gj_old[i]*areaMatrix[i]/rho_old[i]) + (epsilon_old[i-1]*rho_l_old[i-1]*rho_g_old[i-1]*Dhfg[i-1]*V_gj_old[i-1]*areaMatrix[i-1]/rho_old[i-1])
            DT1 = -(self.pressureList[self.timeCount][i%self.nFaces] * self.areaMatrixList[self.timeCount][i%self.nFaces] - P_old[i] * areaMatrix[i])*(self.Dz/self.dt) + (self.rhoList[self.timeCount][i%self.nFaces] * self.enthalpyList[self.timeCount][i%self.nFaces] * areaMatrix[i] * (self.Dz / self.dt))
            VAR_VFM_Class.set_ADi(i, ci =  - rho_old[i-1] * U_old[i-1] * areaMatrix[i-1] + rho_old[i] * areaMatrix[i] * (self.Dz / self.dt),
                ai = rho_old[i] * U_old[i] * areaMatrix[i],
                bi = 0,
                di =  self.q__[self.timeCount][i-1] * self.DV * (self.poro[i]) + DI + DI2 + DT1)
        
        self.FVM = VAR_VFM_Class

    #Create the transient matrix for the velocity and pressure coupling resolution equation system
    def createSystemVelocityPressureTransient(self):

        U_old = self.U[-1]
        P_old = self.P[-1]
        H_old = self.H[-1]
        epsilon_old = self.voidFraction[-1]
        rho_old = self.rho[-1]
        rho_g_old = self.rhoG[-1]
        rho_l_old = self.rhoL[-1]
        areaMatrix = self.areaMatrix
        areaMatrix_old_1 = self.areaMatrix_1[-1]
        areaMatrix_old_2 = self.areaMatrix_2[-1]
        Dhfg = self.Dhfg[-1]
        x_th_old = self.xTh[-1]
        f = self.f[-1]
        V_gj_old = self.Vgj[-1]
        Vgj_prime = self.VgjPrime[-1]
        C0 = self.C0[-1]

        VAR_old = self.mergeVar(U_old,P_old)
        rho_old = self.mergeVar(rho_old, rho_old)
        rho_g_old = self.mergeVar(rho_g_old, rho_g_old)
        rho_l_old = self.mergeVar(rho_l_old, rho_l_old)
        epsilon_old = self.mergeVar(epsilon_old, epsilon_old)
        areaMatrix = self.mergeVar(areaMatrix, areaMatrix)
        areaMatrix_old_1 = self.mergeVar(areaMatrix_old_1, areaMatrix_old_1)
        areaMatrix_old_2 = self.mergeVar(areaMatrix_old_2, areaMatrix_old_2)
        V_gj_old = self.mergeVar(V_gj_old, V_gj_old)
        Vgj_prime = self.mergeVar(Vgj_prime, Vgj_prime)
        Dhfg = self.mergeVar(Dhfg, Dhfg)
        C0 = self.mergeVar(C0, C0)
        x_th_old = self.mergeVar(x_th_old, x_th_old)

        VAR_VFM_Class = FVM(A00 = 1, A01 = 0, Am0 = 0, Am1 = 1, D0 = self.uInlet, Dm1 = self.pOutlet, N_vol = 2*self.nFaces, H = self.height)
        VAR_VFM_Class.boundaryFilling()
        for i in range(1, 2*self.nFaces-1):
            #Inside the velocity submatrix
            if i < self.nFaces-1:
                VAR_VFM_Class.set_ADi(i, ci = - rho_old[i-1]*areaMatrix[i-1],
                ai = rho_old[i]*areaMatrix[i],
                bi = 0,
                di =  ( self.rhoList[self.timeCount][i%self.nFaces] *areaMatrix[i] - rho_old[i] *areaMatrix[i] ) * (self.Dz / self.dt))
            elif i == self.nFaces-1:
                VAR_VFM_Class.set_ADi(i, 
                ci = - rho_old[i-1]*areaMatrix[i-1],
                ai = rho_old[i]*areaMatrix[i],
                bi = 0,
                di =  ( self.rhoList[self.timeCount][i%self.nFaces] * areaMatrix[i] - rho_old[i] *areaMatrix[i] ) * (self.Dz / self.dt))

            #Inside the pressure submatrix
            elif i == self.nFaces:
                DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 * areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 * areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
                VAR_VFM_Class.set_ADi(self.nFaces, 
                ci = 0,
                ai = - areaMatrix[i],
                bi = areaMatrix[i+1],
                di = - ((rho_old[i+1]+ rho_old[i])* self.g * self.DV * ((self.poro[i%self.nFaces]+ self.poro[(i+1)%self.nFaces])/2) / 2) + DI + (self.rhoList[self.timeCount][i%self.nFaces] * areaMatrix[i] * self.velocityList[self.timeCount][i%self.nFaces] * (self.Dz / self.dt)))
            
                VAR_VFM_Class.fillingOutsideBoundary(i, i-self.nFaces,
                ai = - rho_old[i]*VAR_old[i-self.nFaces]*areaMatrix_old_2[i] + rho_old[i]*areaMatrix[i]*(self.Dz/self.dt),
                bi = rho_old[i+1]*VAR_old[i-self.nFaces+1]*areaMatrix_old_1[i+1])

            elif i > self.nFaces and i < 2*self.nFaces-1:
                DI = -((epsilon_old[i+1] * rho_g_old[i+1] * rho_l_old[i+1] * V_gj_old[i+1]**2 * areaMatrix[i+1] )/ ((1 - epsilon_old[i+1])*rho_old[i+1]) )  + ((epsilon_old[i] * rho_g_old[i] * rho_l_old[i] * V_gj_old[i]**2 * areaMatrix[i] )/ ((1 - epsilon_old[i])*rho_old[i]) )     
                VAR_VFM_Class.set_ADi(i, ci = 0,
                ai = - areaMatrix[i],
                bi = areaMatrix[i+1],
                di = - ((rho_old[i+1]+ rho_old[i])* self.g * self.DV * ((self.poro[i%self.nFaces]+ self.poro[(i+1)%self.nFaces])/2)/ 2) + DI + (self.rhoList[self.timeCount][i%self.nFaces] * areaMatrix[i] * self.velocityList[self.timeCount][i%self.nFaces] * (self.Dz / self.dt)))
            
                VAR_VFM_Class.fillingOutsideBoundary(i, i-self.nFaces,
                ai = - rho_old[i]*VAR_old[i-self.nFaces]*areaMatrix_old_2[i] + rho_old[i]*areaMatrix[i]*(self.Dz/self.dt),
                bi = rho_old[i+1]*VAR_old[i+1-self.nFaces]*areaMatrix_old_1[i+1])

        self.FVM = VAR_VFM_Class

    #Calculate the residuals
    def calculateResiduals(self):#change les residus
        self.EPSresiduals.append(np.linalg.norm(self.voidFraction[-1] - self.voidFraction[-2]))
        self.rhoResiduals.append(np.linalg.norm((self.rho[-1] - self.rho[-2])/self.rho[-1]))
        self.xThResiduals.append(np.linalg.norm(self.xTh[-1] - self.xTh[-2]))

    #Checking for convergence
    def testConvergence(self, k):#change rien et return un boolean
        print(f'Convergence test number {k}, RES: errEPS: {self.EPSresiduals[-1]}, errRHO: {self.rhoResiduals[-1]}, errQua: {self.xThResiduals[-1]}')
        if self.EPSresiduals[-1] < self.epsOuterIteration and self.rhoResiduals[-1] < self.epsOuterIteration: #and self.xThResiduals[-1] < 1e-3 :
            return True
        else:
            return False

    #Creting the initial fields for the transient resolution
    def setInitialFieldsTransient(self):
        
        self.velocityList = np.zeros((len(self.timeList), self.nFaces))
        self.pressureList = np.zeros((len(self.timeList), self.nFaces))
        self.enthalpyList = np.zeros((len(self.timeList), self.nFaces))
        self.voidFractionList = np.zeros((len(self.timeList), self.nFaces))
        self.rhoList = np.zeros((len(self.timeList), self.nFaces))
        self.rhoGList = np.zeros((len(self.timeList), self.nFaces))
        self.rhoLList = np.zeros((len(self.timeList), self.nFaces))
        self.xThList = np.zeros((len(self.timeList), self.nFaces))
        self.DhfgList = np.zeros((len(self.timeList), self.nFaces))
        self.fList = np.zeros((len(self.timeList), self.nFaces))
        self.areaMatrix_1List = np.zeros((len(self.timeList), self.nFaces))
        self.areaMatrix_2List = np.zeros((len(self.timeList), self.nFaces))
        self.areaMatrixList = np.zeros((len(self.timeList), self.nFaces))
        self.VgjList = np.zeros((len(self.timeList), self.nFaces))
        self.C0List = np.zeros((len(self.timeList), self.nFaces))
        self.VgjPrimeList = np.zeros((len(self.timeList), self.nFaces))

        self.velocityList[:,0] = self.uInlet
        self.pressureList[:,-1] = self.pOutlet
        self.enthalpyList[:,0] = self.hInlet
        self.velocityList[0,:] = self.uInlet
        self.pressureList[0,:] = self.pOutlet
        self.enthalpyList[0,:] = self.hInlet

    #Function to visualise the live evolution of the resuaduels
    def residualsVisu(self):
        # Mise à jour des données de la ligne
        self.line.set_xdata(self.I)
        self.line.set_ydata(self.rhoResiduals)

        # Ajuste les limites des axes si nécessaire
        self.ax.relim()         # Recalcule les limites des données
        self.ax.autoscale_view()  # Réajuste la vue automatiquement

        # Dessine les modifications
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
    
    #The pressure inlet change so the inlet velocity and enthalpy need to be updated between the pressure-velocity and the enthalpy resolution
    def updateInlet(self):
        #Update uInlet
        self.rhoInlet = IAPWS97(T = self.tInlet, P = self.P[-1][0]*10**(-6)).rho #kg/m3
        self.uInlet = self.qFlow / (self.flowArea * self.rhoInlet) #m/s
        #Update hInlet
        #print(f"uInlet:{self.uInlet}")
        self.hInlet = IAPWS97(T = self.tInlet, P = self.P[-1][0]*10**(-6)).h*1000 #J/kg

    #Main function to solve the drift flux model
    def resolveDFM(self):

        #Steady state
        if self.dt == 0:

            self.setInitialFields()
            # Active le mode interactif
            #plt.ion()
            # Crée la figure et l'axe
            #self.fig, self.ax = plt.subplots()
            # Initialisation de la ligne qui sera mise à jour
            #self.line, = self.ax.plot(self.I, self.rhoResiduals, 'r-', marker='o')  # 'r-' pour une ligne rouge avec des marqueurs

            #Loop for the outer iterations (velocity-pressure and enthalpy)
            for k in range(self.maxOuterIteration):
                
                self.createSystemVelocityPressure()
                resolveSystem = numericalResolution(self.FVM,self.mergeVar(self.U[-1], self.P[-1]), self.epsInnerIteration, self.maxInnerIteration, self.numericalMethod)
                Utemp, Ptemp = self.splitVar(resolveSystem.x)
                
                self.U.append(Utemp)
                self.P.append(Ptemp)
                
                self.updateInlet()
                
                self.createSystemEnthalpy()
                resolveSystem = numericalResolution(self.FVM, self.H[-1], self.epsInnerIteration, self.maxInnerIteration, self.numericalMethod)
                
                Htemp = resolveSystem.x

                self.H.append(Htemp)
                updateVariables = statesVariables(self.U[-1], self.P[-1], self.H[-1], self.voidFraction[-1], self.D_h, self.areaMatrix, self.DV, self.voidFractionCorrel, self.frfaccorel, self.P2Pcorel, self.Dz, self.q__, self.qFlow, self.fuelRadius, self.pitch/2)
                updateVariables.updateFields()

                self.xTh.append(updateVariables.xThTEMP)
                self.rhoL.append(updateVariables.rholTEMP)
                self.rhoG.append(updateVariables.rhogTEMP)
                self.rho.append(updateVariables.rhoTEMP)
                self.voidFraction.append(updateVariables.voidFractionTEMP)
                self.Dhfg.append(updateVariables.DhfgTEMP)
                self.f.append(updateVariables.fTEMP)
                self.areaMatrix_1.append(updateVariables.areaMatrix_1TEMP)
                self.areaMatrix_2.append(updateVariables.areaMatrix_2TEMP)
                self.Vgj.append(updateVariables.VgjTEMP)
                self.C0.append(updateVariables.C0TEMP)
                self.VgjPrime.append(updateVariables.VgjPrimeTEMP)

                self.sousRelaxation()
                self.calculateResiduals()
                self.I.append(k)
                #self.residualsVisu()

                convergence = self.testConvergence(k)

                self.updateInlet()

                if convergence == True:
                    print(f'Convergence reached at iteration number: {k}')
                    break

                elif k == self.maxOuterIteration - 1:
                    raise ValueError('Convergence not reached in the resolution of the drift flux model, not enough iterations. k = ', k)


            self.Ul = updateVariables.Ul
            self.Ug = updateVariables.Ug
            self.Rel = updateVariables.Rel
            #plt.ioff()
            #plt.show()
        
        #Transient
        elif self.dt != 0:

            self.setInitialFieldsTransient()
            # Active le mode interactif
            plt.ion()
            # Crée la figure et l'axe
            self.fig, self.ax = plt.subplots()
            self.figh, self.axh = plt.subplots()   
            # Initialisation de la ligne qui sera mise à jour
            self.line, = self.ax.plot(self.I, self.rhoResiduals, 'r-', marker='o')  # 'r-' pour une ligne rouge avec des marqueurs
            
            for t in range(0, len(self.timeList)-1):
                self.timeCount = t
                self.setInitialFields()

                for k in range(self.maxOuterIteration):
                    self.createSystemVelocityPressureTransient()
                    resolveSystem = numericalResolution(self.FVM,self.mergeVar(self.U[-1], self.P[-1]), self.epsInnerIteration, self.maxInnerIteration, self.numericalMethod)
                    Utemp, Ptemp = self.splitVar(resolveSystem.x)
                    
                    self.U.append(Utemp)
                    self.P.append(Ptemp)
                    self.updateInlet()
                    
                    self.createSystemEnthalpyTransient()
                    resolveSystem = numericalResolution(self.FVM, self.H[-1], self.epsInnerIteration, self.maxInnerIteration, self.numericalMethod)
                    
                    Htemp = resolveSystem.x

                    self.H.append(Htemp)
                    updateVariables = statesVariables(self.U[-1], self.P[-1], self.H[-1], self.voidFraction[-1], self.D_h, self.areaMatrix, self.DV, self.voidFractionCorrel, self.frfaccorel, self.P2Pcorel, self.Dz, self.q__, self.fuelRadius, self.pitch/2)
                    updateVariables.updateFields()

                    self.xTh.append(updateVariables.xThTEMP)
                    self.rhoL.append(updateVariables.rholTEMP)
                    self.rhoG.append(updateVariables.rhogTEMP)
                    self.rho.append(updateVariables.rhoTEMP)
                    self.voidFraction.append(updateVariables.voidFractionTEMP)
                    self.Dhfg.append(updateVariables.DhfgTEMP)
                    self.f.append(updateVariables.fTEMP)
                    self.areaMatrix_1.append(updateVariables.areaMatrix_1TEMP)
                    self.areaMatrix_2.append(updateVariables.areaMatrix_2TEMP)
                    self.Vgj.append(updateVariables.VgjTEMP)
                    self.C0.append(updateVariables.C0TEMP)
                    self.VgjPrime.append(updateVariables.VgjPrimeTEMP)
                    #(f'rho: {self.rho[-1]}, U: {self.U[-1]}, P: {self.P[-1]}, H: {self.H[-1]}')
                    self.sousRelaxation()
                    self.calculateResiduals()
                    if self.I == []:
                        self.I.append(0)
                    else:
                        self.I.append(1+self.I[-1])
                    self.residualsVisu()

                    convergence = self.testConvergence(k)

                    self.updateInlet()

                    if convergence == True:
                        self.velocityList[self.timeCount+1] = self.U[-1]
                        self.pressureList[self.timeCount+1] = self.P[-1]
                        self.enthalpyList[self.timeCount+1] = self.H[-1]
                        self.voidFractionList[self.timeCount+1] = self.voidFraction[-1]
                        self.xThList[self.timeCount+1] = self.xTh[-1]
                        self.rhoList[self.timeCount+1] = self.rho[-1]
                        self.rhoGList[self.timeCount+1] = self.rhoG[-1]
                        self.rhoLList[self.timeCount+1] = self.rhoL[-1]
                        self.DhfgList[self.timeCount+1] = self.Dhfg[-1]
                        self.fList[self.timeCount+1] = self.f[-1]
                        self.areaMatrix_1List[self.timeCount+1] = self.areaMatrix_1[-1]
                        self.areaMatrix_2List[self.timeCount+1] = self.areaMatrix_2[-1]
                        self.areaMatrixList[self.timeCount+1] = self.areaMatrix
                        self.VgjList[self.timeCount+1] = self.Vgj[-1]
                        self.C0List[self.timeCount+1] = self.C0[-1]
                        self.VgjPrimeList[self.timeCount+1] = self.VgjPrime[-1]

                        print(f'Convergence reached at iteration number: {k} for time step: {t}')
                        break

                    elif k == self.maxOuterIteration - 1:
                        raise ValueError('Convergence not reached')
                    
                """ champs = [Utemp, Ptemp, Htemp, self.voidFraction[-1], self.rho[-1], self.rhoG[-1]]
                titres = ['U', 'P', 'H', 'epsilon', 'rho', 'rho_g']
                fig, axs = plt.subplots(3, 2, figsize=(12, 8))
                for i, ax in enumerate(axs.flat):
                    ax.plot(self.z_mesh, champs[i])
                    ax.set_title(titres[i])
                    ax.set_xlabel('x')
                plt.tight_layout()
                plt.title(f'Champs pour le pas de temps {t}')
                plt.legend()
                plt.show() """
                    
                #print(f'U: {self.velocityList}, P: {self.pressureList}, H: {self.enthalpyList}')

            plt.ioff()
            plt.show()

            #besoin d'interpoler pour avoir le bon nombre de valeur (passer de node centered a cell centered)
        
        self.T_water = np.zeros(self.nFaces)
        for i in range(self.nFaces):
            self.T_water[i] = IAPWS97(P=self.P[-1][i]*10**-6, h=self.H[-1][i]*10**-3).T

        self.interpolate()

    #Function to interpolate the values from the nodes to the cells
    def interpolate(self):
        Ptemp, Utemp, Htemp, voidFractionTemp, rhoTemp, rhoGTemp, rhoLTemp, xThTemp, DhfgTemp, fTemp, areaMatrix_1Temp, areaMatrix_2Temp, areaMatrixTemp, VgjTemp, C0Temp, VgjPrimeTemp = self.P[-1], self.U[-1], self.H[-1], self.voidFraction[-1], self.rho[-1], self.rhoG[-1], self.rhoL[-1], self.xTh[-1], self.Dhfg[-1], self.f[-1], self.areaMatrix_1[-1], self.areaMatrix_2[-1], self.areaMatrix, self.Vgj[-1], self.C0[-1], self.VgjPrime[-1]
        T_water_temp = self.T_water
        self.Ultemp, self.Ugtemp, self.Reltemp = self.Ul, self.Ug, self.Rel
        self.Ul, self.Ug, self.Rel = [], [], []
        self.P[-1], self.U[-1], self.H[-1], self.voidFraction[-1], self.rho[-1], self.rhoG[-1], self.rhoL[-1], self.xTh[-1], self.Dhfg[-1], self.f[-1], self.areaMatrix_1[-1], self.areaMatrix_2[-1], self.areaMatrix, self.Vgj[-1], self.C0[-1], self.VgjPrime[-1] = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
        self.T_water = []
        for i in range(self.nCells):
            self.P[-1].append((Ptemp[i] + Ptemp[i+1])/2)
            self.U[-1].append((Utemp[i] + Utemp[i+1])/2)
            self.H[-1].append((Htemp[i] + Htemp[i+1])/2)
            self.voidFraction[-1].append((voidFractionTemp[i] + voidFractionTemp[i+1])/2)
            self.rho[-1].append((rhoTemp[i] + rhoTemp[i+1])/2)
            self.rhoG[-1].append((rhoGTemp[i] + rhoGTemp[i+1])/2)
            self.rhoL[-1].append((rhoLTemp[i] + rhoLTemp[i+1])/2)
            self.xTh[-1].append((xThTemp[i] + xThTemp[i+1])/2)
            self.Dhfg[-1].append((DhfgTemp[i] + DhfgTemp[i+1])/2)
            self.f[-1].append((fTemp[i] + fTemp[i+1])/2)
            self.areaMatrix_1[-1].append((areaMatrix_1Temp[i] + areaMatrix_1Temp[i+1])/2)
            self.areaMatrix_2[-1].append((areaMatrix_2Temp[i] + areaMatrix_2Temp[i+1])/2)
            self.Vgj[-1].append((VgjTemp[i] + VgjTemp[i+1])/2)
            self.C0[-1].append((C0Temp[i] + C0Temp[i+1])/2)
            self.VgjPrime[-1].append((VgjPrimeTemp[i] + VgjPrimeTemp[i+1])/2)
            self.T_water.append((T_water_temp[i] + T_water_temp[i+1])/2)
            self.Ul.append((self.Ultemp[i] + self.Ultemp[i+1])/2)
            self.Ug.append((self.Ugtemp[i] + self.Ugtemp[i+1])/2)
            self.Rel.append((self.Reltemp[i] + self.Reltemp[i+1])/2)

        self.z_mesh = np.linspace(self.Dz/2, self.height-self.Dz/2, self.nCells)

    #Compute the surface temperature
    def compute_T_surf(self):
        self.Pfin = self.P[-1]
        self.h_z = self.H[-1]
        self.T_surf = np.zeros(self.nCells)
        self.Hc = np.zeros(self.nCells)
        for i in range(self.nCells):
            #print(f'At axial slice = {i}, Pfin = {self.Pfin[i]}, h_z = {self.h_z[i]}')
            Pr_number = IAPWS97(P=self.Pfin[i]*10**-6, h=self.h_z[i]*10**-3).Liquid.Prandt
            Re_number = self.getReynoldsNumber(i)
            k_fluid = IAPWS97(P=self.Pfin[i]*10**-6, h=self.h_z[i]*10**-3).Liquid.k
            #print(f"At axial slice = {i}, computed Reynold # = {Re_number}, computed Prandt # = {Pr_number}, k_fluid = {k_fluid}")
            self.Hc[i] = (0.023)*(Pr_number)**0.4*(Re_number)**0.8*k_fluid/self.D_h[i]
            #print(f'self.Hc[i]: {self.Hc[i]}, \n self.q__[i]: {self.q__[i]} ,\n 2*np.pi*self.cladRadius: {2*np.pi*self.cladRadius}')
            self.T_surf[i] = ((self.q__[i]*self.flowArea)/(2*np.pi*self.cladRadius)/self.Hc[i]+self.T_water[i])
    
        return self.T_surf

    #Function to use the sous relaxation
    def sousRelaxation(self):

        for i in range(self.nFaces):
            self.voidFraction[-1][i] = self.voidFraction[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.voidFraction[-2][i]
            self.rho[-1][i] = self.rho[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.rho[-2][i]
            self.rhoG[-1][i] = self.rhoG[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.rhoG[-2][i]
            self.rhoL[-1][i] = self.rhoL[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.rhoL[-2][i]
            self.xTh[-1][i] = self.xTh[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.xTh[-2][i]
            self.Vgj[-1][i] = self.Vgj[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.Vgj[-2][i]
            self.C0[-1][i] = self.C0[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.C0[-2][i]
            self.VgjPrime[-1][i] = self.VgjPrime[-1][i] * self.sousRelaxFactor + (1-self.sousRelaxFactor)*self.VgjPrime[-2][i]
    
    #Function to merge 2 lists into one
    def mergeVar(self, U, P):
        VAR = np.concatenate((U, P))
        return VAR
    
    #Function to split a list into 2
    def splitVar(self, VAR):
        U = VAR[:self.nFaces]
        P = VAR[self.nFaces:]
        return U, P
    
    #Function to get the phases velocity
    def getPhasesVelocity(self):
        water = statesVariables(self.U[-1], self.P[-1], self.H[-1], self.voidFraction[-1], self.D_h, self.areaMatrix, self.DV, self.voidFractionCorrel, self.frfaccorel, self.P2Pcorel, self.Dz, self.q__, self.qFlow, self.fuelRadius, self.pitch/2)
        Ul = [water.getUl(i) for i in range(self.nCells)]
        Ug = [water.getUg(i) for i in range(self.nCells)]
        return Ul, Ug
    
    #Function to get the phases saturation enthalpy
    def createBoundaryEnthalpy(self):
        for i in range(self.nFaces):
            self.hlSat.append(self.getPhasesEnthalpy(i)[0])
            self.hgSat.append(self.getPhasesEnthalpy(i)[1]) 
 
    #get the Reynolds number
    def getReynoldsNumber(self, i):
        return (self.U[-1][i] * self.D_h[i] * self.rho[-1][i]) / IAPWS97(P=self.P[-1][i]*10**-6, x=0).Liquid.mu
     
