# This file contains the classes and functions to create fields and give the values of the variables in the drift flux model reoslution file
# This class models the thermohydraulic behavior of a fluid flow in a system divided into multiple cells. Each cell is characterized by properties such as velocity (U), pressure (P), enthalpy (H), void fraction, and various geometric parameters. 
#The class employs different correlations to compute key thermophysical properties, and performs iterative updates to simulate two-phase flow dynamics.
# The file is used in the BWR/driftFluxModel/thermalHydraulicsTransitoire/THM_convection.py file
#Authors: Clément Huet
#Date: 2024-09-16

import numpy as np
from iapws import IAPWS97


class statesVariables():

    """
    Attributes:
    - nCells (int): Number of cells in the system.
    - U (array): Velocity values for each cell.
    - P (array): Pressure values for each cell.
    - H (array): Enthalpy values for each cell.
    - voidFraction (array): Initial void fraction values for each cell.
    - voidFractionCorrel (str): Correlation model used for void fraction calculations (e.g., 'modBestion', 'HEM1').
    - frfaccorel (str): Correlation model used for friction factor calculations (e.g., 'base', 'blasius', 'Churchill').
    - P2Pcorel (str): Correlation model used for two-phase pressure multiplier (e.g., 'base', 'HEM1').
    - D_h (float): Hydraulic diameter of the system.
    - flowArea (float): Flow area of the system.
    - DV (float): Differential volume between cells.
    - g (float): Gravitational acceleration (9.81 m/s²).
    - K_loss (float): Loss coefficient in the system.
    - Dz (float): Distance between cells.
    - Poro (float): Porosity of the system.

    Methods:
    - createFields(): Initializes temporary arrays for cell-specific properties such as densities, void fractions, friction factors, and geometric areas. Uses specific methods to populate these arrays.
    - updateFields(): Updates key properties (e.g., void fraction, densities) based on the selected void fraction correlation model. Calls corresponding methods for each model.
    - modBestion(): Updates void fraction and related properties using the "modBestion" correlation with iterative convergence.
    - HEM1(): Updates void fraction and related properties using the "HEM1" correlation with iterative convergence.
    - GEramp(): Updates void fraction and related properties using the "GEramp" correlation with iterative convergence.
    - EPRIvoidModel(): Updates void fraction and related properties using the "EPRIvoidModel" correlation with iterative convergence.
    - getDensity(i): Retrieves the liquid and vapor densities for a given cell using IAPWS97 thermodynamic models.
    - getQuality(i): Calculates the quality (phase fraction) based on enthalpy values.
    - getVoidFraction(i): Computes void fraction using different correlations, depending on the selected model.
    - getVgj(i): Calculates the drift velocity for a given cell based on the selected void fraction correlation.
    - getC0(i): Computes the slip ratio for a given cell based on the selected correlation model.
    - getFrictionFactor(i): Calculates the friction factor based on Reynolds number and the selected friction factor correlation.
    - getPhi2Phi(i): Computes the two-phase pressure multiplier for a given cell.
    - getAreas(i): Calculates the positive and negative flow areas for a given cell.
    - getPhasesEnthalpy(i): Retrieves the enthalpy values for the liquid and vapor phases at a given pressure.
    - getReynoldsNumber(i): Computes the Reynolds number for flow in a given cell.
    """

    def __init__(self, U, P, H, voidFraction, D_h, areaMatrix, DV, voidFractionCorrel, frfaccorel, P2Pcorel, Dz, q__, qFlow, rf, rw):
        
        self.nCells = len(U)
        self.U = U
        self.P = P
        self.H = H
        self.voidFraction = voidFraction
        self.voidFractionCorrel = voidFractionCorrel
        self.frfaccorel = frfaccorel
        self.P2Pcorel = P2Pcorel
        self.g = 9.81
        self.D_h = D_h
        self.areaMatrix = areaMatrix
        self.K_loss = 0#0.1
        self.Dz = Dz
        self.DV = DV
        self.q__ = q__
        self.qFlow = qFlow
        self.rf = rf
        self.rw = rw
        self.height = self.Dz * self.nCells

    #Create the array to store the temporary values of the variables
    def createFields(self):
        self.areaMatrixTEMP = np.ones(self.nCells)
        self.rholTEMP, self.rhogTEMP, self.rhoTEMP, self.voidFractionTEMP, self.DhfgTEMP, self.fTEMP, self.areaMatrix_1TEMP, self.areaMatrix_2TEMP, self.areaMatrix_2TEMP, self.VgjTEMP, self.C0TEMP, self.VgjPrimeTEMP = np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells)
        self.voidFractionTEMP = self.voidFraction
        self.xThTEMP = np.ones(self.nCells)
        for i in range(self.nCells):
            self.xThTEMP[i] = self.getQuality(i)
            self.areaMatrixTEMP[i] = self.areaMatrix[i]
            self.rholTEMP[i], self.rhogTEMP[i], self.rhoTEMP[i] = self.getDensity(i)
            self.C0TEMP[i] = self.getC0(i)
            self.VgjTEMP[i] = self.getVgj(i)
            self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
            self.DhfgTEMP[i] = self.getHfg(i)
            self.fTEMP[i] = self.getFrictionFactor(i)
            self.areaMatrix_1TEMP[i], self.areaMatrix_2TEMP[i] = self.getAreas(i)

    #Update the values of the variables using the selected void fraction correlation model
    def updateFields(self):

        self.xThTEMP = np.ones(self.nCells)
        for i in range(self.nCells):
            self.xThTEMP[i] = self.getQuality(i)
        
        if self.voidFractionCorrel == 'modBestion':
            self.modBestion()

        elif self.voidFractionCorrel == 'HEM1':
            self.HEM1()

        elif self.voidFractionCorrel == 'GEramp':
            self.GEramp()

        elif self.voidFractionCorrel == 'EPRIvoidModel':
            self.EPRIvoidModel()
        else:
            raise ValueError('Invalid void fraction correlation model')

    #Update the values of the variables using the selected void fraction correlation model: modBestion
    def modBestion(self):
        self.rholTEMP, self.rhogTEMP, self.rhoTEMP, self.voidFractionTEMP, self.DhfgTEMP, self.fTEMP, self.areaMatrix_1TEMP, self.areaMatrix_2TEMP, self.areaMatrix_2TEMP, self.VgjTEMP, self.C0TEMP, self.VgjPrimeTEMP = np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells)
        self.voidFractionOld = self.voidFraction
        self.Ul = np.ones(self.nCells)
        self.Ug = np.ones(self.nCells)
        self.Rel = np.ones(self.nCells)
        for i in range(self.nCells):
            self.rholTEMP[i], self.rhogTEMP[i], self.rhoTEMP[i] = self.getDensity(i)
            self.C0TEMP[i] = self.getC0(i)
            self.VgjTEMP[i] = self.getVgj(i)
            self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
            self.DhfgTEMP[i] = self.getHfg(i)
            voidFractionNew = self.getVoidFraction(i)
            self.voidFractionTEMP[i] = voidFractionNew
            self.rhoTEMP[i] = self.getDensity(i)[2]
            self.voidFractionTEMP[i] = voidFractionNew
            self.rhoTEMP[i] = self.getDensity(i)[2]
            self.fTEMP[i] = self.getFrictionFactor(i)
            self.areaMatrix_1TEMP[i], self.areaMatrix_2TEMP[i] = self.getAreas(i)
            self.Ul[i] = self.getUl(i)
            self.Ug[i] = self.getUg(i)
            self.Rel[i] = self.getReynoldsNumberLiquid(i)
    
    #Update the values of the variables using the selected void fraction correlation model: HEM1
    def HEM1(self):
        self.rholTEMP, self.rhogTEMP, self.rhoTEMP, self.voidFractionTEMP, self.DhfgTEMP, self.fTEMP, self.areaMatrix_1TEMP, self.areaMatrix_2TEMP, self.areaMatrix_2TEMP, self.VgjTEMP, self.C0TEMP, self.VgjPrimeTEMP = np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells)
        self.voidFractionOld = self.voidFraction
        self.Ul = np.ones(self.nCells)
        self.Ug = np.ones(self.nCells)
        self.Rel = np.ones(self.nCells)
        for i in range(self.nCells):
            self.rholTEMP[i], self.rhogTEMP[i], self.rhoTEMP[i] = self.getDensity(i)
            self.C0TEMP[i] = self.getC0(i)
            self.VgjTEMP[i] = self.getVgj(i)
            self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
            self.DhfgTEMP[i] = self.getHfg(i)
            voidFractionNew = self.getVoidFraction(i)
            self.voidFractionTEMP[i] = voidFractionNew
            self.rhoTEMP[i] = self.getDensity(i)[2]
            self.voidFractionTEMP[i] = voidFractionNew
            self.rhoTEMP[i] = self.getDensity(i)[2]
            self.fTEMP[i] = self.getFrictionFactor(i)
            self.areaMatrix_1TEMP[i], self.areaMatrix_2TEMP[i] = self.getAreas(i)
            self.Ul[i] = self.getUl(i)
            self.Ug[i] = self.getUg(i)
            self.Rel[i] = self.getReynoldsNumberLiquid(i)

    #Update the values of the variables using the selected void fraction correlation model: GEramp
    def GEramp(self):
        self.rholTEMP, self.rhogTEMP, self.rhoTEMP, self.voidFractionTEMP, self.DhfgTEMP, self.fTEMP, self.areaMatrix_1TEMP, self.areaMatrix_2TEMP, self.areaMatrix_2TEMP, self.VgjTEMP, self.C0TEMP, self.VgjPrimeTEMP = np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells)
        self.voidFractionOld = self.voidFraction
        self.Ul = np.ones(self.nCells)
        self.Ug = np.ones(self.nCells)
        self.Rel = np.ones(self.nCells)
        for i in range(self.nCells):
            self.rholTEMP[i], self.rhogTEMP[i], self.rhoTEMP[i] = self.getDensity(i)
            self.C0TEMP[i] = self.getC0(i)
            self.VgjTEMP[i] = self.getVgj(i)
            self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
            self.DhfgTEMP[i] = self.getHfg(i)
            for j in range(1000):
                voidFractionNew = self.getVoidFraction(i)
                if np.linalg.norm(voidFractionNew - self.voidFractionTEMP[i]) < 1e-3:
                    self.voidFractionTEMP[i] = voidFractionNew
                    self.rhoTEMP[i] = self.getDensity(i)[2]
                    self.C0TEMP[i] = self.getC0(i)
                    self.VgjTEMP[i] = self.getVgj(i)
                    self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
                    break
                elif j == 999:
                    raise ValueError('Convergence in update fields not reached')
                    break
                else:
                    self.voidFractionTEMP[i] = voidFractionNew
                    self.rhoTEMP[i] = self.getDensity(i)[2]
                    self.C0TEMP[i] = self.getC0(i)
                    self.VgjTEMP[i] = self.getVgj(i)
                    self.VgjPrimeTEMP[i] = self.getVgj_prime(i)

            self.fTEMP[i] = self.getFrictionFactor(i)
            self.areaMatrix_1TEMP[i], self.areaMatrix_2TEMP[i] = self.getAreas(i)
            self.Ul[i] = self.getUl(i)
            self.Ug[i] = self.getUg(i)
            self.Rel[i] = self.getReynoldsNumberLiquid(i)

    #Update the values of the variables using the selected void fraction correlation model: EPRIvoidModel
    def EPRIvoidModel(self):
        self.rholTEMP, self.rhogTEMP, self.rhoTEMP, self.voidFractionTEMP, self.DhfgTEMP, self.fTEMP, self.areaMatrix_1TEMP, self.areaMatrix_2TEMP, self.areaMatrix_2TEMP, self.VgjTEMP, self.C0TEMP, self.VgjPrimeTEMP = np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells), np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells),np.ones(self.nCells)
        self.voidFractionOld = self.voidFraction
        self.Ul = np.ones(self.nCells)
        self.Ug = np.ones(self.nCells)
        self.Rel = np.ones(self.nCells)
        for i in range(self.nCells):
            self.rholTEMP[i], self.rhogTEMP[i], self.rhoTEMP[i] = self.getDensity(i)
            self.C0TEMP[i] = self.getC0(i)
            self.VgjTEMP[i] = self.getVgj(i)
            self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
            self.DhfgTEMP[i] = self.getHfg(i)
            for j in range(1000):
                voidFractionNew = self.getVoidFraction(i)
                if np.linalg.norm(voidFractionNew - self.voidFractionTEMP[i]) < 1e-3:
                    self.voidFractionTEMP[i] = voidFractionNew
                    self.rhoTEMP[i] = self.getDensity(i)[2]
                    self.C0TEMP[i] = self.getC0(i)
                    self.VgjTEMP[i] = self.getVgj(i)
                    self.VgjPrimeTEMP[i] = self.getVgj_prime(i)
                    break
                elif j == 999:
                    raise ValueError('Convergence in update fields not reached')
                    break
                else:
                    self.voidFractionTEMP[i] = voidFractionNew
                    self.rhoTEMP[i] = self.getDensity(i)[2]
                    self.C0TEMP[i] = self.getC0(i)
                    self.VgjTEMP[i] = self.getVgj(i)
                    self.VgjPrimeTEMP[i] = self.getVgj_prime(i)

            self.fTEMP[i] = self.getFrictionFactor(i)
            self.areaMatrix_1TEMP[i], self.areaMatrix_2TEMP[i] = self.getAreas(i)
            self.Ul[i] = self.getUl(i)
            self.Ug[i] = self.getUg(i)
            self.Rel[i] = self.getReynoldsNumberLiquid(i)

    #Get the density of the liquid and vapor phases for a given cell
    def getDensity(self, i):
        vapor = IAPWS97(P = self.P[i]*(10**(-6)), x = 1)
        liquid = IAPWS97(P = self.P[i]*(10**(-6)), x = 0)
        rho_g = vapor.rho
        rho_l = liquid.rho
        rho = rho_l * (1 - self.voidFractionTEMP[i]) + rho_g * self.voidFractionTEMP[i]
        return rho_l, rho_g, rho
    
    #Get the quality  based on enthalpy values
    def getQuality(self, i):
        correl = "simple"
        if correl == 'simple':
            hl, hg = self.getPhasesEnthalpy(i)
            H = self.H[i]
            hfg = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).h
            if H*0.001 < hl:
                x = (H*0.001 - hl)/(hg - hl)
                return 0
            elif H*0.001 > hg:
                return 1
            
            elif H*0.001 <= hg and H*0.001 >= hl:
                #print(f'Quality: {(H*0.001 - hl)/(hg - hl)}')
                return (H*0.001 - hl)/(hg - hl)
        
        elif correl == 'EPRI':
            epriCorrel = '1'
            Xs = 0.05
            Xh = Xs / 2
            hl, hg = self.getPhasesEnthalpy(i)
            hl = hl
            hg = hg
            H = self.H[i]
            p = self.P[i]
            xeq = (hl - H*0.001) / (hl - hg)
            
            if epriCorrel == "1":
                if xeq >= Xs:
                    QUALITY = xeq
                else:
                    rhol = self.rholTEMP[i]
                    rhog = self.rhogTEMP[i]
                    #print(f'i: {i}, rho_l: {rhol}, rho_g: {rhog}, hauteur: {i*self.Dz}')
                    u = self.U[i]
                    muf = IAPWS97(P = p*(10**(-6)), x = 1).mu
                    Re = rhol * abs(u) * self.D_h[i] / muf

                    Cpf = IAPWS97(P = p*(10**(-6)), x = 1).cp
                    k_f = IAPWS97(P = p*(10**(-6)), x = 1).k
                    Pr = Cpf * muf / k_f

                
                    #qdp = (p * self.D_h[i])  /(2 * np.pi * self.rf * self.D_h[i])
                    qdp = self.q__[i] * self.DV / (2 * np.pi * self.rw * self.height)

                    # Calculate heat transfer coefficients
                    hb = np.exp(p / 4.35e6) / (22.7)**2 * 1000.0  # W/(m^2·K)
                    Chn = 0.2 / 4.0 * self.D_h[i] / self.rf
                    hhn = Chn * Re**0.662 * Pr * k_f / (self.D_h[i])
                    Cdb = (0.033 *self.areaMatrix[i] / (self.areaMatrix[i] +  np.pi * self.rf**2 + np.pi *self.rw**2) + 0.013)
                    hdb = Cdb * Re**0.8 * Pr**0.4 * k_f / (self.D_h[i])

                    # Intermediate calculations
                    tmp1 = 4.0 * hb * (hdb + hhn)**2
                    tmp2 = 2.0 * hdb**2 * (hhn + hdb / 2.0) + 8.0 * qdp * hb * (hdb + hhn)
                    tmp3 = qdp * (4.0 * hb * qdp + hdb**2)

                    # Calculate characteristic quality xd  ///////////////PROBLEM HERE
                    delta_h = hg - hl
                    xd = -Cpf / (delta_h) * (
                    (-tmp2 + np.sqrt(tmp2**2 - 4.0 * tmp1 * tmp3)) / (2.0 * tmp1)
                    )

                    # Determine quality based on xeq and xd
                    if xeq <= 0.0:
                        if xeq <= -xd:
                            QUALITY = 0.0       
                        else:
                            tmp1 = 1 + xeq / xd
                            tmp2 = tmp1**2
                            QUALITY = xd * tmp2 * (0.1 + 0.087 * tmp1 + 0.05 * tmp2)
                    elif Xh > xd:
                        if xeq >= 2 * xd:
                            QUALITY = xeq
                        else:
                            tmp1 = xeq / xd
                            QUALITY = xd * (0.237 + tmp1 * (0.661 + tmp1 * (0.153 + tmp1 * (-0.01725 - tmp1 * 0.0020625))))
                    else:
                        tmp1 = xeq / Xh
                        tmp2 = xd / Xh
                        tmp3 = 0.237 * tmp2
                        QUALITY = Xh * (tmp3 + tmp1 * (0.661 + tmp1 * (0.5085 - 0.3555 * tmp2 + tmp1 * (tmp3 - 0.25425 + tmp1 * (0.042375 - 0.0444375 * tmp2)))))

                if QUALITY >= 1.0:
                    return 0.99
                else:
                    return QUALITY

    #Get the void fraction for a given cell using different correlations for drift flulx model
    def getVoidFraction(self, i):
        correl = 'paths'
        if correl == 'simple':
            x_th = self.xThTEMP[i]
            rho_l = self.rholTEMP[i]
            rho_g = self.rhogTEMP[i]
            if x_th == 0:
                return 0
            elif x_th == 1:
                return 0.99
            else:
                return (x_th * rho_l)/(x_th * rho_l + (1 - x_th) * rho_g)
        elif correl == 'paths':
            x_th = self.xThTEMP[i]
            rho_l = self.rholTEMP[i]
            rho_g = self.rhogTEMP[i]
            u = self.U[i]
            V_gj = self.VgjTEMP[i]
            C0 = self.C0TEMP[i]
            if x_th == 0:
                return 0
            elif x_th == 1:
                return 0.99
            else:
                return x_th / (C0 * (x_th + (rho_g / rho_l) * (1 - x_th)) + (rho_g * V_gj) / (rho_l * u))
    
    #Get the drift velocity for a given cell based on the selected void fraction correlation
    def getVgj(self, i):
        if self.voidFractionCorrel == 'GEramp':
            if self.rhogTEMP[i] == 0:
                return 0
            if self.rholTEMP[i] == 0:
                return 0
            
            sigma = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).sigma
            if sigma == 0:
                return 0
            
            Vgj0 = ((self.g * sigma * (self.rholTEMP[i] - self.rhogTEMP[i]) / self.rholTEMP[i]**2)**0.25)

            if self.voidFractionTEMP[i] <= 0.65:
                return 2.9 * Vgj0
            elif self.voidFractionTEMP[i] > 0.65:
                return (2.9/0.35)*(1-self.voidFractionTEMP[i]) * Vgj0
        
        if self.voidFractionCorrel == 'modBestion':
            if self.rhogTEMP[i] == 0:
                return 0
            if self.rholTEMP[i] == 0:
                return 0
            return 0.188 * np.sqrt(((self.rholTEMP[i] - self.rhogTEMP[i]) * self.g * self.D_h[i] ) / self.rhogTEMP[i] )
        
        if self.voidFractionCorrel == 'EPRIvoidModel':
            if self.rhogTEMP[i] == 0:
                return 0
            if self.rholTEMP[i] == 0:
                return 0
            sigma = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).sigma
            Vgj = (np.sqrt(2)*(self.g * sigma * (self.rholTEMP[i] - self.rhogTEMP[i]) / self.rholTEMP[i]**2)**0.25) * (1 + self.voidFractionTEMP[i])**(3/2)
            return Vgj
        
        if self.voidFractionCorrel == 'HEM1':
            return 0
            
    #Get the distribution parameter for a given cell based on the selected correlation model
    def getC0(self, i):
        if self.voidFractionCorrel == 'GEramp':
            rho_g = self.rhogTEMP[i]
            rho_l = self.rholTEMP[i]
            if rho_g == 0:
                return 0
            if rho_l == 0:
                return 0
            if self.voidFractionTEMP[i] <= 0.65:
                return 1.1
            elif self.voidFractionTEMP[i] > 0.65:
                return 1 + (0.1/0.35)*(1-self.voidFractionTEMP[i])
        
        if self.voidFractionCorrel == 'modBestion':
            rho_g = self.rhogTEMP[i]
            rho_l = self.rholTEMP[i]
            if rho_g == 0:
                return 0
            if rho_l == 0:
                return 0
            return 1.2 - 0.2*np.sqrt(rho_g / rho_l)
        
        if self.voidFractionCorrel == 'EPRIvoidModel':
            rho_g = self.rhogTEMP[i]
            rho_l = self.rholTEMP[i]
            Pc = 22060000
            P = self.P[i]
            Re = self.getReynoldsNumber(i)
            C1 = (4 * Pc**2)/(P*(Pc - P))
            k1 = min(0.8, 1/(1 + np.exp(-Re /60000)))
            k0 = k1 + (1-k1) * (rho_g / rho_l)**2
            r = (1+1.57*(rho_g/rho_l))/(1-k1)
            C0 = (((k0 + (1 - k0) * (self.voidFractionTEMP[i]**r))**(-1)) * ((1 - np.exp(-C1 * self.voidFractionTEMP[i]))/(1 - np.exp(-C1))))
            return C0

        if self.voidFractionCorrel == 'HEM1':
            return 1
            
    #Get the drift velocity for a given cell based on the selected void fraction correlation
    def getVgj_prime(self, i):
        U = self.U[i]
        C0 = self.C0TEMP[i]
        Vgj = self.VgjTEMP[i]
        Vgj_prime = Vgj + (C0-1) * U
        return Vgj_prime
    
    def getHfg(self, i):
        vapor = IAPWS97(P = self.P[i]*(10**(-6)), x = 1)
        liquid = IAPWS97(P = self.P[i]*(10**(-6)), x = 0)
        return (vapor.h - liquid.h)
    
    def getFrictionFactor(self, i):
        U = self.U[i]
        P = self.P[i]
        Re = self.getReynoldsNumberLiquid(i)
        #print('voidFraction list: ', self.voidFractionTEMP)
        if (self.voidFractionTEMP[i]<0.002):
            #print(f"Void fraction: {self.voidFractionTEMP[i]}")
            return 0.316 * Re**(-0.25)
        if self.frfaccorel == 'base': #Validated
            return 0.003
        elif self.frfaccorel == "null": #Validated
            return 0
        elif self.frfaccorel == 'blasius': #Validated
            #return 0.316 * Re**(-0.25)
            return 0.079 * Re**(-0.25)
            #return 0.316 * Re**(-0.25)
        elif self.frfaccorel == 'Churchill': #Validated
            #old 0.4 When Ra increased pressure drop increase
            #Ra = 0.8 * (10**(-6)) #Roughness
            #Ra = 1.0 * (10**(-6)) #Roughness
            Ra = 0.4 * (10**(-6)) #Roughness
            R = Ra / self.D_h[i]
            frict=8*(((8.0/Re)**12)+((2.475*np.log(((7/Re)**0.9)+0.27*R))**16+(37530/Re)**16)**(-1.5))**(1/12)
            return frict
        elif self.frfaccorel == 'Churchill_notOK':
            Re = self.getReynoldsNumberLiquid(i)
            B = (37530/Re)**16
            A = (2.475*np.log(1/(((7/Re)**0.9)+0.27*(0.4/self.D_h[i]))))**16
            f  = 8*(((8/Re)**12) + (1/(A+B)**1.5))**(1/12)
            return f
        else:
            raise ValueError('Invalid friction factor correlation model')

    #Get the two-phase pressure multiplier for a given cell
    def getPhi2Phi(self, i):
        x_th = self.xThTEMP[i]
        rho_l = self.rholTEMP[i]
        rho_g = self.rhogTEMP[i]
        rho = self.rhoTEMP[i]
        P = self.P[i]
        epsilon = self.voidFractionTEMP[i]
        if epsilon <= 0.002:
            return 1
        if self.P2Pcorel == 'base': #Validated
            phi2phi = 1 + 3*epsilon
        elif self.P2Pcorel == 'lockhartMartinelli':
            #return np.sqrt(self.lockhartMartinelli(i))
            return self.lockhartMartinelli(i)
        elif self.P2Pcorel == 'HEM1': #Validated
            phi2phi = (rho/rho_l)*((rho_l/rho_g)*x_th + +1)
        elif self.P2Pcorel == 'HEM2': #Validated    
            m = IAPWS97(P = P*(10**(-6)), x = 0).mu / IAPWS97(P = P*(10**(-6)), x = 1).mu
            phi2phi = (rho/rho_l)*((m-1)*x_th + 1)*((rho_l/rho_g)*x_th + +1)**(0.25)
        elif self.P2Pcorel == 'MNmodel': #Validated
            phi2phi = (1.2 * (rho_l/rho_g -1)*x_th**(0.824) + 1)*(rho/rho_l)
            #print(f'Phi2phi : {phi2phi}')
        else:
            raise ValueError('Invalid two-phase pressure multiplier correlation model')
        return phi2phi
    
    #Get the positive and negative flow areas for a given cell (take into account the two-phase pressure multiplier and friction factor)
    def getAreas(self, i):
        if self.voidFractionTEMP[i] > - 0.001:
            A_chap_pos = self.areaMatrix[i-1] +  (self.getPhi2Phi(i-1)/4) * ((self.fTEMP[i-1] / self.D_h[i-1]) + (self.K_loss / self.Dz)) * self.DV
            A_chap_neg = self.areaMatrix[i] - (self.getPhi2Phi(i)/4) * ((self.fTEMP[i] / self.D_h[i]) + (self.K_loss / self.Dz)) * self.DV
            
        return A_chap_pos, A_chap_neg

    #Get the enthalpy values for the liquid and vapor phases at a given pressure
    def getPhasesEnthalpy(self, i):
        P = self.P[i]
        vapor = IAPWS97(P = P*(10**(-6)), x = 1)
        liquid = IAPWS97(P = P*(10**(-6)), x = 0)
        return liquid.h, vapor.h
    
    #Get the Reynolds number for flow in a given cell
    def getReynoldsNumber(self, i):
        U = self.U[i]
        rho = self.rhoTEMP[i]
        P = self.P[i]
        alpha = self.voidFractionTEMP[i]
        ml = IAPWS97(P = P*(10**(-6)), x = 0).mu
        mv = IAPWS97(P = P*(10**(-6)), x = 1).mu
        m = (mv * ml) / ( ml * (1 - alpha) + mv * alpha )
        
        return rho * abs(U) * self.D_h[i] / m
    
    #Get the Reynolds number for the liquid phase in a given cell
    def getReynoldsNumberLiquid(self, i):
        Ul = self.getUl(i)
        rho = self.rholTEMP[i]
        P = self.P[i]
        m = IAPWS97(P = P*(10**(-6)), x = 0).mu
        return rho * abs(Ul) * self.D_h[i] / m

    #Get the Reynolds number for the vapor phase in a given cell
    def getReynoldsNumberVapor(self, i):
        Ug = self.getUg(i)
        rho = self.rhogTEMP[i]
        P = self.P[i]
        m = IAPWS97(P = P*(10**(-6)), x = 1).mu
        return rho * abs(Ug) * self.D_h[i] / m
    
    #Get the liquid velocity in a given cell
    def getUl(self, i):
        return self.U[i] - (self.voidFractionTEMP[i] / ( 1 - self.voidFractionTEMP[i])) * (self.rhogTEMP[i] / self.rhoTEMP[i]) * self.VgjPrimeTEMP[i]
    
    #Get the vapor velocity in a given cell
    def getUg(self, i):
        return self.U[i] + (self.rholTEMP[i] / self.rhoTEMP[i]) * self.VgjPrimeTEMP[i]

    #Get the two-phase pressure multiplier for a given cell
    def lockhartMartinelli(self, i):
        Ul = self.getUl(i)
        Ug = self.getUg(i)
        Um = self.U[i]
        rhom = self.rhoTEMP[i]
        rho_l = self.rholTEMP[i]
        rho_g = self.rhogTEMP[i]
        X = np.sqrt(rho_l/rho_g)*(Ul/Ug)
        return 1 + 1/np.sqrt(X) + 20/X
    
        """ mul = IAPWS97(P = self.P[i]*(10**(-6)), x = 0).Liquid.mu
        mulg = IAPWS97(P = self.P[i]*(10**(-6)), x = 1).Vapor.mu

        Rel = self.getReynoldsNumberLiquid(i)
        Reg = self.getReynoldsNumberVapor(i)
        fliq = 0.079 * (Rel)**(-0.25)
        fgas = 0.079 * (Reg)**(-0.25)

        return ((1.2*(rho_l/rho_g - 1)*self.xThTEMP[i]**0.824 + 1)*(rhom/rho_l)*(rho_l/rho_g))**2 #*self.xThTEMP[i] + 1)**0.25
        X = (fliq*rho_l*Ul**2)/(fgas*rho_g*Ug**2)
        return 1+20/X """
    
    #get the velocity values for a given cell
    def getVelocity(self):
        Ul = []
        Ug = []
        for i in range(self.nCells):
            Ul.append(self.getUl(i))
            Ug.append(self.getUg(i))
        return Ul, Ug