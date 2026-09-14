#This class implements several numerical methods for solving linear systems of the form Ax = b. 
#It supports different solvers and preconditioning techniques for improving convergence, and it is designed to handle large systems.
# The file is used in the BWR/driftFluxModel/thermalHydraulicsTransitoire/THM_convection.py file
#Authors: Clément Huet
#Date: 2024-09-16
#uses : - newton method
#       - Biconjhugate gradient stabilized method
#       - Gauss Siedel method
#       - Biconjugate gradient method
#       - matrix inversion method
#       - preconditionner method ILU and SPAI
#       - scipy.sparse library 


import numpy as np
from iapws import IAPWS97


class numericalResolution():

    """"
    Attributes:
    - `A`: Matrix A of the system Ax = b.
    - `b`: Right-hand side vector of the system Ax = b.
    - `x0`: Initial guess for the solution vector.
    - `tol`: Tolerance for convergence criteria.
    - `maxIter`: Maximum number of iterations allowed for iterative solvers.
    - `numericalMethod`: String indicating the chosen numerical method ('FVM', 'BiCStab', 'GaussSiedel', 'BiCG').

    Methods:
    - `resolve`: Chooses the appropriate solver based on the `numericalMethod` attribute.
    - `resolveInversion`: Directly solves the system using matrix inversion with preconditioning (LU factorization).
    - `resolveGaussSiedel`: Implements the Gauss-Seidel iterative solver with preconditioning.
    - `resolveBiConjugateGradient`: Implements the BiConjugate Gradient method (BiCG) with preconditioning.
    - `resolveBiCGStab`: Implements the BiCGStab method, an improved version of BiCG with stabilization for better convergence.
    - `preconditionner`: Constructs the preconditioner matrix using either ILU (Incomplete LU) or SPAI (Sparse Approximate Inverse) methods.
    - `scalBIC`: Computes the scalar product of two vectors using the conjugate transpose of the first vector.
    """

    """
    
    Numerical Methods:
    1. `FVM`: Solves the system using matrix inversion with preconditioning.
    2. `GaussSiedel`: Applies the Gauss-Seidel iterative solver.
    3. `BiCG`: Uses the BiConjugate Gradient method for solving non-symmetric or indefinite matrices.
    4. `BiCStab`: Applies the BiCGStab (BiConjugate Gradient Stabilized) method to ensure faster and more stable convergence.

    Preconditioning:
    - `preconditionner`: A function to improve convergence by applying either ILU or SPAI techniques to approximate the inverse of matrix A.
    - ILU: Incomplete LU factorization (used by default).
    - SPAI: Sparse Approximate Inverse (optional for special cases).

    Usage:
    - Create an instance of `numericalResolution` by providing the system's matrix `A`, vector `b`, initial guess `x0`, tolerance `tol`, maximum iterations `maxIter`, and the desired numerical method.
    - Call the `resolve` method to solve the system using the chosen solver.
    - Preconditioning is automatically applied to improve convergence for iterative methods.
    """

    def __init__(self, obj, x0, tol, maxIter, numericalMethod):
        self.A = obj.A
        self.b = obj.D
        self.x0 = x0
        self.tol = tol
        self.maxIter = maxIter
        self.numericalMethod = numericalMethod
        self.n = self.A.shape[0]

        self.resolve()

    #Choose the appropriate solver based on the numericalMethod attribute
    def resolve(self):

        if self.numericalMethod == 'FVM':
            self.x = self.resolveInversion()
        elif self.numericalMethod == 'BiCGStab':
            self.x = self.resolveBiCGStab()
        elif self.numericalMethod == 'GaussSiedel':
            self.x = self.resolveGaussSiedel()
        elif self.numericalMethod == 'BiCG':
            self.x = self.resolveBiConjugateGradient()
        else:
            raise ValueError('Numerical method not recognized')

    #Directly solve the system using matrix inversion with preconditioning (LU factorization)
    def resolveInversion(self):
        L, U = self.preconditionner(self.A)
        M = np.dot(L,U)
        VAR = np.linalg.solve(self.A,self.b)
        return VAR

    #Implement the Gauss-Seidel iterative solver with preconditioning
    def resolveGaussSiedel(self):
        L, U = self.preconditionner(self.A)
        M = L @ U
        MStar = np.linalg.inv(M)

        A = np.dot(MStar,self.A)
        b = np.dot(MStar, self.b)

        x = self.x0
        n = len(x)
        err0 = 0
        Ax = np.zeros(n)
        for i in range(n):
            Ax[i] = np.dot(A[i, :], x)

        for i in range(n):
            err = b[i] - Ax[i]
            err0 += err**2
        err0 = np.sqrt(err0)

        for m in range(1,1000):
            esum = 0
            for i in range(n):
                x_old = x[i]
                sum = 0
                for j in range(n):
                    if j != i:
                        sum += A[i, j] * x[j]
                x[i] = (b[i] - sum) / A[i, i]
                esum += (x[i] - x_old)**2
            
            erout = np.sqrt(esum)
            if np.sqrt(erout) <= self.tol:
                break

        return x
    
    #Construct the preconditioner matrix using either ILU or SPAI techniques
    def preconditionner(self, A):
        m = 50 
        ILU = True
        SPAI = False

        if SPAI == True:
            """Perform m step of the SPAI iteration."""
            from scipy.sparse import identity
            from scipy.sparse import diags
            from scipy.sparse.linalg import onenormest
            from scipy.sparse import csr_array

            A = csr_array(A)
            ident = identity(n, format='csr')
            alpha = 2 / onenormest(A @ A.T)
            M = alpha * A
                
            for index in range(m):
                C = A @ M
                G = ident - C
                AG = A @ G
                trace = (G.T @ AG).diagonal().sum()
                alpha = trace / np.linalg.norm(AG.data)**2
                M = M + alpha * G
            return M.todense()
        
        if ILU == True:
            #Initialize L and U as copies of A
            L = np.eye(self.n)
            U = np.copy(A)

            #Perform ILU factorization
            for i in range(1,self.n):
                for k in range(i):
                    if U[k,k] != 0:
                        L[i,k] = U[i,k] / U[k,k]
                        U[i,k:] = U[i,k:] - L[i,k] * U[k,k:]

            return L,U

    #Implement the BiConjugate Gradient method (BiCG) with preconditioning
    def resolveBiConjugateGradient(self): #Notworking

        L,U = self.preconditionner(self.A)
        M = L @ U
        self.condNUMBERB = np.linalg.cond(self.A)
        self.condNUMBER = np.linalg.cond(np.dot(np.linalg.inv(M), self.A))
        
        MStar = np.linalg.inv(M)
        AStar = np.transpose(self.A)
        x0 = self.x0
        r0 = self.b - np.dot(self.A,x0)
        r0Star = np.transpose(self.b) - np.dot(np.transpose(x0),AStar)
        p0 = np.dot(MStar,r0)
        p0Star = np.dot(r0Star, MStar)
        x0Star = np.transpose(x0)
        for k in range(100000):
            alpha = np.dot(np.dot(r0Star,MStar), r0) / np.dot(p0Star,np.dot(self.A, p0))
            alphaBar = np.conjugate(alpha)
            x = x0 + alpha * p0
            xStar  = x0Star + alphaBar * p0Star
            r = r0 - alpha * np.dot(self.A, p0)
            rStar = r0Star - alphaBar * np.dot(p0Star, AStar)
            if np.linalg.norm(r) < self.tol:
                break
            if k == 99999:
                raise ValueError('BiConjugateGradient did not converge')
            beta = np.dot(np.dot(rStar, MStar),r) / np.dot(np.dot(r0Star, MStar), r0)
            betaBar = np.conjugate(beta)
            p = np.dot(MStar,r) + beta * p0
            pStar = np.dot(rStar, MStar) + betaBar * p0Star

            r0 = r
            r0Star = rStar
            x0 = x
            x0Star = xStar
            p0 = p
            p0Star = pStar

        return x

    #Compute the scalar product of two vectors using the conjugate transpose of the first vector
    def scalBIC(self, a, b):
        return np.dot(np.transpose(a), b)
    
    #Implement the BiCGStab method, an improved version of BiCG with stabilization for better convergence
    def resolveBiCGStab(self): #Not working
        K1, K2 = self.preconditionner(self.A)
        K1Star = np.linalg.inv(self.A)
        K2Star = K1Star.copy()
        x0 = self.x0
        r0 = self.b - np.dot(self.A,self.x0)
        r0Star = r0
        rho0 = self.scalBIC(r0Star, r0)
        p0 = r0

        for k in range(100000):
            y = K2Star @ K1Star @ p0
            v = self.A @ y
            alpha = rho0 / self.scalBIC(r0Star, v)
            h = x0 + alpha * y
            s = r0 - alpha * v
            r = self.A @ h - self.b
            if np.linalg.norm(r) < self.tol:
                x = h
                break
            z = K2Star @ K1Star @ s
            t = self.A @ z
            omega = self.scalBIC(K1Star @ t, K1Star @ s) / self.scalBIC(K1Star @ t, K1Star @ t)
            x = h + omega * z
            r = s - omega * t
            res = self.A @ x - self.b
            if np.linalg.norm(res) < self.tol:
                break
            rho = self.scalBIC(r0Star, r)
            beta = (rho / rho0) * (alpha / omega)
            p = r + beta * (p0 - omega * v)
            
            if k == 99999:
                raise ValueError('BiCGStab did not converge')
            rho0 = rho
            r0 = r
            x0 = x
            p0 = p

        #print(f'end of resolveBiCGStab with k = {k}')
        return x
 
#Class to fill the matrix A and D and solve the differential equation of enthalpy
class FVM:
    def __init__(self, A00, A01, Am0, Am1, D0, Dm1, N_vol, H):
        self.A00 = A00
        self.A01 = A01
        self.Am0 = Am0
        self.Am1 = Am1
        self.D0 = D0
        self.Dm1 = Dm1
        self.N_vol = N_vol
        self.H = H
        self.dz = H / N_vol
        self.z = np.linspace(0, H, self.N_vol)
        self.A, self.D = np.eye(self.N_vol), np.zeros(self.N_vol)

    #function to set the matrix A and D
    def set_ADi(self, i, ci, ai, bi, di):
        self.A[i, i-1:i+2] = [ci, ai, bi]
        self.D[i] = di
        return
    
    #function to set the boundary conditions
    def set_CL(self, A0, Am1, D0, Dm1):
        self.A[0], self.A[-1] = A0, Am1
        self.D[0], self.D[-1] = D0, Dm1
        return
    
    #function to solve the system of equations
    def resoudre_h(self):
        return np.linalg.solve(self.A, self.D)
    
    #function to set the transient parameters
    """ def set_transitoire(self, t_tot, Tini, dt):
        self.t_tot, self.dt = t_tot, dt           
        self.N_temps = round(self.t_tot / self.dt) # pas de temps (timesteps), il faut etre un nombre entier
        self.T = np.zeros((self.N_temps, self.N_vol)) # tableau 2D de temperature. 
        self.T[0] = Tini # Tini est une liste de temperature initiale """

    #function to set the boundary conditions
    def boundaryFilling(self):
        A0, Am1 = np.zeros(self.N_vol), np.zeros(self.N_vol)
        A0[:2] = [self.A00, self.A01]
        Am1[-2:] = [self.Am0, self.Am1]
        D0 = self.D0
        Dm1 = self.Dm1
        self.set_CL(A0, Am1, D0, Dm1)

    def fillingOutside(self, i, j, ci, ai, bi):
        self.A[i, j-1:j+2] = [ci, ai, bi]
        return 
    
    def fillingOutsideBoundary(self, i, j, ai, bi):
        self.A[i, j:j+2] = [ai, bi]
        return 
    
    #function to solve the differential equation of enthalppie
    def differential(self):
        self.boundaryFilling()
        self.h = self.resoudre_h()
        return self.h
    
    #function to calculate the temperature of the surface of the fuel using the fluid parameters, the distribution of enthalpie and the heat flux
    def verticalResolution(self):
        self.differential()
