### PyGan script used as python driver for THM calls and coupled neutronics-TH calculations
# Goal : 1) initialize call to THM
#        2) Call THM as long as convergence is not reached


import lifo
import lcm
import cle2000


def initializeModel(name_compo, Iz, fuel_mass, height, pitch, powi, guessPow, pOutlet, tInlet, massFlow, area, fuelRadius, gapRadius, cladRadius, THM_correlation, THM_relax_f):
    """
    Initial call to THM: with guess axial power profile,
    Initial calls to GEO:, USPLIT:, RESINI: etc 

    returns : 
        - TH data to be used in the neutronics solution,
        - FuelMap to be used in coupled calculations,
        - Material indexation to be used in coupled calculations,
        
    """
    # Initialize stack with PARAMETERS first : 
    # Fmap Matex Cpo Track THData guessPow
    ipLifoIni = lifo.new()
    ipLifoIni.pushEmpty("Fmap", "LCM")
    ipLifoIni.pushEmpty("Matex", "LCM")
    ipLifoIni.pushEmpty("Cpo", "LCM")
    ipLifoIni.pushEmpty("Track", "LCM")
    ipLifoIni.pushEmpty("Thm", "LCM")
    ipLifoIni.push(guessPow) # Guess power profile for first THM: run
    # Add parameters
    ipLifoIni.push(name_compo) # Name of the composition
    ipLifoIni.push(Iz) # Number of axial nodes
    ipLifoIni.push(fuel_mass) # Fuel mass
    ipLifoIni.push(height) # Active core height
    ipLifoIni.push(pitch) # Pitch of the fuel assembly / channel
    ipLifoIni.push(powi) # Total reactor power (MW)
    ipLifoIni.push(pOutlet) # Outlet pressure (Pa)
    ipLifoIni.push(tInlet) # Inlet temperature (K)
    ipLifoIni.push(massFlow) # Mass flow rate (kg/s)
    ipLifoIni.push(area) # Area of the channel (m^2)
    ipLifoIni.push(fuelRadius) # Fuel radius (m)
    ipLifoIni.push(gapRadius) # Gap radius (m)
    ipLifoIni.push(cladRadius) # Cladding radius (m)
    ipLifoIni.push(THM_correlation) # THM correlation to be used
    ipLifoIni.push(THM_relax_f) # THM under-relaxation factor

    # Call cle2000 procedure to initialize THM related LCM data structures
    Initialization = cle2000.new('IniDONJON',ipLifoIni,1)
    Initialization.exec()

    # Get the results
    Fmap = ipLifoIni.node("Fmap")
    Matex = ipLifoIni.node("Matex")
    Cpo = ipLifoIni.node("Cpo")
    Track = ipLifoIni.node("Track")
    Thm = ipLifoIni.node("Thm")
    
    # empty the Lifo stack for IniDONJON
    while ipLifoIni.getMax() > 0:
        ipLifoIni.pop();

    return Fmap, Matex, Cpo, Track, Thm


def solNeutronicsDONJON(iter, powi, Nz, Fmap, Matex, Cpo, Track, Flux):
    """
    Run DONJON5 neutronics calculation using TRIVAC5 Finite Element Solver
    Input parameters :
        - iter (integer) : iteration number
        - powi (double) : total reactor power (MW)
        - Nz (integer) : number of axial nodes
    
    LCM Inputs : 
        - FuelMap (Fmap)
        - Material indexation (Matex)
        - Cpo (MULTICOMPO to read XS from)
        - Track (FEM numbering and mesh to be used in TRIVAC:)
        - Flux (Neutron flux distribution)
        - Power (Power distribution)

    returns : 
        - Power distribution (in W/m^3) to be used in THM calculation
        - FuelMap to be used in coupled calculations, 
    """


    """
    PARAMETER  Fmap Matex Flux Power Cpo Track THData  :: ; ; 

    INTEGER iter ;
    :: >>iter<< ;
    DOUBLE Dpowi  ;
    :: >>Dpowi<< ;
    INTEGER Nz ;
    :: >>Nz<< ;
    """
    ipLifoN = lifo.new()
    ipLifoN.push(Fmap) # Fuel map
    ipLifoN.push(Matex) # Material indexation
    if iter > 1:
        ipLifoN.push(Flux) # Neutron flux distribution
    else:
        ipLifoN.pushEmpty("Flux", "LCM")
    ipLifoN.pushEmpty("Power", "LCM") # Power distribution : to store output Power distribution
    ipLifoN.push(Cpo) # MULTICOMPO to read XS from
    ipLifoN.push(Track) # FEM numbering and mesh to be used in TRIVAC:
    ipLifoN.push(iter) # Iteration number
    ipLifoN.push(powi) # Total reactor power (MW)
    ipLifoN.push(Nz) # Number of axial nodes


    # Call cle2000 procedure to initialize THM related LCM data structures
    Neutronics = cle2000.new('solNeutro', ipLifoN, 1)
    Neutronics.exec()
    
    # Get the results
    Fmap = ipLifoN.node("Fmap")
    Matex = ipLifoN.node("Matex")
    Track = ipLifoN.node("Track")
    Power = ipLifoN.node("Power")
    Flux = ipLifoN.node("Flux")

    # empty the Lifo stack for solNeutro
    while ipLifoN.getMax() > 0:
        ipLifoN.pop();

    return Power, Flux, Fmap, Matex, Track

def solTHM(iter, Fmap, Thm, pOutlet, tInlet, massFlowRate, area, fuelRadius, gapRadius, cladRadius, correlation, relaxation_factor):
    """
    Simplified Thermal Hydraulics calculation using THM:
    LCM Inputs :
        - THM data (Thm)
        - FuelMap (Fmap)
        
    Input parameters :
        - pOutlet (real) : outlet pressure (Pa)
        - tInlet (real) : inlet temperature (K)
        - massFlowRate (real) : mass flow rate (kg/s)
        - area (real) : area of the channel (m^2)
        - fuelRadius (real) : fuel radius (m)
        - gapRadius (real) : gap radius (m)
        - cladRadius (real) : cladding radius (m)
        - iter (integer) : iteration number
        - correlation (str) : THM correlation to be used for sub-cooled boiling treatment
        - relaxation_factor (real) : THM under-relaxation factor
    returns :
        - Thm (THM data to be used in coupled calculations)
        - Fmap (FuelMap to be used in coupled calculations)
        
    
    """

    ipLifoTHM = lifo.new()
    ipLifoTHM.push(Thm)
    ipLifoTHM.push(Fmap)
    ipLifoTHM.push(pOutlet)
    ipLifoTHM.push(tInlet)
    ipLifoTHM.push(massFlowRate)
    ipLifoTHM.push(area)
    ipLifoTHM.push(fuelRadius)
    ipLifoTHM.push(gapRadius)
    ipLifoTHM.push(cladRadius)
    ipLifoTHM.push(iter)
    ipLifoTHM.push(correlation)  # THM sub-cooled boiling correlation to be used
    ipLifoTHM.push(relaxation_factor)  # THM under-relaxation factor
    # Call cle2000 procedure to initialize THM related LCM data structures
    THM = cle2000.new('solTHM', ipLifoTHM, 1)
    THM.exec()
    # Get the results
    Thm = ipLifoTHM.node("Thm")
    Fmap = ipLifoTHM.node("Fmap")
    

    # empty the Lifo stack for solTHM
    while ipLifoTHM.getMax() > 0:
        ipLifoTHM.pop();

    return Fmap, Thm



    





