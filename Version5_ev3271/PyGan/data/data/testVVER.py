##############################################################################
#                                                                            #
# Description : PyGan scritp for VVER pin cell simulation with DRAGON5       #
# Author      : L. Fede                                                      #
# Date        : 2023                                                         #
#                                                                            #
##############################################################################
#
import lifo
import lcm
import cle2000
import numpy as np
from vverGEO import *
from vverMIX import *
from vverTRACK import *
from vverBU import *
from assertS import *

# Parameters of calculus
assembly='V1_S1'
Library='DLIB_J3p1p1'
Multicompo=1

Lifo3=lifo.new()

#################################
# ----------------------- vverGEO
#################################
pyGEOM = vverGEO("GEOM")

#################################
# ----------------------- vverMIX
#################################
pyMIX = vverMIX("LIBRARY",Library)

###################################
# ----------------------- vverTRACK
###################################
[pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC] = vverTRACK("UOX_TBH","FIG1.ps","TRACK","TF_EXC",pyGEOM)

#################################
# ------------------ BU iteration
#################################
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array([0.0,    15.0,    30.0,    50.0], dtype='f')
StepList['ListAutop'] = np.array([15.0,   30.0,    50.0], dtype='f')
StepList['ListCompo'] = np.array([0.0 ,   30.0,    50.0], dtype='f')
StepList.close() # close without erasing
pyCOMPO = vverBU("COMPO",pyMIX,pyTRACK,pyTF_EXC,StepList,Multicompo)

Keff_conv = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['K-INFINITY']
print("testVVER: step 1 converged K-effective=", Keff_conv)
Microlib = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]
Microlib.lib()
assertS(Microlib,'K-EFFECTIVE', 0, 1.2699327)

Keff_conv = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][2]['K-INFINITY']
print("testVVER: step 3 converged K-effective=", Keff_conv)
Microlib = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][2]
Microlib.lib()
assertS(Microlib,'K-EFFECTIVE', 0, 1.2717869)
print("test testVVER completed")
