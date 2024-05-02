#
# Equivalence_SPH_APEX: incorporation of SPH factors in an APEX file
#
import lifo
import lcm
import cle2000
from assertS import *
import numpy as np

# recover the macro-geometry
ipLifo=lifo.new()
ipLifo.pushEmpty("Geom", "LCM")
ipLifo.push('UOX_5x5_TG6_sym8_multiDom.h5')
ipLifo.push(int) # ncals
GeomProc = cle2000.new('Umzone5x5',ipLifo,1)
GeomProc.exec()
ipGeom = ipLifo.node("Geom")
ncals = ipLifo.node(2)
print("test Equivalence_SPH_APEX: number of calculations=",ncals)
ipGeom.lib()

# perform transport-diffusion SPH equivalence
SPH_RT_Proc = cle2000.new('SPH_RT_Proc',ipLifo,1)
SPH_RT_Proc.exec()

# perform transport-SN SPH equivalence
SPH_RT_Proc = cle2000.new('SPH_SN_Proc',ipLifo,1)
SPH_RT_Proc.exec()
print("test Equivalence_SPH_APEX completed")
