#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, adapted by R. Guasch                       #
# Date        : 2023/2024                                           #
# Purpose     : Depletion calculus for pin cell                     #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_C(namCOMPO,pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(pyTRACK)
  myLifo.push(pyTF_EXC)
  myLifo.push(pyTRACK_SS)
  myLifo.push(pyTF_EXC_SS)
  myLifo.push(StepList)
  myLifo.push(Multicompo)
  myLifo.push(name_geom)
  myLifo.push(name_compo)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('BU_C',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
