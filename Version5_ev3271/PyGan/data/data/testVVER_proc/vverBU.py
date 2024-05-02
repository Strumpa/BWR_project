##############################################################################
#                                                                            #
# Description : PyGan scritp for VVER pin cell simulation with DRAGON5       #
# Author      : L. Fede                                                      #
# Date        : 2023                                                         #
#                                                                            #
##############################################################################
#
import lifo
import cle2000

def vverBU(namCOMPO,pyMIX,pyTRACK,pyTF_EXC,StepList,Multicompo):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(pyTRACK)
  myLifo.push(pyTF_EXC)
  myLifo.push(StepList)
  myLifo.push(Multicompo)
  myLifo.lib()

  # Execution 
  calcVVER = cle2000.new('vverBU',myLifo,1)
  calcVVER.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
