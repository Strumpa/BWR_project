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

def vverTRACK(namUOX_TBH,namFIG,namTRACK,namTF_EXC,pyGEOM):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namUOX_TBH, "ASCII")
  myLifo.pushEmpty(namFIG, "ASCII")
  myLifo.pushEmpty(namTRACK, "LCM")
  myLifo.pushEmpty(namTF_EXC, "BINARY")
  myLifo.push(pyGEOM)

  # Execution
  trackVVER = cle2000.new('vverTRACK',myLifo,1)
  trackVVER.exec()

  # Recover
  myLifo.lib()
  pyUOX_TBH = myLifo.node(namUOX_TBH)
  pyFIG = myLifo.node(namFIG)
  pyTRACK = myLifo.node(namTRACK)
  pyTF_EXC = myLifo.node(namTF_EXC)
  return [pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC]
