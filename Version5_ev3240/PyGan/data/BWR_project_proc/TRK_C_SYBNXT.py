#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2024                                                #
# Purpose     : Tracking for pin cells with SYB+NXT modules         #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_C_SYBNXT(namTRACK,namTF_EXC,namTRACK_SS,pyGEOM,pyGEOM_SS):
  # Lifo
  myLifo=lifo.new()


  myLifo.pushEmpty(namTRACK, "LCM")
  myLifo.pushEmpty(namTF_EXC, "BINARY")
  myLifo.pushEmpty(namTRACK_SS, "LCM")
  #myLifo.pushEmpty(namTF_EXC_SS, "BINARY")

  myLifo.push(pyGEOM)
  myLifo.push(pyGEOM_SS)

  # Execution
  trackBWR = cle2000.new('TRK_C_SYBNXT',myLifo,1)
  trackBWR.exec()

  # Recover
  myLifo.lib()

  pyTRACK = myLifo.node(namTRACK)
  pyTF_EXC = myLifo.node(namTF_EXC)

  pyTRACK_SS = myLifo.node(namTRACK_SS)
  #pyTF_EXC_SS = myLifo.node(namTF_EXC_SS)

  return [pyTRACK,pyTF_EXC,pyTRACK_SS]



