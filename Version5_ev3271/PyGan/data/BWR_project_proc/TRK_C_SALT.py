#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, adapted by R.Guasch                        #
# Date        : 2024                                                #
# Purpose     : Tracking for pin cells with SALT module             #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_C_SALT(namUOX_TBH,namFIG,namTRACK,namTF_EXC,namUOX_TBH_SS,namFIG_SS,namTRACK_SS,namTF_EXC_SS,pyGEOM,pyGEOM_SS):
  # Lifo
  myLifo=lifo.new()

  myLifo.pushEmpty(namUOX_TBH, "ASCII")
  myLifo.pushEmpty(namFIG, "ASCII")
  myLifo.pushEmpty(namTRACK, "LCM")
  myLifo.pushEmpty(namTF_EXC, "BINARY")

  myLifo.pushEmpty(namUOX_TBH_SS, "ASCII")
  myLifo.pushEmpty(namFIG_SS, "ASCII")
  myLifo.pushEmpty(namTRACK_SS, "LCM")
  myLifo.pushEmpty(namTF_EXC_SS, "BINARY")

  myLifo.push(pyGEOM)
  myLifo.push(pyGEOM_SS)

  # Execution
  trackBWR = cle2000.new('TRK_C_SALT',myLifo,1)
  trackBWR.exec()

  # Recover
  myLifo.lib()

  pyUOX_TBH = myLifo.node(namUOX_TBH)
  pyFIG = myLifo.node(namFIG)
  pyTRACK = myLifo.node(namTRACK)
  pyTF_EXC = myLifo.node(namTF_EXC)

  pyUOX_TBH_SS = myLifo.node(namUOX_TBH_SS)
  pyFIG_SS = myLifo.node(namFIG_SS)
  pyTRACK_SS = myLifo.node(namTRACK_SS)
  pyTF_EXC_SS = myLifo.node(namTF_EXC_SS)

  return [pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC,pyUOX_TBH_SS,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS]







