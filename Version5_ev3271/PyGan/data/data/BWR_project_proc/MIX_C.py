#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, used by R.Guasch                           #
# Date        : 2023                                                #
# Purpose     : MIX definition for pin cells                        #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_C(namMIX,Library,name_geom):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.push(name_geom)
  #myLifo.lib()

  # Execution
  mixBWR = cle2000.new('MIX_C',myLifo,1)
  mixBWR.exec()

  # Recover
  #myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
