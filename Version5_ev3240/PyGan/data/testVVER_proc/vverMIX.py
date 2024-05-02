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

def vverMIX(namMIX,Library):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.lib()

  # Execution
  mixVVER = cle2000.new('vverMIX',myLifo,1)
  mixVVER.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
