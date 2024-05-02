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

def vverGEO(namGEOM):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")

  # Execution 
  geoVVER = cle2000.new('vverGEO',myLifo,1)
  geoVVER.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  return pyGEOM
