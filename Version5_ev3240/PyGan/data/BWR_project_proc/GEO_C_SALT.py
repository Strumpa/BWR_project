#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, adapted by R.Guasch                        #
# Date        : 2024                                                #
# Purpose     : Geometry definition for pin cells                   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEO_C_SALT(namGEOM,namGEOM_SS,name_geom):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.push(name_geom)
  myLifo.lib()

  # Execution 
  geoBWR = cle2000.new('GEO_C_SALT',myLifo,1)
  geoBWR.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  return [pyGEOM,pyGEOM_SS]
