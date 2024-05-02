##################################################################################
#                                                                                #
# Description : PyGan scritp for BWR simulation with DRAGON5                     #
# Author      : R. Guasch, adapted from L. Fede                                  #
# Date        : 2024                                                             #
# Purpose     : Perform lattice calculations for BWR cells and assemblies        #
#                                                                                #
##################################################################################
#
########################################################################################################################################################################################
#                                #
#             IMPORTS            #
#                                #
##################################
import sys
#sys.path.append('/usr/local/lib/python3.8/dist-packages')
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib
#matplotlib.use('Agg')

#from matplotlib import pyplot
#import matplotlib.pyplot as plt
from assertS import *
# POST-PROCESSING class
from POSTPROC import *

# --- CELL
# GEOMETRY class
from GEO_C_SALT import *
from GEO_C_NXT import *
# MIXTURES class
from MIX_C import *
# TRACKING class
from TRK_C_SALT import *
from TRK_C_SYBNXT import *
# BU DEPLETION class
from BU_C import *
from BU_C_SYBNXT import *

# --- OTHERS
from getLists import *



########################################################################################################################################################################################
#                                                        #
#   PARAMETERS OF FLUX CALCULATION : TO BE MODIFIED      #
#                                                        #
##########################################################
#
# Library = name of the used Library
Library = 'J311_295'
# 
#
# TYPE = type of calculus
#           - 'ALAMOS'
#           - 'HYBRIDE'
#           - 'NATIVE'
TYPE = 'HYBRIDE'
#
# name_geom = name_geom of the simulated case 
# - cells :
#           - 'AT10_24UOX'
#           - 'AT10_32UOX'
#           - 'AT10_42UOX'
#           - 'AT10_45UOX'
#           - 'AT10_48UOX'
#           - 'AT10_50UOX'
#           - 'AT10_45Gd'
#           - 'AT10_42Gd'
# - assemblies : 
#           - 'ATRIUM-10XM_NL24-2'
#           - 'ATRIUM-10XM_NL24-2_ctrl'
name_geom = 'AT10_UOX' # 'AT10_UOX_Gd2O3'
name_mix = "AT10_24UOX"
tracking_module = "SYBNXT" #"SALT"
#
# Multicompo = 1 if you want to generate a MULTICOMPO object 
#            = 0 either
Multicompo = 1
#
#
# burnup_points = string to use predefine burnup points matching with SERPENT2 results
#        - 'test' : short depletion to test the code
#        - 'UOx'  : used for UOx fuel without Gd poison
#        - 'Gd'   : used for UOx fuel with Gd poison
#        - 'free' : modify the burnup points as you wish
burnup_points = 'UOx4_autop5'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
suffixe = tracking_module+"_"+burnup_points
#
# AUTOP: self-shielding parameter
#        - 'ALL'
#        - 'REGI'
#        - 'NONE'
AUTOP = 'REGI' ;
#
# SALT: tracking parameters for REL2005 calculus
an2d = 18 ;
#
# form = format of the figures generated with vverCOMP.py
#        - 'eps' (recommended)
#        - 'png'
form = 'eps'
#
# Nmin = minimum BU step for DELTA calculation and plot
Nmin = 0
#
# Visualisation options for post-processing functions
# 	= 1 if you want to plot and save the results
# 	= 0 either
# NB : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_DRAGON=1
visu_SERPENT=1
visu_COMP=1
visu_DELTA=1
#
########################################################################################################################################################################################
#                                                 #
#         LOCAL PARAMETERS : DO NOT MODIFY        #
#                                                 #
###################################################
#
# case = 'CELL' / 'ASSEMBLY' - string used to determine wich python class must be called
case = 'CELL'

#
# names for exportation
name_compo = './_COMPO_'+name_geom+'_'+name_mix+suffixe
name_BUvector = name_geom+'_'+name_mix+suffixe
name_fig = name_geom+'_'+name_mix+suffixe+'_2L.ps'
name_fig_SS = name_geom+'_'+name_mix+suffixe+'_SS.ps'
name_fig_1L = name_geom+'_'+name_mix+suffixe+'_1L.ps'
namFIG_SS_NXT=name_geom+'_'+name_mix+suffixe+'_SS_NXT.ps'
namFIG_1L_NXT=name_geom+'_'+name_mix+suffixe+'_1L_NXT.ps'

#
# Visualisation parameters stored in VISU_param
VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]
#
# Creation of results directory
path=os.getcwd()

a=os.path.exists('BWRresults_PyGan_'+name_geom)
if a==False:
	os.mkdir('BWRresults_PyGan_'+name_geom)

SAVE_DIR='BWRresults_PyGan_'+name_geom+'/'+name_mix+'_'+suffixe+'_postprocess'
a=os.path.exists(SAVE_DIR)
if a==False:
	os.mkdir(SAVE_DIR)

#
# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
#
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing
#
# Save BU vector
os.chdir(path+'/'+'BWRresults_PyGan_'+name_geom)
np.savetxt(name_BUvector+'_'+name_mix+'_BUvector.txt',ListeCOMPO)
os.chdir(path)

########################################################################################################################################################################################
#                                                 #
#             EXECUTION OF PROCEDURES             #
#                                                 #
###################################################
#
# --------------------------------------
#                 CELL                 
# --------------------------------------
if case == 'CELL' :
    if tracking_module == 'SALT':

            # ---------- GEOMETRY
            [pyGEOM,pyGEOM_SS] = GEO_C_SALT("GEOM","GEOM_SS",name_geom)

            # ---------- MIXTURES
            pyMIX = MIX_C("LIBRARY",Library,name_mix)

            # ---------- TRACKING
            [pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC,pyUOX_TBH_SS,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS] = TRK_C_SALT("UOX_TBH",name_fig,"TRACK","TF_EXC","UOX_TBH_SS",name_fig_SS,"TRACK_SS","TF_EXC_SS",pyGEOM,pyGEOM_SS)

            # ---------- BU DEPLETION
            pyCOMPO = BU_C("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo)

    elif tracking_module == "SYBNXT":

                        # ---------- GEOMETRY
            [pyGEOM,pyGEOM_SS] = GEO_C_NXT("GEOM","GEOM_SS",name_geom)

            # ---------- MIXTURES
            pyMIX = MIX_C("LIBRARY",Library,name_mix)

            # ---------- TRACKING
            [pyTRACK,pyTF_EXC,pyTRACK_SS] = TRK_C_SYBNXT("TRACK","TF_EXC","TRACK_SS",pyGEOM,pyGEOM_SS)

            # ---------- BU DEPLETION
            pyCOMPO = BU_C_SYBNXT("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo)

"""
# --------------------------------------
#               ASSEMBLIES                
# --------------------------------------
elif case == 'ASSEMBLY':

# --------------------------
#                    REL2005
# --------------------------
	if TYPE =='HYBRIDE' :

		# ---------- GEOMETRY 
		[pyGEOM_1L,pyGEOM_SS,pyGEOM_1L_NXT,pyGEOM_SS_NXT] = GEO_A_H("GEOM_1L","GEOM_SS","GEOM_1L_NXT","GEOM_SS_NXT",name_geom)
		
		# ---------- MIXTURES
		pyMIX = MIX_A_H("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK_2L,pyTF_EXC,pyTRACK_SS,pyFIG_2L,pyFIG_SS,pyTRACK_1L,pyFIG_1L,pyFIG_1L_NXT,pyFIG_SS_NXT] = TRK_A_H("TRACK_2L","TF_EXC","TRACK_SS",name_fig,name_fig_SS,"TRACK_1L",name_fig_1L,pyGEOM_SS,pyGEOM_1L,pyGEOM_SS_NXT,pyGEOM_1L_NXT,namFIG_SS_NXT,namFIG_1L_NXT,name_alamos_SSH_ap3,name_alamos_MOC_ap3,name_alamos_SSH_za,name_alamos_MOC_za,an2d)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_H("COMPO",pyMIX,pyTRACK_2L,pyTF_EXC,pyTRACK_1L,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo,AUTOP)


# --------------------------
#                    ALAMOS
# --------------------------
	elif TYPE =='ALAMOS' :

		# ---------- MIXTURES
		pyMIX = MIX_A_A("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,pyFIG,pyFIG_SS] = TRK_A_A("TRACK","TF_EXC","TRACK_SS","TF_EXC_SS",name_fig,name_fig_SS,name_alamos_MOC_ap3,name_alamos_MOC_za,name_alamos_SSH_ap3,name_alamos_SSH_za)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_A("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo,AUTOP)

# --------------------------
#                    NATIVE
# --------------------------
	elif TYPE =='NATIVE' :

		# ---------- GEOMETRY 
		[pyGEOM_1L,pyGEOM_SS] = GEO_A_N("GEOM_1L","GEOM_SS",name_geom)
		
		# ---------- MIXTURES
		pyMIX = MIX_A_N("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK_1L,pyTRACK_SS,pyFIG_1L,pyFIG_SS] = TRK_A_N("TRACK_1L","TRACK_SS",name_fig,name_fig_SS,pyGEOM_1L,pyGEOM_SS,an2d)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_N("COMPO",pyMIX,pyTRACK_1L,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo)
"""
# -----------------------------------
#    POST-PROCESSING OF GLOBAL VALUES
# -----------------------------------

POSTPROC(pyCOMPO,ListeCOMPO,name_geom,name_mix,suffixe,VISU_param,form,Nmin)
