# PyGan script for burnup evolution simplified pincell problems
# Constant POWER normalization of rates over time stages in EVO: -- > keyword POWER 26.5 [W/gU] 
# Author: R. Guasch
# Date: 2026-08-05
# Hypothesis to test : Difference (D5-S2) in N_Gd155/157 is due to differente energy deposition models between DRAGON5 and SERPENT2
# 

import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib.pyplot as plt

# --- DRAGON5 procedures
from assertS import *
# POST-PROCESSING class
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import multiD5S2_comparisons as multiD5S2
from postproc_cst_pow_evol import D5multiS2_comparisons as D5multiS2 
from postproc_hfactors import plot_microlib_info, plot_self_shielded_microlib_info
#from MULTI_SERP_POSTPROC import *
# --- geometry imports
from GEO_C_SALT import GEO_C_SALT
# --- tracking imports
from TRK_C_SALT import TRK_C_SALT
# --- library imports
from Mix import Mix
# --- burnup evolution imports
from PIN_cpow import cpow_EXTR

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
# run options :
case_name = "UO2_BWR_PINCELL" # "GD_BWR_PINCELL"
void_fraction = 0.0 # 40.0, 80.0

if void_fraction == 0.0:
    N_H, N_O = 4.94546E-02, 2.47298E-02
    void_id = "00"
elif void_fraction == 40.0:
    N_H, N_O = 3.05627063E-02, 1.52813531E-02
    void_id = "40"
elif void_fraction == 80.0:
    N_H, N_O = 1.18604298E-02, 5.93021488E-03
    void_id = "80"
#
# Library = name of the draglib used to build LIB: object
dlib_name = "e81d1mt458" # "e81d1heatr", "e81d1mt18", "e81d1mt458", "e81d1draq"
#"endfb8r1_295" #"e81d1mt18v2" #"e81d1mt458" # "endfb8r1_295" # "eb81dec1" or "eb81dec0" or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

ssh_module = "USS" #"USS", "AUTO"
ssh_method = "PT" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "CORR" # "CORR", "noCORR"
time_intergator = "EXTR"
#
elump_flag = "LUMPE"
glob_opt = "EDP0" #"GLOB" # "NOGL", "EDP0"
#
if case_name == "UO2_BWR_PINCELL":
	burnup_points = "UOX"
elif case_name == "GD_BWR_PINCELL":
	burnup_points = "Gd"
tracked_nuclides = ["U235","U238","Pu239", 
					"Gd157","Gd155", "Xe135","Sm149"]
#burnup_points = "test"

#
# Creation of results directory
path=os.getcwd()

save_dir_D5 = f"{path}/BWR_PINCELL_TESTS_results/{case_name}/D5"
save_dir_S2 = f"{path}/BWR_PINCELL_TESTS_results/{case_name}/S2"
save_dir_comparison = f"{path}/BWR_PINCELL_TESTS_results/{case_name}/Comparison"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
	os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)

# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing

# --- Call to DRAGON5 CLE-2000 procedures :

pyGEOM, pyGEOM_SSH = GEO_C_SALT(case_name)
pyTRACK, pyTF_EXC, pyTRACK_SSH, pyTF_EXC_SSH = TRK_C_SALT(pyGEOM, pyGEOM_SSH)
pyLIB = Mix(dlib_name, case_name, ssh_method, elump_flag, N_H, N_O)

compo_name = f"_CPO_{void_id}_{case_name}_{dlib_name}_{elump_flag}_{ssh_module}_{ssh_method}_{correlation}_{burnup_points}_{glob_opt}_{time_intergator}"
# run DRAGON5 calculation with BU evolution
# namCOMPO, pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo
pyCOMPO = cpow_EXTR(case_name, "COMPO", pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SSH, pyTF_EXC_SSH, StepList, compo_name)
#

D5_NOM_case = D5_case(pyCOMPO = pyCOMPO,
				dlib_name = dlib_name,
				bu_points = burnup_points,
				ssh_opt = "PT",
				correlation = correlation,
				sat = "NODI",
				depl_sol = "RUNG",
				tracked_nuclides = tracked_nuclides,
				BU_lists = BU_lists,
				save_dir = save_dir_D5)
D5_NOM_case.set_BUscheme("predictor", "LE")

evaluation = "endfb8r1"
if dlib_name == "e81d1heatr": 
    edep_mode = 2
else:
	edep_mode = 0
origin_of_S2_data = "pynjoy2012_kerma"
S2_ref_case = S2_case(case_name = case_name,
						lib_name = f"{evaluation}_{origin_of_S2_data}",
						edep_id = edep_mode,
						areQfissSet = False,
						isEcaptSet = False,
						pcc_id = 2,
						specific_power = 26.5,
						tracked_nuclides = tracked_nuclides,
						save_dir = save_dir_S2)

S2_ref_case.read_S2_outputs(f"{case_name}_{void_id}_endfb8r1_edep{edep_mode}_pcc2_mc")

comparisonD5_EDP0_S2 = multiD5S2(f"void fraction {void_id} D5 {ssh_method}-{correlation} {elump_flag}+{glob_opt} - S2  vs edep{edep_mode} pcc 2", [D5_NOM_case], S2_ref_case, burnup_points, tracked_nuclides, f"{save_dir_comparison}/{dlib_name}/{glob_opt}")

comparisonD5_EDP0_S2.compare_keffs()
comparisonD5_EDP0_S2.plot_delta_keff()

comparisonD5_EDP0_S2.compare_Ni()
comparisonD5_EDP0_S2.plot_delta_Ni()