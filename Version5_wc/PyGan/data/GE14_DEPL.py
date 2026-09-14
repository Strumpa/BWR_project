## Depletion calculation for GE14 lattices
# Ref. "BWR Progression Problems" - ORNL report.
# Author : R. Guasch
# Date : 11/08/2026

# Goal : test depletion calculations for different energy deposition modes
# Aim : Validate DRAGON EDP0 implementation vs Serpent2 edepmode 0.


import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib.pyplot as plt
import re

def format_isotope_name(isotope):
    return re.sub(r"([A-Za-z]+)(\d+)", r"\1-\2", isotope)

# --- DRAGON5 procedures
from assertS import *
# POST-PROCESSING class
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import multiD5S2_comparisons as multiD5S2

from MIX_GE14 import Mix 
from TRK_SALT import TRK
from burnupSchemes import burnup2LScheme
from getLists import getLists


exec_case = False # Set to false allows for post treatment of a draglib for a specified case, in case it has already been ran
post_treat_void_fractions = False
post_treat_DOM_VAN_void_fractions = True
# Assembly id / general parameters
assembly_id = "VAN" # "VAN" # "DOM"
void_fraction = 80.0 # 80.0 # 0.0

if void_fraction == 0.0:
    N_H, N_O = 4.926497E-02, 2.463249E-02
    void_id = "00"
elif void_fraction == 20.0:
    N_H, N_O = 3.99138446E-02, 1.99569223E-02
    void_id = "20"
elif void_fraction == 40.0:
    N_H, N_O = 3.05627063E-02, 1.52813531E-02
    void_id = "40"
elif void_fraction == 60.0:
    N_H, N_O = 2.12115680E-02, 1.06057840E-02
    void_id = "60"
elif void_fraction == 80.0:
    N_H, N_O = 1.18604298E-02, 5.93021488E-03
    void_id = "80"

#### LIB parameters :
# Draglib name :
dlib_name = "e81d1mt18v2"
# Temperatures :
TFuel, TCool = 900.0, 600.0
# self-shielding
ssh_method = "PT"
# Anisotropy and transport correction
anisotropy_level = 4
transport_correction = "NONE"
# Include lumped isotopes energy
edpmode_flag = "EDP0" # "LUMP"

### SALT: parameters :
# tdt file name matching for DOM or VAN
tdt_file_name_ssh = f"GE14_{assembly_id}_DIAG_SSH_3x3_IC_TISO_MACRO.dat" 
tdt_file_name_L1 = f"GE14_{assembly_id}_DIAG_FLUXL1_3x3_6x6_IC_TISO_MACRO.dat"
tdt_file_name_L2 = f"GE14_{assembly_id}_DIAG_FLUXL2_MOC_TSPC.dat"

### burnup scheme parameters :
# "BWRPP_BU_autop" : investigate self-shielding for VAN case : more Pu breeding
burnup_points = "BWRPPBU" # "BWRPP_BU_autop" # "BWRPP_BU_autop" #"test" # "BWRPPBU" #"Gd_autop3" # "ATRIUM_10_S2_BU", "CASMO-5"
# specific power from "BWR Progression Problems"
specific_power = 40.0 # W/g 
# reaction rates normalisation option :
rates_normalisation_option = "EDP0" # "QFIS", "GLOB", "NOGL"
# max group to perform SPH correction :
max_sph_group = 19
# U234 U235 U238 Pu239 Pu240 Pu241 Pu242 Am241 Gd155 Gd157 Xe135 Sm149
tracked_nuclides = ["U234", "U235", "U238", #"U236"
                    "Np237",
                    "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", 
                    "Am241", "Am243",
					#"Gd152", "Gd160",
                    "Gd154", "Gd155", "Gd156", "Gd157", "Gd158",
                    #"Mo95", -- > need to un-lump Zr95 so that Zr95 decays into Mo95
                    "Tc99",
                    "Ru101",
                    #"Rh103", -- > need to un-lump Ru103
                    "Ag109",
                    #"Xe131",
                    "Xe133", "Xe135",
                    "Cs133",
                    #"Nd143",
                    "Nd145", "Nd148", 
                    # "Sm150" --> need to un-lump Nd150
                    "Sm147", "Sm149", "Sm151", "Sm152",
                    "Eu153", #"Eu155", - > not in S2 results
                    ]

### Get burnup lists from getLists.py
BU_lists = getLists(burnup_points)
BU_pts = BU_lists[0]
SelfShielding_pts = BU_lists[1]
COMPO_pts = BU_lists[2]
RATES_pts = BU_lists[3]
StepList = lcm.new('LCM','StepList')
StepList['BUList']    = np.array(BU_pts, dtype='f')
StepList['SSHList'] = np.array(SelfShielding_pts, dtype='f')
StepList['CompoList'] = np.array(COMPO_pts, dtype='f')
StepList['RatesList'] = np.array(RATES_pts, dtype='f')
StepList.close() # close without erasing
print(f"Burnup points (MWd/kgU) : {BU_pts}")
print(f"Self-shielding calculation points (MWd/kgU) : {SelfShielding_pts}")
print(f"COMPO points (MWd/kgU) : {COMPO_pts}")
print(f"RATES points (MWd/kgU) : {RATES_pts}")

### output handling 

# Creation of results directory
path=os.getcwd()

save_dir_D5 = f"{path}/GE14_DEPL_results/GE14_{assembly_id}_{void_id}/D5"
save_dir_S2 = f"{path}/GE14_DEPL_results/GE14_{assembly_id}_{void_id}/S2"
save_dir_comparison = f"{path}/GE14_DEPL_results/GE14_{assembly_id}_{void_id}/Comparison"
save_dir_post_treatment = f"{path}/GE14_DEPL_results/GE14_figures_{assembly_id}"
save_dir_case = f"{save_dir_D5}/{dlib_name}/{rates_normalisation_option}/{burnup_points}"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
	os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)
if not os.path.exists(save_dir_case):
    os.makedirs(save_dir_case)
if not os.path.exists(save_dir_post_treatment):
    os.makedirs(save_dir_post_treatment)

compo_name = f"_CPO_{assembly_id}_{void_id}_{dlib_name}_{edpmode_flag}_USS_{ssh_method}_{burnup_points}_{rates_normalisation_option}_EXTR_sph{max_sph_group}"

if exec_case:

    pyLIB = Mix(dlib_name, assembly_id, ssh_method, anisotropy_level, transport_correction, TFuel, TCool, N_H, N_O, edpmode_flag)

    trk_ssh, trkfil_ssh, trk_l1, trkfil_l1, trk_l2, trkfil_l2 = TRK(tdt_file_name_ssh, tdt_file_name_L1, tdt_file_name_L2)

    pyCOMPO = burnup2LScheme(pyLIB,trk_ssh,trkfil_ssh,trk_l1,trkfil_l1,trk_l2,trkfil_l2,StepList,compo_name,assembly_id,specific_power,rates_normalisation_option,max_sph_group)

    # Save the MULTICOMPO to a specific directory
    shutil.copyfile(compo_name, f"{save_dir_case}/{compo_name}")

else:
    path_to_saved_compo = f"{os.environ["PYGAN_RESULTS"]}/GE14_DEPL_results/GE14_{assembly_id}_{void_id}/D5/{dlib_name}/{rates_normalisation_option}/{burnup_points}"
    shutil.copyfile(f"{path_to_saved_compo}/{compo_name}", compo_name)
    print(os.listdir())
    pyCOMPO = lcm.new('LCM_INP', compo_name[1:], impx=0)


if exec_case == True or post_treat_void_fractions == False and post_treat_DOM_VAN_void_fractions == False: 
    D5_NOM_case = D5_case(pyCOMPO = pyCOMPO,
                        dlib_name = dlib_name,
                        bu_points = burnup_points,
                        ssh_opt = "PT",
                        correlation = "noCORR",
                        sat = "NODI",
                        depl_sol = "RUNG",
                        tracked_nuclides = tracked_nuclides,
                        BU_lists = BU_lists,
                        save_dir = save_dir_D5)
    D5_NOM_case.set_BUscheme("predictor", "LE")

    evaluation = "endfb8r1"
    origin_of_S2_data = "pynjoy2012_kerma"
    S2_ref_case = S2_case(case_name = f"GE14_{assembly_id}_{void_id}",
                            lib_name = f"{evaluation}_{origin_of_S2_data}",
                            edep_id = 0,
                            areQfissSet = False,
                            isEcaptSet = False,
                            pcc_id = 2,
                            specific_power = specific_power,
                            tracked_nuclides = tracked_nuclides,
                            save_dir = save_dir_S2)


    S2_ref_case.read_S2_outputs(f"GE14_DEPL/GE14_{assembly_id}_{void_id}_endfb8r1_depl.serp")
    pccid = 2

    comparisonD5_EDP0_S2 = multiD5S2(f"void {void_id} D5 {ssh_method} {edpmode_flag} - S2 : {rates_normalisation_option} vs edep0 pcc {pccid}", [D5_NOM_case], S2_ref_case, burnup_points, tracked_nuclides, save_dir_case)

    comparisonD5_EDP0_S2.compare_keffs()
    comparisonD5_EDP0_S2.plot_delta_keff()

    comparisonD5_EDP0_S2.compare_Ni()
    comparisonD5_EDP0_S2.plot_delta_Ni()

if post_treat_void_fractions:
    ### keff plot with all S2 results for a given assembly id but several void fractions
    # First sublot : keff from S2 for reference
    # second subplot : D5 - S2 delta k in pcm
    figKeff, axKeff = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
    figDeltaKeff, axDeltaKeff = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
    figures_Niso = {}
    figures_delta_Niso = {}
    for iso in tracked_nuclides:
        figures_Niso[iso] = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
        figures_delta_Niso[iso] = plt.subplots(nrows=1, ncols=1, figsize=(10,6))

    for void_fraction_id in ["00", "40", "80"]:
        path_to_saved_compo = f"{os.environ["PYGAN_RESULTS"]}/GE14_DEPL_results/GE14_{assembly_id}_{void_fraction_id}/D5/{dlib_name}/{rates_normalisation_option}/{burnup_points}"
        compo_name = f"_CPO_{assembly_id}_{void_fraction_id}_{dlib_name}_{edpmode_flag}_USS_{ssh_method}_{burnup_points}_{rates_normalisation_option}_EXTR_sph{max_sph_group}"
        shutil.copyfile(f"{path_to_saved_compo}/{compo_name}", compo_name)
        pyCOMPO = lcm.new('LCM_INP', compo_name[1:], impx=1)
        D5_assembly_void_id = D5_case(pyCOMPO = pyCOMPO,
                    dlib_name = dlib_name,
                    bu_points = burnup_points,
                    ssh_opt = "PT",
                    correlation = "noCORR",
                    sat = "NODI",
                    depl_sol = "RUNG",
                    tracked_nuclides = tracked_nuclides,
                    BU_lists = BU_lists,
                    save_dir = save_dir_D5)
        D5_assembly_void_id.set_BUscheme("predictor", "LE")
        D5_assembly_void_id.assembly_id = assembly_id
        D5_assembly_void_id.void_id = void_fraction_id

        evaluation = "endfb8r1"
        origin_of_S2_data = "pynjoy2012_kerma"
        S2_ref_case = S2_case(case_name = f"GE14_{assembly_id}_{void_fraction_id}",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 0,
                                areQfissSet = False,
                                isEcaptSet = False,
                                pcc_id = 2,
                                specific_power = specific_power,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_S2)


        S2_ref_case.read_S2_outputs(f"GE14_DEPL/GE14_{assembly_id}_{void_fraction_id}_endfb8r1_depl.serp") 
        pccid = 2
        comparisonD5_EDP0_S2 = multiD5S2(f"void {void_fraction_id} D5 {ssh_method} {edpmode_flag} - S2 : {rates_normalisation_option} vs edep0 pcc {pccid}", [D5_assembly_void_id], S2_ref_case, burnup_points, tracked_nuclides, save_dir_case)

        comparisonD5_EDP0_S2.compare_keffs()
        comparisonD5_EDP0_S2.compare_Ni()

        axKeff.plot(S2_ref_case.BU, S2_ref_case.keff, label=fr"$k_\infty$ S2 : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")
        case_key = list(comparisonD5_EDP0_S2.delta_keffs.keys())[0]
        axDeltaKeff.plot(comparisonD5_EDP0_S2.D5_BU_common, comparisonD5_EDP0_S2.delta_keffs[case_key], label=fr"$\Delta k_\infty$ : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")

        for iso in tracked_nuclides:
            figures_Niso[iso][1].plot(S2_ref_case.BU, S2_ref_case.Ni[iso], label=fr"N {iso} S2 : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")
            case_key = list(comparisonD5_EDP0_S2.delta_Niso_dict[iso].keys())[0]
            figures_delta_Niso[iso][1].plot(comparisonD5_EDP0_S2.D5_BU_common, comparisonD5_EDP0_S2.delta_Niso_dict[iso][case_key], label=fr"$\Delta$ N {iso} : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")

    axDeltaKeff.axhline(y = 250.0, color = 'r', linestyle = '-')
    axDeltaKeff.axhline(y = -250.0, color = 'r', linestyle = '-')  

    axKeff.set_xlabel("Burnup (MWd/t)")
    axKeff.set_ylabel(r"$k_\infty$ (Serpent2)")

    axDeltaKeff.set_xlabel("Burnup (MWd/t)")
    axDeltaKeff.set_ylabel(r"Difference on $k_\infty$ (pcm)")

    axKeff.grid(True)
    axDeltaKeff.grid(True)

    axKeff.legend()
    axDeltaKeff.legend(loc="upper right")

    figKeff.tight_layout()
    figDeltaKeff.tight_layout()

    figKeff.savefig(f"{save_dir_post_treatment}/keff_{assembly_id}_{burnup_points}_all_void_ids.png")
    figDeltaKeff.savefig(f"{save_dir_post_treatment}/delta_keff_{assembly_id}_{burnup_points}_all_void_ids.png")

    plt.close(figKeff)
    plt.close(figDeltaKeff)

    for iso in tracked_nuclides:
        figures_delta_Niso[iso][1].axhline(y=2.0, color='r', linestyle='-')
        figures_delta_Niso[iso][1].axhline(y=-2.0, color='r', linestyle='-')

        figures_Niso[iso][1].set_xlabel("Burnup (MWd/t)")
        figures_Niso[iso][1].set_ylabel(f"{format_isotope_name(iso)} isotopic density (isotope/b/cm) (Serpent2)")

        figures_delta_Niso[iso][1].set_xlabel("Burnup (MWd/t)")
        figures_delta_Niso[iso][1].set_ylabel(f"Relative difference on {format_isotope_name(iso)} isotopic density (%)")

        figures_Niso[iso][1].grid(True)
        figures_delta_Niso[iso][1].grid(True)

        figures_Niso[iso][1].legend()
        figures_delta_Niso[iso][1].legend()

        figures_Niso[iso][0].tight_layout()
        figures_delta_Niso[iso][0].tight_layout()

        figures_Niso[iso][0].savefig(f"{save_dir_post_treatment}/Ni_{iso}_{assembly_id}_{burnup_points}_all_void_ids.png")
        figures_delta_Niso[iso][0].savefig(f"{save_dir_post_treatment}/delta_Ni_{iso}_{assembly_id}_{burnup_points}_all_void_ids.png")

        plt.close(figures_Niso[iso][0])
        plt.close(figures_delta_Niso[iso][0])

if post_treat_DOM_VAN_void_fractions:
    ### Plot all DOM, VAN void fractions on the same plots :
    save_dir_post_treatment_all_cases = f"{path}/GE14_DEPL_results/GE14_figures_DOM_VAN"
    if not os.path.exists(save_dir_post_treatment_all_cases):
        os.makedirs(save_dir_post_treatment_all_cases)
    figKeff, axKeff = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
    figDeltaKeff, axDeltaKeff = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
    figures_Niso = {}
    figures_delta_Niso = {}
    for iso in tracked_nuclides:
        figures_Niso[iso] = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
        figures_delta_Niso[iso] = plt.subplots(nrows=1, ncols=1, figsize=(10,6))

    for assembly_id in ["DOM", "VAN"]:
        for void_fraction_id in ["00", "40", "80"]:
            print(f"Post-treat {assembly_id} {void_fraction_id}")
            path_to_saved_compo = f"{os.environ["PYGAN_RESULTS"]}/GE14_DEPL_results/GE14_{assembly_id}_{void_fraction_id}/D5/{dlib_name}/{rates_normalisation_option}/{burnup_points}"
            compo_name = f"_CPO_{assembly_id}_{void_fraction_id}_{dlib_name}_{edpmode_flag}_USS_{ssh_method}_{burnup_points}_{rates_normalisation_option}_EXTR_sph{max_sph_group}"
            shutil.copyfile(f"{path_to_saved_compo}/{compo_name}", compo_name)
            pyCOMPO = lcm.new('LCM_INP', compo_name[1:], impx=0)
            D5_assembly_void_id = D5_case(pyCOMPO = pyCOMPO,
                        dlib_name = dlib_name,
                        bu_points = burnup_points,
                        ssh_opt = "PT",
                        correlation = "noCORR",
                        sat = "NODI",
                        depl_sol = "RUNG",
                        tracked_nuclides = tracked_nuclides,
                        BU_lists = BU_lists,
                        save_dir = save_dir_D5)
            D5_assembly_void_id.set_BUscheme("predictor", "LE")
            D5_assembly_void_id.assembly_id = assembly_id
            D5_assembly_void_id.void_id = void_fraction_id

            evaluation = "endfb8r1"
            origin_of_S2_data = "pynjoy2012_kerma"
            S2_ref_case = S2_case(case_name = f"GE14_{assembly_id}_{void_fraction_id}",
                                    lib_name = f"{evaluation}_{origin_of_S2_data}",
                                    edep_id = 0,
                                    areQfissSet = False,
                                    isEcaptSet = False,
                                    pcc_id = 2,
                                    specific_power = specific_power,
                                    tracked_nuclides = tracked_nuclides,
                                    save_dir = save_dir_S2)


            S2_ref_case.read_S2_outputs(f"GE14_DEPL/GE14_{assembly_id}_{void_fraction_id}_endfb8r1_depl.serp") 
            pccid = 2
            comparisonD5_EDP0_S2 = multiD5S2(f"void {void_fraction_id} D5 {ssh_method} {edpmode_flag} - S2 : {rates_normalisation_option} vs edep0 pcc {pccid}", [D5_assembly_void_id], S2_ref_case, burnup_points, tracked_nuclides, save_dir_case)

            comparisonD5_EDP0_S2.compare_keffs()
            comparisonD5_EDP0_S2.compare_Ni()

            axKeff.plot(S2_ref_case.BU/1000, S2_ref_case.keff, label=fr"$k_\infty$ S2 : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")
            case_key = list(comparisonD5_EDP0_S2.delta_keffs.keys())[0]
            axDeltaKeff.plot(comparisonD5_EDP0_S2.D5_BU_common/1000, comparisonD5_EDP0_S2.delta_keffs[case_key], label=fr"$\Delta k_\infty$ : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")

            for iso in tracked_nuclides:
                figures_Niso[iso][1].plot(S2_ref_case.BU/1000, S2_ref_case.Ni[iso], label=fr"N {format_isotope_name(iso)} S2 : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")
                case_key = list(comparisonD5_EDP0_S2.delta_Niso_dict[iso].keys())[0]
                figures_delta_Niso[iso][1].plot(comparisonD5_EDP0_S2.D5_BU_common/1000, comparisonD5_EDP0_S2.delta_Niso_dict[iso][case_key], label=fr"$\Delta$ N {format_isotope_name(iso)} : {assembly_id} {void_fraction_id} %", marker="x", linestyle="--")

    axDeltaKeff.axhline(y = 250.0, color = 'r', linestyle = '-')
    axDeltaKeff.axhline(y = -250.0, color = 'r', linestyle = '-')  

    axKeff.set_xlabel("Burnup (GWd/t)", fontsize=14)
    axKeff.set_ylabel(r"$k_\infty$ (Serpent2)", fontsize=14)

    axDeltaKeff.set_xlabel("Burnup (GWd/t)", fontsize=14)
    axDeltaKeff.set_ylabel(r"Difference on $k_\infty$ (pcm)", fontsize=14)

    axKeff.grid(True)
    axDeltaKeff.grid(True)

    axKeff.legend(loc="best", fontsize=14)
    axDeltaKeff.legend(loc="best", fontsize=14)

    figKeff.tight_layout()
    figDeltaKeff.tight_layout()

    figKeff.savefig(f"{save_dir_post_treatment_all_cases}/keff_DOM_VAN_{burnup_points}_all_void_ids.png")
    figDeltaKeff.savefig(f"{save_dir_post_treatment_all_cases}/delta_keff_DOM_VAN_{burnup_points}_all_void_ids.png")

    plt.close(figKeff)
    plt.close(figDeltaKeff)

    for iso in tracked_nuclides:
        figures_delta_Niso[iso][1].axhline(y=2.0, color='r', linestyle='-')
        figures_delta_Niso[iso][1].axhline(y=-2.0, color='r', linestyle='-')

        figures_Niso[iso][1].set_xlabel("Burnup (GWd/t)", fontsize=14)
        figures_Niso[iso][1].set_ylabel(f"{format_isotope_name(iso)} isotopic density (isotope/b/cm) (Serpent2)", fontsize=14)

        figures_delta_Niso[iso][1].set_xlabel("Burnup (GWd/t)", fontsize=14)
        figures_delta_Niso[iso][1].set_ylabel(f"Relative difference on {format_isotope_name(iso)} isotopic density (%)", fontsize=14)

        figures_Niso[iso][1].grid(True)
        figures_delta_Niso[iso][1].grid(True)

        figures_Niso[iso][1].legend(loc="best", fontsize=14)
        figures_delta_Niso[iso][1].legend(loc="best", fontsize=14)

        figures_Niso[iso][0].tight_layout()
        figures_delta_Niso[iso][0].tight_layout()

        figures_Niso[iso][0].savefig(f"{save_dir_post_treatment_all_cases}/Ni_{iso}_DOM_VAN_{burnup_points}_all_void_ids.png")
        figures_delta_Niso[iso][0].savefig(f"{save_dir_post_treatment_all_cases}/delta_Ni_{iso}_DOM_VAN_{burnup_points}_all_void_ids.png")

        plt.close(figures_Niso[iso][0])
        plt.close(figures_delta_Niso[iso][0])