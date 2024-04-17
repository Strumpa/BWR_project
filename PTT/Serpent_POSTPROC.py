# Python3 script for post treating Serpent2 test cases for JEFF311-PyNjoy2016 study
# Author : R. Guasch
# test 1 : removed all metastable isotopes that caused error ReadACEFile function
# test 2 : changed all metastable isotopes to fundamental states
# test 3 : changed all fundamentals to metastables

import matplotlib.pyplot as plt
import serpentTools
from serpentTools.settings import rc
import numpy as np
import os


path_to_ref = "/home/p117902/Serpent2/Linux_x86_64/"
ref_name = "AT10_24UOX_mc"

test_name='24UOX_test_1_mc' # name of the 

path_to_test1 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/"
isotopes_=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd155','Gd157','Xe135','Sm149']

# --- Chemin du répertoire courant
path=os.getcwd()
#print("$$$ ---------- path = ",path)

# --- Création du répertoire de stockage de résultats 
a=os.path.exists('Serpent2_PyNjoy2016_tests')
if a==False:
    os.mkdir('Serpent2_PyNjoy2016_tests')

SAVE_DIR='Serpent2_PyNjoy2016_tests/'

a=os.path.exists(SAVE_DIR)
if a==False:
    os.mkdir(SAVE_DIR)
form = "png"

def load_serpent2_keffs(path_to_serpent2_res, case_name):
    # --- Keff
    res=serpentTools.read(path_to_serpent2_res+case_name+"_res.m")
    serpent_keff=res.resdata["absKeff"]
    np.savetxt(SAVE_DIR+f'serpent_keff_{case_name}.txt',serpent_keff)
    SERPENT_keff=np.loadtxt(SAVE_DIR+f'serpent_keff_{case_name}.txt',dtype=float)
    SERPENT_Keff=np.zeros(lenBU_SERPENT)
    Ls2=np.shape(SERPENT_keff)
    lenBU_SERPENT=Ls2[0]    
    for k in range(lenBU_SERPENT):
        SERPENT_Keff[k]=SERPENT_keff[k][0] # reshaping array to a 1D list
    return SERPENT_Keff

def load_Serpent_BU(path_to_serpent2_res, case_name):   
    # --- BU
    depFile = path_to_serpent2_res+case_name+'_dep.m'
    dep = serpentTools.read(depFile)
    fuel=dep.materials['total']
    serpent_BU=fuel.burnup
    np.savetxt(SAVE_DIR+f'/serpent_BU_{case_name}.txt',serpent_BU)
    SERPENT_BU=np.loadtxt(SAVE_DIR+f'/serpent_BU_{case_name}.txt',dtype=float)
    Ls1=np.shape(SERPENT_BU)
    lenBU_SERPENT=Ls1[0]
    for k in range(lenBU_SERPENT):
        SERPENT_BU[k]=1000*SERPENT_BU[k]
    return SERPENT_BU

def load_Serpent_IsoDens(path_to_serpent2_res, case_name, isotope):
# --- ISOTOPES DENSITIES
    depFile = path_to_serpent2_res+case_name+'_dep.m'
    dep = serpentTools.read(depFile)
    fuel=dep.materials['total']
    serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotope)
    np.savetxt(SAVE_DIR+f'serpent_ISOTOPESDENS_{case_name}_{isotope}.txt',serpent_ISOTOPESDENS)
    SERPENT_ISOTOPESDENS=np.loadtxt(SAVE_DIR+f'serpent_ISOTOPESDENS_{case_name}_{isotope}.txt',dtype=float)
    SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)
    return SERPENT_ISOTOPESDENS

    

#print('$$$ ---------------- SERPENT_BU =',SERPENT_BU)
#print("$$$ ---------------- SERPENT_Keff = ",SERPENT_Keff)    
#print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS)

def plot_comparisons_keff(listBU, listsKeff, case_names, dpi):
    print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff) ')
    fig,ax = plt.subplots(dpi=dpi)
    for i in range(len(case_names)):
        ax.plot(listBU, listsKeff[i], label=f"Keff for {case_names[i]}")
    ax.legend()
    ax.set_xlabel("BU (MWj/t)")
    ax.set_ylabel("Keff Serpent2")
    ax.set_title("Comparison of Keffs for PyNjoy2016 tests")
    fig.savefig()
def plot_comparisons_isotopeDens(listBU, listsIsoDens, case_names, isotope, dpi):
    print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens) ')
    fig,ax = plt.subplots(dpi=dpi)
    for i in range(len(case_names)):
        ax.plot(listBU, listsIsoDens[i], label=f"Keff for {case_names[i]}")
    ax.legend()
    ax.set_xlabel("BU (MWj/t)")
    ax.set_ylabel("Isotope density")
    ax.set_title("Comparison of Keffs for PyNjoy2016 tests")



