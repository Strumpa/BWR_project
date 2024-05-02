#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, R. Guasch adapted for BWR_project          #
# Date        : 2024                                                #
# Purpose     : Post-processing for global values (Keff,Ni)         #
#                                                                   #
#####################################################################
#
import os
import shutil # move files or folders
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import serpentTools
from serpentTools.settings import rc


def POSTPROC(pyCOMPO,ListeCOMPO,name_geom,name_mix,suffixe,VISU_param,form,Nmin):

    ########################################################
    #  Paramètres et initialisation et paramètres locaux 
    ########################################################

    # Récupération des paramètres de visualisation
    visu_DRAGON=VISU_param[0]
    visu_SERPENT=VISU_param[1]
    visu_COMP=VISU_param[2]
    visu_DELTA=VISU_param[3]

    # --- Paramètres de figures
    SIZE=(6,4)

    # --- Options de tracé des erreurs
    LEGENDE_COMP=['DRAGON5','Serpent2']
    LEGENDE_ERROR=['DRAGON5-Serpent2']

    # --- Chemin du répertoire courant
    path=os.getcwd()
    #print("$$$ ---------- path = ",path)

    # --- Création du répertoire de stockage de résultats 
    a=os.path.exists('BWRresults_PyGan_'+name_geom)
    if a==False:
        os.mkdir('BWRresults_PyGan_'+name_geom)

    SAVE_DIR='BWRresults_PyGan_'+name_geom+'/'+name_mix+'_'+suffixe+'_postprocess/'

    a=os.path.exists(SAVE_DIR)
    if a==False:
        os.mkdir(SAVE_DIR)

    # --- Numéro de MIX pour récupérer résultats homogénéisés et condensés 
    NMIX=0

    # --- Ordre des isotopoes post-processés
    isotopes_SOUHAITES=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd155','Gd157','Xe135','Sm149']

    # Chemin d'accès aux résultats Serpent2
    burnup_points=suffixe.split("_")[1]
    #SERPENT_path=f'/home/p117902/Serpent2/Linux_x86_64/' # path to Serpent2 results with sss_jeff311 library.
    SERPENT_path=f'/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/' #for results with JEFF311_Pynjoy2016 acelib <--- working on fixing that.
    if burnup_points != "UOx":
        serpent_suffix = burnup_points+"_"
    else:
        serpent_suffix = ""




    ################################################################
    #   CREACTION DES MATRICES DE DONNEES : DRAGON - SERPENT - ERROR
    ################################################################

    # ------------------------------
    #   MATRICE DES RESULTATS DRAGON 
    # ------------------------------
    if visu_DRAGON==1 or visu_COMP==1 or visu_DELTA==1 :

        if 'ASS' in name_geom :
            DIR='EDIBU_HOM'
        else :
            DIR='EDIBU'   

        lenBU_DRAGON=np.shape(ListeCOMPO)[0]
         
        ISOTOPES=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        #lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1
        lenISOT_DRAGON=np.shape(ISOTOPES)[0]

        DRAGON_BU=ListeCOMPO
        DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        DRAGON_Keff=np.zeros(lenBU_DRAGON)

        #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

        for k in range(lenBU_DRAGON):
            DRAGON_Keff[k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            #print("$$$ ---------------- ISOTOPES BU step ",k,"/",lenBU_DRAGON," = ",pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'])    
            for j in range(lenISOT_DRAGON):
                #print("$$$ ---------------- ISOTOPES ",j,"/",lenISOT_DRAGON," = ",pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                DRAGON_ISOTOPESDENS[j][k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]

        print('$$$ ---------------- DRAGON_BU =',DRAGON_BU)
        print('$$$ ---------------- LEN_DRAGON_BU =',len(DRAGON_BU))
        #print("$$$ ---------------- DRAGON_Keff = ",DRAGON_Keff)    
        #print("$$$ ---------------- DRAGON_ISOTOPESDENS = ",DRAGON_ISOTOPESDENS)

        # --------- Liste des isotopes recuperes dans la Multicompo
        isotopes2=[]
        isotopes=[]
        for k in range(lenISOT_DRAGON):
            isotopes2=isotopes2+[pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
        for k in range(lenISOT_DRAGON):
            if isotopes2[k][0]=='U':
                isotopes=isotopes+[isotopes2[k][0:4]]
            else:
                isotopes=isotopes+[isotopes2[k][0:5]]

        indices=np.zeros(len(isotopes_SOUHAITES))
        for n in range(len(isotopes_SOUHAITES)):
            for m in range(len(isotopes)):
                if isotopes_SOUHAITES[n]==isotopes[m]:
                    indices[n]=m

        #print("$$$ ---------------- DRAGON isotopes = ",isotopes)
        #print("$$$ ---------------- DRAGON isotopes souhaites = ",isotopes_SOUHAITES)
        #print("$$$ ---------------- indices correspondances = ",indices)
            
        # Store in DRAGON_ALL
        DRAGON_ALL=[
               DRAGON_BU,
               DRAGON_Keff,
               DRAGON_ISOTOPESDENS[int(indices[0]),:],
               DRAGON_ISOTOPESDENS[int(indices[1]),:],
               DRAGON_ISOTOPESDENS[int(indices[2]),:],
               DRAGON_ISOTOPESDENS[int(indices[3]),:],
               DRAGON_ISOTOPESDENS[int(indices[4]),:],
               DRAGON_ISOTOPESDENS[int(indices[5]),:],
               DRAGON_ISOTOPESDENS[int(indices[6]),:],
               DRAGON_ISOTOPESDENS[int(indices[7]),:],
               DRAGON_ISOTOPESDENS[int(indices[8]),:],
               DRAGON_ISOTOPESDENS[int(indices[9]),:],
               DRAGON_ISOTOPESDENS[int(indices[10]),:],
           ]
        #print("$$$ ---------------- DRAGON_ALL",DRAGON_AL)

    # -------------------------------
    #   MATRICE DES RESULTATS SERPENT 
    # -------------------------------
    if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :

        # --- Keff
        #res=serpentTools.read(SERPENT_path+name_mix+"_mc_"+serpent_suffix+"res.m")
        # testing data from test3 : PyNjoy2016 xs + metastables from sss_jeff311u
        res=serpentTools.read(SERPENT_path+"24UOX_test_3_mc_res.m")
        serpent_keff=res.resdata["absKeff"]
        np.savetxt('serpent_keff.txt',serpent_keff)
        SERPENT_keff=np.loadtxt('serpent_keff.txt',dtype=float)
            
        # --- BU
        depFile = SERPENT_path+"24UOX_test_3_mc_dep.m"
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        serpent_BU=fuel.burnup
        np.savetxt('serpent_BU.txt',serpent_BU)
        SERPENT_BU=np.loadtxt('serpent_BU.txt',dtype=float)
        
        # --- ISOTOPES DENSITIES
        serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes_SOUHAITES)
        np.savetxt('serpent_ISOTOPESDENS.txt',serpent_ISOTOPESDENS)
        SERPENT_ISOTOPESDENS=np.loadtxt('serpent_ISOTOPESDENS.txt',dtype=float)
        SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)

        Ls1=np.shape(SERPENT_BU)
        print('$$$ ---------------- SERPENT_BU shape =',Ls1)

        Ls2=np.shape(SERPENT_keff)
        lenISOT_SERPENT2=Ls2[0]
        lenBU_SERPENT2=Ls2[1]
        print('$$$ ---------------- SERPENT_keff shape =',Ls2)

        Ls=np.shape(SERPENT_ISOTOPESDENS)
        lenISOT_SERPENT=Ls[0]
        lenBU_SERPENT=Ls[1]
        print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
        
        # Modification unite BU pour match avec DRAGON
        SERPENT_Keff=np.zeros(lenBU_SERPENT)    
        for k in range(lenBU_SERPENT):
            SERPENT_BU[k]=1000*SERPENT_BU[k]
            SERPENT_Keff[k]=SERPENT_keff[k][0]

        #print('$$$ ---------------- SERPENT_BU =',SERPENT_BU)
        #print("$$$ ---------------- SERPENT_Keff = ",SERPENT_Keff)    
        #print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS)

        SERPENT_ALL=[
           SERPENT_BU,
           SERPENT_Keff,
           SERPENT_ISOTOPESDENS[0,:],
           SERPENT_ISOTOPESDENS[1,:],
           SERPENT_ISOTOPESDENS[2,:],
           SERPENT_ISOTOPESDENS[3,:],
           SERPENT_ISOTOPESDENS[4,:],
           SERPENT_ISOTOPESDENS[5,:],
           SERPENT_ISOTOPESDENS[6,:],
           SERPENT_ISOTOPESDENS[7,:],
           SERPENT_ISOTOPESDENS[8,:],
           SERPENT_ISOTOPESDENS[9,:],
           SERPENT_ISOTOPESDENS[10,:],
           ]
        #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)


    # -------------------------------
    #   MATRICE DES ERREUR : ERROR 
    # -------------------------------
    if visu_DELTA==1 :

        ERROR=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        #LE=np.shape(ERROR)
        #print('$$$ ------------------------ ERROR shape=',LE)

        for k in range(lenISOT_DRAGON):
            for j in range(Nmin,lenBU_DRAGON):
                #print('$$$ ----------------------- k=',k,'    j=',j)
                #print('$$$ ----------------------- SERPENT_ALL[k][j]=',SERPENT_ALL[k][j])
                if k==0: # Vecteur BU
                    ERROR[k][j-Nmin]=SERPENT_ALL[k][j]
                elif k==1:  # Vecteur Keff --> erreur en pcm
                    ERROR[k][j-Nmin]=1.0E+5*(DRAGON_ALL[k][j]-SERPENT_ALL[k][j])
                else: # Vecteur isotopique --> erreur en %
                    if SERPENT_ALL[k][j]==0:
                        ERROR[k][j-Nmin]=0
                    else:
                        ERROR[k][j-Nmin]=100*(DRAGON_ALL[k][j]-SERPENT_ALL[k][j])/SERPENT_ALL[k][j]
        #print("$$$ ---------------- ERROR",ERROR)


    ################################################################
    #             TRACE ET SAUVEGARDE DES FIGURES
    ################################################################

    if visu_DRAGON==1:
        print('$$$ -------- POSTPROC.py : DRAGON5 figures ')

        for k in range(lenISOT_DRAGON):

            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1)
            plt.xlabel('BU (MWj/t)')
            plt.grid()
            plt.legend(['DRAGON5'])

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=name_mix+'_DRAGON5_Keff'
                fig_name=name_mix+' - Keff'
            else : # Erreur sur isotopes
                plt.ylabel('Concentration atomique (a/barn.cm)')
                save_name=name_mix+'_DRAGON5_'+isotopes_SOUHAITES[k-1]
                fig_name=name_mix+' - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            #plt.savefig(save_name+"_test3"+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+"_test3"+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_SERPENT==1:
        print('$$$ -------- POSTPROC.py : Serpent2 figures ')

        for k in range(lenISOT_SERPENT):

            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(SERPENT_ALL[0],SERPENT_ALL[k+1],'2-',linewidth=1)
            plt.xlabel('BU (MWj/t)')
            plt.grid()
            plt.legend(['Serpent2'])

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=name_mix+'_Serpent2_Keff'
                fig_name=name_mix+' - Keff'
            else : # Erreur sur isotopes
                plt.ylabel('Concentration atomique (a/barn.cm)')
                save_name=name_mix+'_Serpent2_'+isotopes_SOUHAITES[k-1]
                fig_name=name_mix+' - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            #plt.savefig(save_name+"_test3"+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+"_test3"+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_COMP==1:
        print('$$$ -------- POSTPROC.py : Comparison DRAGON5 / Serpent2 figures ')

        for k in range(lenISOT_DRAGON):

            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1)
            plt.plot(SERPENT_ALL[0],SERPENT_ALL[k+1],'2-',linewidth=1)
            plt.xlabel('BU (MWj/t)')
            plt.grid()
            plt.legend(LEGENDE_COMP)

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=name_geom+'_COMP_Keff'
                fig_name=name_geom+' - Keff'
            else : # Erreur sur isotopes
                plt.ylabel('Concentration atomique (a/barn.cm)')
                save_name=name_geom+'_COMP_'+isotopes_SOUHAITES[k-1]
                fig_name=name_geom+' - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            #plt.savefig(save_name+"_test3"+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+"_test3"+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        

    if visu_DELTA==1 :
        print('$$$ -------- POSTPROC.py : ERROR DRAGON5-Serpent2 figures ')

        for k in range(lenISOT_DRAGON-1):

            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(ERROR[0],ERROR[k+1],'2-',linewidth=1)
            plt.xlabel('BU (MWj/t)')
            plt.grid()
            plt.legend(LEGENDE_ERROR)

            if k == 0: # Erreur sur Keff
                plt.plot([0,60000],[300,300],'r-.') # limite +300pcm
                plt.plot([0,60000],[-300,-300],'r-.') # limite -300pcm
                plt.ylabel('\u0394 Keff (pcm)')
                save_name=name_mix+'_ERROR_Keff'
                fig_name=name_mix+' - \u0394 Keff'
            else : # Erreur sur isotopes
                plt.plot([0,60000],[2,2],'r-.') # limite +2%
                plt.plot([0,60000],[-2,-2],'r-.') # limite -2%
                plt.ylabel('Erreur relative (%)')
                save_name=name_mix+'_ERROR_'+isotopes_SOUHAITES[k-1]
                fig_name=name_mix+' - \u0394 '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            #plt.savefig(save_name+"_test3"+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+"_test3"+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        



























