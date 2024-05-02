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


class Serpent2_case:
    def __init__(self, path, file_name, case_name, temperature, library, isotopes, save_dir):
        self.path_to_data = path
        self.file_name = file_name
        self.case_name = case_name
        self.Temp = temperature
        self.library_used = library
        self.isotopes = isotopes
        self.savedir = save_dir

        if "test" in self.file_name:
            self.case_is_ref = False
        else:
            self.case_is_ref = True

        self.load_Serpent_BU() # initialize SERPENT_BU attribute
        self.load_serpent2_keffs() # initialize SERPENT_Keff attribute
        self.load_Serpent_IsoDens() # initialize isotopic densities in evolutions attribute, stored in SERPENT_isotopes_dens dictionnary with key = isotope, values = list of atomic densities


    def load_Serpent_BU(self):   
        """
        Initialize SEREPENT_BU attribute = BU steps used for depletion calculations
        """
        depFile = self.path_to_data+self.file_name+'_dep.m'
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        serpent_BU=fuel.burnup
        np.savetxt(self.savedir+f'/serpent_BU_{self.file_name+"_"+str(self.Temp)}.txt',serpent_BU)
        self.SERPENT_BU=np.loadtxt(self.savedir+f'/serpent_BU_{self.file_name+"_"+str(self.Temp)}.txt',dtype=float)
        Ls1=np.shape(self.SERPENT_BU)
        lenBU_SERPENT=Ls1[0]
        for k in range(lenBU_SERPENT):
            self.SERPENT_BU[k]=1000*self.SERPENT_BU[k]
        return

    def load_serpent2_keffs(self):
        """
        Initialize SERPENT_Keff attribute
        """
        res=serpentTools.read(self.path_to_data+self.file_name+"_res.m")
        serpent_keff=res.resdata["absKeff"]
        np.savetxt(self.savedir+f'serpent_keff_{self.file_name+"_"+str(self.Temp)}.txt',serpent_keff)
        SERPENT_keff=np.loadtxt(self.savedir+f'serpent_keff_{self.file_name+"_"+str(self.Temp)}.txt',dtype=float)
        Ls2=np.shape(SERPENT_keff)
        lenBU_SERPENT=Ls2[0]
        self.SERPENT_Keff=np.zeros(lenBU_SERPENT)    
        for k in range(lenBU_SERPENT):
            self.SERPENT_Keff[k]=SERPENT_keff[k][0] # reshaping array to a 1D list
        return

    def load_Serpent_IsoDens(self):
        """
        Initialize isotopic evolution data
        """
        depFile = self.path_to_data+self.file_name+'_dep.m'
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        self.SERPENT_isotopes_dens = {}
        for isotope in self.isotopes:
            serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotope)
            np.savetxt(self.savedir+f'serpent_ISOTOPESDENS_{self.file_name+"_"+str(self.Temp)}_{isotope}.txt',serpent_ISOTOPESDENS)
            SERPENT_ISOTOPESDENS=np.loadtxt(self.savedir+f'serpent_ISOTOPESDENS_{self.file_name+"_"+str(self.Temp)}_{isotope}.txt',dtype=float)
            SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)
            self.SERPENT_isotopes_dens[isotope] = np.array(SERPENT_ISOTOPESDENS)
        return



class Serpent2_post_treatment:
    def __init__(self, comparision_name, comparision_list, save_dir):
        """
        comparison_name = string to describe the comparison to plot
        comparision_list = list with a sequence of reference case followed by its associated test case, 
        save_dir = directory where to save figures
        """
        self.print_name = comparision_name
        self.ref_cases = [case for case in comparision_list if case.case_is_ref]
        print(f"ref_cases = {self.ref_cases}")
        self.test_cases = [case for case in comparision_list if not case.case_is_ref]           
        print(f"TEST CASES = {self.test_cases}")
        print(f"len(TEST CASES) = {len(self.test_cases)}")
        self.savedir = save_dir
        for ref in self.ref_cases:
            self.listBU = ref.SERPENT_BU
        self.save_name = self.print_name.replace(" ","_")
        print(f"SAVE NAME is = {self.save_name}")
        # Compare test cases to reference case according to their temperature attribute
        # Generalize this to allow for the comparison of a test case with a reference case based on the .Temp attribute

        # Keff plots
        self.plot_comparisons_keff() # plot Keffs for reference case + all test cases given as input
        self.error_keffs_list = []
        for i in range(len(self.ref_cases)):
            self.compute_errors_Keff_to_ref(self.ref_cases[i], self.test_cases[i]) # creation of the error_keffs_list attribute
        self.plot_errorK_to_ref() # plot errors on Keff between reference case and test cases considered


        # Isotopic density plots
        self.plot_comparisons_isotopeDens() # plot isotopic densities for reference case + all test cases given as input
        self.error_iso_dens_for_tests = []
        for i in range(len(self.ref_cases)):
            print(f"{self.test_cases[i]}")
            self.error_iso_dens_for_tests.append(self.compute_error_isoDens(self.ref_cases[i], self.test_cases[i])) # creation of the error_iso_dens_for_tests attribute
        print("ERROR ISO DENS FOR TESTS IS :")
        print(self.error_iso_dens_for_tests)
        print(f"length of error list attribute is {len(self.error_iso_dens_for_tests)}")
        self.plot_error_isoDens()

    def plot_comparisons_keff(self): 
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff) ')
        fig,ax = plt.subplots(dpi=250)
        for ref in self.ref_cases:
            ax.plot(self.listBU, ref.SERPENT_Keff, label=ref.case_name, marker = "D", linewidth=1, markersize=2, linestyle = "--")
        for test in self.test_cases:
            print(test)
            ax.plot(self.listBU, test.SERPENT_Keff, label=test.case_name, marker = "x", linewidth=1, markersize=2, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU (MWj/t)")
        ax.set_ylabel("Keff Serpent2")
        ax.set_title(f"Comparison of Keffs for PyNjoy2016 tests in {self.print_name}", y=1.05)
        fig.savefig(self.savedir+f"Keffs_comparison_{self.save_name}")
        plt.close("all")
        return
    
    def compute_errors_Keff_to_ref(self, reference, test):
        tmp_err =[]
        for i in range(len(test.SERPENT_Keff)):
            print("In compute errors")
            print(test.case_name)
            tmp_err.append((1/reference.SERPENT_Keff[i]-1/test.SERPENT_Keff[i])*1e5)
        self.error_keffs_list.append(tmp_err)
        return

    def plot_errorK_to_ref(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff error) ')
        fig,ax = plt.subplots(dpi=250)
        for i in range(len(self.error_keffs_list)):
            print(f"i = {i}")
            print(f"len(self.error_keffs_list) = {len(self.error_keffs_list)}")
            print(f"len(self.test_cases) = {len(self.test_cases)}")
            ax.plot(self.listBU, self.error_keffs_list[i], label=self.test_cases[i].case_name, marker = "x", linewidth=1, markersize=3, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU (MWj/t)")
        ax.set_ylabel(f"Error on Keff (pcm)")
        ax.set_title(f"{self.print_name} Keff errors PyNjoy2016 tests.", y=1.05)
        fig.savefig(self.savedir+f"Keff_error_{self.save_name}")
        plt.close("all")
        return
    
    
    def plot_comparisons_isotopeDens(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens) ')
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(dpi=250)
            for ref in self.ref_cases:
                ax.plot(self.listBU, ref.SERPENT_isotopes_dens[iso], label=ref.case_name, marker = "D", linewidth=1, markersize=2, linestyle = "--")
            for test in self.test_cases:
                ax.plot(self.listBU, test.SERPENT_isotopes_dens[iso], label=test.case_name, marker = "x", linewidth=1, markersize=2, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU (MWj/t)")
            ax.set_ylabel(f"Isotope density of {iso} (a/b*cm)")
            ax.set_title(f"{iso} density for PyNjoy2016 tests, in {self.print_name}", y=1.05)
            fig.savefig(self.savedir+f"{iso}_comparison_{self.save_name}")
            plt.close("all")
        return

    def compute_error_isoDens(self, reference, test):
        # Temporarily suppress the warning for invalid division
        with np.errstate(divide='ignore', invalid='ignore'):
            error_per_isotope = {}
            for iso in reference.SERPENT_isotopes_dens.keys():
                ref_dens = reference.SERPENT_isotopes_dens[iso]
                test_dens = test.SERPENT_isotopes_dens[iso]
                relative_error = np.where(ref_dens != 0.0, (test_dens - ref_dens) * 100 / ref_dens, 0.0)
                error_per_isotope[iso] = relative_error
        return error_per_isotope
    
    def plot_error_isoDens(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens error) ')
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(dpi=250)
            print(f"number of errors associated to tests = {len(self.error_iso_dens_for_tests)}")
            print(f"errors array is = {self.error_iso_dens_for_tests}")
            for i in range(len(self.error_iso_dens_for_tests)):
                ax.plot(self.listBU, self.error_iso_dens_for_tests[i][iso], label=self.test_cases[i].case_name, marker = "x", linewidth=1, markersize=2, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU (MWj/t)")
            ax.set_ylabel(f"Error on isotopic density (%)")
            
            ax.set_title(f"{self.print_name} errors on isotopic density of {iso}, PyNjoy2016 tests", y=1.05)
            fig.savefig(self.savedir+f"Error_{iso}_dens_{self.save_name}")
            plt.close("all")
        return


#print('$$$ ---------------- SERPENT_BU =',SERPENT_BU)
#print("$$$ ---------------- SERPENT_Keff = ",SERPENT_Keff)    
#print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS)



testsssssss = {"bob": 50, "alice": 50}
print(testsssssss.keys())

cell = "AT10_24UOX_750K"
cell_temp = "AT10_24UOX_900K"
path_to_data = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/"
ref_T750_file_name = "AT10_24UOX_mc"
ref_T900_file_name = 'AT10_24UOX_temp_mc'

test1_file_name='24UOX_test_1_mc' # name of the tests' outputs
test3_file_name='24UOX_test_3_mc'
test3_900_file_name='24UOX_test_3_temp_mc'
path_to_tests = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/"

# testing correct implementation of error calculations

#path_to_tests = path_to_ref
#test1_name = ref_name
isotopes_=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd155','Gd157','Xe135','Sm149']


#postTreat_Tests(cell, path_to_data, path_to_data, isotopes_, [ref_T750_case_name, ref_T900_case_name], [test1_name, test3_name, test3_temp_name], SAVE_DIR)

reference_case_sss_jeff311_750 = Serpent2_case(path_to_data, ref_T750_file_name, "Serpent1 jeff311, @750K", 750, "sss_jeff311u.data", isotopes_, SAVE_DIR)
reference_case_sss_jeff311_900 = Serpent2_case(path_to_data, ref_T900_file_name, "Serpent1 jeff311, @900K", 900, "sss_jeff311u.data", isotopes_, SAVE_DIR)

test1_case_Pynjoy_nometastables_750 = Serpent2_case(path_to_tests, test1_file_name, "PyNjoy2016, no metastables, @750K", 750, "JEFF311PynJoy2016", isotopes_, SAVE_DIR)
test3_case_Pynjoy_metastables_from_sss1_750 = Serpent2_case(path_to_tests, test3_file_name, "PyNjoy2016 + meta stables from sss_jeff311u, @750K", 750, "JEFF311PynJoy2016", isotopes_, SAVE_DIR)
test3_case_Pynjoy_metastables_from_sss1_900 = Serpent2_case(path_to_tests, test3_900_file_name, "PyNjoy2016 + meta stables from sss_jeff311u @900K", 900, "JEFF311PynJoy2016", isotopes_, SAVE_DIR)

#compare_references = Serpent2_post_treatment("24UOX references", {reference_case_sss_jeff311_750:[reference_case_sss_jeff311_900]}, SAVE_DIR)

#compare_tests_at_750 = Serpent2_post_treatment("24UOX tests 1 and 3", {reference_case_sss_jeff311_750:[test1_case_Pynjoy_nometastables_750, test3_case_Pynjoy_metastables_from_sss1_750]}, SAVE_DIR)
#compare_tests_at_900 = Serpent2_post_treatment("24UOX test 3 900K", {reference_case_sss_jeff311_900:[test3_case_Pynjoy_metastables_from_sss1_900]}, SAVE_DIR)
comparison_list = [reference_case_sss_jeff311_750, test1_case_Pynjoy_nometastables_750, reference_case_sss_jeff311_750, test3_case_Pynjoy_metastables_from_sss1_750, reference_case_sss_jeff311_900, test3_case_Pynjoy_metastables_from_sss1_900]

compare_all_tests_to_ref = Serpent2_post_treatment("24UOX all tests", comparison_list, SAVE_DIR)

"""
    
def postTreat_Gd():
    cell = "AT10_45Gd_750K"
    cell_temp = "AT10_45Gd_900K"
    path_to_ref = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/"
    ref_name = "AT10_45Gd_mc"
    ref_temp_name = "AT10_45Gd_temp_mc"

    test1_name='45Gd_test_1_mc' # name of the tests' outputs
    test3_name='45Gd_test_3_mc'
    test3_temp_name = '45Gd_test_3_temp_mc'

    path_to_tests = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/"

    # testing correct implementation of error calculations

    #path_to_tests = path_to_ref
    #test1_name = ref_name
    isotopes_=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd155','Gd157','Xe135','Sm149']

    reference_keffs=load_serpent2_keffs(path_to_ref, ref_name)
    test1_keffs=load_serpent2_keffs(path_to_tests, test1_name)
    test3_keffs=load_serpent2_keffs(path_to_tests, test3_name)

    ref_temp_keffs = load_serpent2_keffs(path_to_ref, ref_temp_name)
    test3_temp_keffs = load_serpent2_keffs(path_to_tests, test3_temp_name)
    

    serpBU=load_Serpent_BU(path_to_ref, ref_name)
    print(serpBU)

    isotopes_data_ref = []
    isotopes_data_test1 = []
    isotopes_data_test3 = []
    isotopes_data_ref_temp = []
    isotopes_data_test3_temp = []
    for iso in isotopes_:
        isotopes_data_ref.append(load_Serpent_IsoDens(path_to_ref, ref_name, iso))
        isotopes_data_ref_temp.append(load_Serpent_IsoDens(path_to_ref, ref_temp_name, iso))
        isotopes_data_test1.append(load_Serpent_IsoDens(path_to_tests, test1_name, iso))
        isotopes_data_test3.append(load_Serpent_IsoDens(path_to_tests, test3_name, iso))
        isotopes_data_test3_temp.append(load_Serpent_IsoDens(path_to_tests, test3_temp_name, iso))

    errors_keff = compute_errors_Keff_to_ref(reference_keffs, [test1_keffs,test3_keffs])
    errors_keff_temp = compute_errors_Keff_to_ref(ref_temp_keffs, [test3_temp_keffs])
    errors_isodens_test1 = compute_error_isoDens(isotopes_data_ref, isotopes_data_test1)
    errors_isodens_test3 = compute_error_isoDens(isotopes_data_ref, isotopes_data_test3)
    errors_isodens_test3_temp = compute_error_isoDens(isotopes_data_ref_temp, isotopes_data_test3_temp)
    print(errors_keff)
    print(len(errors_isodens_test1))
    print(len(isotopes_))
    print(len(isotopes_data_test1))
    plot_comparisons_keff(serpBU,[reference_keffs, test1_keffs, test3_keffs, ref_temp_keffs, test3_temp_keffs], ["Serpent1 lib @750K", "PyNjoy2016, no metastables @750K", "PynJoy2016, metastables from Serpent1 lib @750K", "Serpent1 lib @900K", "PynJoy2016, metastables from Serpent1 lib @900K"], cell, dpi=250, save_dir=SAVE_DIR)
    for iso_num in range(len(isotopes_)):
        plot_comparisons_isotopeDens(serpBU, [isotopes_data_ref[iso_num], isotopes_data_test1[iso_num], isotopes_data_test3[iso_num]], ["Serpent1 lib : sss_jeff311u", "Test 1, modified PyNjoy2016, no metastables", "Test 3 : metastables from sss_jeff311u"], isotopes_[iso_num], cell,dpi=250, save_dir=SAVE_DIR)
        plot_comparisons_isotopeDens(serpBU, [isotopes_data_ref_temp[iso_num], isotopes_data_test3_temp[iso_num]], ["Serpent1 lib, T=900K", "Test 3, T=900K"], isotopes_[iso_num], cell,dpi=250, save_dir=SAVE_DIR)
        plot_error_isoDens(serpBU, [errors_isodens_test1[iso_num], errors_isodens_test3[iso_num]], isotopes_[iso_num], [f"Test 1", "Test 3"], cell, save_dir=SAVE_DIR)
        plot_error_isoDens(serpBU, [errors_isodens_test3_temp[iso_num]], isotopes_[iso_num], ["Test 3, T=900K"], cell_temp, save_dir=SAVE_DIR)
    plot_errorK_to_ref(serpBU, errors_keff, ["Test 1 vs sss_jeff311u", "Test 3 vs sss_jeff311u"], cell, save_dir=SAVE_DIR)
    #plot_errorK_to_ref(serpBU, ["Test 1 vs sss_jeff311u", "Test 3 vs sss_jeff311u", "Test 3 vs sss_jeff311u @900K"], cell, save_dir=SAVE_DIR)
    plot_errorK_to_ref(serpBU, errors_keff_temp, ["Test 3 vs sss_jeff311u, T=900K"], cell_temp, save_dir=SAVE_DIR)
"""