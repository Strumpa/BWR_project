#---------------------------------------------------------------------------
#
#  Makefile for executing the Donjon non-regression tests
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
ifeq ($(intel),1)
  fcompiler = intel
else
  ifeq ($(nvidia),1)
    fcompiler = nvidia
  else
    fcompiler = custom
  endif
endif
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rdonjon -c $(fcompiler) -q 1Sap_vs_3Sap.x2m
	./rdonjon -c $(fcompiler) -q AFMtest.x2m
	./rdonjon -c $(fcompiler) -q Candu6.x2m
	./rdonjon -c $(fcompiler) -q candu6_iram.x2m
	./rdonjon -c $(fcompiler) -q Equivalence_SPH_SN_MCO.x2m
	./rdonjon -c $(fcompiler) -q Example1.x2m
	./rdonjon -c $(fcompiler) -q Example2.x2m
	./rdonjon -c $(fcompiler) -q Fessenheim.x2m
	./rdonjon -c $(fcompiler) -q FourChannels.x2m
	./rdonjon -c $(fcompiler) -q Hex_Model3D.x2m
	./rdonjon -c $(fcompiler) -q NCRtest.x2m
	./rdonjon -c $(fcompiler) -q SCR_RES_Interp.x2m
	./rdonjon -c $(fcompiler) -q SCWR64Core.x2m
	./rdonjon -c $(fcompiler) -q TC_EVO.x2m
	./rdonjon -c $(fcompiler) -q TC_XENON.x2m
	./rdonjon -c $(fcompiler) -q TINST.x2m
	./rdonjon -c $(fcompiler) -q Test4x4coupl_nompi.x2m
	./rdonjon -c $(fcompiler) -q Test_microdepl_NCR.x2m
	./rdonjon -c $(fcompiler) -q Test_IDET.x2m
	./rdonjon -c $(fcompiler) -q Test_PMAXS_Macro.x2m
	./rdonjon -c $(fcompiler) -q Test_PMAXS_Map.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH_void_diff_picard.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH_void.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH_void_bfgs.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH_void_lbfgs.x2m
	./rdonjon -c $(fcompiler) -q Test_SPH_void_newton.x2m
	./rdonjon -c $(fcompiler) -q Triga_Model3D.x2m
	./rdonjon -c $(fcompiler) -q VHTR_SPH.x2m
	./rdonjon -c $(fcompiler) -q VHTRc.x2m
	./rdonjon -c $(fcompiler) -q channel_mphy.x2m
	./rdonjon -c $(fcompiler) -q d2p_mco.x2m
	./rdonjon -c $(fcompiler) -q d2p_rfl.x2m
	./rdonjon -c $(fcompiler) -q d2p_sap.x2m
	./rdonjon -c $(fcompiler) -q fbr_core.x2m
	./rdonjon -c $(fcompiler) -q history_mphy.x2m
	./rdonjon -c $(fcompiler) -q pincell_mphy_thm.x2m
	./rdonjon -c $(fcompiler) -q pmox_sap.x2m
	./rdonjon -c $(fcompiler) -q pulseTHM_0d.x2m
	./rdonjon -c $(fcompiler) -q rep1300_reflector.x2m
	./rdonjon -c $(fcompiler) -q rep1300_reflector_sp1.x2m
	./rdonjon -c $(fcompiler) -q rep900.x2m
	./rdonjon -c $(fcompiler) -q rep900_anm.x2m
	./rdonjon -c $(fcompiler) -q rep900cluster.x2m
	./rdonjon -c $(fcompiler) -q rep900EnrichCOMPO.x2m
	./rdonjon -c $(fcompiler) -q rep900EnrichCOMPOhet.x2m
	./rdonjon -c $(fcompiler) -q rep900_iram.x2m
	./rdonjon -c $(fcompiler) -q rep900_mphy.x2m
	./rdonjon -c $(fcompiler) -q rep900_mphy_alt.x2m
	./rdonjon -c $(fcompiler) -q rep900_msap.x2m
	./rdonjon -c $(fcompiler) -q rep900_sim_defrod.x2m
	./rdonjon -c $(fcompiler) -q rep900_sim_recopy.x2m
	./rdonjon -c $(fcompiler) -q rep900_sim.x2m
	./rdonjon -c $(fcompiler) -q ssrcell_thm.x2m
	./rdonjon -c $(fcompiler) -q testNAPhet.x2m
	./rdonjon -c $(fcompiler) -q test_channBWR.x2m
	./rdonjon -c $(fcompiler) -q thm_salt.x2m
	./rdonjon -c $(fcompiler) -q vver_test.x2m
ifeq ($(hdf5),1)
	./rdonjon -c $(fcompiler) -q AFA_180_310_type1_dual.x2m
	./rdonjon -c $(fcompiler) -q Equivalence_SPH_RT_APEX.x2m
	./rdonjon -c $(fcompiler) -q Equivalence_SPH_SN_APEX.x2m
	./rdonjon -c $(fcompiler) -q Equivalence_SPH_RT_MPO.x2m
	./rdonjon -c $(fcompiler) -q SmallCore_rt_mpo.x2m
	./rdonjon -c $(fcompiler) -q SmallCore_anm_mpo_mpo.x2m
	./rdonjon -c $(fcompiler) -q UOX_5x5_TG6_sym8_multiDom_rt2.x2m
	./rdonjon -c $(fcompiler) -q mpo_mco_refl.x2m
	./rdonjon -c $(fcompiler) -q pincell_mpo_concat.x2m
	./rdonjon -c $(fcompiler) -q Reflector_beavrs_DF-NEM.x2m
	./rdonjon -c $(fcompiler) -q Reflector_beavrs_DF-ANM.x2m
	./rdonjon -c $(fcompiler) -q Reflector_beavrs_ERM-NEM.x2m
	./rdonjon -c $(fcompiler) -q Reflector_beavrs_ERM-ANM.x2m
	./rdonjon -c $(fcompiler) -q Reflector_beavrs_Koebke.x2m
	./rdonjon -c $(fcompiler) -q SmallCore_BaffRefl_Benchmark_mpo.x2m
endif
