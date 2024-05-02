#---------------------------------------------------------------------------
#
#  Makefile for executing the Dragon non-regression tests
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
	./rdragon -c $(fcompiler) -q iaea2d.x2m
	./rdragon -c $(fcompiler) -q g2s_prestation.x2m
	./rdragon -c $(fcompiler) -q salmacro.x2m
	./rdragon -c $(fcompiler) -q tmacro.x2m
	./rdragon -c $(fcompiler) -q VanDerGucht.x2m
	./rdragon -c $(fcompiler) -q VanDerGucht-295.x2m
	./rdragon -c $(fcompiler) -q tdraglib.x2m
	./rdragon -c $(fcompiler) -q twimsE.x2m
	./rdragon -c $(fcompiler) -q twlup.x2m
	./rdragon -c $(fcompiler) -q tmatxs2.x2m
	./rdragon -c $(fcompiler) -q trowland_shem295_jeff3.1.x2m
	./rdragon -c $(fcompiler) -q OSC_CASEA_1level_rse.x2m
	./rdragon -c $(fcompiler) -q sens.x2m
	./rdragon -c $(fcompiler) -q testVVER7.x2m
	./rdragon -c $(fcompiler) -q fbr_tone.x2m
	./rdragon -c $(fcompiler) -q testDuo.x2m
	./rdragon -c $(fcompiler) -q testDuo_B1.x2m
	./rdragon -c $(fcompiler) -q uo2_evo_xsm.x2m
	./rdragon -c $(fcompiler) -q C2D20.x2m
	./rdragon -c $(fcompiler) -q CFC-CELL.x2m
	./rdragon -c $(fcompiler) -q ErmBeavrsPwrRefl.x2m
	./rdragon -c $(fcompiler) -q rep900_het_gff_jef2p2.x2m
	./rdragon -c $(fcompiler) -q RegtestCNG_mccg.x2m
	./rdragon -c $(fcompiler) -q RegtestLZC_mccg.x2m
	./rdragon -c $(fcompiler) -q pincell_mco.x2m
	./rdragon -c $(fcompiler) -q pincell_sap.x2m
ifeq ($(hdf5),1)
	./rdragon -c $(fcompiler) -q UOX_5x5_TG6_sym8_multiDom.x2m
	./rdragon -c $(fcompiler) -q uo2_evo_hdf.x2m
	./rdragon -c $(fcompiler) -q pincell_mpo.x2m
	./rdragon -c $(fcompiler) -q ASSBLY_CASEA_1level_mpo.x2m
	./rdragon -c $(fcompiler) -q ASSBLY_VVER_1level_mpo.x2m
endif
