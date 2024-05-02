#---------------------------------------------------------------------------
#
#  Makefile for executing the Trivac non-regression tests
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
	./rtrivac -c $(fcompiler) -q DIFtst.x2m
	./rtrivac -c $(fcompiler) -q Ktests.x2m
	./rtrivac -c $(fcompiler) -q SPNtst.x2m
	./rtrivac -c $(fcompiler) -q iaea3d.x2m
	./rtrivac -c $(fcompiler) -q monju3D_spn.x2m
	./rtrivac -c $(fcompiler) -q multigroup_albedo_2d.x2m
	./rtrivac -c $(fcompiler) -q multigroup_albedo_3d.x2m
	./rtrivac -c $(fcompiler) -q takedaM4_spn.x2m
	./rtrivac -c $(fcompiler) -q NodalTests.x2m
