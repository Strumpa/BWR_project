#---------------------------------------------------------------------------
#
#  Makefile for executing the Ganlib non-regression tests
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
	./rganlib -c $(fcompiler) -q testgan1.x2m
	./rganlib -c $(fcompiler) -q testgan2.x2m
	./rganlib -c $(fcompiler) -q testgan3.x2m
ifeq ($(hdf5),1)
	./rganlib -c $(fcompiler) -q testgan4.x2m
endif
