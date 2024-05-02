
#---------------------------------------------------------------------------
#
#  Makefile for executing the Pygan non-regression tests
#  Author : A. Hebert (2020-7-3)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
ifeq ($(intel),1)
  fcompiler = -c intel
else
  ifeq ($(nvidia),1)
    fcompiler = -c nvidia
  else
    fcompiler = -c custom
  endif
endif

ifeq ($(pip),1)
  MAKE_PY = $(MAKE) -f Makefile_pip
  options = -q -e
else
  MAKE_PY = $(MAKE)
  options = -q
endif

all :
	$(MAKE_PY) donjon -C src
ganlib :
	$(MAKE_PY) ganlib -C src
trivac :
	$(MAKE_PY) trivac -C src
dragon :
	$(MAKE_PY) dragon -C src
donjon :
	$(MAKE_PY) donjon -C src
clean :
	$(MAKE_PY) clean -C src
tests :
	./rpython $(options) $(fcompiler) test_lcm.py
	./rpython $(options) $(fcompiler) test_lifo.py
	./rpython $(options) $(fcompiler) fact.py
	./rpython $(options) $(fcompiler) simplePOW.py
	./rpython $(options) $(fcompiler) Exemple_jdd.py
	./rpython $(options) $(fcompiler) testVVER.py
ifeq ($(hdf5),1)
	./rpython $(options) $(fcompiler) Equivalence_SPH_APEX.py
	./rpython $(options) $(fcompiler) pincell_mpo_concat.py
endif
