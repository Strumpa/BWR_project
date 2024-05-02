
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

ifeq ($(pip),1)
  MAKE_PY = $(MAKE) -f Makefile_pip
  options = -q -p
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
	./rpython $(options) test_lcm.py
	./rpython $(options) test_lifo.py
	./rpython $(options) fact.py
	./rpython $(options) simplePOW.py
	./rpython $(options) Exemple_jdd.py
	./rpython $(options) testVVER.py
ifeq ($(hdf5),1)
	./rpython $(options) Equivalence_SPH_APEX.py
	./rpython $(options) pincell_mpo_concat.py
endif
