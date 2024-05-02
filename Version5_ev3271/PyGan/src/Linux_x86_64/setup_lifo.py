#
# python3 setup_lifo.py install --home=.
#
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

def main():
  import os
  incdir = os.path.join(get_python_lib(plat_specific=1), "numpy/core/include")
  mach = os.path.basename(os.getcwd())
  Compiler = os.environ.get("COMPILER", None) # Compiler selection
  if Compiler == "NVTOOLS":
    libdir="../../lib/"+mach+"_nvidia"
  elif Compiler == "INTELTOOLS":
    libdir="../../lib/"+mach+"_intel"
  else:
    libdir="../../lib/"+mach
  setup(name="Lifo",
          version="5.0",
          description="Python interface for the lifo C library function",
          author="Alain Hebert",
          author_email="alain.hebert@polymtl.ca",
          license="LGPL",
          ext_modules=[Extension("lifo", ["lifomodule.c"],
                      include_dirs=["../../../Ganlib/src",incdir],
                      library_dirs=[libdir],
                      libraries=["Ganlib"] ) ])

if __name__ == "__main__":
    main()
