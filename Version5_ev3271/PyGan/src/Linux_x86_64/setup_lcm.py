#
# python3 setup_lcm.py install --home=.
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
    libNv=os.environ.get("NVTOOLS", None)+"/../lib"
    extralink=["-lnvc","-lnvcpumath"]
  elif Compiler == "INTELTOOLS":
    libdir="../../lib/"+mach+"_intel"
    libNv=" "
    extralink=[ ]
  else:
    libdir="../../lib/"+mach
    libNv=" "
    extralink=[ ]
  setup(name="Lcm",
          version="5.0",
          description="Python interface for the lcm C library API",
          author="Alain Hebert",
          author_email="alain.hebert@polymtl.ca",
          license="LGPL",
          ext_modules=[Extension("lcm", ["lcmmodule.c"],
                      extra_link_args = extralink,
                      include_dirs=["../../../Ganlib/src",incdir],
                      library_dirs=[libdir,libNv],
                      runtime_library_dirs=[libNv],
                      libraries=["Ganlib"] ) ])

if __name__ == "__main__":
    main()
