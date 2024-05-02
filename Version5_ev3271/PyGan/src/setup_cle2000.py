#
# python3 setup_cle2000.py install --home=.
#
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

def main():
  import os
  from distutils.sysconfig import get_config_var
  incdir = os.path.join(get_python_lib(plat_specific=1), "numpy/core/include")
  mach = os.path.basename(os.getcwd())
  Code = os.environ.get("CODE_EMBEDDED", None) # Code selection
  Compiler = os.environ.get("COMPILER", None) # Compiler selection
  FortranLib = os.environ.get(Compiler, None) # directory with libgfortran.a
  HDF5Lib = os.environ.get("HDF5_API", None) # directory with libhdf5.a
  pylib = os.path.basename(get_config_var("LIBDIR")) # get lib or lib64
  print("install Cle2000 binding to", Code, "on directory",mach, "pylib=",pylib, "Compiler=",Compiler)
  if Compiler == "NVTOOLS":
    libdir="../../lib/"+mach+"_nvidia"
    libUtl="../../../Utilib/lib/"+mach+"_nvidia"
    libTri="../../../Trivac/lib/"+mach+"_nvidia"
    libDra="../../../Dragon/lib/"+mach+"_nvidia"
    libDon="../../../Donjon/lib/"+mach+"_nvidia"
    libNv=os.environ.get("NVTOOLS", None)+"/../lib"
    extralink=["-lnvc","-lnvcpumath","-lnvf"]
  elif Compiler == "INTELTOOLS":
    libdir="../../lib/"+mach+"_intel"
    libUtl="../../../Utilib/lib/"+mach+"_intel"
    libTri="../../../Trivac/lib/"+mach+"_intel"
    libDra="../../../Dragon/lib/"+mach+"_intel"
    libDon="../../../Donjon/lib/"+mach+"_intel"
    libNv=" "
    extralink=[ ]
  else:
    libdir="../../lib/"+mach
    libUtl="../../../Utilib/lib/"+mach
    libTri="../../../Trivac/lib/"+mach
    libDra="../../../Dragon/lib/"+mach
    libDon="../../../Donjon/lib/"+mach
    libNv=" "
    extralink=[ ]
  print("debug Compiler=",Compiler,"libdir=",libdir,"Code=",Code)

  if Code == "GANLIB":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with GANLIB",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     extra_link_args = extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Ganlib","gfortran","hdf5"] ) ])
  elif Code == "TRIVAC":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with TRIVAC",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__trivac__', None)],
                     extra_link_args = extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DRAGON":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DRAGON",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__dragon__', None)],
                     extra_link_args = extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libDra,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DONJON":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DONJON",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__donjon__', None)],
                     extra_link_args = extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libDra,libDon,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Donjon","Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "GANLIB_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with GANLIB_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"]+extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Ganlib","gfortran","hdf5"] ) ])
  elif Code == "TRIVAC_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with TRIVAC_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__trivac__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"]+extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DRAGON_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DRAGON_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__dragon__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"]+extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libDra,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  elif Code == "DONJON_OMP":
    setup (name="Cle2000",
       version="5.0",
       description="Python bindings for Cle-2000 with DONJON_OMP",
       author="Alain Hebert",
       author_email="alain.hebert@polymtl.ca",
       license="LGPL",
       ext_modules=[Extension('cle2000',sources=['cle2000module.c'],
                     define_macros=[('__donjon__', None)],
                     extra_f90_compile_args = ["-fopenmp"],
                     extra_link_args = ["-lgomp"]+extralink,
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=[libdir,FortranLib,HDF5Lib,libUtl,libTri,libDra,libDon,libNv],
                     runtime_library_dirs=[HDF5Lib,libNv],
                     libraries=["Donjon","Dragon","Trivac","Utilib","Ganlib","gfortran","hdf5"] ) ])
  else:
    raise ValueError(Code+" is not implemented for distutils bindings")
if __name__ == "__main__":
  main()
