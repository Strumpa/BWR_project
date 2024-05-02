#!/bin/env python3
""" generation of PYTHONPATH environment variable for PyGan using pth file"""

import os, sys
mach=os.environ['MACH']
Compiler = os.environ.get("COMPILER", None) # Compiler selection
if Compiler == "NVTOOLS":
  name=os.getcwd() + "/../PyGan/lib/" + mach+"_nvidia" + "/python"
elif Compiler == "INTELTOOLS":
  name=os.getcwd() + "/../PyGan/lib/" + mach+"_intel" + "/python"
else:
  name=os.getcwd() + "/../PyGan/lib/" + mach + "/python"
if os.path.isfile(name + "/easy-install.pth"):
  pythonpath = ""
  count = 0
  f = open(name + "/easy-install.pth", 'r')
  for line in f:
    count += 1
    location = "{}".format(name + line.strip()[1:])
    pythonpath = pythonpath + location + ":"
  f.close()
else:
  pythonpath = name
print(pythonpath)
