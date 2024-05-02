#!/bin/env python3
""" generation of pylib variable for use in Makefiles """

from distutils.sysconfig import get_config_var
pylib = get_config_var("LIBDIR")
if pylib.find('lib64') >= 0:
    print('lib64')
elif pylib.find('lib') >= 0:
    print('lib')
else:
    raise RuntimeError("find_pylib: doesn't contains lib substring")
