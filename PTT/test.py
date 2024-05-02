import numpy as np
import os
import sys

# Get the directory where NumPy is imported from
numpy_dir = os.path.dirname(np.__file__)

# Get the directory where the numpy module is located
numpy_module_dir = os.path.dirname(os.path.dirname(numpy_dir))

# Get the full path to where NumPy is installed
numpy_install_dir = os.path.dirname(numpy_module_dir)

print(numpy_install_dir)
