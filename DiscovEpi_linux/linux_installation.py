#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscovEpi
# Version 1.0
# A tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 20.11.2023: Product.
#################################################

import subprocess
import sys


# This file checks and installs all required packages from pip.
# This function checks and installs packages.
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    return 0


install("re")
install("time")
install("numpy")
install("pandas")
install("seaborn")
install("pathlib")
install("PySide6")
install("datetime")
install("operator")
install("requests")
install("xlsxwriter")
install("webbrowser")
install("collections")
install("matplotlib.pylab")
install("ctypes.wintypes")

print("\n\nAll required packages installed\nmain_gui.py is ready for use!")