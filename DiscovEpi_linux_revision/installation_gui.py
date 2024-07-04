#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscovEpi
# Version 1.0
# A tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 17.10.2023: Product.
#################################################


import subprocess
import sys
import importlib.metadata as ilm

# Third-party packages
third_party_packages = [
    "numpy",
    "pandas",
    "seaborn",
    "PySide6",
    "requests",
    "xlsxwriter",
    "matplotlib"
]

# Standard_libraries
standard_library_modules = [
    "time",
    "pathlib",
    "datetime",
    "operator",
    "webbrowser",
    "collections",
    "re",
]

# Function – Third-party checking
def install_third_party_package(package):
    try:
        dist = ilm.version(package)
        print(f"{package} is already installed (version {dist})")
    except ilm.PackageNotFoundError:
        print(f"{package} is not installed, installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Function – Standard libraries checking
def check_standard_library_module(module):
    try:
        __import__(module)
        print(f"{module} is a standard library module and is already available")
    except ImportError:
        print(f"{module} is not a standard library module or is missing")

# Install third-party packages
for package in third_party_packages:
    install_third_party_package(package)

# Check standard library modules
for module in standard_library_modules:
    check_standard_library_module(module)

print("\n\nAll required packages are installed or available\nmain_gui.py is ready for use!")
# End of installation_gui.py