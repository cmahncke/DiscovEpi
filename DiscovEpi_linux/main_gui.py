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


import window


try:
    window.open_win()
except Exception as e:
    print("Error occured", e)
