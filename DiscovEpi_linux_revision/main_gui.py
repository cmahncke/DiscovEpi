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
import time


try:
    start = time.time()
    window.open_win()
    end = time.time()
    print("Overall time elapsed: ", end - start)
except Exception as e:
    print("Error occured", e)
