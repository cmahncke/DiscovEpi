#################################################
# Cedric Mahncke
# cedric.mahncke@leibniz-liv.de
# DiscovEpi
# Version 1.1
# Automatically retrieve protein data,
# predict corresponding epitopes and produce
# an epitope map for whole proteomes.
# 08.07.2024: Product.
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
