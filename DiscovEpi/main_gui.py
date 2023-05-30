#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscoEpi
# A Tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 05.05.2023: Product.
#################################################


import window


try:
    window.open_win()
except Exception as e:
    print("Error occured", e)
