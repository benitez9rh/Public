# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:04:08 2022

@author: s2132627
"""

"""
This script uses the vcol of a surface section (A) in order to blind predict another surface section (B) basepoints (based on surface A's spatial continuity) and then krig the gaps.


Inputs:
    vcol of surface section A.
    Upscaled surface section A (x,y) positions.
    Variogram map of the surface section A.

Outputs:
    Blind prediction of surface section B.
"""


import os
import pathlib
from scipy import stats
from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
from scipy.stats import lognorm, norm
from numpy import *
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import math as m
#import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
#import geostatspy.geostats as geostats                            # functions I needed directly to this script
import time
from time import ctime
from sklearn.metrics import r2_score
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
cmap = plt.cm.plasma                    # color map
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging


"""#################################################################################################################################################################################
########################################################################### User's Input ###########################################################################################
#################################################################################################################################################################################"""

''' Save .pngs and .csvs? '''
save = False

''' Name of the surface A DataFrame's value column'''
vcol = 'zr'

''' density of the krigged output and vertical exageration '''
xdens = 1 # Use max density of 0.5. Density of 0.1 causes error: MemoryError: Unable to allocate 64.9 GiB for an array with shape (952300, 9152) and data type float64
ydens = 1
zexag = 0 # vertical exageration (%)

# APERTURE

''' Path to .csv containing Upscaled vcol of surface A '''
#set working directory and filename
wdu = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Residuals\Upscaled').as_posix() #flips the backslashes to plug in open()
inputfilenameu = pathlib.PureWindowsPath(r'GW1Q4_LR_Upscaled1_zr_FractureMap').as_posix()
bs="\\"; wdu=wdu+bs                               # Use this instead in Windows

# Change working directory to save
wdsave = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Residuals\Upscaled').as_posix() #flips the backslashes to plug in open()
bs="\\"; wdsave = wdsave + bs 

''' Variogram model of surface A '''
variogram_model = "spherical"

ang=float(45)
Mrange = 18
mrange = 15
sill = 1
nugget = 0
dic={'sill': sill, 'range': Mrange, 'nugget': nugget}
"""#################################################################################################################################################################################
#################################################################### Functions & System Variables ##################################################################################
#################################################################################################################################################################################"""
ratio = mrange/Mrange
param = [dic["sill"], dic["range"], dic["nugget"]] #sill, range, nugget
cmap = plt.cm.plasma                    # color map

# set the colormap and centre the colorbar. Thanks to http://chris35wills.github.io/matplotlib_diverging_colorbar/
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f' days: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
def isanydigit(n: str) -> bool: # Tests if is a digit because isdigit() doesn't recognise floats
    try:
        float(n)
        n.isdigit()
        return True
    except ValueError:
        return False


"""#################################################################################################################################################################################
############################################################################# Computing ############################################################################################
#################################################################################################################################################################################"""
stop0 = time.time()
print("Importing surfaces...")
# =====================================================================================================================================================
# Upscaled DataFrame
dfu = pd.read_csv(wdu + inputfilenameu + extension_csv, index_col = 0); dfu.reset_index(drop=True, inplace=True); # Sometimes the index is all zeros and it causes problemas in other functions like xymap
dfu[f'N{vcol}'], tvz, tnsz = nscore(dfu, f'{vcol}')
try:
    scale = [i for i in ''.join((ch if ch in '0123456789.-e' else ' ') for ch in inputfilenameu[inputfilenameu.find("Upscaled")+len("Upscaled"):inputfilenameu.find("Upscaled")+len("Upscaled")+2]).split() if isanydigit(i) == True][-1] # Extracts the scale. 2 is hardcoded because at most my scales have 2 digits https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python
except:
    pass
# =====================================================================================================================================================

print("Kriging...")
# =====================================================================================================================================================
gridx = np.arange(dfu['x'].min(), dfu['x'].max(), xdens) # Based on original points' coordinates
gridy = np.arange(dfu['y'].max(), dfu['y'].min(), -1*ydens)
data = np.array([  dfu['x'], dfu['y'], dfu[f'N{vcol}']  ]).transpose()

if variogram_model == "custom":
    OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2],
        variogram_model="custom",
        anisotropy_angle = -ang,
        variogram_parameters = param,
        anisotropy_scaling = ratio,
        verbose=True,
        enable_plotting=True,) # scaling is the ratio between the major and minor directions' ranges
elif variogram_model == "spherical":
    OK = OrdinaryKriging(data[:, 0], data[:, 1], data[:, 2],
        variogram_model="spherical",
        anisotropy_angle = -ang,
        variogram_parameters = dic,
        anisotropy_scaling = ratio,
        verbose=True,
        enable_plotting=True,) # scaling is the ratio between the major and minor directions' ranges

# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See OrdinaryKriging.__doc__ for more information.)
    z, ss = OK.execute("grid", gridx, gridy)
#================================================================================================================================================================================ 
print("Kriging complete.")

print("Plotting...")
#================================================================================================================================================================================ 
kt.write_asc_grid(gridx, gridy, z, filename="output.asc") # Writes the kriged grid to an ASCII grid file and plot it.
plt.imshow(z, cmap = cmap)
plt.show()

xx, yy = np.meshgrid(gridx, gridy)
im = plt.contourf(xx, yy, z, cmap = cmap, vmin = z.min(), vmax = z.max(), levels = np.linspace(z.min(), z.max(), 100),)                                               ####################################################
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(z.min(), z.max(), 10))                                                                              ## contour map, colour bar and colourbar label??? ##
cbar.set_label(f'{vcol} (mm)', rotation=90, labelpad=20)                                                                                                              ####################################################
plt.imshow(z, interpolation = None, extent = [dfu['x'].min(), dfu['x'].max(), dfu['y'].min(), dfu['y'].max()], vmin = z.min(), vmax = z.max(), cmap = cmap)
plt.title(f"Normalised Residuals Ordinary Kriging\n" + f"Scale: {scale}, Density: {str(xdens)}x{str(ydens)}y, Mr: {Mrange}, mr: {mrange}, S: {sill}, N: {nugget}"); plt.xlabel('X (mm)'); plt.ylabel('Y (mm)')
plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
plt.tight_layout()
if save == True:
    if 'scale' in globals():
        plt.savefig(f"{wdsave}NOK_{vcol}" + f"_scale{scale}_Dens{str(xdens).replace('.','p')}x{str(ydens).replace('.','p')}{extension_png}", dpi=300, bbox_inches = "tight")
    else:
        plt.savefig(f"{wdsave}NOK_{vcol}" + f"_scale{scale}_Dens{str(xdens).replace('.','p')}x{str(ydens).replace('.','p')}{extension_png}", dpi=300, bbox_inches = "tight")

plt.show()
if save == True:
    zz = z.data.flatten() # gets the data from a np masked array and flattens it from 2D to 1D array
    krigdf = pd.DataFrame(list(zip(xx.flatten(), yy.flatten(), zz.flatten())), columns = ['x', 'y', vcol])
    if 'scale' in globals():
        krigdf.to_csv(f"{wdsave}NOK_{vcol}" + f"_scale{scale}_Dens{str(xdens).replace('.','p')}x{str(ydens).replace('.','p')}{extension_csv}")
    else:
        krigdf.to_csv(f"{wdsave}NOK_{vcol}" + f"_scale{scale}_Dens{str(xdens).replace('.','p')}x{str(ydens).replace('.','p')}{extension_csv}")
plt.show()
#================================================================================================================================================================================ 
print("Plotting complete.")

duration()

