# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 12:11:54 2022

@author: s2132627
"""
import os
import pathlib
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm

cmap = plt.cm.plasma                    # color map

extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
kwords = ["lagx","lagy","lagz"]
"""
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
##################################                Variogram map read and plot code        ##################################
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
"""
# Aperture Map .csv file location (needed for the accurate plot of x,y)
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename_a = pathlib.PureWindowsPath(r'z_ORIGvsORIG-rot(-0.02, 0.31, 0)_ApertureMap').as_posix()

# Variogam Map .csv file location
wdv = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals').as_posix() #flips the backslashes to plug in open()
bs="//"; wdv=wdv+bs                              # Use this instead in linux
inputfilename_v = pathlib.PureWindowsPath(r'z_ORIGvsORIG-rot(-0.02, 0.31, 0)_ApertureMap_ApertureR_VariogramMap').as_posix()
os.chdir(wdv)                                   # set the working directory

df = pd.read_csv(wd + inputfilename_a + extension_csv, index_col = 0, header = 0)
dfv = pd.read_csv(wdv + inputfilename_v + extension_csv, index_col = 0, header = 0)
arr = dfv.to_numpy()
masked_arr = np.ma.masked_equal(arr, -999.0, copy=False) # we ignore the -999 values for the colour scale
# arr = np.delete(arr, 0 , 0) # delete 0 index of axis 0 (rows)
# arr = np.delete(arr, -1 , 0) # delete -1 (i.e. last) index of axis 0(rows)
# arr = np.delete(arr, 0 , 1) # delete 0 index of axis 1 (columns)
# arr = np.delete(arr, -1 , 1) # delete -1 (i.e. last) index of axis 0(columns)


xmin = df['x'].min(); ymin = df['y'].min(); xmax = df['x'].max(); ymax = df['y'].max(); vmin = np.min(masked_arr)  ; vmax = np.max(masked_arr);


# # Extract from the inputfilename the lag at which the variogram map was calculated
# lags = inputfilename_v.split("_"); lags = [l for l in lags if l.rfind("lag") != -1];
# if len(lags) == 1:
#     lag = int(lags[0][lags[0].rfind("lag")+len("lag"):]); 
#     dxlag, dylag = lag, lag
# else:
#     dxlag, dylag, *dzlag = [int(lag[lag.rfind(word) + len(word):]) for lag in lags for word in kwords if word in lag]
#     for lag in ["dxlag","dylag","dzlag"]:
#         try:
#             exec(f"{lag} = {lag}[0]")
#         except:
#             continue
#### For Testing
lag = 1
dxlag, dylag = lag, lag
print(lag)
print(dxlag)
print(dylag)

nxlag = int(round( (xmax - xmin ) / dxlag )) #number of lags in x-direction from a central (0, 0) point (excluding); As a first pass. use  xtotal/dxlag
nylag = int( round( ( ymax - ymin ) / dylag ))
mesh_xmin = -(dxlag*nxlag)-(dxlag/2) ; mesh_ymin = -(dylag*nylag)-(dylag/2); mesh_xmax = (dxlag*nxlag)+(dxlag/2); mesh_ymax = (dylag*nylag)+(dylag/2);

print(dxlag)
print(dylag)
print(nxlag)
print(nylag)
print(mesh_xmin)
print(mesh_ymin)
print(mesh_xmax)
print(mesh_ymax)


Marrow_dir = 157.5
Marrow_dir = 90 - Marrow_dir # Convert Major arrow direction in convetional degrees (clockwise from North) to mathematical degrees (counterclockwise from East) 
marrow_dir = Marrow_dir + 90
Marrow_r = 40
marrow_r = 2

save = False
plot = True
plt.subplot(111)
xx, yy = np.meshgrid(np.arange(mesh_xmin, mesh_xmax, dxlag), np.arange(mesh_ymax, mesh_ymin, -1 * dylag))
im = plt.contourf(xx, yy, arr, cmap = cmap, vmin = vmin, vmax = vmax, levels = np.linspace(vmin, vmax, 100),)                            
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(vmin, vmax, 10))                                                                                        
cbar.set_label('Variogram Value', rotation=90, labelpad=20)            
plt.imshow(arr, interpolation = None, extent = [mesh_xmin, mesh_xmax, mesh_ymin, mesh_ymax], vmin = vmin, vmax = vmax, cmap = cmap)
plt.title(f'GWQ4 Original Aperture lag{lag} Variogram Map'); plt.xlabel('X (mm)'); plt.ylabel('Y (mm)')
plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)

# MajorContinuityVector = plt.arrow(0, 0, np.cos(np.radians(Marrow_dir))*Marrow_r, np.sin(np.radians(Marrow_dir))*Marrow_r , color = "green", edgecolor = "black", head_width =  Marrow_r*0.35, length_includes_head = True)
# MinorContinuityVector = plt.arrow(0, 0, np.cos(np.radians(marrow_dir)) * marrow_r, np.sin(np.radians(marrow_dir)) * marrow_r , color = "red", edgecolor = "black", head_width = marrow_r*0.35, length_includes_head = True)

plt.tight_layout()
if save == True:
    plt.savefig(f'{wdv}GWQ4 aperture lag{lag}_VariogramMap{extension_png}', dpi=600, bbox_inches = "tight")
if plot == True:
    plt.show()
else:
    plt.close()
    
    
"""
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
"""



"""
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
##################################                Semi-Variogram read and plot code         ##################################
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
"""  


# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\upscale19').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'SemiVariogram_Dir67p5').as_posix()
# bs="\\"; wd=wd+bs
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols = [1,2], header = 0)

# direction = 157.5

# plt.plot(df['lag'],df['gamma'],'x',color = 'black',label = 'Azimuth ' + f"{direction}")
# plt.plot([0,2000],[1.0,1.0],color = 'black')
# plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
# plt.ylabel(r'$\gamma \bf(h)$')
# plt.title(f'GW4 Normalised Upscale19 Aperture {direction}deg Semi-Variogram')
# plt.xlim([0, df['lag'].max()+1])
# plt.ylim([0, df['gamma'].max()+1 ])
# plt.legend(loc='upper left')
# plt.grid(True)
# #plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=4.2, wspace=0.2, hspace=0.3)
# if save == True:
#     plt.savefig(f'{wd}{inputfilename}{extension_png}',dpi=300, bbox_inches = "tight")
# if plot == True:
#     plt.show()
# else:
#     plt.close()
    
    
"""
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
"""