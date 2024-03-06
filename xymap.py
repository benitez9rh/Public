# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 18:37:06 2022

@author: s2132627
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
import math
# import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
# import geostatspy.geostats as geostats                            # functions I needed directly to this script
# import geostatspynscore                                             # Copied function and function dependencies from geostatspy GitHub
import time
from time import ctime
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"

# # =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'GW1_Q4Res').as_posix()
# os.chdir(wd)                                   # set the working directory
# =============================================================================

# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap').as_posix()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows
# 
# vcol = 'aperture'
# =============================================================================
#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT AND NORMALISE DATAFRAME ################################################################
###########################################################################################################################################################'''
# df = pd.read_csv(wd + inputfilename + extension_csv, index_col = 0)                  #This new line is because reading the .csv file with a header and an index column will create some issue


# stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
# vmin = stats['min']['aperture']; vmax = stats['max']['aperture']; #nvmin = stats['min']['Naperture']; nvmax = stats['max']['Naperture']

# =============================================================================
# 
# ''' Variogram map parameters '''
# tmin=-999 #trim values below this
# tmax=999 #trim values above this
# dxlag=1; dxlag=int(dxlag)#lag size in x-direction (i.e., bin size in x-direction); As a first pass. use  2
# dylag=dxlag #lag size in y-direction (i.e., bin size in y-direction); As a first pass. use  2
# nxlag = int(round( (xmax - xmin ) / dxlag )) #number of lags in x-direction from a central (0, 0) point (excluding); As a first pass. use  xtotal/dxlag
# nylag = int( round( ( ymax - ymin ) / dylag )) #number of lags in y-direction from a central (0, 0) point (excluding); As a first pass. use  ytotal/dylag
# 
# =============================================================================
def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f' days: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')

def xymap(df, xcol, ycol, tmin, tmax, origxmin, origxmax, origymin, origymax, nxlag, nylag, dxlag, dylag):
    """
    Parameters
    ----------
    df : Pandas DataFrame
        DESCRIPTION.
    xcol : String
        Name of the df column corresponding to the x-coordinate information.
    ycol : String
        Name of the df column corresponding to the y-coordinate information.
    tmin : Float
        Value below or equal which the vcol data in the df is filtered out.
    tmax : Float
        Value above or equal which the vcol data in the df is filtered out.
    origxmin : Float
        Minimum value of x from the original dataset.
    origxmax : Float
        Maximum value of x from the original dataset.
    origymin : Float
        Minimum value of y from the original dataset.
    origymax : Float
        Maximum value of y from the original dataset.    
    nxlag : Integer
        Number of lags in the x-direction.
    nylag : Integer
        Number of lags in the y-direction.
    dxlag : Float
        Value of the lag distance in the x-direction.
    dylag : Float
        Value of the lag distance in the y-direction.

    Returns
    -------
    xymap : Pandas DataFrame
        Creates a pandas DataFrame with the rows size of the input df with two columns ("xmap", "ymap") which contain the number of the x,y cell coordinates of the variogram map in which the point falls into..
    """
    global xymap
    # df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]    # trim values outside tmin and tmax
    nd = len(df)
    xymap = pd.DataFrame(np.zeros((df.shape[0],2)), columns = ['xmap', 'ymap']).astype('int')
    # stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
    xcoord = np.linspace(origxmin, origxmax, nxlag)
    ycoord = np.linspace(origymin, origymax, nylag)
    for i in range(0,nd):
        if i%1000 ==0:                                  # Track the progress
            print('xymap Loop 1/3: ' + str(i) + "/" + str(nd) + ' ' + ctime())
        try:    # Exception handling in case x-coordinate of points is above xmax
            xi = np.where(xcoord <= df.loc[i]['x'])[-1][-1]  # Get the last x index in which point i falls into in the variogram map indexing
            xymap.loc[i]['xmap'] = xi
        except:
            print(f"\n## Exception: ##\ndf.loc[{i}]['x'] = {df.loc[i]['x']}, origxmax = {origxmax}\n")
        try:    # Exception handling in case y-coordinate of points is above ymax
            yi = np.where(ycoord <= df.loc[i]['y'])[-1][-1]  # Get the last y index in which point i falls into in the variogram map indexing
            xymap.loc[i]['ymap'] = yi
        except:
            print(f"\n## Exception: ##\ndf.loc[{i}]['y'] = {df.loc[i]['y']}, origymax = {origymax}\n")
    return xymap


# stop0 = time.time()
# xymap = xymap(df,'x','y', tmin, tmax, nxlag=nxlag, nylag=nylag, dxlag=dxlag, dylag=dylag)
# exec(f'pd.DataFrame(xymap).to_csv("{wd}{inputfilename}_xymap{extension_csv}")')

# =============================================================================
# # Plot
# cmap = plt.cm.plasma                    # color map
# xmin = -(dxlag*nxlag)-(dxlag/2); ymin = -(dylag*nylag)-(dylag/2); xmax = (dxlag*nxlag)+(dxlag/2); ymax = (dylag*nylag)+(dylag/2);
# plt.subplot(111)
# xx, yy = np.meshgrid(np.linspace(xmin, xmax, nxlag), np.linspace(ymin, ymax, nylag))
# im = plt.contourf(xx, yy, xymap, cmap=cmap, vmin=vmin, vmax=vmax, levels=np.linspace(vmin, vmax, 100),)                                                                          ####################################################
# cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(vmin, vmax, 10))                                                                                        #contour map, colour bar and colourbar label???
# cbar.set_label('Variogram Value', rotation=270, labelpad=20)                                                                                                              ####################################################
# plt.imshow(vmap,interpolation = None, extent = [xmin,xmax,ymin,ymax], vmin = vmin, vmax = vmax, cmap = cmap)
# plt.title(f'{inputfilename} Residual-Z Variogram Map'); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')
# 
# plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
# plt.tight_layout()
# plt.savefig(f'{wd}{inputfilename}_Res_VariogramMap{extension_png}', bbox_inches = 'tight')
# plt.show()
# exec(f'pd.DataFrame(vmap).to_csv("{wd}{inputfilename}_Res_VariogramMap{extension_csv}")')
# exec(f'pd.DataFrame(npmap).to_csv("{wd}{inputfilename}_Res_NPMap{extension_csv}")')
# =============================================================================

# duration()