# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:29:36 2022

@author: s2132627
"""

"""
Difference (or aperture) between original surface 1 and Original surface 2

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
# import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
# import geostatspy.geostats as geostats                            # functions I needed directly to this script
import time
from time import ctime
from sklearn.metrics import r2_score
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
extension_txt = ".txt"
cmap = plt.cm.plasma                    # color map
import numpy as np
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
from Euler3DRotation import rotate_np


class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    # set the colormap and centre the colorbar. Thanks to http://chris35wills.github.io/matplotlib_diverging_colorbar/
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
# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean_Q4').as_posix()
# os.chdir(wd)                                   # set the working directory
# =============================================================================
save = False

xcol = 'x'
ycol = 'y'
vcol = 'aperture'

xdens = 0.5 # density of 0.1 causes error: MemoryError: Unable to allocate 64.9 GiB for an array with shape (952300, 9152) and data type float64
ydens = 0.5
zexag = 0 # vertical exageration (%)

''' Path to save directory ''' # Where to save the aperture files between top and bottom surfaces
savewd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\upscale\upscale1\BlindPrediction').as_posix() #flips the backslashes to plug in open()
bs="\\"; savewd=savewd+bs                               # Use this instead in Windows

''' Path to .csv containing vcol of surface 1 '''
#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Aperture').as_posix() #flips the backslashes to plug in open()
inputfilename = pathlib.PureWindowsPath(r'GW1vsGW2_ApertureMap').as_posix()
bs="\\"; wd=wd+bs                               # Use this instead in Windows

''' Use second surface? '''
surf2_switch = 1 #1 is on, 0 is off
""" Auto-centre df2 with df1 so they are above one another """
auto_centre = False

''' Use surfaces offsets (mm)? '''
x1off = 0
y1off = 0
z1off = 0
x2off = 0
y2off = 0
z2off = 0

"""Apply rotations"""  # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/
s1xr = 0 # Surface 1 rotations in the x-direction
s1yr = 0 # Surface 1 rotations in the y-direction
s1zr = 0 # Surface 1 rotations in the z-direction
s2xr = 0 # Surface 2 rotations in the x-direction
s2yr = 0 # Surface 2 rotations in the y-direction
s2zr = 0 # Surface 2 rotations in the z-direction

''' Path to .csv containing vcol of surface 2 '''
#set working directory and filename
wd2 = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\upscale\upscale1').as_posix() #flips the backslashes to plug in open()
inputfilename2 = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_LR_Upscaled1_aperture_FractureMap').as_posix()
bs="\\"; wd2=wd2+bs                               # Use this instead in Windows

''' Calculate Aperture between surfaces? '''
Aperture_switch = 1 #1 is on, 0 is off

#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT SURFACE 1 ################################################################
###########################################################################################################################################################'''
df1 = pd.read_csv(wd + inputfilename , index_col = 0);
df1 = pd.read_csv(wd + inputfilename , index_col = 0,  header = 0); #index_col = False, skiprows = 2, usecols=[0,1,2], names = ['x', 'y', vcol]

# Apply offsets for surface 1
df1['x'] = df1['x'] + x1off
df1['y'] = df1['y'] + y1off
df1[vcol] = (df1[vcol] + z1off) * (1+zexag)

stats=df1.describe().transpose(); xmin = stats['min'][xcol]; xmax = stats['max'][xcol]; ymin = stats['min'][ycol]; ymax = stats['max'][ycol];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]

df1 = df1[[xcol, ycol, vcol]] # Filter DataFrame to only the necessary columns

#================================================================================================================================================================================ 
'''###########################################################################################################################################################
###################################################################### SECOND SURFACE ########################################################################
###########################################################################################################################################################'''



    
# df2 = pd.read_csv(wd2 + inputfilename2 + extension_csv, index_col = 0);    
df2 = pd.read_csv(wd2 + inputfilename2 + extension_csv, index_col = 0,  header = 0); # , usecols=[0,1,2],  , skiprows = , index_col = False,

if auto_centre == True:
    df_xctr = df['x'].mean()                  # Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
    df_yctr = df['y'].mean()
    df2_xctr = df2['x'].mean()            # Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
    df2_yctr = df2['y'].mean()
    ctr_xdiff = (df_xctr - df2_xctr )
    ctr_ydiff = (df_yctr - df2_yctr )  
    
    df2['x'] = df2['x'] + ctr_xdiff
    df2['y'] = df2['y'] + ctr_ydiff




# Apply offsets for surface 1
df2[xcol] = df2[xcol] + x2off
df2[ycol] = df2[ycol] + y2off
df2[vcol] = (df2[vcol] + z2off) * (1+zexag)

df2 = df2[[xcol, ycol,vcol]] # Filter DataFrame to only the necessary columns

stats=df2.describe().transpose(); xmin = stats['min'][xcol]; xmax = stats['max'][xcol]; ymin = stats['min'][ycol]; ymax = stats['max'][ycol];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]

''' Using the same as surface 1'''
# =============================================================================
# dxlag=1; dxlag=int(dxlag)#lag size in x-direction (i.e., bin size in x-direction); As a first pass. use  2
# dylag=dxlag #lag size in y-direction (i.e., bin size in y-direction); As a first pass. use  2
# nxlag = int(round( (xmax - xmin ) / dxlag )) #number of lags in x-direction from a central (0, 0) point (excluding); As a first pass. use  xtotal/dxlag
# nylag = int( round( ( ymax - ymin ) / dylag )) #number of lags in y-direction from a central (0, 0) point (excluding); As a first pass. use  ytotal/dylag
# minnp=1 ; minnp=int(minnp)#minimum number of points in lag bin need to be minnp+1
# isill=1; isill=int(isill) #sill
# # step=1
# 
# ''' Semi-Variograms' parameters '''
# # tmin = -9999.; tmax = 9999.                             # This has to be set before the semi-variograms down below or it will mess up the variogram map
# lag_dist = dxlag; lag_tol = (dxlag/2); nlag = nxlag;            # maximum lag is 700m and tolerance > 1/2 lag distance for smoothing
# =============================================================================
xcol2 = xcol
ycol2 = ycol
vcol2 = vcol
gridx = np.arange(xmin, xmax, xdens)
gridy = np.arange(ymax, ymin, -1*ydens)


if Aperture_switch == 1:
    arr1 = df1.to_numpy(); #arr1 = arr1[:,1:] # gets rid of first column
    arr2 = df2.to_numpy(); # arr2  = arr2[:,1:] # gets rid of first column
    
    """Apply rotations"""  # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/
    # Centralise
    arr1xavg = arr1[:,0].mean();
    arr1yavg = arr1[:,1].mean();
    arr2xavg = arr2[:,0].mean();
    arr2yavg = arr2[:,1].mean();
    print(arr1xavg, arr1yavg, arr2xavg, arr2yavg)
    arr1[:,0] = arr1[:,0] - arr1xavg; arr1[:,1] = arr1[:,1] - arr1yavg;
    arr2[:,0] = arr2[:,0] - arr2xavg; arr2[:,1] = arr2[:,1] - arr2yavg;
    #rotate
    arr1 = rotate_np(arr1, s1xr, s1yr, s1zr)
    arr2 = rotate_np(arr2, s2xr, s2yr, s2zr)
    # De-centralise
    arr1[:,0] = arr1[:,0] + arr1xavg; arr1[:,1] = arr1[:,1] + arr1yavg;
    arr2[:,0] = arr2[:,0] + arr2xavg; arr2[:,1] = arr2[:,1] + arr2yavg;
    
    arr1 = arr1[~np.isnan(arr1).any(axis=1)]   # Removes all rows with NaNs but keeps the shape of the initial array https://stackoverflow.com/questions/11620914/how-do-i-remove-nan-values-from-a-numpy-array
    arr2 = arr2[~np.isnan(arr2).any(axis=1)]

    """ Creating the PolyData """
    df1cloud = pv.PolyData(arr1)
    df2cloud = pv.PolyData(arr2)
    
    """ Create triangulated surfaces """
    Surf_1 = df1cloud.delaunay_2d()
    Surf_2 = df2cloud.delaunay_2d()
    
    
    diff = np.zeros(Surf_2.n_points)
    Surf_2_verticals = np.zeros((Surf_2.n_points,3))
    Surf_1_verticals = np.zeros((Surf_2.n_points,3))
    for i in range(Surf_2.n_points):
        if i%10_000 ==0:                                  # Track the progress
            print('Surf_2.n_points loop: ' + str(i) + "/" + str(Surf_2.n_points) + ' ' + ctime())
        p = Surf_2.points[i]
        intersect, cell = Surf_1.ray_trace(p , [p[0],p[1], 9999], first_point=True)
        if intersect.size != 0:
            if intersect.size == 3:
                intersectz = intersect[2]
                pz = p[2]
            dist = np.sqrt(np.sum((intersectz - pz) ** 2))
            Surf_1_verticals[i] = intersect
            Surf_2_verticals[i] = p
            
        else:
            intersect, cell = Surf_1.ray_trace(p , [p[0],p[1], -9999], first_point=True)
            if intersect.size != 0:
                if intersect.size == 3:
                    intersectz = intersect[2]
                    pz = p[2]
                dist = np.sqrt(np.sum((intersectz - pz) ** 2))
                Surf_1_verticals[i] = intersect
                Surf_2_verticals[i] = p
            elif intersect.size == 0:
                dist = np.nan
                Surf_1_verticals[i] = np.nan
                Surf_2_verticals[i] = np.nan
        diff[i] = dist
    
    x = Surf_2.points[:,0].T; x = x[~np.isnan(diff)]        # Remove indexes where diff contains NaNs, i.e., all non-intersections
    y = Surf_2.points[:,1].T; y = y[~np.isnan(diff)]        # Remove indexes where diff contains NaNs, i.e., all non-intersections
    Surf_1_verticals = Surf_1_verticals[~np.isnan(diff)]    # Remove indexes where diff contains NaNs, i.e., all non-intersections
    Surf_2_verticals = Surf_2_verticals[~np.isnan(diff)]    # Remove indexes where diff contains NaNs, i.e., all non-intersections
    diff = diff[~np.isnan(diff)]                            # Remove indexes containing NaNs, i.e., all non-intersections

    r2 = 1 - np.sum(diff) / np.sum(np.array([(i-Surf_1_verticals[:,2].mean())**2 for i in Surf_1_verticals[:,2]]) ) # r2 = 1 - ( SSres / SStotal):           SumSquares of residuals = sum of (observed - predicted)**2.        SumSquares total = sum of (observed - observed mean)**2.
    
    #Create df of aperture
    ApertureMap = pd.DataFrame(np.vstack((np.array(x).flatten(), np.array(y).flatten(), diff)).transpose(), columns = [xcol, ycol, "aperture"]).drop_duplicates()
    # ApertureMap = ApertureMap.drop(101535) # There is an outlier point
    # Filter out points that are outside the xy boundaries of the upscaled pointset to avoid artefacts at the edges
    # ApertureMap = ApertureMap.loc[(ApertureMap[xcol] >= dfu2xmin) & (ApertureMap[xcol] <= dfu2xmax) & (ApertureMap[ycol] >= dfu2ymin) & (ApertureMap[ycol] <= dfu2ymax)]
    # Save the aperture df in .csv
    if save == True:
        # ApertureMap.to_csv(f"{savewd}{inputfilename}_{vcol}_ORIGvsORIG-r({s2xr}, {s2yr}, {s2zr})-offset({x1off}, {y1off}, {z1off})_ApertureMap{extension_csv}")
        ApertureMap.to_csv(f"{savewd}{inputfilename2} vs BPred {vcol}_ORIGvsORIG-r({s2xr}, {s2yr}, {s2zr})-offset({x1off}, {y1off}, {z1off})_ErrorMap{extension_csv}")
        
    """ Fracture Map """
    if amin(ApertureMap['aperture'])<0:
        cmap = plt.cm.seismic                    # color map
        plt.scatter(ApertureMap['x'], ApertureMap['y'], c=ApertureMap['aperture'], cmap=cmap, norm=MidpointNormalize(midpoint=0,vmin=amin(ApertureMap['aperture']), vmax=amax(ApertureMap['aperture'])), marker = ".", alpha=1) #alpha is transparency
    else:
        cmap = plt.cm.plasma                    # color map
        plt.scatter(ApertureMap['x'], ApertureMap['y'], c=ApertureMap['aperture'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
    plt.annotate(f'r^2 = {r2:.4f}', (float(x.min()-10), float(y.min())-40), annotation_clip = False)
    plt.colorbar(mappable = None, label = 'Error (mm)', orientation="vertical", ticks=np.linspace(amin(ApertureMap['aperture']), amax(ApertureMap['aperture']), 10))                                                                                                                                                                         #contour map, colour bar and colourbar label???
    plt.xlabel(r'$X (mm)$', fontsize=15)
    plt.ylabel(r'$Y (mm)$', fontsize=15)
    plt.title(f'{inputfilename}\n{vcol}_ORIGvsORIG-r:({s2xr}, {s2yr}, {s2zr})_ApertureMap')
    # plt.title(f'{inputfilename2}\nvs Blind Prediction {vcol}_ORIGvsORIG-r:({s2xr}, {s2yr}, {s2zr})_ErrorMap')
    plt.tight_layout()
    if save == True:
        plt.savefig(f'{savewd}{inputfilename}_{vcol}_ORIGvsORIG-r({s2xr}, {s2yr}, {s2zr})-offset({x1off}, {y1off}, {z1off})ApertureMap{extension_png}', dpi=1000, bbox_inches = 'tight')
        # plt.savefig(f'{savewd}{inputfilename2}vs BPred {vcol}_ORIGvsORIG-r({s2xr}, {s2yr}, {s2zr})-offset({x1off}, {y1off}, {z1off})_ErrorMap{extension_png}', bbox_inches = 'tight')
    plt.show()
    
    plt.hist(ApertureMap['aperture'], bins=np.linspace(ApertureMap['aperture'].min(), ApertureMap['aperture'].max(), 1000))
    plt.xlabel(f"Aperture", fontsize=15)
    plt.ylabel(f'Frequency', fontsize=15)
    plt.title(f'Aperture Histogram')       
    plt.annotate( r"$\mu$" + f": {ApertureMap['aperture'].mean():.3f}\n" + r"$\sigma:$" + f" {ApertureMap['aperture'].std():.3f}\n",(0.85,0.8), xycoords='axes fraction')
    if save == True:
        plt.savefig(f'{savewd}{inputfilename}_{vcol}_ApertureHist{extension_png}', dpi=1000, bbox_inches = 'tight')
        # plt.savefig(f'{savewd}{inputfilename2}vs BPred {vcol}_ORIGvsORIG-r({s2xr}, {s2yr}, {s2zr})-offset({x1off}, {y1off}, {z1off})_ErrorMap{extension_png}', bbox_inches = 'tight')
    plt.show()      
# =============================================================================
# # df doesn't have rotations applied
# cmap = plt.cm.plasma                    # color map
# plt.scatter(df['x'], df['y'], c=df['z'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
# plt.scatter(df2['x'], df2['y'], c=df2['z'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
# plt.colorbar(mappable = None, label = 'Aperture (mm)', orientation="vertical", ticks=np.linspace(amin(ApertureMap['aperture']), amax(ApertureMap['aperture']), 10))                                                                                                                                                                         #contour map, colour bar and colourbar label???
# plt.xlabel(r'$X (mm)$', fontsize=15)
# plt.ylabel(r'$Y (mm)$', fontsize=15)
# plt.title(f'{inputfilename} {vcol}_ORIGvsORIG-r:({s2xr}, {s2yr}, {s2zr})_ApertureMap')
# plt.tight_layout()
# plt.show()        
# # arr has rotations applied
# cmap = plt.cm.plasma                    # color map
# plt.scatter(arr1[:,0], arr1[:,1], c=arr1[:,2], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
# plt.scatter(arr2[:,0], arr2[:,1], c=arr2[:,2], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
# plt.colorbar(mappable = None, label = 'Aperture (mm)', orientation="vertical", ticks=np.linspace(amin(ApertureMap['aperture']), amax(ApertureMap['aperture']), 10))                                                                                                                                                                         #contour map, colour bar and colourbar label???
# plt.xlabel(r'$X (mm)$', fontsize=15)
# plt.ylabel(r'$Y (mm)$', fontsize=15)
# plt.title(f'{inputfilename} {vcol}_ORIGvsORIG-r:({s2xr}, {s2yr}, {s2zr})_ApertureMap')
# plt.tight_layout()
# plt.show()
# =============================================================================

# Create Aperture Delaunay Surface
ApertureNParr = np.vstack((np.array(x).flatten(), np.array(y).flatten(), diff)).transpose() # Create Aperture np.array
dfacloud = pv.PolyData(ApertureNParr)
Surf_a = dfacloud.delaunay_2d() # Create Aperture surface




# =============================================================================
# # Rotations in PyVista surfaces (not the data itself)
# Surf_1 = Surf_1.rotate_x(0, inplace=False)
# Surf_1 = Surf_1.rotate_y(0, inplace=False)
# Surf_1 = Surf_1.rotate_z(0, inplace=False)
# Surf_2 = Surf_1.rotate_x(0, inplace=False)
# Surf_2 = Surf_1.rotate_y(0, inplace=False)
# Surf_2 = Surf_1.rotate_z(0, inplace=False)
# =============================================================================

# Visualise
p = pv.Plotter();
p.add_mesh(Surf_1, color="r", label = "Orig Top Surface", opacity=0.5);
p.add_mesh(Surf_2, color="b", label = "Orig Bottom Surface");
p.add_axes(); p.add_legend();
p.add_camera_orientation_widget();
if save == True:
    p.export_html(f'{savewd}{inputfilename}FractureTopBottom_VExag{zexag}.html')  
p.show()

# Visualise Aperture surface
p = pv.Plotter();
p.add_mesh(Surf_a, scalars = Surf_a.points[:,2], show_edges=False, point_size = 50, label = "Aperture Surface", opacity=0.5); # color="r",scalars = Surf_1.points[:,2], show_edges=True
p.add_axes(); p.add_legend();
p.add_camera_orientation_widget();
if save == True:
    p.export_html(f'{savewd}{inputfilename}FractureAperture_VExag{zexag}.html')  
p.show()
