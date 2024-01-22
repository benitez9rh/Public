# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:33:16 2023

@author: s2132627
"""


"""
This script uses the spatial continuity data of surface A to calculate the areal extension to which it's spatial continuity profile has correlation around surface A's data points.
It then creates a grid of points up to those edges (hereafter surface B) and kriges those grid points using surface A's spatial continuity.

    B B B B B B B B
    B a a a a a a B
    B a a a a a a B
    B a a a a a a B
    B B B B B B B B

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
import numpy as np
#import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
# from Euler3DRotation import rotate_np
# from xymap import xymap
from convexhull import convexhull
# import convexhull
import polygon_inclusion
import shapely
from shapely import Polygon, Point
"""#################################################################################################################################################################################
########################################################################### User's Input ###########################################################################################
#################################################################################################################################################################################"""

''' Save .pngs and .csvs? '''
save = True # True
plot = True
''' Name of the surface A DataFrame's value column'''
vcol = 'aperture'
style="points"  # Style of the kriging algorithm

''' density of the krigged output and vertical exageration '''
xdens = 1 #  Density of 0.1 causes error: MemoryError: Unable to allocate 64.9 GiB for an array with shape (952300, 9152) and data type float64
ydens = 1
# ZR
''' Path to .csv containing vcol of surface A '''
#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\upscale\upscale1').as_posix() #flips the backslashes to plug in open()
inputfilenameu = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_LR_Upscaled1_aperture_FractureMap').as_posix()
bs="\\"; wd=wd+bs                               # Use this instead in Windows

# Change working directory to save
wdsave = pathlib.PureWindowsPath(wd+r'\BlindPrediction').as_posix() #flips the backslashes to plug in open()
bs="//"; wdsave=wdsave+bs

''' Variogram model of surface A '''
rangeM = 6 # Major direction range
rangem = 4 # minor direction range
angle = 67.5 # in degrees Counted clockwise from North (up)
variogram_model = "Spherical"
sill = 1
nugget = 0
# =============================================================================
# # APERTURE
# ''' Path to .csv containing vcol of surface A '''
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap').as_posix()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows
# 
# ''' Path to .csv containing Upscaled vcol of surface A '''
# #set working directory and filename
# wdu = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\upscale').as_posix() #flips the backslashes to plug in open()
# inputfilenameu = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_LR_Upscaled22_aperture_FractureMap').as_posix()
# bs="\\"; wdu=wdu+bs                               # Use this instead in Windows
# 
# ''' Path to .csv containing the vcol Variogram map of surface section A '''
# wdv = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
# inputfilenamevmap = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_Naperture_lag1_VariogramMap').as_posix()
# bs="\\"; wdv=wdv+bs # Use this instead in Windows
# 
# # Change working directory to save
# wdsave = pathlib.PureWindowsPath(wd+r'\BlindPrediction').as_posix() #flips the backslashes to plug in open()
# bs="//"; wdsave=wdsave+bs 
# 
# ''' Variogram model of surface A '''
# variogram_model = "spherical"
# dic={'sill': 1, 'range': 100, 'nugget': 0}
# ang=float(67.5)
# ratio = 75/100
# param = [1, 100, 0] #sill, range, nugget
# =============================================================================

"""#################################################################################################################################################################################
#################################################################### Functions & System Variables ##################################################################################
#################################################################################################################################################################################"""
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
    print(f'days: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
def isanydigit(n: str) -> bool: # Tests if is a digit because isdigit() doesn't recognise floats
    try:
        float(n)
        n.isdigit()
        return True
    except ValueError:
        return False
def inside_convex_polygon(point, vertices):
    """
    Check if point is inside convex polygon # https://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates
    """
    previous_side = None
    n_vertices = len(vertices)
    for n in range(n_vertices):
        a, b = vertices[n], vertices[(n+1)%n_vertices]
        affine_segment = v_sub(b, a)
        affine_point = v_sub(point, a)
        current_side = get_side(affine_segment, affine_point)
        if current_side is None:
            return False #outside or over an edge
        elif previous_side is None: #first segment
            previous_side = current_side
        elif previous_side != current_side:
            return False
    return True
def get_side(a, b):
    x = cosine_sign(a, b)
    if x < 0:
        return "LEFT"
    elif x > 0: 
        return "RIGHT"
    else:
        return None
def v_sub(a, b):
    return (a[0]-b[0], a[1]-b[1])
def cosine_sign(a, b):
    return a[0]*b[1]-a[1]*b[0]
def quad(degree):
    full_rotations = degree//360
    degree = degree%360
    quadrant = math.floor(degree/90 % 4 + 1)
    #print(f"{degree} degree angle ({round(np.radians(degree), 3)} radians) after {full_rotations} full rotations falls in Quadrant {quadrant}")
    return quadrant, full_rotations

dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
angOK=float(-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html
angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
ratio = rangem/rangeM
param = [sill, rangeM, nugget] #sill, range, nugget
cmap = plt.cm.plasma                    # color map
#scale = [i for i in ''.join((ch if ch in '0123456789.-e' else ' ') for ch in inputfilenameu[inputfilenameu.find("Upscaled")+len("Upscaled"):inputfilenameu.find("Upscaled")+len("Upscaled")+2]).split() if isanydigit(i) == True][-1] # Extracts the scale. 2 is hardcoded because at most my scales have 2 digits https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python


"""#################################################################################################################################################################################
############################################################################# Computing ############################################################################################
#################################################################################################################################################################################"""
stop0 = time.time()
print("Importing surfaces...")
# =====================================================================================================================================================
# Original DataFrame
df = pd.read_csv(wd + inputfilenameu + extension_csv, usecols=[1,2,3,4,5,6]); # df = df[[xcol, ycol, vcol]] # Filter DataFrame to only the necessary columns
stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol];
xdiff = xmax-xmin; ydiff = ymax-ymin

df_xctr = (df['x'].max() - df['x'].min()) /2# Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
df_yctr = (df['y'].max() - df['y'].min()) /2


vertices_indexes = convexhull(df, xcol = "x", ycol = "y")

dfvs = df[df.index.isin(vertices_indexes)] # Create a df with only the points of vindex (but it changes the order to ascending!)
dfvs = dfvs.reindex(vertices_indexes) # reorder back
vslst = list(map(list, zip( dfvs['x'], dfvs['y'] ) )) # Create a list of lists containing the xy values pairs for the inside_convex_polygon because it doesn't accept dataframes
vslst.append(vslst[0])

# Create 4 points for each convexhull vertice located at the end of each spatial continuity elipse direction vector
expanded = []
for i, vi in enumerate(vertices_indexes):
    point1 = [dfvs.loc[vi]['x'] + rangeM * np.cos(np.radians(angmath)),        dfvs.loc[vi]['y'] + rangeM * np.sin(np.radians(angmath)) ]
    point2 = [dfvs.loc[vi]['x'] + rangem * np.cos(np.radians(angmath + 90)),     dfvs.loc[vi]['y'] + rangem * np.sin(np.radians(angmath + 90)) ]
    point3 = [dfvs.loc[vi]['x'] + rangeM * np.cos(np.radians(angmath + 180)),    dfvs.loc[vi]['y'] + rangeM * np.sin(np.radians(angmath + 180)) ]
    point4 = [dfvs.loc[vi]['x'] + rangem * np.cos(np.radians(angmath + 270)),     dfvs.loc[vi]['y'] + rangem * np.sin(np.radians(angmath + 270)) ]
    expanded.append(point1)
    expanded.append(point2)
    expanded.append(point3)
    expanded.append(point4)
# expanded_round = expanded; expanded_round.append(expanded_round[0])
expandeddf = pd.DataFrame(expanded, columns = ['x','y'])
expanded_indexes = convexhull(expandeddf, xcol = "x", ycol = "y"); 
expandeddf = expandeddf[expandeddf.index.isin(expanded_indexes)] # Create a df with only the points of expanded_indexes (but it changes the order to ascending!)
expandeddf = expandeddf.reindex(expanded_indexes)
expandedlst = expandeddf.values.tolist()
expandedlst.append(expandedlst[0])

shapely_vslst = [(x,y) for x,y in vslst] # Format vslst into the shapely accepted format [(x1,y1), (x2,y2), ..., (xn,yn)]
shapely_expanded = [(x,y) for x,y in expandedlst] # Format expanded into the shapely accepted format [(x1,y1), (x2,y2), ..., (xn,yn)]

polygon_vslst = Polygon(shapely_vslst)
polygon_expanded = Polygon(shapely_expanded)

X,Y = np.meshgrid(np.arange(expandeddf['x'].min(), expandeddf['x'].max(), xdens), np.arange(expandeddf['y'].min(), expandeddf['y'].max(), ydens))
points = list( map( list,zip(X.flatten(),Y.flatten())))
points_shapely = [(x,y) for x,y in points]

pointsdf = pd.DataFrame(points_shapely, columns = ['x','y']) # Add two zeros columns to points and convert to pd df
pointsdf['outVSLSTinEXPANDED'] = [(polygon_expanded.contains(Point(point)) == True) & (Point(point).within(polygon_vslst) == False) for point in points_shapely]
pointsdf = pointsdf[pointsdf["outVSLSTinEXPANDED"] == True]

# Add the first point of the polygon to the end so the polygons close in matplotlib plotting.
expandeddf = pd.concat(  [expandeddf, pd.DataFrame([np.array(expandeddf.iloc[0])], columns=expandeddf.columns ).set_index([pd.Index(  [expandeddf.index.values[0]] )]  ) ]   )
dfvs = pd.concat(  [dfvs, pd.DataFrame([np.array(dfvs.iloc[0])], columns=dfvs.columns ).set_index([pd.Index(  [dfvs.index.values[0]] )]  ) ]   )

if plot == True:
    plt.plot(pointsdf['x'], pointsdf['y'], 'o', color="green", markersize = 1)
    plt.plot(expandeddf['x'], expandeddf['y'],  color="blue")
    plt.plot(dfvs['x'], dfvs['y'],  color="red")
    # plt.plot(dfvs['x'], dfvs['y'], 'r--', lw=1)
    plt.show()

# =====================================================================================================================================================
print("Kriging...")
# =====================================================================================================================================================
data = np.array(
    [
        df['x'],
        df['y'],
        df[vcol]
    ]
    ).transpose()

if variogram_model == "custom":
    OK = OrdinaryKriging(
        data[:, 0],
        data[:, 1],
        data[:, 2],
        variogram_model="custom",
        anisotropy_angle = angOK,
        variogram_parameters = dic,
        anisotropy_scaling = ratio,
        exact_values  = True,
        verbose=True,
        enable_plotting=True,
    ) # scaling is the ratio between the major and minor directions' ranges
elif variogram_model == "Spherical":
    OK = OrdinaryKriging(
        data[:, 0],
        data[:, 1],
        data[:, 2],
        variogram_model=variogram_model.lower(),
        anisotropy_angle = angOK,
        variogram_parameters = dic,
        anisotropy_scaling = ratio,
        exact_values  = True,
        verbose=True,
        enable_plotting=True,
    ) # scaling is the ratio between the major and minor directions' ranges

###############################################################################
# Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
# grid of points, on a masked rectangular grid of points, or with arbitrary points.
# (See OrdinaryKriging.__doc__ for more information.)
# style (str) – Specifies how to treat input kriging points. Specifying ‘grid’ treats xpoints and ypoints as two arrays of x and y coordinates that define a rectangular grid.
# Specifying ‘points’ treats xpoints and ypoints as two arrays that provide coordinate pairs at which to solve the kriging system.
# Specifying ‘masked’ treats xpoints and ypoints as two arrays of x and y coordinates

if style == "grid":
    gridx = np.arange(xmin-xdiff, xmax+xdiff, xdens) # Based on original points' coordinates
    gridy = np.arange(ymin-ydiff, ymax+ydiff, ydens)
    xx, yy = np.meshgrid(gridx, gridy)

zz, ss = OK.execute(style, np.array(pointsdf['x']), np.array(pointsdf['y']), backend = "vectorized")
# z is the value, ss should stand for sigma squared which is the variance. Backend options are vectorized or loop: vectorized faster but memory intensive whereas loop is slower but also less memory-intensive.
###############################################################################

# kt.write_asc_grid(gridx, gridy, z, filename="output.asc") # Writes the kriged grid to an ASCII grid file and plot it.
# plt.imshow(z, cmap = cmap)
# plt.show()
if style=="points":
    plt.scatter(pointsdf['x'], pointsdf['y'], c=zz, cmap=cmap, s=2, marker = ".", alpha=1 )                                              
    plt.scatter(df['x'], df['y'], c=df[vcol], cmap=cmap, s=2, marker = ".", alpha=1 )  
    plt.colorbar(mappable = None, label = f'{vcol} (mm)', orientation="vertical", ticks=np.linspace(amin(df[vcol]), amax(df[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
    plt.xlabel(r'$X (mm)$', fontsize=15)
    plt.ylabel(r'$Y (mm)$', fontsize=15)
    plt.title(f'{vcol.capitalize()} Ordinary Kriging Extrapolation from Scale 1')
    plt.plot(expandeddf['x'], expandeddf['y'],  color="blue", linewidth = 1)
    plt.plot(dfvs['x'], dfvs['y'],  color="red", linewidth = 0.5)
    if save == True:
        plt.savefig(f"{wdsave}{inputfilenameu}KrigBPAperMap{extension_png}", dpi=1000, bbox_inches = "tight")
        krigdf = pd.DataFrame(list(zip(pointsdf['x'], pointsdf['y'], zz.flatten())), columns = ['x', 'y', vcol])
        krigdf.to_csv(f"{wdsave}{inputfilenameu}KrigBPAperMap{extension_csv}")
    
    
elif style == "grid":
    im = plt.contourf(xx, yy, zz, cmap = cmap, vmin = zz.min(), vmax = zz.max(), levels = np.linspace(zz.min(), zz.max(), 100),)                                               ####################################################
    plt.show()
    cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(zz.min(), zz.max(), 10))                                                                              ## contour map, colour bar and colourbar label??? ##
    cbar.set_label(f'Blind Prediction {vcol} Value', rotation=270, labelpad=20)                                                                                                              ####################################################
    plt.imshow(zz, vmin = zz.min(), vmax = zz.max(), extent = [xx.min(),xx.max(),yy.min(),yy.max()], cmap = cmap) # 
    plt.title("Blind Prediction Around map " + f"\nBased on {inputfilenameu}\nand Upscale{scale}"); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
    plt.tight_layout()
    if save == True:
        plt.savefig(f"{wdsave}{inputfilenameu}KrigBPAperMap{extension_png}", dpi=1000, bbox_inches = "tight")
        krigdf = pd.DataFrame(list(zip(xx.flatten(), yy.flatten(), zz.flatten())), columns = ['x', 'y', vcol])
        krigdf.to_csv(f"{wdsave}{inputfilenameu}KrigBPAperMap{extension_csv}")
plt.show()

#================================================================================================================================================================================ 
print("Kriging complete.")


duration()

