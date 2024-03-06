# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:04:08 2022

@author: s2132627
"""

"""
Difference (or aperture) between original surface 1 and Ordinary Krigged surface 2

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
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
from Euler3DRotation import rotate_np

# set the colormap and centre the colorbar. Thanks to http://chris35wills.github.io/matplotlib_diverging_colorbar/
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work their way either side from a prescribed midpoint value)

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
# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean_Q4').as_posix()
# os.chdir(wd)                                   # set the working directory
# =============================================================================


xcol = 'x'
ycol = 'y'
vcol = 'z'

xdens = 0.5 # density of 0.1 causes error: MemoryError: Unable to allocate 64.9 GiB for an array with shape (952300, 9152) and data type float64
ydens = 0.5
zexag = 0 # vertical exageration (%)

''' Path to .csv containing vcol of surface 1 '''
#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q3').as_posix() #flips the backslashes to plug in open()
inputfilename = pathlib.PureWindowsPath(r'GW1Q3').as_posix()
bs="\\"; wd=wd+bs                               # Use this instead in Windows


''' Use second surface? '''
surf2_switch = 1 #1 is on, 0 is off

xcol2 = 'x'
ycol2 = 'y'
vcol2 = 'zr'


''' Use surfaces offsets (mm)? '''
x1off = 0
y1off = 0
z1off = 0
x2off = 0
y2off = 0
z2off = 0

"""Apply rotations""" # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/
s1xr = 0 # Surface 1 Rotation in the x-direction
s1yr = 0 # Surface 1 Rotation in the y-direction
s1zr = 0 # Surface 1 Rotation in the z-direction
s2xr = 0 # Surface 2 Rotation in the x-direction
s2yr = 0 # Surface 2 Rotation in the y-direction
s2zr = 0 # Surface 2 Rotation in the z-direction

''' Path to .csv containing vcol of surface 2 '''
#set working directory and filename
wd2 = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
inputfilename2 = pathlib.PureWindowsPath(r'GW1Q4').as_posix()
bs="\\"; wd2=wd2+bs                               # Use this instead in Windows

''' Path to .csv containing Upscaled vcol of surface 2 '''
#set working directory and filename
wdu2 = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
inputfilenameu2 = pathlib.PureWindowsPath(r'GW1Q4_SampUps2Res').as_posix()
bs="\\"; wdu2=wdu2+bs                               # Use this instead in Windows

''' Surface 2 input_C ''' # input_C is the plane's trend parameters
input_C2 = pathlib.PureWindowsPath(r'GW1Q4_C').as_posix()

''' Variogram model of surface 2 '''
variogram_model2 = "spherical"
dic2={'sill': 1, 'range': 15, 'nugget': 0}
ang2=float(45)
ratio2 = 15/18
param2 = [1, 40, 0] #sill, range, nugget


''' Apply trend surfaces back? '''
TrendApply_switch = 1 #1 is on, 0 is off

''' Calculate Aperture between krigged surfaces? '''
Aperture_switch = 1 #1 is on, 0 is off

#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT SURFACE 1 ################################################################
###########################################################################################################################################################'''
df = pd.read_csv(wd + inputfilename + extension_csv, index_col = 0);
pseudovcol =  "pseudo" + vcol2
C2 = pd.read_csv(wd2 + input_C2 + extension_csv, index_col = 0); C2 = C2.to_numpy(); C2 = [C2[0][0], C2[1][0], C2[2][0]] # Import surf2 trend parameters
df[pseudovcol] = df[vcol] - (C2[0]*df[xcol] + C2[1]*df[ycol] + C2[2]) # subtracting surf2 trend from surf1 - Because we're blind to surf1 and we don't know what the actual trend is.
vcol = pseudovcol
df = df[[xcol, ycol, vcol]] # Filter DataFrame to only the necessary columns

# Apply offsets for surface 1
df[xcol] = df[xcol] + x1off
df[ycol] = df[ycol] + y1off
df[vcol] = (df[vcol] + z1off) * (1+zexag)


stats=df.describe().transpose(); xmin = stats['min'][xcol]; xmax = stats['max'][xcol]; ymin = stats['min'][ycol]; ymax = stats['max'][ycol];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]

gridx = np.arange(xmin, xmax, xdens)
gridy = np.arange(ymax, ymin, -1*ydens)


#================================================================================================================================================================================ 
'''###########################################################################################################################################################
###################################################################### SECOND SURFACE ########################################################################
###########################################################################################################################################################'''

if surf2_switch == 1:

    
    df2 = pd.read_csv(wd2 + inputfilename2 + extension_csv, index_col = 0);
    df2 = df2[[xcol2, ycol2, vcol2]] # Filter DataFrame to only the necessary columns
    
    
    
    dfu2 = pd.read_csv(wdu2 + inputfilenameu2 + extension_csv, index_col = 0)#; ind = list(dfu2.columns);  # Print out the columns so you know what the np columns are
    
    df_xctr = df['x'].mean()                  # Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
    df_yctr = df['y'].mean()
    dfu2_xctr = dfu2['x'].mean()            # Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
    dfu2_yctr = dfu2['y'].mean()
    ctr_xdiff = (df_xctr - dfu2_xctr )
    ctr_ydiff = (df_yctr - dfu2_yctr )  
    
    dfu2['x'] = dfu2['x'] + ctr_xdiff
    dfu2['y'] = dfu2['y'] + ctr_ydiff
    
    dfu2xmin = dfu2['x'].min(); dfu2xmax = dfu2['x'].max(); dfu2ymin = dfu2['y'].min(); dfu2ymax = dfu2['y'].max()
    
    stats=dfu2.describe().transpose(); xmin = stats['min'][xcol]; xmax = stats['max'][xcol]; ymin = stats['min'][ycol]; ymax = stats['max'][ycol];
    vmin = stats['min'][vcol2]; vmax = stats['max'][vcol2]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]
    
    xdens=0.5
    ydens=0.5
    gridx = np.arange(xmin, xmax, xdens)
    gridy = np.arange(ymax, ymin, -1*ydens)
    
    data = np.array(
        [
            dfu2[xcol2],
            dfu2[ycol2],
            dfu2[vcol2]
        ]
        ).transpose()
    
    if variogram_model2 == "custom":
        OK = OrdinaryKriging(
            data[:, 0],
            data[:, 1],
            data[:, 2],
            variogram_model="custom",
            anisotropy_angle = -ang2,
            variogram_parameters = param2,
            anisotropy_scaling = 15/18,
            verbose=True,
            enable_plotting=True,
        ) # scaling is the ratio between the major and minor directions' ranges
    elif variogram_model2 == "spherical":
        OK = OrdinaryKriging(
            data[:, 0],
            data[:, 1],
            data[:, 2],
            variogram_model="spherical",
            anisotropy_angle = -ang2,
            variogram_parameters = dic2,
            anisotropy_scaling = ratio2,
            verbose=True,
            enable_plotting=True,
        ) # scaling is the ratio between the major and minor directions' ranges
    
    ###############################################################################
    # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
    # grid of points, on a masked rectangular grid of points, or with arbitrary points.
    # (See OrdinaryKriging.__doc__ for more information.)
        z, ss = OK.execute("grid", gridx, gridy)
    ###############################################################################
    
    kt.write_asc_grid(gridx, gridy, z, filename="output.asc") # Writes the kriged grid to an ASCII grid file and plot it.
    # plt.imshow(z, cmap = cmap)
    # plt.show()
    
    xx, yy = np.meshgrid(gridx, gridy)
    im = plt.contourf(xx, yy, z, cmap = cmap, vmin = z.min(), vmax = z.max(), levels = np.linspace(z.min(), z.max(), 100),)                                                                          ####################################################
    cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(z.min(), z.max(), 10))                                                                                        #contour map, colour bar and colourbar label???
    cbar.set_label(vcol+' Value', rotation=270, labelpad=20)                                                                                                              ####################################################
    plt.imshow(z, interpolation = None, extent = [xmin,xmax,ymin,ymax], vmin = z.min(), vmax = z.max(), cmap = cmap)
    basepts = inputfilenameu2.strip('GW1Q4_'); basepts = basepts.strip('_zr_FractureMap')
    plt.title(f"{inputfilename} {vcol} via {inputfilename2} OK {basepts} ang:{ang2} r:{dic2['range']} s:{dic2['sill']} n:{dic2['nugget']} FracMap"); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)', fontsize = 10)
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
    plt.tight_layout()
    plt.show()
    plt.savefig(f"{wd}{inputfilename}_{vcol}_via_{inputfilename2}_OK_{basepts}_ang{ang2}_r{dic2['range']}_s{dic2['sill']}_n{dic2['nugget']}_FracMap{extension_png}")
    
    zz = z.data.flatten() # gets the data from a np masked array and flattens it from 2D to 1D array
    
    # Save as df
    krigdf2 = pd.DataFrame(list(zip(xx.flatten(), yy.flatten(), zz.flatten())), columns = [xcol2, ycol2, vcol2])
    #krigdf2.to_csv(f"{wd2}{inputfilename2}_{vcol}_OK_{basepts}_ang:{ang2:.1f}_r{dic2['range']}_s{dic2['sill']}_n{dic2['nugget']}_FracMap{extension_csv}")
    
    # Apply offsets for surface 2
    krigdf2[xcol2] = krigdf2[xcol2] + x2off
    krigdf2[ycol2] = krigdf2[ycol2] + y2off
    krigdf2[vcol2] = (krigdf2[vcol2] + z2off) * (1 + zexag)


    if TrendApply_switch == 1:
        krigdf2['z'] = krigdf2[vcol2] + (C2[0]*krigdf2['x'] + C2[1]*krigdf2['y'] + C2[2])
        krigdf2 = krigdf2[[xcol2, ycol2, 'z']]

        plt.scatter(krigdf2[xcol2], krigdf2[ycol2], c = krigdf2['z'], cmap=cmap, marker = '.', alpha=1)
        plt.colorbar(mappable = None, label = vcol, orientation="vertical", ticks=np.linspace(amin(krigdf2['z']), amax(krigdf2['z']), 10))                                                                                      #contour map, colour bar and colourbar label???
        plt.xlabel(r'$X (mm)$', fontsize=15)
        plt.ylabel(r'$Y (mm)$', fontsize=15)
        plt.title(f"{inputfilename} OKRecon via {inputfilename2} {basepts} ang:{ang2} r:{dic2['range']} s:{dic2['sill']} n:{dic2['nugget']} FracMap"); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')
        plt.tight_layout()
        plt.savefig(f'{wd}{inputfilename}_Recon_{basepts}_FracMap{extension_png}')
        plt.show()
        
    if Aperture_switch == 1:
        arr1 = df.to_numpy();
        krigarr2 = krigdf2.to_numpy(); 
        
        """Apply rotations""" # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/
        # Centralise
        arr1xavg = arr1[:,0].mean();
        arr1yavg = arr1[:,1].mean();
        arr2xavg = krigarr2[:,0].mean();
        arr2yavg = krigarr2[:,1].mean();
        print(arr1xavg, arr1yavg, arr2xavg, arr2yavg)
        arr1[:,0] = arr1[:,0] - arr1xavg; arr1[:,1] = arr1[:,1] - arr1yavg;
        krigarr2[:,0] = krigarr2[:,0] - arr2xavg; krigarr2[:,1] = krigarr2[:,1] - arr2yavg;
        #rotate
        arr1 = rotate_np(arr1, s1xr, s1yr, s1zr)
        krigarr2 = rotate_np(krigarr2, s2xr, s2yr, s2zr)
        # De-centralise
        arr1[:,0] = arr1[:,0] + arr1xavg; arr1[:,1] = arr1[:,1] + arr1yavg;
        krigarr2[:,0] = krigarr2[:,0] + arr2xavg; krigarr2[:,1] = krigarr2[:,1] + arr2yavg;

        """ Creating the PolyData """
        df1cloud = pv.PolyData(arr1)
        krigdf2cloud = pv.PolyData(krigarr2)
        
        """ Create triangulated surfaces """
        Surf_1 = df1cloud.delaunay_2d()
        Surf_2 = krigdf2cloud.delaunay_2d()
        
        diff = np.zeros(Surf_2.n_points)
        Surf_1_verticals = np.zeros((Surf_2.n_points,3))
        Surf_2_verticals = np.zeros((Surf_2.n_points,3))
        for i in range(Surf_2.n_points):
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
                else:
                    dist = np.nan
                    Surf_1_verticals[i] = np.nan
                    Surf_2_verticals[i] = np.nan
            diff[i] = dist
        
        x = Surf_2.points[:,0].T; x = x[~np.isnan(diff)] # Remove indexes where Surf_2_verticals contains NaNs, i.e., all non-intersections
        y = Surf_2.points[:,1].T; y = y[~np.isnan(diff)] # Remove indexes where Surf_2_verticals contains NaNs, i.e., all non-intersections
        Surf_1_verticals = Surf_1_verticals[~np.isnan(diff)]    # Remove indexes where diff contains NaNs, i.e., all non-intersections
        Surf_2_verticals = Surf_2_verticals[~np.isnan(diff)]    # Remove indexes where diff contains NaNs, i.e., all non-intersections
        diff = diff[~np.isnan(diff)]                            # Remove indexes containing NaNs, i.e., all non-intersections

        r2 = 1 - np.sum(diff) / np.sum(np.array([(i-Surf_1_verticals[:,2].mean())**2 for i in Surf_1_verticals[:,2]]) ) # r2 = 1 - ( SSres / SStotal):           SumSquares of residuals = sum of (observed - predicted)**2.        SumSquares total = sum of (observed - observed mean)**2.
        
        #Create df of apertures
        ApertureMap = pd.DataFrame(np.vstack((np.array(x).flatten(), np.array(y).flatten(), diff)).transpose(), columns = [xcol, ycol, "aperture"]).drop_duplicates()
        # Filter out points that are outside the xy boundaries of the upscaled pointset to avoid artefacts at the edges
        ApertureMap = ApertureMap.loc[(ApertureMap[xcol] >= dfu2xmin) & (ApertureMap[xcol] <= dfu2xmax) & (ApertureMap[ycol] >= dfu2ymin) & (ApertureMap[ycol] <= dfu2ymax)]
        ApertureMap.to_csv(f'{wd}{inputfilename}_ExtrapFrom{inputfilename}_OK-rot({s2xr},{s2yr},{s2zr})_ApertureMap{extension_csv}')
        
        """ Fracture Map """
        if amin(diff)<0:
            cmap = plt.cm.seismic                    # color map
            plt.scatter(ApertureMap['x'], ApertureMap['y'], c=ApertureMap['aperture'], cmap=cmap, norm=MidpointNormalize(midpoint=0,vmin=amin(ApertureMap['aperture']), vmax=amax(ApertureMap['aperture'])), marker = ".", alpha=1) #alpha is transparency
        else:
            cmap = plt.cm.plasma                    # color map
            plt.scatter(ApertureMap['x'], ApertureMap['y'], c=ApertureMap['aperture'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
        plt.colorbar(mappable = None, label = 'Aperture (mm)', orientation="vertical", ticks=np.linspace(amin(ApertureMap['aperture']), amax(ApertureMap['aperture']), 10))                                                                                      #contour map, colour bar and colourbar label???
        plt.xlabel(r'$X (mm)$', fontsize=15)
        plt.ylabel(r'$Y (mm)$', fontsize=15)
        plt.title(f'{inputfilename}_{vcol}_ORIGvsOK-rot:({s2xr}, {s2yr}, {s2zr})_ApertureMap')
        # plt.tight_layout()
        plt.savefig(f'{wd}{inputfilename}_{vcol}_ORIGvsOK-rot({s2xr}, {s2yr}, {s2zr})_ApertureMap{extension_png}')
        plt.show()
                    

# Create Aperture Delaunay Surface
ApertureNParr = np.vstack((np.array(x).flatten(), np.array(y).flatten(), Surf_2_verticals[:,2])).transpose() # Create Aperture np.array
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

#Add the trend to the Residual upscaled points, convert to np.array and create pyvista polydata
dfu2b = dfu2[[xcol2, ycol2, vcol2]]; dfu2b[vcol2] = dfu2b[vcol2] + (C2[0]*dfu2b['x'] + C2[1]*dfu2b['y'] + C2[2])
arru2 = dfu2b.to_numpy()
dfu2cloud = pv.PolyData(arru2)

# Visualise kriged surface and upscaled points
p = pv.Plotter();
p.add_mesh(Surf_2, color="b", label = "Kriged Surface", opacity = 0.25); # show_edges=True, 
p.add_mesh(dfu2cloud, color="w", label = "Upscaled points", opacity = 0.25);
p.add_axes(); p.add_legend();
p.add_camera_orientation_widget();
p.show()


# Visualise original vs kriged surfaces
p = pv.Plotter();
p.add_mesh(Surf_1, color="r", show_edges=True, label = "Orig Surface", opacity=1); # show_edges=True, 
p.add_mesh(dfu2cloud, color="w", label = "Upscaled points", opacity = 0.25);
p.add_mesh(Surf_2, color="b", show_edges=True, label = "Kriged Surface", opacity = 0.75); # show_edges=True, 
p.add_axes(); p.add_legend();
p.add_camera_orientation_widget();
p.show()




# Visualise Aperture surface
p = pv.Plotter();
p.add_mesh(Surf_a, scalars = Surf_a.points[:,2], show_edges=True, point_size = 50, label = "Aperture Surface", opacity=0.5); # color="r",scalars = Surf_1.points[:,2], show_edges=True
p.add_axes(); p.add_legend();
p.add_camera_orientation_widget();
p.show()



