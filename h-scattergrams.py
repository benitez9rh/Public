# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 15:05:39 2021

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
import geostatspy.geostats as geostats                            # functions I needed directly to this script
import time
from time import ctime
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"

#set working directory and filename
wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans/Q4Aperture/upscale').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_LR_Upscaled1_aperture_FractureMap').as_posix() 

wdsave = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans/Q4Aperture/upscale/upscale9/h-scattergram').as_posix() #flips the backslashes to plug in open()
bs="//"; wdsave=wdsave+bs  
os.chdir(wdsave)                                  # set the working directory
# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Res').as_posix() #flips the backslashes to plug in open()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
# inputfilename = pathlib.PureWindowsPath(r'GW1_Q4Res').as_posix()
# os.chdir(wd)                                   # set the working directory
# =============================================================================

vcol = 'aperture'
#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT AND NORMALISE DATAFRAME ################################################################
###########################################################################################################################################################'''
#df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
df = pd.read_csv(wd + inputfilename + extension_csv,  index_col = 0)                              #This new line is because reading the .csv file with a header and an index column will create some issue
# ResVarMap = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                              # Read csv's of the variogram maps


df[f'Z{vcol}'] = zscore(df[vcol]) #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
df[f'N{vcol}'], tvz, tnsz = geostats.nscore(df, vcol) # nscore transform for all facies porosity #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba





vr = f'N{vcol}'


stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
vmin = stats['min'][vr]; vmax = stats['max'][vr]; nvmin = stats['min'][vr]; nvmax = stats['max'][vr]


''' Variogram map parameters '''
tmin=-999 #trim values below this
tmax=999 #trim values above this
dxlag=1 #lag size in x-direction (i.e., bin size in x-direction); As a first pass. use  2
dylag=dxlag #lag size in y-direction (i.e., bin size in y-direction); As a first pass. use  2
nxlag = int(round( (xmax - xmin ) / dxlag ) )#number of lags in x-direction from a central (0, 0) point (excluding); As a first pass. use  xtotal/dxlag
nylag = int(round( ( ymax - ymin ) / dylag ) )#number of lags in y-direction from a central (0, 0) point (excluding); As a first pass. use  ytotal/dylag
minnp=1 #minimum number of points in lag bin need to be minnp+1
isill=1 #sill
# step=1

''' Semi-Variograms' parameters '''
tmin = -9999.; tmax = 9999.                             # This has to be set before the semi-variograms down below or it will mess up the variogram map
lag_dist = dxlag; lag_tol = (dxlag/2); nlag = nxlag;            # maximum lag is 700m and tolerance > 1/2 lag distance for smoothing
bandh = 9999.9; atol = 22.5                             # no bandwidth, directional variograms
isill = 1                                               # standardize sill
azi_mat = [0, 22.5, 45, 67.5, 90, 112.5, 135]           # directions in azimuth to consider





#================================================================================================================================================================================ 
'''###########################################################################################################################################################
###################################################################### H-SCATTERGRAMS ########################################################################
###########################################################################################################################################################'''
stop0 = time.time()

def hscattergram(df, x, y, vr, xlag, ylag, xltol, nlag, azm, atol, bandwh):
    """Calculate the h-scattergram by looping over combinatorial of data pairs.
    :param x: x values
    :param y: y values
    :param vr: property values
    :param xlag: lag distance
    :param xltol: lag distance tolerance
    :param nlag: number of lags to calculate
    :param azm: azimuth
    :param atol: azimuth tolerance
    :param bandwh: horizontal bandwidth / maximum distance offset orthogonal to
                   azimuth
    :return: TODO
    """
    # Load the data
    df_extract = df.loc[(df[vr] >= tmin) & (df[vr] <= tmax)]    # trim values outside tmin and tmax
    nd = len(df_extract)                                            #  number of total points
    x = df_extract[x].values
    y = df_extract[y].values
    vr = df_extract[vr].values
        
    # Summary statistics for the data after trimming
    avg = vr.mean()                                                 # Values average
    stdev = vr.std()                                                # Values standard deviation
    sills = stdev**2.0
    ssq = sills
    vrmin = vr.min()                                                # values minimum
    vrmax = vr.max()                                                # values maximum
        
    comb = int(math.factorial(nd)/(math.factorial(2)*math.factorial(nd-2)))
    
    
        
    # Allocate the needed memory
    nvarg = 1
    mxdlv = nlag + 2  # in gamv the npp etc. arrays go to nlag + 2
    dis = np.zeros(mxdlv)
    lag = np.zeros(mxdlv)  # TODO: not used
    vario = np.zeros(mxdlv)
    hm = np.zeros(mxdlv)
    tm = np.zeros(mxdlv)
    hv = np.zeros(mxdlv)  # TODO: not used
    npp = np.zeros(mxdlv)
    ivtail = np.zeros(nvarg + 2)
    ivhead = np.zeros(nvarg + 2)
    ivtype = np.ones(nvarg + 2)
    ivtail[0] = 0
    ivhead[0] = 0
    ivtype[0] = 0

    EPSLON = 1.0e-20
    nd = len(x)
    
    azmuth = (90.0 - azm) * math.pi / 180.0 # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
    uvxazm = math.cos(azmuth) # u-vector cosine (x component)
    uvyazm = math.sin(azmuth) # u-vector sine (y component)
    if atol <= 0.0:
        csatol = math.cos(45.0 * math.pi / 180.0) # if angle tolerance is not set by the user, i..e == 0, angle tolerance is set to 45 degrees and converted to mathematical azimuth in radians
    else:
        csatol = math.cos(atol * math.pi / 180.0) # cosine of angle tolerance and convertion to mathematical azimuth in radians

    # Initialize the arrays for each direction, variogram, and lag
    nsiz = nlag + 2  # TODO: not used
    dismxs = ((float(nlag) + 0.5 - EPSLON) * xlag) ** 2 # maximum distance in x-direction squared
    # Initialise DataFrame for scattergrams
    hscat = pd.DataFrame(columns = ["direction", "lag", "tail", "head"])
# =============================================================================
#     for ilag in range(1, nlag + 1):
#         exec(f"scatx{lagend} = []")
#         exec(f"scaty{lagend} = []")
# =============================================================================
    # Main loop over all pairs
    for i in range(0, nd):
        if i % 1_000 == 0:
            print(f'Direction: {iazi+1}/{len(azi_mat)}: ' + f'{i:,}/{nd:,}' + ' ' + ctime())
        for j in range(0, nd):
            
            # Definition of the lag corresponding to the current pair
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dxs = dx * dx   # distance in x-direaction squared
            dys = dy * dy   # distance in y-direaction squared
            hs = dxs + dys  # head distance squared
            if hs <= dismxs: #  if head distance squared is within the limit/tolerance
                if hs < 0.0:
                    hs = 0.0
                h = np.sqrt(hs)     # Calculate head distance

                # Determine which lag this is and skip if outside the defined distance tolerance
                if h <= EPSLON:     # If head distance is greater than "infinitesimally" small, lag = 0. h is always positive from now on because we squared the distances and squarerooted it
                    lagbeg = 0
                    lagend = 0
                else:
                    lagbeg = -1
                    lagend = -1
                    for ilag in range(1, nlag + 1):
                        # reduced to -1
                        if (    
                            (xlag * float(ilag - 1) - xltol)
                            <= h
                            <= (xlag * float(ilag - 1) + xltol)
                        ):                                                  # If head distance is within -lag tolerance and + lag tolerance
                            if lagbeg < 0:                                  # If it was larger than infinitesimally small, assign the lag that it corresponds to
                                lagbeg = ilag                               
                            lagend = ilag                                   # If  head distance is within -lag tolerance and + lag tolerance, assign corresponding lag
                if lagend >= 0:                                             # If lag is not infinitesimally small and within the range
                    
                    # Definition of the direction corresponding to the current pair. All directions are considered (overlapping of direction tolerance cones is allowed)

                    # Check for an acceptable azimuth angle
                    dxy = np.sqrt(max((dxs + dys), 0.0))                    # Calculate real distance (couldn't have just used h?)
                    if dxy < EPSLON:
                        dcazm = 1.0     # If distance is infinitesimally small,  = 1
                    else:
                        dcazm = (dx * uvxazm + dy * uvyazm) / dxy    # If non-infiitesimally small, find the magnitude of the u-vector (scaled to the unit circle)

                    # Check the horizontal bandwidth criteria (maximum deviation perpendicular to the specified direction azimuth)
                    band = uvxazm * dy - uvyazm * dx    # Calculate the orthogonal distance of the point pair from the azimuth

                    # Apply all the previous checks at once to avoid a lot of nested if statements
                    if (abs(dcazm) >= csatol) and (abs(band) <= bandwh):    # If magnitude u-vector is => magnitue of the tolerance angle AND orthogonal distance of the u-vector from the azimuth is within the band threshold
                        # Check whether or not an omni-directional variogram is being computed
                        omni = False
                        if atol >= 90.0:
                            omni = True

                        # For this variogram, sort out which is the tail and the head value
                        iv = 0  # hardcoded just one variogram
                        it = ivtype[iv]  # TODO: not used
                        if dcazm >= 0.0:    # if dcazm is positive j is closer to the origin hence it's the tail value, otherwise i is the tail value
                            vrh = vr[i]
                            vrt = vr[j]
                            if omni:
                                vrtpr = vr[i]
                                vrhpr = vr[j]
                        else:
                            vrh = vr[j]
                            vrt = vr[i]
                            if omni:
                                vrtpr = vr[j]
                                vrhpr = vr[i]

                        # Reject this pair on the basis of missing values Data was trimmed at the beginning
                        
                        # assign value of tail and head to corresponding lag
                        hscat.loc[-1] = [azm, lagend, vrt, vrh]
                        hscat.index = hscat.index + 1
                        
    hscat = hscat.sort_index()  # sorting by index
# =============================================================================
#     hscat = pd.DataFrame(list(zip(scatx,scaty)), columns = ["t", "h"] )    # Creates an Dataframe size (comb, 2) with columns tail and head. The replace is for situations when the angle has decimals in which case "." can't be used in the name hence replaced by "p".
# =============================================================================
    for dire in hscat['direction'].unique():
        for lag in exec(f"hscat{azm}['lag'].unique()"):
            hscat_extract = hscat.loc[(hscat['direction'] == dire) & (hscat['lag'] == lag)]
            exec(f'hscat_extract.to_csv("hscat_extract{dire}.csv")') # and save it in a .csv
            plt.scatter( hscat_extract['tail'], hscat_extract['head'], s=1 )
            plt.axline((0, 0), slope=1, color="black", linestyle='--', linewidth = 0.5)
            plt.xlabel('Tail values (Z)'); plt.ylabel('Head values (Z)'); plt.title(f'h-scattergram {vcol} azimuth: {dire}, lag: {int(lag+1)}')
            plt.grid(True)
            plt.show()   
    
    
    
    
stop0 = time.time()
for iazi in range(0, len(azi_mat)):
    #hscattergram(df['x'], df['y'], df[vr], dxlag, dylag, lag_tol, nlag, azi_mat[iazi], atol, bandh)         
    hscattergram(df, 'x', 'y', vr, dxlag, dylag, lag_tol, nlag, azi_mat[iazi], atol, bandh)         

    exec( f'hscattergram{str(azi_mat[iazi]).replace(".","p")} = zeros ' )  # Create a DF with the correct name
#del hscat # Manage memory

duration()