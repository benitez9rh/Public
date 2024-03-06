# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:49:16 2021

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
import geostatspynscore                                             # Copied function and function dependencies from geostatspy GitHub
import time
from time import ctime
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"

def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f'\n**Duration:**\ndays: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
def isanydigit(n: str) -> bool: # Tests if is a digit because isdigit() doesn't recognise floats
    try:
        float(n)
        n.isdigit()
        return True
    except ValueError:
        return False

# =============================================================================
# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Freiberg Gneiss/aperture').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'FG_bottom_ctr_trim vs BPred z_ORIGvsORIG-r(0, 0, 0)-offset(0, 0, 0)_ErrorMap').as_posix()
# wdsave = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Freiberg Gneiss/aperture/dx3dy1').as_posix() 
# os.chdir(wd)                                   # set the working directory
# =============================================================================
# =============================================================================

#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename_a = pathlib.PureWindowsPath(r'Extracted M CNL_Aperture').as_posix()                             # Use this instead in Windows

vcol = 'ApertureR'
save = True
plot = True
#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT AND NORMALISE DATAFRAME ################################################################
###########################################################################################################################################################'''
#df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
df = pd.read_csv(wd + inputfilename + extension_csv, index_col=0)                  #This new line is because reading the .csv file with a header and an index column will create some issue
# ResVarMap = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                              # Read csv's of the variogram maps


# =============================================================================
# df[f'Z{vcol}'] = zscore(df[f'{vcol}']) #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
# # df['Nzr'], tvz, tnsz = geostats.nscore(df, 'zr') # nscore transform for all facies porosity #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
#df[f'N{vcol}'], tvz, tnsz = geostatspynscore.nscore(df, f'{vcol}')
df[f'N{vcol}'], tvz, tnsz = geostats.nscore(df, f'{vcol}')
# =============================================================================

stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
vmin = stats['min'][f'{vcol}']; vmax = stats['max'][f'{vcol}']; nvmin = stats['min'][f'N{vcol}']; nvmax = stats['max'][f'N{vcol}']


# =============================================================================
# scale = [i for i in ''.join((ch if ch in '0123456789.-e' else ' ') for ch in inputfilename[inputfilename.find("Upscaled")+len("Upscaled"):inputfilename.find("Upscaled")+len("Upscaled")+2]).split() if isanydigit(i) == True][-1] # Extracts the scale. 2 is hardcoded because at most my scales have 2 digits https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python
# inputfilename = inputfilename[inputfilename.find("GW"):inputfilename.find("GW")+2] + inputfilename[inputfilename.find("GW")+3:inputfilename.find("GW")+5] + f'_Upscaled{scale}_apertureR_FractureMap'
# wdsave = pathlib.PureWindowsPath(wd + r'\Upscale'+f'{scale}').as_posix() #flips the backslashes to plug in open()
# bs="\\"; wdsave=wdsave+bs                               # Use this instead in Windows
# =============================================================================
 

''' Variogram map parameters '''
tmin=-999; tmin=int(tmin) #trim values below this
tmax=999 ; tmax=int(tmax)#trim values above this
dxlag=1; dxlag=int(dxlag)#lag size in x-direction (i.e., bin size in x-direction); As a first pass. use  2
dylag=1; dylag=int(dylag) #lag size in y-direction (i.e., bin size in y-direction); As a first pass. use  2
nxlag = int(round( (xmax - xmin ) / dxlag )) #number of lags in x-direction from a central (0, 0) point (excluding); As a first pass. use  xtotal/dxlag
nylag = int( round( ( ymax - ymin ) / dylag )) #number of lags in y-direction from a central (0, 0) point (excluding); As a first pass. use  ytotal/dylag
minnp=1 ; minnp=int(minnp)#minimum number of points in lag bin need to be minnp+1
isill=1; isill=int(isill) #sill
# step=1

''' Semi-Variograms' parameters '''
# tmin = -9999.; tmax = 9999.                             # This has to be set before the semi-variograms down below or it will mess up the variogram map
lag_dist = dxlag; lag_tol = (dxlag/2); nlag = nxlag;            # maximum lag is 700m and tolerance > 1/2 lag distance for smoothing
bandh = 9999.9; atol = 22.5                             # no bandwidth, directional variograms
isill=1; isill=int(isill)                                           # standardize sill
azi_mat = [67.5, 157.5, 0.0, 22.5, 45.0, 90.0, 112.5, 135]           # directions in azimuth to consider





#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################ ORIGINAL AND NORMALISED HISTOGRAMS ##############################################################
###########################################################################################################################################################'''
stop0 = time.time()
# Distribution
plt.subplot(121)                                        # plot original and normalised histograms
plt.hist(df[f'{vcol}'], facecolor='red',bins=np.linspace(df[f'{vcol}'].min(),df[f'{vcol}'].max(),1000), histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black',label='Original')
# plt.xlim([0.05,0.25
# plt.ylim([0,5.0])
plt.xlabel(f'{vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} {vcol}')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(122)  
plt.hist(df[f'N{vcol}'], facecolor='blue',bins=np.linspace(-3.0,3.0,1000),histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black',label = 'Normal-Transformed')
plt.xlim([-3.0,3.0]); plt.ylim([0,1.0])
plt.xlabel(f'Nscore {vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} Nscore {vcol}')
plt.legend(loc='upper left')
plt.grid(True)
plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=2.2, wspace=0.2, hspace=0.3)
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_{vcol}_Hist{extension_png}',dpi=1000, bbox_inches = "tight")
if plot == True:
    plt.show()
else:
    plt.close()

# With fitted lognormal and normal dists
# =============================================================================
# mu, sigma = df['zr'].mean(), df['zr'].std()
# x = np.linspace(0, abs(max(df['zr']))+abs(min(df['zr'])), 10000)
# pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
#        / (x * sigma * np.sqrt(2 * np.pi)))
# #x = x - abs(min(df['zr']))
# x = x - 1.15
# 
# stop0 = time.time()
# # Distribution
# plt.subplot(121)                                        # plot original and normalised histograms
# plt.hist(df['zr'], facecolor='red',bins=np.linspace(df['zr'].min(),df['zr'].max(),1000), histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black',label='Original')
# plt.plot(x, pdf, linewidth=2, color='r', label = f'Lognormal Distribution with mu: {mu:.2f}, sigma: {sigma:.2f}')
# # plt.xlim([0.05,0.25
# # plt.ylim([0,5.0])
# plt.xlabel('Residual-Z'); plt.ylabel('Frequency'); plt.title('Residual-Z')
# plt.legend(loc='upper left')
# plt.grid(True)
# 
# mu, sigma = df['Nzr'].mean(), df['Nzr'].std()
# mu = 0
# variance = 1
# sigma = math.sqrt(variance)
# x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
# 
# 
# plt.subplot(122)  
# plt.hist(df['Nzr'], facecolor='blue',bins=np.linspace(-3.0,3.0,1000),histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black',label = 'Normal-Transformed')
# plt.plot(x, stats.norm.pdf(x, mu, sigma), linewidth = 2, color = 'r', label = f'Normal Distribution with mu: {mu:.2f}, sigma: {sigma:.2f}')
# plt.xlim([-3.0,3.0]); plt.ylim([0,1.0])
# plt.xlabel('Nscore Residual-Z'); plt.ylabel('Frequency'); plt.title('Nscore Residual-Z')
# plt.legend(loc='upper left')
# plt.grid(True)
# plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=2.2, wspace=0.2, hspace=0.3)
# plt.savefig(f'{wd}{inputfilename}_Res_FitHist{extension_png}',dpi=600, bbox_inches = "tight")
# plt.show()
# =============================================================================


# Cumulative
plt.subplot(121)                                        # plot original and normalised histograms
plt.hist(df[f'{vcol}'], facecolor='red',bins=np.linspace(df[f'{vcol}'].min(),df[f'{vcol}'].max(),1000), histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label='Original')
# plt.xlim([0.05,0.25
# plt.ylim([0,5.0])
plt.xlabel('Original Cumulative {vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} Cumulative {vcol}')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(122)  
plt.hist(df[f'N{vcol}'], facecolor='blue',bins=np.linspace(-3.0,3.0,1000),histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label = 'Normal-Transformed')
plt.xlim([-3.0,3.0]); plt.ylim([0,1.0])
plt.xlabel(f'Nscore Cumulative {vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} Nscore Cumulative {vcol}')
plt.legend(loc='upper left')
plt.grid(True)
plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=2.2, wspace=0.2, hspace=0.3)
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_N{vcol}_Cum_Hist{extension_png}',dpi=1000, bbox_inches = "tight")
if plot == True:
    plt.show()
else:
    plt.close()








""" Fracture Map """
cmap = plt.cm.plasma                    # color map
# Original
plt.scatter(df['x'], df['y'], c=df[f'{vcol}'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
plt.colorbar(mappable = None, label = f'{vcol} Value', orientation="vertical", ticks=np.linspace(amin(df[f'{vcol}']), amax(df[f'{vcol}']), 10))                                                                                      #contour map, colour bar and colourbar label???
plt.xlabel(r'$X (mm)$', fontsize=15)
plt.ylabel(r'$Y (mm)$', fontsize=15)
plt.title(f'{inputfilename} Original {vcol} Fracture map')
plt.axis('equal')
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_{vcol}FractureMapOrig{extension_png}',dpi=1000, bbox_inches = "tight")
plt.tight_layout()
if plot == True:
    plt.show()
else:
    plt.close()
#Normalised
plt.scatter(df['x'], df['y'], c=df[f'N{vcol}'], cmap=cmap, marker = ".", alpha=1) #alpha is transparency
plt.colorbar(mappable = None, label = f'{vcol} Value', orientation="vertical", ticks=np.linspace(amin(df[f'N{vcol}']), amax(df[f'N{vcol}']), 10))                                                                                      #contour map, colour bar and colourbar label???
plt.xlabel(r'$X (mm)$', fontsize=15)
plt.ylabel(r'$Y (mm)$', fontsize=15)
plt.title(f'{inputfilename} Normalised {vcol} Fracture map')
plt.axis('equal')
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_N{vcol}FractureMap{extension_png}',dpi=1000, bbox_inches = "tight")
plt.tight_layout()
if plot == True:
    plt.show()
else:
    plt.close()


#================================================================================================================================================================================                     
'''###########################################################################################################################################################
########################################################################## VARIOGRAM MAP######################################################################
###########################################################################################################################################################'''
# GSLIB's VARMAP program (Deutsch and Journel, 1998) converted from the original Fortran to Python 
# by Michael Pyrcz, the University of Texas at Austin (Jan, 2019)
# Note simplified for 2D, irrelgular data only

def varmapv(df,xcol,ycol,vcol,tmin,tmax,nxlag,nylag,dxlag,dylag,minnp,isill): 
    global gamf, nppf
    # Parameters - consistent with original GSLIB    
    # df - DataFrame with the spatial data, xcol, ycol, vcol coordinates and property columns
    # tmin, tmax - property trimming limits
    # xlag, xltol - lag distance and lag distance tolerance
    # nlag - number of lags to calculate
    # azm, atol - azimuth and azimuth tolerance
    # bandwh - horizontal bandwidth / maximum distance offset orthogonal to azimuth
    # isill - 1 for standardize sill
    
    # Load the data

    df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]    # trim values outside tmin and tmax
    nd = len(df_extract)
    x = df_extract[xcol].values
    y = df_extract[ycol].values
    vr = df_extract[vcol].values
    
    # Summary statistics for the data after trimming
    avg = vr.mean()
    stdev = vr.std()
    sills = stdev**2.0
    ssq = sills
    vrmin = vr.min()
    vrmax = vr.max()   
    
    # Initialize the summation arrays
    npp = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    gam = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    nppf = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    gamf = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    hm = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    tm = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    hv = np.zeros((int(nylag*2+1),int(nxlag*2+1)))
    tv = np.zeros((int(nylag*2+1),int(nxlag*2+1)))

    # First fix the location of a seed point: 
    for i in range(0,nd):
        if i%1000 ==0:                                  # Track the progress
            print('Variogram map Loop 1/3: ' + str(i) + "/" + str(nd) + ' ' + ctime())
            # Second loop over the data: 
        for j in range(0,nd): 
            # The lag:
            ydis = y[j] - y[i]
            iyl = int(nylag) + int(ydis/dylag)
            if iyl < 0 or iyl > nylag*2: # acocunting for 0,...,n-1 array indexing
                continue
            xdis = x[j] - x[i]
            ixl = nxlag + int(xdis/dxlag)
            if ixl < 0 or ixl > nxlag*2: # acocunting for 0,...,n-1 array indexing
                continue
                
            # We have an acceptable pair, therefore accumulate all the statistics
            # that are required for the variogram:
            npp[iyl,ixl] = npp[iyl,ixl] + 1 # our ndarrays read from the base to top, so we flip
            tm[iyl,ixl] = tm[iyl,ixl] + vr[i]
            hm[iyl,ixl] = hm[iyl,ixl] + vr[j]
            tv[iyl,ixl] = tm[iyl,ixl] + vr[i]*vr[i]
            hv[iyl,ixl] = hm[iyl,ixl] + vr[j]*vr[j]
            gam[iyl,ixl] = gam[iyl,ixl] + ((vr[i]-vr[j])*(vr[i]-vr[j]))

    # Get average values for gam, hm, tm, hv, and tv, then compute
    # the correct "variogram" measure:
    for iy in range(0,nylag*2+1):
        if iy%1000 ==0:                                  # Track the progress
            print('Variogram map Loop 2/3: ' + str(iy) + "/" + str(nylag*2+1) + ' ' + ctime())
            
        for ix in range(0,nxlag*2+1): 
            if npp[iy,ix] <= minnp:
                gam[iy,ix] = -999.
                hm[iy,ix]  = -999.
                tm[iy,ix]  = -999.
                hv[iy,ix]  = -999.
                tv[iy,ix]  = -999.
            else:
                rnum = npp[iy,ix]
                gam[iy,ix] = gam[iy,ix] / (2*rnum) # semivariogram
                hm[iy,ix] = hm[iy,ix] / rnum
                tm[iy,ix] = tm[iy,ix] / rnum
                hv[iy,ix] = hv[iy,ix] / rnum - hm[iy,ix]*hm[iy,ix]
                tv[iy,ix] = tv[iy,ix] / rnum - tm[iy,ix]*tm[iy,ix]
                
                # Attempt to standardize:
            if isill > 0:
                gamf[iy,ix] = gamf[iy,ix]/sills

    for iy in range(0,nylag*2+1):
        if iy%1000 ==0:                                  # Track the progress
            print('Variogram map Loop 3/3: ' + str(iy) + "/" + str(nylag*2+1) + ' ' + ctime())
            
        for ix in range(0,nxlag*2+1):             
            gamf[iy,ix] = gam[nylag*2-iy,ix]
            nppf[iy,ix] = npp[nylag*2-iy,ix]
            
    return gamf, nppf
#================================================================================================================================================================================ 
stop0 = time.time()
vmap, npmap = varmapv(df,'x','y',f'N{vcol}',tmin=tmin,tmax=tmax,nxlag=nxlag,nylag=nylag,dxlag=dxlag,dylag=dylag,minnp=minnp,isill=isill)

xmin = -(dxlag*nxlag)-(dxlag/2); ymin = -(dylag*nylag)-(dylag/2); xmax = (dxlag*nxlag)+(dxlag/2); ymax = (dylag*nylag)+(dylag/2); vmin = np.amin(gamf[gamf!=tmin])  ; vmax = np.amax(gamf[gamf!=tmax])
plt.subplot(111)
xx, yy = np.meshgrid(np.arange(xmin, xmax, dxlag), np.arange(ymax, ymin, -1 * dylag))
im = plt.contourf(xx,yy,vmap,cmap=cmap,vmin=vmin,vmax=vmax,levels=np.linspace(vmin, vmax, 100),)                                                                          ####################################################
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(vmin, vmax, 10))                                                                                        #contour map, colour bar and colourbar label???
cbar.set_label('Variogram Value', rotation=270, labelpad=20)                                                                                                              ####################################################
plt.imshow(vmap,interpolation = None, extent = [xmin,xmax,ymin,ymax], vmin = vmin, vmax = vmax, cmap = cmap)
plt.title(f'{inputfilename} {vcol} Variogram Map'); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')

plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
plt.tight_layout()
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_{vcol}_VariogramMap{extension_png}', bbox_inches = 'tight')
if plot == True:
    plt.show()
else:
    plt.close()
if save == True:
    exec(f'pd.DataFrame(vmap).to_csv("{wdsave}{inputfilename}_{vcol}_VariogramMap{extension_csv}")')
    exec(f'pd.DataFrame(npmap).to_csv("{wdsave}{inputfilename}_{vcol}_NPMap{extension_csv}")')
duration()
print('The shape of the output is ' + str(vmap.shape))
#================================================================================================================================================================================ 
'''###########################################################################################################################################################
########################################################################## SEMI-VARIOGRAMS ###################################################################
###########################################################################################################################################################'''
# Define gamv to avoid importing geostatspy module because it is always giving trouble. Thanks anyway Geostatsguy!
def gamv(
    df,
    xcol,
    ycol,
    vcol,
    tmin,
    tmax,
    xlag,
    xltol,
    nlag,
    azm,
    atol,
    bandwh,
    isill,
):
    """GSLIB's GAMV program (Deutsch and Journel, 1998) converted from the
    original Fortran to Python by Michael Pyrcz, the University of Texas at
    Austin (Jan, 2019).
    Note simplified for 2D, semivariogram only and one direction at a time.
    :param df: pandas DataFrame with the spatial data
    :param xcol: name of the x coordinate column
    :param ycol: name of the y coordinate column
    :param vcol: name of the property column
    :param tmin: property trimming limit
    :param tmax: property trimming limit
    :param xlag: lag distance
    :param xltol: lag distance tolerance
    :param nlag: number of lags to calculate
    :param azm: azimuth
    :param atol: azimuth tolerance
    :param bandwh: horizontal bandwidth / maximum distance offset orthogonal to
                   azimuth
    :param isill: 1 for standardize sill
    :return: TODO
    """
    # Load the data
    # Trim values outside tmin and tmax
    df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]
    nd = len(df_extract)  # TODO: not used
    x = df_extract[xcol].values
    y = df_extract[ycol].values
    vr = df_extract[vcol].values

    # Summary statistics for the data after trimming
    avg = vr.mean()  # TODO: not used
    stdev = vr.std()
    sills = stdev ** 2.0
    ssq = sills  # TODO: not used
    vrmin = vr.min()  # TODO: not used
    vrmax = vr.max()  # TODO: not used

    # Define the distance tolerance to half lag distance, if it isn't already
    if xltol < 0.0:
        xltol = 0.5 * xlag

    # Loop over combinatorial of data pairs to calculate the variogram
    dis, vario, npp = variogram_loop(
        x, y, vr, xlag, xltol, nlag, azm, atol, bandwh
    )

    # Standardize sill to one by dividing all variogram values by the variance
    for il in range(0, nlag + 2):
        if il%10 == 0:                                  # Track the progress
            print('gamv Loop: ' + str(il+1)+"/" + str(nlag + 2) + ' ' + ctime())
            
        if isill == 1:
            vario[il] = vario[il] / sills

        # Apply 1/2 factor to go from variogram to semivariogram
        vario[il] = 0.5 * vario[il]

    return dis, vario, npp

def variogram_loop(x, y, vr, xlag, xltol, nlag, azm, atol, bandwh):
    """Calculate the variogram by looping over combinatorial of data pairs.
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
    # The mathematical azimuth is measured counterclockwise from EW and
    # not clockwise from NS as the conventional azimuth is
    azmuth = (90.0 - azm) * math.pi / 180.0 # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    uvxazm = math.cos(azmuth) # unit vector x-component (in the unit circle)
    uvyazm = math.sin(azmuth) # unit vector y-component (in the unit circle)
    if atol <= 0.0:
        csatol = math.cos(45.0 * math.pi / 180.0) # if angle tolerance is negative or 0, angle tolerance is set to 45degrees
    else:
        csatol = math.cos(atol * math.pi / 180.0) # cosine of angle tolerance in rads

    # Initialize the arrays for each direction, variogram, and lag
    nsiz = nlag + 2  # TODO: not used
    dismxs = ((float(nlag) + 0.5 - EPSLON) * xlag) ** 2 # maximum distance squared, (nlag*xlag)**2

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

                        # Reject this pair on the basis of missing values

                        # Data was trimmed at the beginning

                        # The Semivariogram (all other types of measures are
                        # removed for now)
                        for il in range(lagbeg, lagend + 1):
                            npp[il] = npp[il] + 1
                            dis[il] = dis[il] + h
                            tm[il] = tm[il] + vrt
                            hm[il] = hm[il] + vrh
                            vario[il] = vario[il] + ((vrh - vrt) * (vrh - vrt))
                            if omni:
                                npp[il] = npp[il] + 1.0
                                dis[il] = dis[il] + h
                                tm[il] = tm[il] + vrtpr
                                hm[il] = hm[il] + vrhpr
                                vario[il] = vario[il] + (
                                    (vrhpr - vrtpr) * (vrhpr - vrtpr)
                                )

    # Get average values for gam, hm, tm, hv, and tv, then compute the correct
    # "variogram" measure
    for il in range(0, nlag + 2):
        if il%10_000 ==0:                                  # Track the progress
            print('Direction: ' + str(iazi+1) + '/' + str(len(azi_mat)) + '.  Variogram Loop 2/2: ' + str(il)+"/" + str(nlag + 2) + ' ' + ctime())
            
        i = il
        if npp[i] > 0:
            rnum = npp[i]
            dis[i] = dis[i] / rnum
            vario[i] = vario[i] / rnum
            hm[i] = hm[i] / rnum
            tm[i] = tm[i] / rnum

    return dis, vario, npp


tmin = -9999.; tmax = 9999.                             # no trimming. This needs to be set down here or it will mess up the vaiogram map above.
lag = np.zeros((len(azi_mat),int(nlag)+2)); gamma = np.zeros((len(azi_mat),int(nlag)+2)); npp = np.zeros((len(azi_mat),int(nlag+2)));
stop0 = time.time()
for iazi in range(0,len(azi_mat)):                      # Loop over all directions
    print('Semi-variogram direction ' + str(iazi+1)+"/"+str(len(azi_mat)) + ' ' + ctime())                                   # Track the progress   
    lag[iazi,:], gamma[iazi,:], npp[iazi,:] = gamv(df,"x","y", vcol, tmin,tmax,lag_dist,lag_tol,nlag,azi_mat[iazi],atol,bandh,isill)
    iazistr = f'{azi_mat[iazi]}';
    if "." in iazistr: iazistr = iazistr.replace(".","p"); # A variable name cannot have "." in it so we replace it with "p"
    if save == True:
        exec(f'SemiVariogram_Dir{iazistr} = pd.DataFrame(list(zip(lag[iazi, :],gamma[iazi,:])), columns=["lag","gamma"] )') # Createa pandas DataFrame with the semi-variogram values for that direction
        exec(f'SemiVariogram_Dir{iazistr}.to_csv("{wdsave}lag{dxlag}_SemiVariogram_Dir{iazistr}.csv")') # and save it in a .csv
    #plt.subplot(4,2,iazi+1)
    plt.plot(lag[iazi,:],gamma[iazi,:], '.', markersize = 4, color = 'black', label = f"Azimuth {azi_mat[iazi]} Experimental variogram")
    plt.plot([0,lag[iazi,:].max()],[1.0,1.0],color = 'black')
    plt.title(f'{inputfilename}\nNormalised {vcol} {dxlag}_{int(azi_mat[iazi])}deg Semi-Variogram')
    plt.xticks(range(0, int(max(lag[iazi,:])), 10))
    plt.grid(which="major", color='#CCCCCC', linestyle='--')
    plt.grid(which='minor', color='#CCCCCC', linestyle=':')
    plt.minorticks_on()
    plt.annotate(f"x-Lag: {dxlag}\ny-Lag: {dylag}" , (lag[iazi,:].max()-15, 0.3))
    # plt.xlim([0, lag_dist*(nlag+1) ])
    # plt.ylim([0, np.amax(gamf[gamf!=tmax])+0.25 ])
    plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
    plt.ylabel(r'$\gamma \bf(h)$')
    plt.legend(loc="lower right")
    #plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=4.2, wspace=0.2, hspace=0.3)
    if save == True:
        plt.savefig(f'{wdsave}{inputfilename}_Res_Semi-variograms_lag{dxlag}_{int(azi_mat[iazi])}deg{extension_png}',dpi=1000, bbox_inches = "tight")
    if plot == True:
        plt.show()
    else:
        plt.close()

#================================================================================================================================================================================ 

duration()





'''
https://www.youtube.com/watch?v=bryRCrtf3hk&list=PLG19vXLQHvSB-D4XKYieEku9GQMQyAzjJ&index=37
https://github.com/GeostatsGuy/PythonNumericalDemos/blob/master/GeostatsPy_spatial_continuity_directions.ipynb


plt.subplot(221)                                        # plot original sand and shale porosity histograms
plt.hist(df['Porosity'], facecolor='red',bins=np.linspace(0.0,0.25,1000),histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label='Original')
plt.xlim([0.05,0.25]); plt.ylim([0,1.0])
plt.xlabel('Porosity (fraction)'); plt.ylabel('Frequency'); plt.title('Porosity')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(222)  
plt.hist(df['NPor'], facecolor='blue',bins=np.linspace(-3.0,3.0,1000),histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label = 'Trans')
plt.xlim([-3.0,3.0]); plt.ylim([0,1.0])
plt.xlabel('Porosity (fraction)'); plt.ylabel('Frequency'); plt.title('Nscore Porosity')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(223)                                        # plot nscore transformed sand and shale histograms
plt.hist(df['Perm'], facecolor='red',bins=np.linspace(0.0,1000.0,100000),histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label='Original')
plt.xlim([0.0,1000.0]); plt.ylim([0,1.0])
plt.xlabel('Porosity (fraction)'); plt.ylabel('Frequency'); plt.title('Permeability')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplot(224)                                        # plot nscore transformed sand and shale histograms
plt.hist(df['NPerm'], facecolor='blue',bins=np.linspace(-3.0,3.0,100000),histtype="stepfilled",alpha=0.2,density=True,cumulative=True,edgecolor='black',label = 'Trans')
plt.xlim([-3.0,3.0]); plt.ylim([0,1.0])
plt.xlabel('Permeability (mD)'); plt.ylabel('Frequency'); plt.title('Nscore Permeability')
plt.legend(loc='upper left')
plt.grid(True)

plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=2.2, wspace=0.2, hspace=0.3)
plt.show()





plt.subplot(131)
GSLIB.locmap_st(df,'X','Y','NPor',0,1000,0,1000,-3,3,'Nscore Porosity - All Facies','X (m)','Y (m)','Nscore Porosity',cmap)

plt.subplot(132)
GSLIB.locmap_st(df,'X','Y','NPerm',0,1000,0,1000,-3,3,'Nscore Permeability - All Facies','X (m)','Y (m)','Nscore Permeability',cmap)

plt.subplot(133)
facies = df['Facies'].values +0.01
plt.scatter(df['NPor'],df['NPerm'],c = facies,edgecolor = 'black',cmap = plt.cm.inferno)
#plt.plot([-3,3],[-3,3],color = 'black')
plt.xlabel(r'Nscore Porosity')
plt.ylabel(r'Nscore Permeability')
plt.title('Nscore Permeability vs. Porosity')
plt.xlim([-3,3])
plt.ylim([-3,3])
plt.grid(True)

plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=0.8, wspace=0.5, hspace=0.3)
plt.show()


# GSLIB's VARMAP program (Deutsch and Journel, 1998) converted from the original Fortran to Python 
# by Michael Pyrcz, the University of Texas at Austin (Jan, 2019)
# Note simplified for 2D, irrelgular data only

def varmapv(df,xcol,ycol,vcol,tmin,tmax,nxlag,nylag,dxlag,dylag,minnp,isill): 

# Parameters - consistent with original GSLIB    
# df - DataFrame with the spatial data, xcol, ycol, vcol coordinates and property columns
# tmin, tmax - property trimming limits
# xlag, xltol - lag distance and lag distance tolerance
# nlag - number of lags to calculate
# azm, atol - azimuth and azimuth tolerance
# bandwh - horizontal bandwidth / maximum distance offset orthogonal to azimuth
# isill - 1 for standardize sill

# Load the data

    df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]    # trim values outside tmin and tmax
    nd = len(df_extract)
    x = df_extract[xcol].values
    y = df_extract[ycol].values
    vr = df_extract[vcol].values
    
# Summary statistics for the data after trimming
    avg = vr.mean()
    stdev = vr.std()
    sills = stdev**2.0
    ssq = sills
    vrmin = vr.min()
    vrmax = vr.max()   
    
# Initialize the summation arrays
    npp = np.zeros((nylag*2+1,nxlag*2+1))
    gam = np.zeros((nylag*2+1,nxlag*2+1))
    nppf = np.zeros((nylag*2+1,nxlag*2+1))
    gamf = np.zeros((nylag*2+1,nxlag*2+1))
    hm = np.zeros((nylag*2+1,nxlag*2+1))
    tm = np.zeros((nylag*2+1,nxlag*2+1))
    hv = np.zeros((nylag*2+1,nxlag*2+1))
    tv = np.zeros((nylag*2+1,nxlag*2+1))

# First fix the location of a seed point: 
    for i in range(0,nd):     
# Second loop over the data: 
        for j in range(0,nd): 
# The lag:
            ydis = y[j] - y[i]
            iyl = nylag + int(ydis/dylag)
            if iyl < 0 or iyl > nylag*2: # acocunting for 0,...,n-1 array indexing
                continue
            xdis = x[j] - x[i]
            ixl = nxlag + int(xdis/dxlag)
            if ixl < 0 or ixl > nxlag*2: # acocunting for 0,...,n-1 array indexing
                continue
                
# We have an acceptable pair, therefore accumulate all the statistics
# that are required for the variogram:
            npp[iyl,ixl] = npp[iyl,ixl] + 1 # our ndarrays read from the base to top, so we flip
            tm[iyl,ixl] = tm[iyl,ixl] + vr[i]
            hm[iyl,ixl] = hm[iyl,ixl] + vr[j]
            tv[iyl,ixl] = tm[iyl,ixl] + vr[i]*vr[i]
            hv[iyl,ixl] = hm[iyl,ixl] + vr[j]*vr[j]
            gam[iyl,ixl] = gam[iyl,ixl] + ((vr[i]-vr[j])*(vr[i]-vr[j]))

# Get average values for gam, hm, tm, hv, and tv, then compute
# the correct "variogram" measure:
    for iy in range(0,nylag*2+1): 
        for ix in range(0,nxlag*2+1): 
            if npp[iy,ix] <= minnp:
                gam[iy,ix] = -999.
                hm[iy,ix]  = -999.
                tm[iy,ix]  = -999.
                hv[iy,ix]  = -999.
                tv[iy,ix]  = -999.
            else:
                rnum = npp[iy,ix]
                gam[iy,ix] = gam[iy,ix] / (2*rnum) # semivariogram
                hm[iy,ix] = hm[iy,ix] / rnum
                tm[iy,ix] = tm[iy,ix] / rnum
                hv[iy,ix] = hv[iy,ix] / rnum - hm[iy,ix]*hm[iy,ix]
                tv[iy,ix] = tv[iy,ix] / rnum - tm[iy,ix]*tm[iy,ix]
                
# Attempt to standardize:
            if isill > 0:
                gamf[iy,ix] = gamf[iy,ix]/sills

    for iy in range(0,nylag*2+1): 
        for ix in range(0,nxlag*2+1):             
            gamf[iy,ix] = gam[nylag*2-iy,ix]
            nppf[iy,ix] = npp[nylag*2-iy,ix]
            
    return gamf, nppf



vmap, npmap = varmapv(df,'X','Y','NPor',tmin=-999,tmax=999,nxlag=11,nylag=11,dxlag=50,dylag=50,minnp=1,isill=1)

xmin = -575;ymin = -575; xmax = 575; ymax = 575; step = 50; vmin = 0.0; vmax = 1.6 
plt.subplot(111)
xx, yy = np.meshgrid(np.arange(xmin, xmax, step), np.arange(ymax, ymin, -1 * step))
#im = plt.contourf(xx,yy,vmap,cmap=cmap,vmin=vmin,vmax=vmax,levels=np.linspace(vmin, vmax, 100),)
plt.imshow(vmap,interpolation = None,extent = [xmin,xmax,ymin,ymax], vmin = vmin, vmax = vmax,cmap = cmap)
plt.title('Variogram Map'); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(vmin, vmax, 10))
cbar.set_label('Variogram Value', rotation=270, labelpad=20)

plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
plt.show()

print('The shape of the output is ' + str(vmap.shape))



tmin = -9999.; tmax = 9999.                             # no trimming 
lag_dist = 100.0; lag_tol = 50.0; nlag = 7;            # maximum lag is 700m and tolerance > 1/2 lag distance for smoothing
bandh = 9999.9; atol = 22.5                             # no bandwidth, directional variograms
isill = 1                                               # standardize sill
azi_mat = [0,22.5,45,67.5,90,112.5,135,157.5]           # directions in azimuth to consider


# Arrays to store the results
lag = np.zeros((len(azi_mat),nlag+2)); gamma = np.zeros((len(azi_mat),nlag+2)); npp = np.zeros((len(azi_mat),nlag+2));

for iazi in range(0,len(azi_mat)):                      # Loop over all directions
    lag[iazi,:], gamma[iazi,:], npp[iazi,:] = geostats.gamv(df,"X","Y","NPor",tmin,tmax,lag_dist,lag_tol,nlag,azi_mat[iazi],atol,bandh,isill)
    plt.subplot(4,2,iazi+1)
    plt.plot(lag[iazi,:],gamma[iazi,:],'x',color = 'black',label = 'Azimuth ' +str(azi_mat[iazi]))
    plt.plot([0,2000],[1.0,1.0],color = 'black')
    plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
    plt.ylabel(r'$\gamma \bf(h)$')
    plt.title('Directional NSCORE Porosity Variogram')
    plt.xlim([0,700])
    plt.ylim([0,1.8])
    plt.legend(loc='upper left')
    plt.grid(True)

plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=4.2, wspace=0.2, hspace=0.3)
plt.show()








'''














