# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 13:42:35 2022

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
import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
import geostatspy.geostats as geostats                            # functions I needed directly to this script
import time
from time import ctime
from sklearn.metrics import r2_score
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
cmap = plt.cm.plasma                    # color map

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

#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Res\Upscaled').as_posix() #flips the backslashes to plug in open()
bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
inputfilename = pathlib.PureWindowsPath(r'GW1Q4_LR_Upscaled1_zr_FractureMap').as_posix()
os.chdir(wd)                                   # set the working directory

xcol = 'x'
ycol = 'y'
vcol = 'zr'
#================================================================================================================================================================================ 
'''###########################################################################################################################################################
############################################################## IMPORT AND NORMALISE DATAFRAME ################################################################
###########################################################################################################################################################'''
df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2,3,4,5,6], index_col = 0); ind = list(df.columns);  # Print out the columns so you know what the np columns are
arr = np.genfromtxt(wd + inputfilename + extension_csv, delimiter = ',', skip_header = 1).astype(np.float32) # import the data from the .csv without the first row (all NaNs) as float32 datatype (supposedly faster and less memory hungry at the cost of some precision).
arr = arr[:,1:] # gets rid of first column
permutations = df.shape[0]**2 # permutations w repetition = n**r,  **2 Because there are 2 loops 
#permutations = int((math.factorial(df.shape[0]) / (math.factorial(df.shape[0]-math.factorial(2)))))+1 # Permutations wo repetition = n! / (n!-r!)
# a[:, 0] get the column 0
# a[0, :] get line 0


# Get mins and maxs of x, y, vcol, and Nvcol (whichever vcol is)
# =============================================================================
# xmin = arr[:, ind.index(xcol)].min(); xmax = arr[:, ind.index(xcol)].max(); ymin = arr[:, ind.index(ycol)].min(); ymax = arr[:, ind.index(ycol)].max();
# vmin = arr[:, ind.index(vcol)].min(); vmax = arr[:, ind.index(vcol)].max(); nvmin = arr[:, ind.index('N'+vcol)].min(); nvmax = arr[:, ind.index('N'+vcol)].max()
# =============================================================================
stats=df.describe().transpose(); xmin = stats['min'][xcol]; xmax = stats['max'][xcol]; ymin = stats['min'][ycol]; ymax = stats['max'][ycol];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]

 


''' Variogram map parameters '''
tmin=-999; tmin=int(tmin) #trim values below this
tmax=999 ; tmax=int(tmax)#trim values above this
dxlag=1; dxlag=int(dxlag)#lag size in x-direction (i.e., bin size in x-direction); As a first pass. use  2
dylag=dxlag #lag size in y-direction (i.e., bin size in y-direction); As a first pass. use  2
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
azi_mat = [0.0, 22.5, 45.0, 67.5, 90.0, 112.5, 135.0, 157.5]           # directions in azimuth to consider

# For testing

# arr2 = arr[:1000,0:4]
# azi_mat = [135.0, 45.0, 90.0]           # directions in azimuth to consider
# bandwh = bandh
# azm = 135.0
# xlag=dxlag
# xltol = lag_tol


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
    x = df_extract[xcol].values; #x = x[:10000]
    y = df_extract[ycol].values; #y = y[:10000]
    vr = df_extract[vcol].values; #vr = vr[:10000]
    arr_stk = np.vstack((x,y,vr)).transpose();
    # Summary statistics for the data after trimming
    avg = vr.mean()
    stdev = vr.std()
    sills = stdev**2.0
    ssq = sills
    vrmin = vr.min()
    vrmax = vr.max()   
    
    # Initialize the summation arrays
    npp = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    gam = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    nppf = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    gamf = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    hm = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    tm = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    hv = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))
    tv = np.zeros((int(nylag)*2+1,int(nxlag)*2+1))

    # First fix the location of a seed point: 
    def sub_variomap_loop(i, j, tail, head):
        #for i in range(0,nd):
        if i%1000 ==0 and j == 0:                                  # Track the progress
            print('Variogram map Loop 1/3: ' + str(i) + "/" + str(shape(arr_stk)[0]) + ' ' + ctime())
        # Second loop over the data: 
        # The lag:
        ydis = head[1] - tail[1]                                                 #distance in y direction
        iyl = int(nylag) + int(ydis/dylag)
        if iyl < 0 or iyl > nylag*2:                                        # check the pair has an acceptable distance in x-direction. accounting for 0,...,n-1 array indexing
            pass                                                        # Continue here means it disregards this pair and "continues"/moves onto the next
        xdis = head[0] - tail[0]                                                  # distance in x-direction
        ixl = int(nxlag) + int(xdis/dxlag)
        if ixl < 0 or ixl > nxlag*2:                                        # check the pair has an acceptable distance in x-direction. accounting for 0,...,n-1 array indexing
            pass
        
        # We have an acceptable pair, therefore accumulate all the statistics
        # that are required for the variogram:
        npp[iyl,ixl] = npp[iyl,ixl] + 1 # our ndarrays read from the base to top, so we flip
        tm[iyl,ixl] = tm[iyl,ixl] + tail[2]
        hm[iyl,ixl] = hm[iyl,ixl] + head[2]
        tv[iyl,ixl] = tm[iyl,ixl] + tail[2]*tail[2]
        hv[iyl,ixl] = hm[iyl,ixl] + head[2]*head[2]
        gam[iyl,ixl] = gam[iyl,ixl] + ((tail[2]-head[2])*(tail[2]-head[2]))
        
    
    
    np.array( [sub_variomap_loop(i, j, tail, head) for i, tail in enumerate(arr_stk) for j, head in enumerate(arr_stk)])
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

stop0 = time.time()
vmap, npmap = varmapv(df,xcol, ycol, vcol, tmin=tmin, tmax=tmax, nxlag=nxlag, nylag=nylag, dxlag=dxlag, dylag=dylag, minnp=minnp, isill=isill)

xmin = -(dxlag*nxlag)-(dxlag/2); ymin = -(dylag*nylag)-(dylag/2); xmax = (dxlag*nxlag)+(dxlag/2); ymax = (dylag*nylag)+(dylag/2); vmin = np.amin(gamf[gamf!=tmin])  ; vmax = np.amax(gamf[gamf!=tmax])
plt.subplot(111)
xx, yy = np.meshgrid(np.arange(xmin, xmax, dxlag), np.arange(ymax, ymin, -1 * dylag))
im = plt.contourf(xx,yy,vmap,cmap=cmap,vmin=vmin,vmax=vmax,levels=np.linspace(vmin, vmax, 100),)                                                                          ####################################################
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(vmin, vmax, 10))                                                                                        #contour map, colour bar and colourbar label???
cbar.set_label('Variogram Value', rotation=270, labelpad=20)                                                                                                              ####################################################
plt.imshow(vmap,interpolation = None, extent = [xmin,xmax,ymin,ymax], vmin = vmin, vmax = vmax, cmap = cmap)
plt.title(inputfilename + ' ' + vcol + ' Variogram Map'); plt.xlabel('X Offset (mm)'); plt.ylabel('Y Offset (mm)')

plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
plt.tight_layout()
#plt.savefig(f'{wd}{inputfilename}_VariogramMap{extension_png}')
plt.show()
#exec(f'pd.DataFrame(vmap).to_csv("{wd}{inputfilename}_VariogramMap{extension_csv}")')
duration()
print(f'The shape of the output is {vmap.shape}')

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
    global hscat, c
    #main = np.zeros((0,4)).astype(float32)
    hscat = np.zeros((permutations,3)).astype(float32) # Columns equivalent to ["direction", "lag", "tail", "head"]
    # Load the data
    # Trim values outside tmin and tmax
    df_extract = df.loc[(df[vcol] >= tmin) & (df[vcol] <= tmax)]
    nd = len(df_extract)  # TODO: not used
    x = df_extract[xcol].values; #x = x[:100]
    y = df_extract[ycol].values; #y = y[:100]
    vr = df_extract[vcol].values; #vr = vr[:100]

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
    #exec(f'hscat{int(iazi)} = hscat') # to do this we need to have a way to make each hscat global through a loop. not sure how to do that atm but at least we can save it into a csv below
    print(f'Saving hscat{int(azi_mat[iazi])} into a .csv it might take a while...')
    exec(f'pd.DataFrame(hscat, columns = ["lagend", "vrt", "vrh"]).to_csv("hscat{int(azi_mat[iazi])}.csv")') # turn np array into pd.DF with headers and save it in a .csv
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
    global hscat, c
    arr_stk = np.vstack((x,y,vr)).transpose();
    #arr_stk[:100,:]
    
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
    nd = shape(arr_stk)[0]*shape(arr_stk)[1]
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
    
    def sub_variogram_loop(i, j, tail, head):
        global hscat, c
        lagend, vrt, vrh = None, None, None
        if i%1000 == 0 and j == 0:                                  # Track the progress
            print('Direction: ' + str(iazi+1) + '/' + str(len(azi_mat)) + '.  Variogram Loop 1/2: ' + str(i)+"/" + str(shape(arr_stk)[0]) + ' ' + ctime())
            
        # if j%10_000 == 0:
        #     print(str(j+1)+"/"+str(nd+1))
        
        
        # Definition of the lag corresponding to the current pair
        dx = head[0] - tail[0]
        dy = head[1] - tail[1]
        dxs = dx * dx
        dys = dy * dy
        hs = dxs + dys
        if hs <= dismxs:
            if hs < 0.0:
                hs = 0.0
            h = np.sqrt(hs)

            # Determine which lag this is and skip if outside the defined
            # distance tolerance
            if h <= EPSLON:
                lagbeg = 0
                lagend = 0
            else:
                lagbeg = -1
                lagend = -1
                for ilag in range(1, nlag + 1):             # Figures out which lag this pair corresponds to
                    # reduced to -1
                    if (
                        (xlag * float(ilag - 1) - xltol)
                        <= h
                        <= (xlag * float(ilag - 1) + xltol)
                    ):
                        if lagbeg < 0:
                            lagbeg = ilag
                        lagend = ilag
            if lagend >= 0:
                # Definition of the direction corresponding to the current
                # pair. All directions are considered (overlapping of
                # direction tolerance cones is allowed)

                # Check for an acceptable azimuth angle
                dxy = np.sqrt(max((dxs + dys), 0.0))
                if dxy < EPSLON:
                    dcazm = 1.0
                else:
                    dcazm = (dx * uvxazm + dy * uvyazm) / dxy

                # Check the horizontal bandwidth criteria (maximum deviation
                # perpendicular to the specified direction azimuth)
                band = uvxazm * dy - uvyazm * dx

                # Apply all the previous checks at once to avoid a lot of
                # nested if statements
                if (abs(dcazm) >= csatol) and (abs(band) <= bandwh): # Check if dcazm is below angle tolerance and band is below bandwidth
                    # Check whether or not an omni-directional variogram is
                    # being computed
                    omni = False
                    if atol >= 90.0:
                        omni = True

                    # For this variogram, sort out which is the tail and
                    # the head value
                    iv = 0  # hardcoded just one variogram
                    it = ivtype[iv]  # TODO: not used
                    if dcazm >= 0.0:
                        vrh = tail[2]
                        vrt = head[2]
                        if omni:
                            vrtpr = tail[2]
                            vrhpr = head[2]
                    else:
                        vrh = head[2]
                        vrt = tail[2]
                        if omni:
                            vrtpr = head[2]
                            vrhpr = tail[2]

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
                    n = i*arr_stk.shape[0] + j # 
                    hscat[n,0] = float32(lagend); hscat[n,1] = float32(vrt); hscat[n,2] = float32(vrh)
                    #hscat = np.vstack([hscat, [float32(azm), float32(lagend), float32(vrt), float32(vrh)]]) # inserts new entry with the data pair at the end of the np.array
        return lagend, vrt, vrh, hscat

    # Main loop over all pairs
    c = np.array( [sub_variogram_loop(i, j, tail, head) for i, tail in enumerate(arr_stk) for j, head in enumerate(arr_stk)]) # list comprehensions quicker than nested for loops and numpy than pandas, apparently.
    #hscat = c.transpose()
    #main = np.vstack((main, c))
    hscat = hscat[~np.all(hscat == 0, axis=1)] # Because not all data pairs will go in the dataframe (because not all of them fall within the seach cone), this will clean all zero rows                    
                        
    # Get average values for gam, hm, tm, hv, and tv, then compute the correct "variogram" measure
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



azi_mat = [ 45.0, 135.0, 90.0]           # directions in azimuth to consider
tmin = -9999.; tmax = 9999.                             # no trimming. This needs to be set down here or it will mess up the vaiogram map above.
lag = np.zeros((len(azi_mat),int(nlag)+2)); gamma = np.zeros((len(azi_mat),int(nlag)+2)); npp = np.zeros((len(azi_mat),int(nlag+2)));
stop0 = time.time()

for iazi in range(0,len(azi_mat)):                      # Loop over all directions
    print('Semi-variogram direction ' + str(iazi+1)+"/"+str(len(azi_mat)) + ' ' + ctime())                                   # Track the progress   
    
    lag[iazi,:], gamma[iazi,:], npp[iazi,:] = gamv(df, xcol, ycol, vcol, tmin, tmax, lag_dist, lag_tol, nlag, azi_mat[iazi], atol,bandh, isill)
    exec(f'{inputfilename}_Res_SemiVariogram_Dir{int(azi_mat[iazi])} = pd.DataFrame(list(zip(lag[iazi,:],gamma[iazi,:])), columns=["lag","gamma"] )') # Createa pandas DataFrame with the semi-variogram values for that direction
    exec(f'{inputfilename}_Res_SemiVariogram_Dir{int(azi_mat[iazi])}.to_csv("{inputfilename}_Res_SemiVariogram_lag{dxlag}_Dir{int(azi_mat[iazi])}.csv")') # and save it in a .csv
    #plt.subplot(4,2,iazi+1)
    plt.plot(lag[iazi,:],gamma[iazi,:],'x',color = 'black', label = f'Azimuth {azi_mat[iazi]}')
    plt.plot([0,2000],[1.0,1.0],color = 'black')
    plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
    plt.ylabel(r'$\gamma \bf(h)$')
    plt.title(f'{inputfilename} {vcol} Semi-Variogram')
    plt.xlim([0, lag_dist*(nlag+1) ])
    plt.ylim([0, gamma[iazi,:].max()+0.25 ])
    plt.legend(loc='upper left')
    plt.grid(True)
    #plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=4.2, wspace=0.2, hspace=0.3)
    plt.savefig(f'{wd}{inputfilename}_Res_Semi-variograms_lag{dxlag}_{int(azi_mat[iazi])}deg{extension_png}',dpi=300, bbox_inches = "tight")
    plt.show()
    
    
    #hscat_extract = exec(f"hscat{int(azm)}.loc[(hscat{int(azm)}['direction'] == dire) & (hscat{int(azm)}['lag'] == lag)]")
    for lg in unique(hscat[:, 0]):
        if hscat[hscat[:,0]==lg, 1].shape[0]>0: # Sometimes there's no values, i.e. it's an empty arraay and causes an error which makes the loop stop
            plt.scatter( hscat[hscat[:,0]==lg, 1], hscat[hscat[:,0]==lg, 2], s=0.1, label = f"{vcol} lag: {int(lg)}" ) # Filter hscat by lag, plot vrt in x-axis [1], vrh in y-axis [2]
            r_squared_value = r2_score(hscat[hscat[:,0]==lg, 1], hscat[hscat[:,0]==lg, 2])
            plt.axline((0, 0), slope=1, color="black", linestyle='--', linewidth = 0.5, label = f"R^2={r_squared_value}")        
            plt.xlabel(f'Tail values ({vcol})'); plt.ylabel(f'Head values ({vcol})'); plt.title(f'h-scattergram {vcol} azimuth: {int(azi_mat[iazi])}, lag: {int(lg)}')
            plt.grid(True)
            plt.legend(loc = 1, fontsize='small') # Location Code:'best'0, 'upper right'1, 'upper left'2, 'lower left'3, 'lower right'4, 'right'5, 'center left'6, 'center right'7, 'lower center'8, 'upper center'9, 'center'10. fontsizeint or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
            plt.savefig(f'{wd}{inputfilename}_hscat_{vcol}_{int(azi_mat[iazi])}_lag_{int(lg)}{extension_png}')
            plt.show() 
        
    
# exec(f"dires = [int(l) for l in np.unique(hscat{int(azm)}[:,0:1], axis=1)]") # Creates an array of just column 0 (directions) and creates a list of all the unique values
# exec(f"lags = [int(l) for l in np.unique(hscat{int(azm)}[:,1:2], axis=1)]") # Creates an array of just column 1 (lags) and creates a list of all the unique values
# for dire in dires:
#     for lag in lags:
#         hscat_extract = exec(f"hscat{int(azm)}.loc[(hscat{int(azm)}['direction'] == dire) & (hscat{int(azm)}['lag'] == lag)]")
#         exec(f'hscat_extract.to_csv("hscat_extract{dire}.csv")') # and save it in a .csv
#         plt.scatter( hscat_extract['tail'], hscat_extract['head'], s=1 )
#         plt.axline((0, 0), slope=1, color="black", linestyle='--', linewidth = 0.5)
#         plt.xlabel('Tail values (Z)'); plt.ylabel('Head values (Z)'); plt.title(f'h-scattergram {vcol} azimuth: {dire}, lag: {int(lag+1)}')
#         plt.grid(True)
#         plt.show() 
#================================================================================================================================================================================ 
duration()
    

