# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 15:06:31 2021

@author: s2132627
"""
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:38:16 2021

@author: s2132627
"""



''' Import Libraries'''
import os
import pathlib
from numpy import *
import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import csv
import pandas as pd
from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
import geostatspy
from time import time, ctime
import time
from scipy.spatial import KDTree
from matplotlib import cm

''' Set variables '''
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
cmap = plt.cm.plasma                    # color map
# =============================================================================
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean_Q4').as_posix()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows
# =============================================================================

# =============================================================================
# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q3Aperture').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'GW1Q3_z_ORIGvsORIG-r(-0.02, 0.31, 0)-offset(0, 0, -0.0658)_ApertureMap').as_posix()
# os.chdir(wd)          
# =============================================================================

# set the working directory
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap').as_posix()
os.chdir(wd)  
                      
# =============================================================================

#set working directory and filename
#wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Res').as_posix() #flips the backslashes to plug in open()
#bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
#inputfilename = pathlib.PureWindowsPath(r'GW1_Q4Res').as_posix()
#os.chdir(wd)                                   # set the working directory
# =============================================================================
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1\DipDirPlot').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# =============================================================================

vcol="aperture"
res = 'low'
# harmonic = "DipDir" #  options are vcol, Dip, DipDir 
E = 30
A = -30

# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2,3,4,5], index_col = 0)                  #This new line is because reading the .csv file with a header and an index column will create some issue
df = pd.read_csv(wd + inputfilename + extension_csv, index_col = 0)                  
# ResVarMap = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                              # Read csv's of the variogram maps

wd = pathlib.PureWindowsPath(wd+r'\DipDirPlot').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs 

#df['Zzr'] = zscore(df['zr']) #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
#df['Nzr'], tvz, tnsz = geostats.nscore(df, 'zr') # nscore transform for all facies porosity #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol];

def planedipdir0(C): # http://www.tjscientific.com/2017/08/16/python-script-to-calculate-strike-and-dip/
    a,b,c = C[0], C[1], C[2]
    ptA, ptB, ptC = (0, 0, 0*a + 0*b + c), (0, 10, 0*a + 10*b + c), (10, 10, 10*a + 10*b + c)
    x1, y1, z1 = float(ptA[0]), float(ptA[1]), float(ptA[2])
    x2, y2, z2 = float(ptB[0]), float(ptB[1]), float(ptB[2])
    x3, y3, z3 = float(ptC[0]), float(ptC[1]), float(ptC[2])
    
    u1 = float(((y1-y2)*(z3-z2)-(y3-y2)*(z1-z2)))
    u2 = float((-((x1-x2)*(z3-z2)-(x3-x2)*(z1-z2))))
    u3 = float(((x1-x2)*(y3-y2)-(x3-x2)*(y1-y2)))
    
    if u3 < 0:
        easting = u2
    else:
        easting = -u2
    
    if u3 > 0:
        northing = u1
    else:
        northing = -u1
    
    if easting >= 0:
        partA_strike = math.pow(easting, 2) + math.pow(northing, 2)
        strike = math.degrees(math.acos(northing / math.sqrt(partA_strike)))
    else:
        partA_strike = northing / math.sqrt(math.pow(easting, 2) + math.pow(northing, 2))
        strike = math.degrees(2 * math.pi - math.acos(partA_strike))
       
    part1_dip = math.sqrt(math.pow(u2, 2) + math.pow(u1, 2))
    part2_dip = math.sqrt(math.pow(u1,2) + math.pow(u2,2) + math.pow(u3,2))
    dip = math.degrees(math.asin(part1_dip / part2_dip))

    return dip, strike


def planedipdir(C, deg = True): # Based on Nicholson, W. K. (2006). Linear Algebra with Applications. Toronto, ON: McGraw-Hill Ryerson edited from http://www.tjscientific.com/2017/08/16/python-script-to-calculate-strike-and-dip/
    a,b,c = C[0], C[1], C[2]
    ptA, ptB, ptC = (0, 0, 0*a + 0*b + c), (0, 10, 0*a + 10*b + c), (10, 10, 10*a + 10*b + c)
    x1, y1, z1 = float(ptA[0]), float(ptA[1]), float(ptA[2])
    x2, y2, z2 = float(ptB[0]), float(ptB[1]), float(ptB[2])
    x3, y3, z3 = float(ptC[0]), float(ptC[1]), float(ptC[2])
    
    # cross product of AB and AC gives us the coordinates u1, u2, u3 fof the normal vector to ABC plane
    u1 = float(((y2-y1)*(z3-z1)-(y3-y1)*(z2-z1))) 
    u2 = float((-((x2-x1)*(z3-z1)-(x3-x1)*(z2-z1)))) # notice u2 has a negative sign
    u3 = float(((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)))
    
    # u3 corresponds to z-coordinate of the normal vector of the plane. 
    # In geology, the convention is that dip direction sits to the right of the strike angle measured from north (or positive y-axis).
    # If positive, dip direction points to the right of the strike, otherwise to the left. Therefore, if the plane is turned upside down, i.e. u3<0, the northing must be made symmetrical. If the plane is upside up, easting remains negative (because East is counted towards the positive side of the x-axis).
    if u3 < 0:
        easting = u2 # i.e. remains negative
        northing = -1*u1
    else:
        easting = -1*u2 # made positive
        northing = u1
    
    if easting >= 0:
        strike = math.degrees(math.acos(northing / math.sqrt(  math.pow(easting, 2) + math.pow(northing, 2)  ))) # i.e. arccos (N / sqrt(E^2+N^2))
    else:
        strike = math.degrees(2 * math.pi - math.acos(  northing / math.sqrt(math.pow(easting, 2) + math.pow(northing, 2))  ))

    dip = math.degrees(math.asin( (math.sqrt(math.pow(u2, 2) + math.pow(u1, 2))) / (math.sqrt(math.pow(u1,2) + math.pow(u2,2) + math.pow(u3,2))) )) # i.e. arcsin ( sqrt(u1^2+u2^2) / sqrt(u1^2+u2^2+u3^2) )
    dip_direction = strike + 90 #  Because dip direction, as per the convention, is to the right of the strike direction, add 90 degrees to get it
    if dip_direction > 360: dip_direction = dip_direction-360
    
    if deg == False:
        dip = math.radians(dip)
        dip_direction = math.radians(dip_direction)
        strike = math.radians(strike)
    # print("dip, dip direction, strike")
    return dip, dip_direction, strike

def grid_fit(dataframe):
    """Calculate the average minimum distance of a 2D dataset of points using K-D Tree.
    :param dataframe: 
    :return: avg_min_dist, average minimum distance
    :return: df2, dataframe with the distances in order of the points given in the input dataframe
    """
    nd = len(df)                                                                                                #  number of total points
    # Find the nearest point from a given point to a large list of points. Useful for finding minimum box size.    
    dist_df = df[['x', 'y']].copy(); ar = dist_df.to_numpy(); dist_df['dist_to_closest_pt'] = 0;                            # Create a copy of the xy values of all points. Memory allocation. Initialise needed list to make the tree
    tree = KDTree(ar)                                                                                           # Create the K-D Tree used for a (very!) quick computation without the need for nested loops (how does this even work so fast??)
    ####  Average distance between points ###
    for i in range(0, nd):
        if i % 1_000 == 0:                                                                                      # Track the progress
            print('Calculating average distance between points: ' + str(i) + "/" + str(nd) + ' ' + ctime())
        pt = (df['x'][i], df['y'][i])                                                                           # records the point's coordinates
        d, j = tree.query(pt, k=2)                                                                            # d = distance, i = index. k = 2 to get the first and second nearest points because the first nearest point will be the point itself, hence dist = 0.
        dist_df.loc[i,'dist_to_closest_pt'] = d[1]                                                                     # Assign the value of the distance between each point in the loop and it's closest point (apart from itself)
    avg_min_dist = dist_df['dist_to_closest_pt'].mean()
    print( "The average minimum distance between points is: " + str(avg_min_dist))  
    return avg_min_dist, dist_df

file_format = "long"
# file_format = "wide"            
Db = True                                                                       #difference between highest and lowest intensity in a box
Dm = False
Davg = False                                                                    #average between highest and lowest intensity in a box
min_box_size = grid_fit(df)[0] * 1.5                                            # minimum box size should be based on the average minimum distance (use grid_fit() method to calculate). However if you simply use avg_min_dist,the chances of points falling within the bins are small, thus using at least 1.5*avg_min_dist is recommended
n_samples = 20                                                                                                                       
norm = False                                                                    #uses the normalised dataset given the max/min parameters
# =============================================================================
# 
# def normvecctr(C, l = 1, xctr = 0, yctr = 0, zctr = 0):
#     """
#     Normal vector to a plane with centre at (xctr, yctr, zctr) and length l.
#     Parameters
#     ----------
#     C : List, tuple
#         List or tuple containing the linear (total of 3) plane equation constants. (Ax + By + Cz = D)
#     l : 
#         Length of unit-vector. Standard = 1
#     xctr  : x-coordinate of vector's tail ( or central point of the plane). Standard value is x = 0.
#     yctr  : y-coordinate of vector's tail ( or central point of the plane). Standard value is y = 0.
#     zctr  : x-coordinate of vector's tail ( or central point of the plane). Standard value is z = 0.
#     
#     Returns
#     -------
#     nv :    np.array()
#             Vector normal to plane.
#     unv :   np.array()
#             Univt Normal Vector to a plane.
# 
#     """    
#     
#     i = C[0]+xctr
#     j = C[1]+yctr
#     k = C[2]+zctr
#     
#     nv = np.array( [i, j, k] ) #  a(x-x0)+b(y-y0)+c(z-z0)=0
#     
#     unv = (nv/(sqrt(nv[0]**2 + nv[1]**2 + nv[2]**2)))*l
#     # unv = (l/abs(nv))*nv
#     
#     return nv, unv
# 
# def vang(n, axis):
#     """
#     Calculates the angle between a vector n and an axis x, y or z. If angle is returned negative, n is pointing up into the positive Z realm, i.e. above horizontal plane XY.
#     
#     Parameters
#     ----------
#     n : np.array or list.
#         np.array or list containg 3 coordinates xyz.
#     axis : String
#         x, y or z
# 
#     Returns
#     -------
#     rads : float
#         angle in radians.
#     degs : float
#         angle in degrees.
# 
#     """
#     xaxis = np.array([1, 0, 0])
#     yaxis = np.array([0, 1, 0])
#     zaxis = np.array([0, 0, -1]) # Because geometric bodies dip towards -Z
#     
#     if axis == "x":
#         v = xaxis
#         rads = np.arccos(np.dot(n,v) / (( sqrt(n[0]**2 + n[1]**2 + n[2]**2) )*( sqrt(v[0]**2 + v[1]**2 + v[2]**2) )))
#     if axis == "y":
#         v = yaxis
#         rads = np.arccos(np.dot(n,v) / (( sqrt(n[0]**2 + n[1]**2 + n[2]**2) )*( sqrt(v[0]**2 + v[1]**2 + v[2]**2) )))
#     if axis == "z":
#         v = zaxis
#         rads = math.pi/180 - np.arccos(np.dot(n,v) / (( sqrt(n[0]**2 + n[1]**2 + n[2]**2) )*( sqrt(v[0]**2 + v[1]**2 + v[2]**2) ))) 
#     
#     
#     degs = rads*(180/math.pi)
#     
#     return rads, degs
# 
# def planedipdir0(C, degrees = True):
#     """
#     Calculates dip and dip directions of a plane from the plane's coefficients (3)
# 
#     Parameters
#     ----------
#     C : tuple or np.array with 3 
#         tuple or np.array with 3 coefficients for the plane's definition formula.
#     
#     degrees: Boolean
#              Flag to return results in degrees (True) or radians (False)
#     
#     Returns
#     -------
#     dip : float
#         Angle between the vector and the negative z-axis.
#     dire : float
#         Angle between the vector and the positive y-axis.
# 
#     """
#     if degrees == True:
#         _,dip = vang(normvecctr(C)[0], "z")
#         if dip < 0:
#             dip = -1*dip
#         elif dip > 90:
#             dip = 90 - dip
#             
#         _,dire = vang(normvecctr(C)[0], "y")
#     else:
#         dip,_ = vang(normvecctr(C)[0], "z")
#         if dip < 0:
#             dip = -1*dip
#         elif dip > 90:
#             dip = 90 - dip
#         
#         dire,_ = vang(normvecctr(C)[0], "y")
#     
#     return dip, dire
# =============================================================================

def surf_fit(dataframe, vcol):
    global C, mn, mx
    points = list(map(list, zip(dataframe['x'], dataframe['y'], dataframe[vcol])))     # Creates a list of lists (instead of a list of tuples if you only use list(zip(df...)) ) with (x ,y ,z) coordinates of the initial data, i.e. by matching the x, y and z of the first row in the first tuple and then adding it to the list and moving on to the next (unzipping).
    data = np.array(points)
    order = 1   # 1: linear, 2: quadratic, 3: cubic
    mn = floor(np.min(data, axis=0)) #mn = np.min(data, axis=0)                                       ##########################################################
    mx = ceil(np.max(data, axis=0)) #mx = np.max(data, axis=0)                                        # This section creates a linearly (regularly) spaced mesh#
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 100), np.linspace(mn[1], mx[1], 100))                 # in which the surface is drawn on                       #
    XX = X.flatten(); YY = Y.flatten()                                                                
    # X = data[:,0]; Y = data[:,1];                                                                       # Use this section instead simply uses the x and y values of the actual data set as the "mesh" for the surface which will allow a much easier computation of the residuals.
    # XX = X.flatten(); YY = Y.flatten()                                                                  #
    
    ''' Create the matrix and get coefficients depending on the desired polynomial order'''  
    if order == 1:                                                                                      # best-fit linear plane
        A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])                                                      # coefficients
        Z = C[0]*X + C[1]*Y + C[2]                                                                      # evaluate it on grid. Or expressed using matrix/vector product: # Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)
        
    elif order == 2:    # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]       # Creates an array A which has (1st) column full of ones, (2nd and 3rd) first 2 columns of data (x and y), (4th) is product of x*y, column is x**2, (5th) column is the product of x*y,  
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])  
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)           # evaluate it on a grid
    
    elif order == 3:    # Best fit cubic curve
        A = np.c_[np.ones(data.shape[0]), data[:,:2], data[:,0]**2, np.prod(data[:,:2], axis=1), \
                  data[:,1]**2, data[:,0]**3, np.prod(np.c_[data[:,0]**2,data[:,1]],axis=1), \
                  np.prod(np.c_[data[:,0],data[:,1]**2],axis=1), data[:,2]**3]                          # M = [ones(size(x)), x, y, x.^2, x.*y, y.^2, x.^3, x.^2.*y, x.*y.^2, y.^3]     Creates an array A which has (1st) column full of ones, (2nd and 3rd) first 2 columns of data (x and y), (4th) column is x**2, (5th) column is the product of x*y, (6th) is y**2, (7th) x**3, (8th) product of (x**2) * y (9th) product of x * (y**2), (10th) y**3
        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])  
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX**2, XX*YY, YY**2, XX**3, XX**2*YY, XX*YY**2, \
                         YY**3], C).reshape(X.shape)                                                    # evaluate it on a grid
    return C

def upscale(dataframe, vcol):
    '''Read .csv straight to np array'''
# =============================================================================
#     with open (wd+inputfilename+extension_csv, 'r') as f:
#         d_array = np.genfromtxt(f,delimiter=',')                                    #data here is in a numpy array. 
#     for i in range(d_array.shape[1]-1, 2, -1):
#         d_array = np.delete(d_array, i, 1)                                          #Keep only the first 3 columns (x,y,z supposedly)
# =============================================================================
    #exec(open(r"C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Python Scripts\temp_functions.py").read())

    
    ''' Basic Statistics'''
    xtotal = dataframe.shape[0]     # number of rows(number of points in long format)
    ytotal = dataframe.shape[1]     #number of columns.
    if file_format == "long":
        total = xtotal 
    else:
        total = xtotal * ytotal #Number of points
    
    ''''determine the scales to measure on'''
    xdiff=max(dataframe['x']) - min(dataframe['x']); ydiff=max(dataframe['y']) - min(dataframe['y']); max_box_size = int(np.floor(np.log2(min(xdiff/2, ydiff/2))))   
                       
    scales = np.floor(np.logspace(max_box_size,min_box_size, num = n_samples, base =2 ))    #default max size is the largest power of 2 that fits in the smallest dimension of the array. Note that this
                                                                                            #function has been edited from a true 3D which takes np.min(array.shape) to nsmall(grid.shape, 2) which finds the second smallest array shape. This is because our dataset
                                                                                            #is a z-slice of a 3D (#x, #y. 1). Hence, if min is used, max_box_size will return 0.
    scales = np.unique(scales).astype('int32')                                              #remove duplicates that could occur as a result of the floor
    scales = scales[ scales < min(xdiff, ydiff)/2]  
    
    Ns=[]
    # =============================================================================
    # a=[]    #uncomment for QC
    # a2=[]   #uncomment for QC
    # a3=[]   #uncomment for QC
    # =============================================================================
    mainpd = pd.DataFrame(columns = ['scale','xctr','yctr', vcol, 'dip', 'dipdir'])
    for count, scale in enumerate(scales):
        bin_edges = [np.arange(start, end, scale) for start in (int(min(dataframe['x'])),int(min(dataframe['y']))) for end in (int(max(dataframe['x'])), int(max(dataframe['y']))) ]; bin_edges = bin_edges[0::3]                    #creates a np.array with the edges of the grid based on the scale. bin edges start at minimums x and y
        #array = np.zeros(((len(bin_edges[0])-1), (len(bin_edges[1])-1) ))  # Not being used: I was testing on creating a grid plot but opted with the scatter
        df_c = df[['x','y', vcol]].copy(); df_c[vcol][:] = np.nan     # Creates a copy of df and fills the 'z' column with NaNs. This is to create a upscaled surface with high resolution, i.e., all the points within the box get the same value corresponding to the z value of the centre of the surface fitted to the box
        exec(f"upscale{scale} = pd.DataFrame( columns= ['x','y', vcol, 'dip', 'dipdir'])")     # Creates the blank DataFrame containing the upscaled points. This will be a low resolution upscaled surface, i.e., each box will render only the point at the centre of the fitted surface, as opposed to replacing the z values of all the points within the box
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        
        for i, (x, x2) in enumerate(zip(bin_edges[0][:], bin_edges[0][1:])):                            #depending on the bin_edges for each particular scale, grab the (x)-value and (x+1)-value
            print(f"Scale: {count+1}/{len(scales)}. Box: {i+1}/{(len(bin_edges[0])-1)}")                       # Reporting
            for j, (y, y2) in enumerate(zip(bin_edges[1][:], bin_edges[1][1:])):                        #depending on the bin_edges for each particular scale, grab the (y)-value and (y+1)-value
                if (x2-x == scale and y2-y == scale):                               #checks that the box is of the correct size
    # =============================================================================
    #                 # a.append([(x,y), (x,y2), (x2,y), (x2,y2)])  #Uncomment for QC
    #                 # print(f'({x},{y}), ({x},{y2}), ({x2},{y}), ({x2},{y2})' )   #Uncomment for QC    
    # =============================================================================
            
                    box = dataframe[(dataframe["x"] > x) & (dataframe["x"] < x2) & (dataframe["y"] > y) & (dataframe["y"] < y2)]     #filtering the dataframe only to the points within the box
                    if box.empty == False:
                        C = surf_fit(box, vcol)   # Fit the surface to the box's points
                        dip, dipdir, strike = planedipdir(C)
                        bcx = x+(x2-x)/2; bcy = y+(y2-y)/2; bcz = bcx*C[0] + bcy*C[1] + C[2] # Find the box's centre in x and y coordinates: Box Centre x and Box Centre y
                        
                        start = [bcx, bcy, bcz]          # The z is calculated from the plane function so that it lies on the plane rather than getting the average z of the data itself which most likely won't fall on the plane.
                        X = linspace(x, x2, 10)          # linspace doesn't need to be very dense because it is just for display purposes.
                        Y = linspace(y, y2, 10)
                        # print(f"{scale}, C: {C}, bcx: {bcx}, bcy: {bcy}, bcz: {bcz}")
                        XX, YY = np.meshgrid(X, Y)
                        Z = C[0]*XX + C[1]*YY + C[2] # Z = A*x + B*y + D

                        surf = ax.plot_surface(XX, YY, Z, cmap = cm.coolwarm)                      # plot the scattered data. c is colour and s is size , cmap = cm.coolwarm

                        # nv, unv = normvecctr(C, l = 10, xctr = bcx, yctr = bcy, zctr = bcz)
                        # end = [start[0] + unv[0], start[1] + unv[1], start[2] + unv[2]]
                        # plt.quiver(start[0], start[1], start[2], end[0], end[1], end[2], color = 'green')
                        
                        exec(f"upscale{scale} = upscale{scale}.append( pd.DataFrame( [[bcx, bcy, bcz, dip, dipdir ]], columns= ['x','y', vcol, 'dip', 'dipdir']) )")   # Append line (point) of the centre of the box to the DataFrame for that scale.
                        #df_c.loc[(df_c['x'] > x) & (df_c['x'] < x2) & (df_c['y'] > y) & (df_c['y'] < y2), vcol] = (C[0]*bcx + C[1]*bcy + C[2]) # Make all the values within the box the same z as the centre point of the fitted surface
        plt.xlabel('X'); plt.ylabel('Y'); ax.set_zlabel('Z')                                            #set labels
        plt.title(f'{inputfilename}\n{vcol} Scale{scale} Dip/DipDirection', wrap = True)
        ax.axis('auto'); ax.axis('tight')
        ax.view_init(elev=E, azim=A)                                                                    # rotate view depending on Elevation and Azimuth
        plt.savefig(f'{wd}{inputfilename}_Scale{scale}_Trends{extension_png}', bbox_inches = 'tight')
        plt.show()
              
        exec(f"mainpd = mainpd.append(upscale{scale}).fillna({scale})")

        
        if res == "low":
            exec(f"plt.scatter(upscale{scale}['x'], upscale{scale}['y'], c = upscale{scale}[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the low resolution upscaled. alpha is transparency
            plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(dataframe[vcol]), amax(dataframe[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
            plt.xlabel(r'$X (mm)$', fontsize=15)
            plt.ylabel(r'$Y (mm)$', fontsize=15)
            plt.title(f'{inputfilename} LR Upscaled{scale} {vcol} Fracture map', wrap = True)
            plt.tight_layout()
            plt.savefig(f'{wd}{inputfilename}_LR_Upscaled{scale}_{vcol}_FracMap{extension_png}', bbox_inches = 'tight')
            plt.show()
            #saves upscaled to csv
            exec(f'upscale{scale}.to_csv("{wd}{inputfilename}_LR_Upscaled{scale}_{vcol}_FractureMap{extension_csv}")') # and save it in a .csv
        
        
    exec(f'mainpd.to_csv("{wd}{inputfilename}_upscaled_{vcol}_Fracs{extension_csv}")')
    globals().update(locals())  # Make all local variables global so that we get the data outside of the function and can work with it (i.e., all the upscale{scale} DataFrames)


upscale(df, vcol)

for count, scale in enumerate(scales):
    exec(f"mainpd = mainpd.append(upscale{scale}).fillna({scale})")
#plt Dip
plt.scatter(mainpd['scale'], mainpd["dip"], s=0.5)
plt.xlabel(r'$Scale$', fontsize=15)
plt.ylabel(f'$Dip (degrees)$', fontsize=15)
plt.title(f'{inputfilename} {vcol} Dip Spread', wrap = True)
plt.tight_layout()
plt.savefig(f'{wd}{inputfilename}_{vcol}_DipSpread{extension_png}', bbox_inches = 'tight')
plt.show()
mainpd[['scale','dip']].to_csv(f"{wd}{inputfilename}_{vcol}_DipSpread{extension_csv}")
for i, scale in enumerate(scales):
    a = mainpd.loc[mainpd["scale"] == scale, "dip"]
    #plt.subplot(3,5,i+1)                                       # plot original and normalised histograms
    plt.hist(a, facecolor='red',bins=np.linspace(a.min(),a.max(),300), alpha=0.2,density=True,cumulative=False,edgecolor='black', label=f'Histogram Dip')
    plt.xlabel(f'{vcol} Dip'); plt.ylabel('Frequency (%)'); plt.title(f'{inputfilename} Histogram upscaled Dip at scale {scale}', wrap = True)
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(f'{wd}{inputfilename}_Upscaled{scale}_Dip_{vcol}_Hist{extension_png}', bbox_inches = 'tight')
    plt.show()
    exec(f'a.to_csv("{wd}{inputfilename}_Upscaled{scale}_Dip_{vcol}_Hist{extension_csv}")') # and save it in a .csv
#plt DipDirection    
plt.scatter(mainpd['scale'], mainpd["dipdir"], s=0.5)
plt.xlabel(r'$Scale$', fontsize=15)
plt.ylabel(f'$DipDirection (degrees)$', fontsize=15)
plt.title(f'{inputfilename} {vcol} DipDirection Spread', wrap = True)
plt.tight_layout()
plt.savefig(f'{wd}{inputfilename}_{vcol}_DipDirSpread{extension_png}')
plt.show()
mainpd[['scale','dipdir']].to_csv(f"{wd}{inputfilename}_{vcol}_DipDirSpread{extension_csv}")
for i, scale in enumerate(scales):
    a = mainpd.loc[mainpd["scale"] == scale, "dipdir"]
    #plt.subplot(3,5,i+1)                                       # plot original and normalised histograms
    plt.hist(a, facecolor='red',bins=np.linspace(a.min(),a.max(),300), histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black', label=f'Histogram DipDirection')
    plt.xlabel(f'{vcol} DipDirection (degrees)'); plt.ylabel('Frequency (%)'); plt.title(f'{inputfilename} Histogram upscaled DipDirection at scale {scale}', wrap = True)
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig(f'{wd}{inputfilename}_Upscaled{scale}_DipDir_{vcol}_Hist{extension_png}', bbox_inches = 'tight')
    plt.show()
    exec(f'a.to_csv("{wd}{inputfilename}_Upscaled{scale}_DipDir_{vcol}_Hist{extension_csv}")') # and save it in a .csv

# =============================================================================
#                         # print(f'({i},{j})' )      #uncomment for QC
#                         # a2.append(grid[i][j][0])  #uncomment for QC
#                         # a3.append([i,j])          #uncomment for QC
# =============================================================================

def plotter_trends(E, A):                                                                                  # E=elevation A=azimuth for the viewing angle
    fig = plt.figure()
    ax = plt.axes(projection='3d')                                                                   #gca stands for "get current axis"
    for Q in Qs:
        exec(f"global Q_xctr; Q_xctr = {Q}['x'].mean()")                  # Translating s1 upscaled points to XY area of s2 by finding s2 centre of gravity and add/subtract centre of gravity XY coordinates to the s1 upscaled points
        exec(f"global Q_yctr; Q_yctr = {Q}['y'].mean()")
        exec(f"global Q_zctr; Q_zctr = {Q}[vcol].mean()")
        exec(f"global C; C = surf_fit({Q}, vcol)")
        start = [Q_xctr, Q_yctr, C[0]*Q_xctr + C[1]*Q_yctr + C[2]]          # The z is calculated from the plane function so that it lies on the plane rather than getting the average z of the data itself which most likely won't fall on the plane.
        exec(f"X = linspace({Q}['x'].min(), {Q}['x'].max(), 100)")          # linspace doesn't need to be very dense because it is just for display purposes.
        exec(f"Y = linspace({Q}['y'].min(), {Q}['y'].max(), 100)")
        print(f"{Q}, C: {C}, Q_xctr: {Q_xctr}, Q_yctr: {Q_yctr}, Q_zctr: {Q_zctr}")
        XX, YY = np.meshgrid(X, Y)
        Z = C[0]*XX + C[1]*YY + C[2] # Z = A*x + B*y + D

        surf = ax.plot_surface(XX, YY, Z, cmap = cm.coolwarm)                      # plot the scattered data. c is colour and s is size , cmap = cm.coolwarm

        nv, unv = normvecctr(C, l = 10, xctr = Q_xctr, yctr = Q_yctr, zctr = Q_zctr)
        end = [start[0] + unv[0], start[1] + unv[1], start[2] + unv[2]]
        plt.quiver(start[0], start[1], start[2], end[0], end[1], end[2], color = 'green')
    plt.xlabel('X'); plt.ylabel('Y'); ax.set_zlabel('Z')                                            #set labels
    plt.title('GW1 Trends', wrap = True)
    ax.axis('auto'); ax.axis('tight')
    ax.view_init(elev=E, azim=A)                                                                    # rotate view depending on Elevation and Azimuth
    # plt.savefig(f'{basewd}Qs_Trends_3D{extension_png}', bbox_inches = 'tight')
    plt.show()

plotter_trends(0, -30)



'''
###  Recognition  ###

https://imagej.nih.gov/ij/plugins/fraclac/FLHelp/Glossary.htm#graydb
Citation: Karperien, A., FracLac for ImageJ. http://rsb.info.nih.gov/ij/plugins/fraclac/FLHelp/Introduction.htm. 1999-2013. 
'''



'''
#####################################################################################################################################
#                                                                                                                                   #
#   Code above based on code below from https://github.com/ChatzigeorgiouGroup/FractalDimension/blob/master/FractalDimension.py     #
#                                                                                                                                   #
#####################################################################################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 09:47:15 2019
@author: daniel
"""
import numpy as np
import matplotlib.pyplot as plt

def fractal_dimension(array, max_box_size = None, min_box_size = 1, n_samples = 20, n_offsets = 0, plot = False):
    """Calculates the fractal dimension of a 3D numpy array.
    
    Args:
        array (np.ndarray): The array to calculate the fractal dimension of.
        max_box_size (int): The largest box size, given as the power of 2 so that
                            2**max_box_size gives the sidelength of the largest box.                     
        min_box_size (int): The smallest box size, given as the power of 2 so that
                            2**min_box_size gives the sidelength of the smallest box.
                            Default value 1.
        n_samples (int): number of scales to measure over.
        n_offsets (int): number of offsets to search over to find the smallest set N(s) to
                       cover  all voxels>0.
        plot (bool): set to true to see the analytical plot of a calculation.
                            
        
    """
    #determine the scales to measure on
    if max_box_size == None:
        #default max size is the largest power of 2 that fits in the smallest dimension of the array:
        max_box_size = int(np.floor(np.log2(np.min(array.shape))))
    scales = np.floor(np.logspace(max_box_size,min_box_size, num = n_samples, base =2 ))
    scales = np.unique(scales) #remove duplicates that could occur as a result of the floor
    
    #get the locations of all non-zero pixels
    locs = np.where(array > 0)
    voxels = np.array([(x,y,z) for x,y,z in zip(*locs)])
    
    #count the minimum amount of boxes touched
    Ns = []
    #loop over all scales
    for scale in scales:
        touched = []
        if n_offsets == 0:
            offsets = [0]
        else:
            offsets = np.linspace(0, scale, n_offsets)
        #search over all offsets
        for offset in offsets:
            bin_edges = [np.arange(0, i, scale) for i in array.shape]
            bin_edges = [np.hstack([0-offset,x + offset]) for x in bin_edges]
            H1, e = np.histogramdd(voxels, bins = bin_edges)
            touched.append(np.sum(H1>0))
        Ns.append(touched)
    Ns = np.array(Ns)
    
    #From all sets N found, keep the smallest one at each scale
    Ns = Ns.min(axis=1)
   
    
    
    #Only keep scales at which Ns changed
    scales  = np.array([np.min(scales[Ns == x]) for x in np.unique(Ns)])
    
    
    Ns = np.unique(Ns)
    Ns = Ns[Ns > 0]
    scales = scales[:len(Ns)]
    #perform fit
    coeffs = np.polyfit(np.log(1/scales), np.log(Ns),1)
    
    #make plot
    if plot:
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(np.log(1/scales), np.log(np.unique(Ns)), c = "teal", label = "Measured ratios")
        ax.set_ylabel("$\log N(\epsilon)$")
        ax.set_xlabel("$\log 1/ \epsilon$")
        fitted_y_vals = np.polyval(coeffs, np.log(1/scales))
        ax.plot(np.log(1/scales), fitted_y_vals, "k--", label = f"Fit: {np.round(coeffs[0],3)}X+{coeffs[1]}")
        ax.legend();
    return(coeffs[0])
'''