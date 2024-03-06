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
# import geostatspy
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
#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename = pathlib.PureWindowsPath(r'Extracted M CNL_Aperture').as_posix()
os.chdir(wd)                                   # set the working directory

wdsave = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\aperture\upscaletest').as_posix() #flips the backslashes to plug in open()
bs="//"; wdsave=wdsave+bs  
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap').as_posix()
# os.chdir(wd)     

# =============================================================================

#set working directory and filename
#wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\GW1_Q4\Res').as_posix() #flips the backslashes to plug in open()
#bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
#inputfilename = pathlib.PureWindowsPath(r'GW1_Q4Res').as_posix()
#os.chdir(wd)                                   # set the working directory


vcol="aperture"
res = 'low'
harmonic = vcol #  options are vcol, std, avg, var
save = True
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2,3,4,5], index_col = 0)                  #This new line is because reading the .csv file with a header and an index column will create some issue
df = pd.read_csv(wd + inputfilename + extension_csv, index_col = 0, header = 0)                  
# ResVarMap = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                              # Read csv's of the variogram maps
# df = df.reset_index(drop=True) # If needed to reset the index

#df['Zzr'] = zscore(df['zr']) #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
# df['N'+vcol], tvz, tnsz = geostats.nscore(df, vcol) # nscore transform for all facies porosity #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
df['N'+vcol], tvz, tnsz = nscore(df, vcol) # if you run geostatspynscore.py from the Python Scripts folder

stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
vmin = stats['min'][vcol]; vmax = stats['max'][vcol]; nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]

plt.hist(df[vcol], 50, histtype = "bar", cumulative = False, density=True, facecolor='g', ec='black', alpha=0.75)
mu = df[vcol].mean()
std = df[vcol].std()
plt.axvline(mu, c= 'k', ls='-', label = "$\mu = $"+f"{mu:2f}")
for i in [j for j in range(-3,4,1) if j != 0]:
    plt.axvline(mu+i*std, c= 'b', ls='--', label = f"{i}"+"$\sigma$")
plt.legend(loc = "upper left")
plt.title(f"{inputfilename}")
plt.show()

def grid_fit(dataframe):
    """Calculate the average minimum distance of a 2D dataset of points using K-D Tree.
    :param dataframe: 
    :return: avg_min_dist, average minimum distance
    :return: df2, dataframe with the distances in order of the points given in the input dataframe
    """
    nd = len(dataframe)                                                                                                #  number of total points
    # Find the nearest point from a given point to a large list of points. Useful for finding minimum box size.    
    dist_df = dataframe[['x', 'y']].copy(); ar = dist_df.to_numpy(); dist_df['dist_to_closest_pt'] = 0;                            # Create a copy of the xy values of all points. Memory allocation. Initialise needed list to make the tree
    tree = KDTree(ar)                                                                                           # Create the K-D Tree used for a (very!) quick computation without the need for nested loops (how does this even work so fast??)
    ####  Average distance between points ###
    for i in range(0, nd):
        if i % 1_000 == 0:                                                                                      # Track the progress
            print('Calculating average distance between points: ' + str(i) + "/" + str(nd) + ' ' + ctime())
        pt = (dataframe['x'][i], dataframe['y'][i])                                                                           # records the point's coordinates
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
# min_box_size = 0.17340891795287267
min_box_size = grid_fit(df)[0] * 1.5                                            # minimum box size should be based on the average minimum distance (use grid_fit() method to calculate). However if you simply use avg_min_dist,the chances of points falling within the bins are small, thus using at least 1.5*avg_min_dist is recommended
n_samples = 20                                                                                                                       
norm = False                                                                    #uses the normalised dataset given the max/min parameters


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
    mainpd = pd.DataFrame(columns = ['scale','x','y', vcol, 'avg', 'var', 'std'])
    for count, scale in enumerate(scales):
        bin_edges = [np.arange(start, end, scale) for start in (int(min(dataframe['x'])),int(min(dataframe['y']))) for end in (int(max(dataframe['x'])), int(max(dataframe['y']))) ]; bin_edges = bin_edges[0::3]                    #creates a np.array with the edges of the grid based on the scale. bin edges start at minimums x and y
        #array = np.zeros(((len(bin_edges[0])-1), (len(bin_edges[1])-1) ))  # Not being used: I was testing on creating a grid plot but opted with the scatter
        df_c = dataframe[['x','y', vcol]].copy(); df_c[vcol][:] = np.nan     # Creates a copy of df and fills the 'z' column with NaNs. This is to create a upscaled surface with high resolution, i.e., all the points within the box get the same value corresponding to the z value of the centre of the surface fitted to the box
        exec(f"upscale{scale} = pd.DataFrame( columns= ['scale', 'x','y', vcol, 'avg', 'var', 'std'])")     # Creates the blank DataFrame containing the upscaled points. This will be a low resolution upscaled surface, i.e., each box will render only the point at the centre of the fitted surface, as opposed to replacing the z values of all the points within the box
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
                        surf_fit(box, vcol)   # Fit the surface to the box's points
                        bcx = x+(x2-x)/2; bcy = y+(y2-y)/2; # Find the box's centre in x and y coordinates: Box Centre x and Box Centre y
                        boxstd = box[vcol].std()
                        boxavg = box[vcol].mean()
                        boxvar = box[vcol].var()
                        exec(f"upscale{scale} = pd.concat( [upscale{scale}, pd.DataFrame( [[scale, bcx, bcy, C[0]*bcx + C[1]*bcy + C[2], boxavg, boxvar, boxstd ]], columns= ['scale', 'x','y', vcol, 'avg', 'var', 'std']) ], ignore_index = True )")   # Append line (point) of the centre of the box to the DataFrame for that scale.
                        #array[i][j] = (C[0]*bcx + C[1]*bcy + C[2])     # Not being used: I was testing on creating a grid plot but opted with the scatter
                        df_c.loc[(df_c['x'] > x) & (df_c['x'] < x2) & (df_c['y'] > y) & (df_c['y'] < y2), vcol] = (C[0]*bcx + C[1]*bcy + C[2]) # Make all the values within the box the same z as the centre point of the fitted surface
        exec(f"mainpd = pd.concat(   [mainpd, upscale{scale}.fillna({scale})]   , ignore_index = True )")
        
        if res == "low":
            exec(f"plt.scatter(upscale{scale}['x'], upscale{scale}['y'], c = upscale{scale}[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the low resolution upscaled. alpha is transparency
            plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(df_c[vcol]), amax(df_c[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
            plt.xlabel(r'$X (mm)$', fontsize=15)
            plt.ylabel(r'$Y (mm)$', fontsize=15)
            plt.title(f'{inputfilename}\nUpscaled{scale} {vcol} Fracture map')
            plt.tight_layout()
            if save == True:
                plt.savefig(f'{wdsave}{inputfilename}_LR_Upscaled{scale}_{vcol}_FracMap{extension_png}', bbox_inches = "tight")
                exec(f'upscale{scale}.to_csv("{wdsave}{inputfilename}_LR_Upscaled{scale}_{vcol}_FractureMap{extension_csv}")') # and save it in a .csv
            plt.show()
            
# =============================================================================
#         elif res == "high":
#             exec(f"plt.scatter(df_c['x'], df_c['y'], c = df_c[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the high resolution upscaled. alpha is transparency
#             plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(dataframe[vcol]), amax(dataframe[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
#             plt.xlabel(r'$X (mm)$', fontsize=15)
#             plt.ylabel(r'$Y (mm)$', fontsize=15)
#             plt.title(f'{inputfilename} HR Upscaled{scale} {vcol} Fracture map')
#             plt.tight_layout()
#             if save == True:
#                 plt.savefig(f'{wdsave}{inputfilename}_HR_Upscaled{scale}_{vcol}_FracMap{extension_png}', bbox_inches = "tight")
#             plt.show()
#         
#         elif res == "both":
#             exec(f"plt.scatter(df_c['x'], df_c['y'], c = df_c[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the high resolution upscaled. alpha is transparency
#             plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(dataframe[vcol]), amax(dataframe[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
#             plt.xlabel(r'$X (mm)$', fontsize=15)
#             plt.ylabel(r'$Y (mm)$', fontsize=15)
#             plt.title(f'{inputfilename} HR Upscaled{scale} {vcol} Fracture map')
#             plt.tight_layout()
#             if save == True:
#                 plt.savefig(f'{wdsave}{inputfilename}_HR_Upscaled{scale}_{vcol}_FracMap{extension_png}', bbox_inches = "tight")
#             plt.show()
#             exec(f"plt.scatter(upscale{scale}['x'], upscale{scale}['y'], c = upscale{scale}[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the low resolution upscaled. alpha is transparency
#             plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(dataframe[vcol]), amax(dataframe[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
#             plt.xlabel(r'$X (mm)$', fontsize=15)
#             plt.ylabel(r'$Y (mm)$', fontsize=15)
#             plt.title(f'{inputfilename}\nUpscaled{scale} {vcol} Fracture map')
#             plt.tight_layout()
#             if save == True:
#                 plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_FractureMap{extension_png}', bbox_inches = "tight")
#             plt.show()
# =============================================================================
        
    if save == True:
        ################################################################################################
        # For some reason this has to be run manually to work or it will save an empty .csv file #######
        ################################################################################################
        exec("mainpd.to_csv(f'{wdsave}{inputfilename}_upscaled_{vcol}_Fracs{extension_csv}')")  # 
        
    globals().update(locals())  # Make all local variables global so that we get the data outside of the function and can work with it (i.e., all the upscale{scale} DataFrames)


upscale(df, vcol)

for count, scale in enumerate(scales):
    exec(f"mainpd = pd.concat( [mainpd, upscale{scale}.fillna({scale})] , ignore_index = True)")
#plt harmonic            
if harmonic == vcol:
    plt.scatter(mainpd['scale'], mainpd[vcol], s=0.5)
    plt.xlabel(r'$Scale$', fontsize=15)
    plt.ylabel(f'${vcol} (mm)$', fontsize=15)
    plt.title(f'{inputfilename}\n{vcol} Harmonic')
    plt.tight_layout()
    if save == True:
        plt.savefig(f'{wdsave}{inputfilename}_{vcol}_Harmonic{extension_png}')
    plt.show()
    for i, scale in enumerate(scales):
        a = mainpd.loc[mainpd["scale"] == scale, vcol]
        #plt.subplot(3,5,i+1)                                       # histograms for each scale
        plt.hist(a, facecolor='red', bins=np.linspace(a.min(),a.max(), 300), histtype="stepfilled", edgecolor='black', alpha=0.4,density=True,cumulative=False, label=f'Histogram {vcol}')
        plt.xlabel(f'{vcol} (mm)'); plt.ylabel('Frequency (%)'); plt.title(f'Histogram upscaled {vcol} at scale {scale}')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.tight_layout()
        if save == True:
            plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_Hist{extension_png}', bbox_inches = "tight")
            exec(f'a.to_csv("{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_Hist{extension_csv}")') # and save it in a .csv
        plt.show()
elif harmonic == "std":
    plt.scatter(mainpd['scale'], mainpd[harmonic], s=0.5)
    plt.xlabel(r'$Scale$', fontsize=15)
    plt.ylabel(f'$Standard Deviation$', fontsize=15)
    plt.title(f'{inputfilename} {vcol} Std Harmonic')
    plt.tight_layout()
    if save == True:
        plt.savefig(f'{wd}{inputfilename}_{vcol}_Std_Harmonic{extension_png}', bbox_inches = "tight")
    plt.show()
    for i, scale in enumerate(scales):
        a = mainpd.loc[mainpd["scale"] == scale, harmonic]
        #plt.subplot(3,5,i+1)                                       # plot original and normalised histograms
        plt.hist(a, facecolor='red',bins=np.linspace(a.min(),a.max(),300), alpha=0.2,density=True,cumulative=False,edgecolor='black', label=f'Histogram {harmonic}')
        plt.xlabel(f'{vcol} Standard Deviation'); plt.ylabel('Frequency (%)'); plt.title(f'Histogram upscaled {harmonic} at scale {scale}')
        plt.legend(loc='upper right')
        plt.grid(True)
        if save == True:
            plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_Hist{extension_png}', bbox_inches = "tight")
            exec(f'a.to_csv("{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_Hist{extension_csv}")') # and save it in a .csv
        if save == True:
            plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_png}', bbox_inches = "tight")
            exec(f'a.to_csv("{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_csv}")') # and save it in a .csv
        plt.show()
elif harmonic == "avg":        
    plt.scatter(mainpd['scale'], mainpd['avg'], s=0.5)
    plt.xlabel(r'$Scale$', fontsize=15)
    plt.ylabel(f'$Average$', fontsize=15)
    plt.title(f'{inputfilename} {vcol} Avg Harmonic')
    plt.tight_layout()
    if save == True:
        plt.savefig(f'{wdsave}{inputfilename}_{vcol}_Avg_Harmonic{extension_png}', bbox_inches = "tight")
    plt.show()
    for i, scale in enumerate(scales):
        a = mainpd.loc[mainpd["scale"] == scale, harmonic]
        #plt.subplot(3,5,i+1)                                       # plot original and normalised histograms
        plt.hist(a, facecolor='red',bins=np.linspace(a.min(),a.max(),300), histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black', label=f'Histogram {harmonic}')
        plt.xlabel(f'{vcol} (mm)'); plt.ylabel('Frequency (%)'); plt.title(f'Histogram upscaled {harmonic} at scale {scale}')
        plt.legend(loc='upper right')
        plt.grid(True)
        if save == True:
            plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_png}', bbox_inches = "tight")
            exec(f'a.to_csv("{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_csv}")') # and save it in a .csv
        plt.show()
elif harmonic == "var":        
    plt.scatter(mainpd['scale'], mainpd['var'], s=0.5)
    plt.xlabel(r'$Scale$', fontsize=15)
    plt.ylabel(f'$Variance$', fontsize=15)
    plt.title(f'{inputfilename} {vcol} Var Harmonic')
    plt.tight_layout()
    if save == True:
        plt.savefig(f'{wdsave}{inputfilename}_{vcol}_Var_Harmonic{extension_png}', bbox_inches = "tight")
    plt.show()
    for i, scale in enumerate(scales):
        a = mainpd.loc[mainpd["scale"] == scale, harmonic]
        #plt.subplot(3,5,i+1)                                       # plot original and normalised histograms
        plt.hist(a, facecolor='red',bins=np.linspace(a.min(),a.max(),300), histtype="stepfilled", alpha=0.2,density=True,cumulative=False,edgecolor='black', label=f'Histogram {harmonic}')
        plt.xlabel(f'{vcol} Variance'); plt.ylabel('Frequency (%)'); plt.title(f'Histogram upscaled {harmonic} at scale {scale}')
        plt.legend(loc='upper right')
        plt.grid(True)
        if save == True:
            plt.savefig(f'{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_png}', bbox_inches = "tight")
            exec(f'a.to_csv("{wdsave}{inputfilename}_Upscaled{scale}_{vcol}_{harmonic}_Hist{extension_csv}")') # and save it in a .csv
        plt.show()
# =============================================================================
#                         # print(f'({i},{j})' )      #uncomment for QC
#                         # a2.append(grid[i][j][0])  #uncomment for QC
#                         # a3.append([i,j])          #uncomment for QC
# =============================================================================





# =============================================================================
# def upscale(dataframe, vcol):
#     global df_c
#     df_c = df.copy();
#     #df_c[vcol][:] = np.nan
#     exec(f"plt.scatter(df_c['x'], df_c['y'], c = df_c[vcol], cmap=cmap, marker = '.', alpha=1)") # Creates the scatter plot of the high resolution upscaled. alpha is transparency
#     plt.colorbar(mappable = None, label = 'Z Value', orientation="vertical", ticks=np.linspace(amin(dataframe['z']), amax(dataframe[vcol]), 10))                                                                                      #contour map, colour bar and colourbar label???
#     plt.xlabel(r'$X (mm)$', fontsize=15)
#     plt.ylabel(r'$Y (mm)$', fontsize=15)
#     plt.title(f'{inputfilename} HR Upscaled{scale} Residuals Fracture map')
#     plt.tight_layout()
#     plt.savefig(f'{wd}{inputfilename}_HR_Upscaled{scale}_ResFractureMap{extension_png}')
#     plt.show()
# scale = 1
# upscale(df, 'zr')
# =============================================================================

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