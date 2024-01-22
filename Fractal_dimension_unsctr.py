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
import matplotlib
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
from time import time, ctime
import time
from scipy.spatial import KDTree
# from upscale import grid_fit

''' Set variables '''
# plt.rcParams['agg.path.chunksize'] = 10000                                                                                                              #set chunksize of matplotlib to higher
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals').as_posix() #flips the backslashes to plug in open()
bs = "\\"
wd=wd+bs                                                                                                                                                      #set a working directiory
inputfilename = pathlib.PureWindowsPath(r'z_ORIGvsORIG-rot(-0.02, 0.31, 0)_ApertureMap').as_posix()
os.chdir(wd) 
wdsave = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals\Upscale\Fractal').as_posix() #flips the backslashes to plug in open()
wdsave=wdsave+bs 

save = False
norm = False  
scale0switch = False  # Include the raw data (scale 0) in the plot
file_format = "long"          
method = "Nayak_Mishra" # Options are
# Db                                                                      # sum of all intensity differences between highest and lowest intensity in a box over each scale (epsilon)
# Dm                                                                     # Mean of all intensity differences between highest and lowest intensity in a box over each scale (epsilon)
# Davg                                                                   # average between highest and lowest intensity in a box
# Gangepain_RoquesCarmes
# Sarkar_Chaudhuri
# Hausdorf_Besicovitch
# Nayak_Mishra

# min_box_size = 0.17340891795287267 * 1.5                                        # minimum box size should be based on the average minimum distance (use grid_fit() method to calculate). However if you simply use avg_min_dist,the chances of points falling within the bins are small, thus using at least 1.5*avg_min_dist is recommended
n_samples = 20                                                                                                                                                                                         #uses the normalised dataset given the max/min parameters
vcol="ApertureR"


'''methods'''
def nsmall(a, n):
    return np.partition(a, n)[n]                                                #finds the n-th smallest value
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
        dist_df.loc[i,'dist_to_closest_pt'] = d[1]                                                                     # Absign the value of the distance between each point in the loop and it's closest point (apart from itself)
    avg_min_dist = dist_df['dist_to_closest_pt'].mean()
    print( "The average minimum distance between points is: " + str(avg_min_dist))  
    return avg_min_dist, dist_df
'''Read .csv straight to np array'''
with open (wd+inputfilename+extension_csv, 'r') as f:
    d_array = np.genfromtxt(f,delimiter=',')                                    #data here is in a numpy array. 
for i in range(d_array.shape[1]-1, 2, -1):
    d_array = np.delete(d_array, i, 1)                                          #Keep only the first 3 columns (x,y,z supposedly)
df = pd.read_csv(wd + inputfilename + extension_csv, index_col=0, header=0)
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])
if norm == True:
    df[f'N{vcol}'], tvz, tnsz = nscore(df, f'{vcol}')
    if df[f"N{vcol}"].min() < 0:
        df[f"N{vcol}"] = df[f"N{vcol}"] + abs(df[f"N{vcol}"].min())
    vcol = f"N{vcol}"

''' Basic Statistics'''
xtotal = d_array.shape[0]     # number of rows(number of points in long format)
ytotal = d_array.shape[1]     #number of columns.
if file_format == "long":
    total = xtotal 
else:
    total = xtotal * ytotal #Number of points

''''determine the scales to measure on'''
min_box_size, _ = grid_fit(df) ; min_box_size = min_box_size* 1.5
xdiff=max(df['x']) - min(df['x']); ydiff=max(df['y']) - min(df['y']); max_box_size = int(np.floor(np.log2(min(xdiff/2, ydiff/2))))   
Am = xdiff*ydiff
                   
scales = np.floor(np.logspace(max_box_size, min_box_size, num = n_samples, base =2 ))    #default max size is the largest power of 2 that fits in the smallest dimension of the array. Note that this
                                                                                        #function has been edited from a true 3D which takes np.min(array.shape) to nsmall(grid.shape, 2) which finds the second smallest array shape. This is because our dataset
                                                                                        #is a z-slice of a 3D (#x, #y. 1). Hence, if min is used, max_box_size will return 0.
scales = np.unique(scales).astype('int32')                                              #remove duplicates that could occur as a result of the floor
scales = scales[ scales < min(xdiff, ydiff)/2]  


Ns = []
# =============================================================================
# a=[]    #uncomment for QC
# a2=[]   #uncomment for QC
# a3=[]   #uncomment for QC
# =============================================================================
for i, scale in enumerate(scales):
    print(f"Scale Loop: {i+1}/{len(scales)}")
    bin_edges = [np.arange(start, end, scale) for start in (int(min(df['x'])),int(min(df['y']))) for end in (int(max(df['x'])), int(max(df['y']))) ]; bin_edges = bin_edges[0::3]                    #creates a np.array with the edges of the grid based on the scale. bin edges start at minimums x and y
    for x,x2 in zip(bin_edges[0][:], bin_edges[0][1:]):                            #depending on the bin_edges for each particular scale, grab the (x)-value and (x+1)-value
        for y,y2 in zip(bin_edges[1][:], bin_edges[1][1:]):                        #depending on the bin_edges for each particular scale, grab the (y)-value and (y+1)-value
            As = (x2-x)*(y2-y)# Area of s (s is the current box)
            r = As/Am
            if (x2-x == scale and y2-y == scale):                               #checks that the box is of the correct size
# =============================================================================
#                 # a.append([(x,y), (x,y2), (x2,y), (x2,y2)])  #Uncomment for QC
#                 # print(f'({x},{y}), ({x},{y2}), ({x2},{y}), ({x2},{y2})' )   #Uncomment for QC    
# =============================================================================

                box = df[(df["x"] > x) & (df["x"] < x2) & (df["y"] > y) & (df["y"] < y2)]     #filtering the df only to the points within the box
                if box.empty == False:
                    Intensity_range = max(box[vcol]) - min(box[vcol]) #sqrt( max(box[vcol])**2 - min(box[vcol])**2)                            #You can add a "1+" to prevent problems from being 0 values in later calculations with the logs. This is because we are inferring something about geometry from pixel intensity rather than directly measuring the binary state of the pixel.
                    Mean = box[vcol].mean()
                    Ns.append([scale, As, r, Intensity_range, Mean, box.shape[0]])
# =============================================================================
#                         # print(f'({i},{j})' )      #uncomment for QC
#                         # a2.append(grid[i][j][0])  #uncomment for QC
#                         # a3.append([i,j])          #uncomment for QC
# =============================================================================

if scale0switch == True:
    # scales = np.insert(scales, 0.1, 0, axis=0) # Insert scale 0 (corresponds to the initial raw data) in position 0 of the np array in axis=0
    df['Scale'] = 0.1 # Add an empty column for the corresponding scale in df. This is a very small number because the scale at which the raw data is compared to scale 1 (min_box_size*min_box_size Area) is infinitesimally small because it collapses on a point. This is important because the first point will impact the slope of the fitting line. EDIT: actually it doesn't seem to matter
    df['As'] = 1
    df['r'] = 1
    df['Mean'] = df[vcol]*1000 # because we want the values in meters
    df['Count'] = 1
    scale0 = list(map(list, zip(df['Scale'], df['As'], df['r'], df[vcol], df['Mean'], df['Count']))) # Create a list with similar structure to Ns with the scale 0 data for later merging
    Ns = scale0 + Ns
# Ns = pd.DataFrame( Ns, columns = ["Scale", "IntensityRange", "Mean", "Count"] )
# Ns2 = pd.concat( (pd.DataFrame(scale0) , Ns)) 
''' Calculating data for ploting'''
scale = list(zip(*Ns)); scale = list(scale[0])                                  #extract the scales, i.e. transpose a list of lists and extract the first row ( previously first column).
I_range = list(zip(*Ns)); I_range = list(I_range[3])                            #extract the intensity range, i.e. transpose a list of lists and extract the second row ( previously second column).
S = float(len(scale))                                                           #sets the maximum numbers to use in the plot
scale_vals = np.linspace(0, S, int(S)+1)                                            #creates a list of values linearly spaced between 0 and S with S number of values, i.e. the values of scale to use in the x axis (S should be in the units of distance, i.e. if unit of distance is m but you use mm, S should be in float)
i_scale_vals = np.searchsorted(scale, scale_vals)                               #BINARY SEARCH: Find the indices into a sorted array scale such that, if the corresponding elements in scale_vals were inserted before the indices, the order of  would be preserved. 
                                                                                #--> In other words, if you try to insert scale_vals into scale and preserve the order, i_scale_vals is the list containing the indices where each element in scale_vals would fall.
                                                                                #This i_scale_vals will basically be a map of the indexes in I_range where each scale (say scale==2 or scale==4) start and stop. That way you can sum all the values with scale==2 or scale==4 for instance without having to search the whole dataset.
Abs = list(zip(*Ns)); Abs = list(Abs[1])                                        # Area of the box for each scale
Abs = np.unique(Abs)
Fractal=[]
# =============================================================================
# # with open('C:/Users/s2132627/Documents/test234.txt', 'w') as f:             #uncomment for QC
# #     f.write(f"i, i_scale_val, N, start, stop")                              #uncomment for QC
# =============================================================================

#############################################################################  Indent this whole block if you wish to use the QC method above  ###########################################################
j=-1
for i, i_scale_val in enumerate(i_scale_vals[1:]-1):    # i starts at the 1st i_scale_vals item (not 0th)
    start, stop = i_scale_vals[i], i_scale_val        #start = index i-1 (first iteration will be 0) and stop = ihval (the index in zh, zsh and zz that represents the last item of a certain h)
    N = stop-start                          #amount of values to sum, so that you can calculate the average
    # f.write(f"\n {i}, {i_scale_val}, {N}, {start}, {stop}")             #uncomment for QC
    F = 0
    if N>0:
        j+=1
        if method == "Db":  # Box counting dimension for greyscale images: sum of all intensities of all boxes for a particular scale
            F =  sum (1 + intensity[0] for intensity in zip(I_range[start:stop]))  #calculate the sum of  all intensity ranges for a particular scale (box size), and add it to the Fractal list. EDIT: For some reason, F =( np.sum(float(intensity[0])) for intensity in zip(I_range[start:stop]) ) was not working. so I played around with the parenthesis.
            Fractal.append([ (scale[i_scale_val] ), F, 1/Abs[j]/Am, N])
        elif method == "Dm": # Box counting dimension for greyscale images: mean of all intensities of all boxes for a particular scale
            F =  sum (1 + intensity[0] for intensity in zip(I_range[start:stop])) / N    
            Fractal.append([ (scale[i_scale_val]), F, 1/Abs[j]/Am, N])
        elif method == "Davg":
            F =  sum(intensity[0] for intensity in zip(I_range[start:stop]))      
            Fractal.append([ (scale[i_scale_val])  / len(I_range[start:stop]), F, 1/Abs[j]/Am,  N])
        elif method == "Gangepain_RoquesCarmes":
            F =  sum(intensity[0] for intensity in zip(I_range[start:stop]))      
            Fractal.append([ (scale[i_scale_val]), F, 1/(Abs[j]*N/Am), N])
        elif method == "Sarkar_Chaudhuri":
            F =  sum(1 + intensity[0] for intensity in zip(I_range[start:stop]));
            Fractal.append([ scale[i_scale_val], F, 1/Abs[j]/Am, N]) # r = As/Am Abs[i] is the As for the current scale. As is the area of the box of that particular scale and Am is the Area of the surface/image. i-1 is because we start at 1 for other reasons.
        elif method == "Nayak_Mishra":
            for I in I_range[start:stop]:
                if I == 0:
                    f = 1
                else:
                    f = math.ceil((1 + I) / (N*Am/Abs[j]))
                F+=f
            Fractal.append([ scale[i_scale_val], F, 1/Abs[j]/Am, N]) # r = As/Am Abs[i] is the As for the current scale. As is the area of the box of that particular scale and Am is the Area of the surface/image. i-1 is because we start at 1 for other reasons.
        
        print(f"i: {i}, start: {start}, stop: {stop}, N: {N}, As: {Abs[j]}, F: {F}")
# if method == "Davg":
#     Fractal =  Fractal[:,1].mean()
############################################################################################################################################################################################################


Fractal = np.array(Fractal)                                                    #transforms result list into np.array
if 0 in Fractal[:, 0]: # If 0 exists in the scale column, it needs to be replaced with any number > 0 because the plot is log-log and log(0)=indefined
    Fractal[ np.where(Fractal[:, 0] == 0)[0][0], 0] = 1E-1 # Replace any scale = 0 with another value for plotting

Fractal_df  = pd.DataFrame(Fractal, columns=["Scale", "SdI", "1/r", "N"])
if save == True:
    exec(f"Fractal_df.to_csv('{wdsave}{inputfilename}_{vcol}_Fractal.csv')") # and save it in a .csv

''' Actual Ploting'''
plotting="1/r" # Scale , 1/r
if plotting == "Scale":
    xs = "Scale"
    ys = "SdI" # [ list(Fractal[:,0]), list(1/Fractal[:,1])]
elif plotting == "1/r":
    xs = "1/r"
    ys = "SdI"    #  [Fractal[:,2], Fractal[:,1]]
fig = plt.figure(); ax = fig.gca()                                              #create an empty figure. axes
fig = plt.loglog(Fractal_df[xs], Fractal_df[ys], 'ro', label = "Fractal")           #Create the fractal scatter plot in loglog 'b-' uses a blue line, 'ro' uses red circles
m, c = np.polyfit(log(Fractal_df[xs]), log(Fractal_df[ys]), 1)                      #Extract the slope and intersect of the data in loglog form - so that it looks like a straight line on a loglog plot. Degree of polynomial == 1, thus linear.
y_fit = np.exp(m*log(Fractal_df[xs]) + c)                                         # calculate the fitted values of y
fig = plt.loglog( Fractal_df[xs], y_fit, 'b-')                                    #plot the original x-values with the fitted y-values in the same loglog plot as the original data.
# =============================================================================
# # print('y = {:.2f}x + {:.2f}'.format(m,c))                                   #Uncomment for QC
# =============================================================================

#############################################################################  Plot format  ################################################################################################################
plt.title(f'{method} Fractal Dimension - GW1Q4 {vcol}')

if method == "Db":
    plt.ylabel(f"{plotting}", fontsize=10 ,color='black')    #r'$Scale (\varepsilon)(mm)$'       #visit https://matplotlib.org/3.1.0/tutorials/text/mathtext.html#symbols for more info on how to display mathematical formulas on the plot axis
    plt.xlabel(r'$I\varepsilon = \sum_{i,j=0}^{N_\varepsilon} 1+ \delta I_{i,j,\varepsilon}$', fontsize=10, color='black')
elif method == "Dm":
    plt.xlabel(f"{plotting}", fontsize=10 ,color='black')
    plt.ylabel(r'$I\varepsilon = \sum_{i,j=0}^{N_\varepsilon} 1+ \delta I_{i,j,\varepsilon}$', fontsize=10, color='black')
elif method == "Davg":
    plt.xlabel(f"{plotting}", fontsize=10 ,color='black')
    plt.ylabel(r'$I\varepsilon = \sum_{i,j=0}^{N_\varepsilon} \delta I_{i,j,\varepsilon}$', fontsize=10, color='black')
else:
    plt.xlabel(f"{plotting}", fontsize=10 ,color='black')
    plt.ylabel(r'$I\varepsilon = \sum_{i,j=0}^{N_\varepsilon} 1+ \delta I_{i,j,\varepsilon}$', fontsize=10, color='black')
ax.xaxis.grid(True, which='both'); ax.yaxis.grid(True, which='both')
ax.set_xticks([i for i in Fractal_df[xs]]);
ax.set_yticks([10**i for i in range(0, len(str(math.ceil(max(Fractal.transpose()[1]))))+1)]); ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()); ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
k = len(str(math.ceil(max(Fractal.transpose()[1]))));
plt.text( max(Fractal_df[xs]), 10**(k)-10**(k-1), 'y = {:.2f}x + {:.2f}             '.format(m, c), horizontalalignment="right", verticalalignment="top")
plt.grid(True)
if save == True:
    plt.savefig(f'{wdsave}\{inputfilename}Davg_Fractal{extension_png}',dpi=1000, bbox_inches = "tight")
############################################################################################################################################################################################################


        

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