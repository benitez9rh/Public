# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 17:20:34 2021

@author: s2132627
"""

import os
import pathlib
from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
import math
from time import time, ctime
import time
from scipy.spatial import KDTree
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"

#set working directory and filename
# wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean').as_posix()

# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans').as_posix() #flips the backslashes to plug in open()
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean_Q4').as_posix()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows

# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\synfrac_Test_64x64').as_posix() #flips the backslashes to plug in open()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows
# inputfilename = pathlib.PureWindowsPath(r'synfrac_Test_64x64.top').as_posix()
# os.chdir(wd)                                   # set the working directory

def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f' days: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')

'''###########################################################################################################################################################
############################################################## IMPORT AND NORMALISE DATAFRAME ################################################################
###########################################################################################################################################################'''
# df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[0,1,2], names=['x', 'y', 'z'])                     # read a .csv file in as a DataFrame. The input file needs to be in a table format, i.e. columns of 'x', 'y', 'z' with #rows = #points. As is, it only reads columns 0,1,2 and ignores the rest.
# df['Nz'] = zscore(df['z']) #  This creates a Z-score of the line(from scipy.stats)    Or in numpy: line_Z= pd.DataFrame(np.abs(stats.zscore(line))).  https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba
# stats=df.describe().transpose(); xmin = stats['min']['x']; xmax = stats['max']['x']; ymin = stats['min']['y']; ymax = stats['max']['y'];
# vmin = stats['min']['z']; vmax = stats['max']['z']; nvmin = stats['min']['Nz']; nvmax = stats['max']['Nz']

'''###########################################################################################################################################################
############################################################ CALCULATE AVERAGE MINIMUM DISTANCE ##############################################################
###########################################################################################################################################################'''
stop0 = time.time()
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
            print(f"Calculating average distance between points: {i} / {nd} {ctime()}")
        pt = (dataframe['x'][i], dataframe['y'][i])                                                                           # records the point's coordinates
        d, j = tree.query(pt, k=[2])                                                                            # d = distance, i = index. k = 2 to get the first and second nearest points because the first nearest point will be the point itself, hence dist = 0.
        dist_df.loc[i,'dist_to_closest_pt'] = d                                                                     # Assign the value of the distance between each point in the loop and it's closest point (apart from itself)
    avg_min_dist = dist_df['dist_to_closest_pt'].mean()
    min_dist = dist_df['dist_to_closest_pt'].min()
    max_dist = dist_df['dist_to_closest_pt'].max()
    print( f"The average minimum distance between closest points is {avg_min_dist}, the minimum distance between closest points is {min_dist} and the maximum distance between closest points  is {max_dist}.") 
    return avg_min_dist, min_dist, max_dist, dist_df
avg_min_dist, min_dist, max_dist, dist_df = grid_fit(df)
duration()