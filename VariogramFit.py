# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:55:17 2022

@author: s2132627
"""


'''###########################################################################################################################################################
########################################################################## VARIOGRAM FIT##########################################################################'''

import os
import pathlib
from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
from numpy import *
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import colors
from matplotlib import cm
import math
# import geostatspy.GSLIB as GSLIB                                  # Geostatspy is always giving trouble importing and impoting numba so I simply copied the
# import geostatspy.geostats as geostats                            # functions I needed directly to this script
import time
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
pltfont = {'fontname':'Arial'}

# =============================================================================
# #set working directory and filename
# wd = pathlib.PureWindowsPath(r'/home/s2132627/Documents/Step 1 - Variograms/Greywacke scans\GW1_Q4').as_posix() #flips the backslashes to plug in open()
# bs="//"; wd=wd+bs                              # Use this instead in linux
# inputfilename = pathlib.PureWindowsPath(r'Greywacke1_matched_clean_Q4').as_posix()
# os.chdir(wd)                                   # set the working directory
# =============================================================================

#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture').as_posix() #flips the backslashes to plug in open()
bs="\\"; wd=wd+bs                               # Use this instead in Windows                             # Use this instead in linux
inputfilename = pathlib.PureWindowsPath(r'lag1_SemiVariogram_Dir45').as_posix()
wdsave = wd
os.chdir(wd)                                   # set the working directory

save = False
lagcol="lag"
varcol="gamma"
# function = "Spherical"
exp_var = pd.read_csv(wd + inputfilename + extension_csv, index_col = 0)



def Find_DF_Row_Index_Between_Which_Value_Falls_In_Column(df, col, value):
    ranges = list(df[col])
    ranges.insert(0, 0) # assume nug = 0
    for b, a in zip(ranges[:len(ranges)-1], ranges[1:]):
        if value >= b and value < a:
    # for b, a in zip(df[col][:len(df[col])-1], df[col][1:]):
    #     if value >= b and value < a:
    #         index = df[col].loc[lambda x: x==b].index[0]
            index = df[col].loc[lambda x: x==b].index[0]
    return index
def reverse_enumerate(iterable, start="end"):
    if start == "end":
        n = len(iterable) - 1
    else:
        n = start
    for elem in iterable:
        yield n, elem
        n -= 1
def insert_position(position, list1, list2):
    return list1[:position] + list2 + list1[position:]

if exp_var.values[-1][0] == 0.0:                                                        #if last row has lag 0 (duplicate)
    exp_var.drop(exp_var.tail(1).index,inplace=True)      # drop last n rows
x = np.array(pd.concat([exp_var[lagcol], pd.DataFrame(arange(0,exp_var[lagcol].max(),0.1))]).sort_values(by=[0]))  # creates a 1-column DF of the lags in the input DF and merges it with a range from 0 to input-DF maximum ever 0.1 then sorts it.
#x = x[np.argsort(x)]                                                                               # Sorting the lag values (x-values)
y = np.zeros((x.shape[0]), dtype=x.dtype).transpose()                                               # adds a column with zeros

def semivar(dataframe, function, params, lagcol="lag", varcol="gamma"): 
    nug, r, isill = params
    if dataframe.values[-1][0] == 0.0:                                                        #if last row has lag 0 (duplicate)
        dataframe.drop(exp_var.tail(1).index,inplace=True)      # drop last n rows
    x = np.array(pd.concat([dataframe[lagcol], pd.DataFrame(arange(0,dataframe[lagcol].max(),0.1))]).sort_values(by=[0]))  # creates a 1-column DF of the lags in the input DF and merges it with a range from 0 to input-DF maximum ever 0.1 then sorts it.
    #x = x[np.argsort(x)]                                                                               # Sorting the lag values (x-values)
    y = np.zeros((x.shape[0]), dtype=x.dtype).transpose()                                               # adds a column with zeros
    if function == "Spherical":
        eq=r'$\gamma(h) = C(0) + 1.5*h + 0.5*h^r$'
        for i in range(len(x)):
            if x[i] == 0:
                y[i] = 0
            elif x[i] > 0 and x[i] <= r:
                y[i] = nug + (isill - nug) * ((3*x[i])/(2*r) - 0.5*(x[i]/r)**3)
            elif x[i] > r:
                y[i] = isill
    elif function == "Power":
        g=0.01
        beta=2
        for i in range(len(x)):
            if x[i] <= r:
                y[i] = nug + g*(x[i]**beta)
                if y[i] > isill:
                    y[i] = isill
            else:
                y[i] = isill
    elif function == "Exponential":
        alpha = 1
        OnW_FP = 1.1   # Oliver & Webster fitting parameter
        for i in range(len(x)):
            if x[i] <= r:
                if OnW_FP: # if defined
                    y[i] =  nug + (isill + nug) * (alpha - math.e**(-x[i]/ (OnW_FP * r) )) # Oliver & Webster 2015
                else:
                    y[i] =  nug + (isill + nug) * (alpha - math.e**(-x[i]/r)) # Oliver & Webster 2015
                # y[i] =   math.e**(-r/x[i]) # Chiles and Delfiner 2012 
            else:
                y[i] = isill
    elif function == "Gaussian":
        for i in range(len(x)):
            if x[i] <= r:
                y[i] =  nug + (isill + nug) * (1 - (math.e)**(-(x[i]**2)/r**2))
            else:
                y[i] = isill
    elif function == "Cauchy":
        beta = 0.5 # 3/2
        for i in range(len(x)):
            if x[i] <= r:
                y[i] =  1 + ( (r**2)/(x[i]**2) )**(-beta/2)
            else:
                y[i] = isill
    elif function == "Nested Spherical":                        # Oliver & Webster 2015
        eq=r'$\gamma(h) = C(0) + 1.5*h + 0.5*h^r$'
        for i in range(len(x)):
            if x[i] == 0:
                y[i] = 0
            elif x[i] > 0 and x[i] <= int_r_1:
                y[i] = nug + psill_1*((3/2)*(x[i]/int_r_1) - (1/2)*(x[i]/int_r_1)**3) + psill_2*((3/2)*(x[i]/r) - (1/2)*(x[i]/r)**3)
            elif x[i] > int_r_1 and x[i] <= r:
                    y[i] = nug + psill_1 + psill_2 * ((3/2)*(x[i]/r) - (1/2)*(x[i]/r)**3)
            elif x[i] > r:
                y[i] = isill
    return x, y


def My_Custom_Function(params, dists):
    """
    Parameters
    ----------
    dataframe : Pandas dataframe
        Contains the experimental variogram with the lags and the corresponding semi-variances columns.
    params : List of Lists
        Provide a list of lists containing, in this order, the nugget (only!), the variogram type, variogram partial sill and variogram range.
        The format is params = [[nugget], ["Variogram 1 type",  variogram 1 partial sill, variogram 1 range], ["Variogram 2 type",  variogram 2 partial sill, variogram 2 range], etc...]
                             
        Example:
            params = [
                [0],                                # Nugget
                ["Spherical", 0.3, 2.5],            # Spherical model uses only partial sill and range
                ["Spherical", 0.7, 14],             # Another speherical model
                ["Power", 0.8, 16, 0.01, 2],        # Power Model uses partial sill, range, Beta (represents the curvature) and g (intensity of varioation or gradient) parameters
                ["Exponential", 0.9, 18]    # Exponential Model uses partial sill and range.  For practical purposes it is usual to assign an effective range, aʹ, which is approximately equal to 3a (Oliver & Webster 2015).
                ]
    dists : Numpy 1-D array.
        DESCRIPTION. Sequence of points to simulate a continuous function. 
    
    Returns
    -------
    Numpy 1-D array Sequence of points to simulate a continuous function.

    """
    def modelvariogrampoint(lag, nug, vmodels, sills, ranges, *vargs):
        variance = 0
        for i, (r1, r2) in enumerate(zip(ranges[:len(ranges)-1], ranges[1:])):
            if lag > r1 and lag <= r2:
                for n in range(i+1):
                    variance += sills[n]
                # ####### Specific to each variogram model #######
                for n in range(i+1, len(ranges)):
                    if vmodels[i] == "Spherical":
                        # print(ranges, sills)
                        # print(n, r)
                        if r != 0:      # Divide by zero exception catch all
                            variance += sills[n] * ( 1.5*(lag/ranges[n]) - 0.5*(lag/ranges[n])**3 )
                    elif vmodel[i] == "Power":
                        beta=2
                        g=0.01
                        variance = nug + g * (lag**beta)
                    elif vmodels[i] == "Exponential":
                        alpha = 1
                        FP = 1.1   # Oliver & Webster fitting parameter
                        variance =  (nug + sills[i]) * (alpha - math.e**(-lag/ (FP * r2) )) # Oliver & Webster 2015           
                    elif vmodels[i] == "Gaussian":
                        variance =  nug + (sills[i] + nug) * (1 - (math.e)**(-(lag**2)/r2**2))
                    elif vmodels[i] == "Cauchy":
                        beta = 0.5 
                        variance =  1 + ( (r2**2)/(lag**2) )**(-beta/2)
                    else:
                        print("Variogram model not recognised.")
                    
            # ####### Common for all variogram models #######        
            elif lag <= 0:                            # if lag below or equal to 0
                variance = 0
            elif lag > r2 and i == len(ranges)-2:     # if lag beyond the last range
                variance = sum(sills)
                # print(f"i: {i},    r1: {r1},    r2: {r2},    lag: {lag},    vmodels: {vmodels[i]},   isills: {sills[i]}")
        # print(variance)
        return variance
    ###### Checks ######
    try:
        for i, v in enumerate(params):
            if i == 0:
                if type(v[0]).__name__ == "int" or type(v[0]).__name__ == "float":
                    nug = params[0][0]
                    print(f"Variogram model {i} successfully checked. Nugget with value {v[0    ]}.")
            else:
                try:
                    if len(v) >= 3:
                        print(f"Variogram model parameters are Complete.")
                except:
                    print(f"Variogram model parameters are incomplete. Only {len(params)} parameters provided: 3 minimum are required (nugget, sill and range) in list format for a simple Spherical or Exponential Model. Follow the format described below for other simp    le models\n \
                         For a Nested model, please provide a list of lists containing in the first list the nugget value and any subsequent lists the necessary parameters for that particular variogram model:\n \
                         Spherical(3) = nugget, sill, range,\
                         Exponential(3) = nugget, sill, range,\
                         Gaussian(3) = nugget, sill, range \
                        Power(5) = nugget, sill, range, Beta (represents the curvature)     and g (intensity of varioation or gradient) \
                        Cauchy(4) = nugget, sill, range, Beta \
                             \ ")
                if type(v[0]).__name__ == "str":
                    if v[0] == "Spherical":
                        if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                    elif v[0] == "Power":
                        if len(v) == 5 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Exponential":
                    if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Cauchy":
                    if len(v) == 4 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
    except:
        print(f"Variogram model {i} with format {v} does not conform with the required format. Either there are too many parameters, not enough parameters or a parameters is not of the required type.")
    finally:
    
        modelled_lags = dists
        modelled_gammas = np.zeros((modelled_lags.shape[0]), dtype=modelled_lags.dtype).transpose()                                               # adds a column with zeros
        if 'nug' in locals():
            params_f = params[1:]
        vmodels, sills, ranges, *vargs = list(map(list, zip(*params_f)))
        sills.insert(0, nug)
        ranges.insert(0, 0)
        if nug != 0.0:
            ranges.insert(0, 0)
        print(vmodels, sills, ranges, *vargs)
        if isinstance(modelled_lags, list) or isinstance(modelled_lags, np.ndarray): # Checks it dists are of list or np.array type
            for i, lag in enumerate(modelled_lags):
               modelled_gammas[i] = modelvariogrampoint(lag, nug, vmodels, sills, ranges, *vargs)
    return modelled_lags, modelled_gammas




def semivar_v2(dataframe, params, lagcol="lag", varcol="gamma"): 
    """
    Parameters
    ----------
    dataframe : Pandas dataframe
        Contains the experimental variogram with the lags and the corresponding semi-variances columns.
    params : List of Lists
        Provide a list of lists containing, in this order, the nugget (only!), the variogram type, variogram partial sill and variogram range.
        The format is params = [[nugget], ["Variogram 1 type",  variogram 1 partial sill, variogram 1 range], ["Variogram 2 type",  variogram 2 partial sill, variogram 2 range], etc...]
                             
        Example:
            params = [
                [0],                                # Nugget
                ["Spherical", 0.3, 2.5],            # Spherical model uses only partial sill and range
                ["Spherical", 0.7, 14],             # Another speherical model
                ["Power", 0.8, 16, 0.01, 2],        # Power Model uses partial sill, range, Beta (represents the curvature) and g (intensity of varioation or gradient) parameters
                ["Exponential", 0.9, 18]    # Exponential Model uses partial sill and range.  For practical purposes it is usual to assign an effective range, aʹ, which is approximately equal to 3a (Oliver & Webster 2015).
                ]
    dists : Numpy 1-D array.
        DESCRIPTION. Sequence of points to simulate a continuous function. 
    

    Returns
    -------
    Numpy 1-D array Sequence of points to simulate a continuous function.

    """
    ###### Checks ######
    if dataframe.values[-1][0] == 0.0:                                                        #if last row has lag 0 (duplicate)
        dataframe.drop(dataframe.tail(1).index,inplace=True)      # drop last n rows
    try:
        for i, v in enumerate(params):
            if i == 0:
                if type(v[0]).__name__ == "int" or type(v[0]).__name__ == "float":
                    nug = params[0][0]
                    print(f"Variogram model {i} successfully checked. Nugget with value {v[0    ]}.")
            else:
                try:
                    if len(v) >= 3:
                        print(f"Variogram model parameters are Complete.")
                except:
                    print(f"Variogram model parameters are incomplete. Only {len(params)} parameters provided: 3 minimum are required (nugget, sill and range) in list format for a simple Spherical or Exponential Model. Follow the format described below for other simp    le models\n \
                         For a Nested model, please provide a list of lists containing in the first list the nugget value and any subsequent lists the necessary parameters for that particular variogram model:\n \
                         Spherical(3) = nugget, sill, range,\
                         Exponential(3) = nugget, sill, range,\
                         Gaussian(3) = nugget, sill, range \
                        Power(5) = nugget, sill, range, Beta (represents the curvature)     and g (intensity of varioation or gradient) \
                        Cauchy(4) = nugget, sill, range, Beta \
                             \ ")
                if type(v[0]).__name__ == "str":
                    if v[0] == "Spherical":
                        if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                    elif v[0] == "Power":
                        if len(v) == 5 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Exponential":
                    if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Cauchy":
                    if len(v) == 4 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
    except:
        print(f"Variogram model {i} with format {v} does not conform with the required format. Either there are too many parameters, not enough parameters or a parameters is not of the required type.")
    finally:
        modelled_lags = np.sort( np.concatenate( (np.array(exp_var[lagcol]), arange(0,exp_var[lagcol].max(),0.1) ), axis=None) )
        #x = x[np.argsort(x)]                                                                               # Sorting the lag values (x-values)
        modelled_gammas = np.zeros((modelled_lags.shape[0]), dtype=modelled_lags.dtype).transpose()                                               # adds a column with zeros
        
        if 'nug' in locals():
            params_f = params[1:]
            # print(nug)
        vmodels, sills, ranges, *vargs = list(map(list, zip(*params_f)))
        sills.insert(0, nug)
        ranges.insert(0, 0)
        if nug != 0.0:
            ranges.insert(0, 0)
        print(vmodels, sills, ranges, *vargs)
                    
        if isinstance(modelled_lags, list) or isinstance(modelled_lags, np.ndarray):
            for i, lag in enumerate(modelled_lags):
                print(lag, nug, vmodels, sills, ranges, *vargs)
                modelled_gammas[i] = modelvariogrampoint(lag, nug, vmodels, sills, ranges, *vargs)
        
    return modelled_lags, modelled_gammas, params

def modelvariogrampoint(lag, nug, vmodels, sills, ranges, *vargs):
    variance = 0
    for i, (r1, r2) in enumerate(zip(ranges[:len(ranges)-1], ranges[1:])):
        if lag > r1 and lag <= r2:
            for n in range(i+1):
                variance += sills[n]
            # ####### Specific to each variogram model #######
            for n in range(i+1, len(ranges)):
                if vmodels[i] == "Spherical":
                    # print(ranges, sills)
                    # print(n, r)
                    if r2 != 0:      # Divide by zero exception catch all
                        variance += sills[n] * ( 1.5*(lag/ranges[n]) - 0.5*(lag/ranges[n])**3 )
                elif vmodel[i] == "Power":
                    beta=2
                    g=0.01
                    variance = nug + g * (lag**beta)
                elif vmodels[i] == "Exponential":
                    alpha = 1
                    FP = 1.1   # Oliver & Webster fitting parameter
                    variance =  (nug + sills[i]) * (alpha - math.e**(-lag/ (FP * r2) )) # Oliver & Webster 2015           
                elif vmodels[i] == "Gaussian":
                    variance =  nug + (sills[i] + nug) * (1 - (math.e)**(-(lag**2)/r2**2))
                elif vmodels[i] == "Cauchy":
                    beta = 0.5 
                    variance =  1 + ( (r2**2)/(lag**2) )**(-beta/2)
                else:
                    print("Variogram model not recognised.")
                
        # ####### Common for all variogram models #######        
        elif lag <= 0:                            # if lag below or equal to 0
            variance = 0
        elif lag > r2 and i == len(ranges)-2:     # if lag beyond the last range
            variance = sum(sills)
            # print(f"i: {i},    r1: {r1},    r2: {r2},    lag: {lag},    vmodels: {vmodels[i]},   isills: {sills[i]}")
    # print(variance)
    return variance

def paramslabelstrconstr(params):
    ###### Checks ######
    try:
        for i, v in enumerate(params):
            if i == 0:
                if type(v[0]).__name__ == "int" or type(v[0]).__name__ == "float":
                    nug = params[0][0]
                    print(f"Variogram model {i} successfully checked. Nugget with value {v[0    ]}.")
            else:
                try:
                    if len(v) >= 3:
                        print(f"Variogram model parameters are Complete.")
                except:
                    print(f"Variogram model parameters are incomplete. Only {len(params)} parameters provided: 3 minimum are required (nugget, sill and range) in list format for a simple Spherical or Exponential Model. Follow the format described below for other simp    le models\n \
                         For a Nested model, please provide a list of lists containing in the first list the nugget value and any subsequent lists the necessary parameters for that particular variogram model:\n \
                         Spherical(3) = nugget, sill, range,\
                         Exponential(3) = nugget, sill, range,\
                         Gaussian(3) = nugget, sill, range \
                        Power(5) = nugget, sill, range, Beta (represents the curvature)     and g (intensity of varioation or gradient) \
                        Cauchy(4) = nugget, sill, range, Beta \
                             \ ")
                if type(v[0]).__name__ == "str":
                    if v[0] == "Spherical":
                        if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                    elif v[0] == "Power":
                        if len(v) == 5 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                            print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Exponential":
                    if len(v) == 3 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
                elif v[0] == "Cauchy":
                    if len(v) == 4 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
                        print(f"Variogram model {i} successfully checked. {v}")
    except:
        print(f"Variogram model {i} with format {v} does not conform with the required format. Either there are too many parameters, not enough parameters or a parameters is not of the required type.")
    finally:
        if 'nug' in locals():
            params = params[1:]
        print(params)
        vmodels, sills, ranges, *vargs = list(map(list, zip(*params)))
        if len(params) == 1:
            string = f"{vmodels[0]} model nug: {nug}, r: {ranges[0]}, sill: {sills[0]}"
        else:
            string = "Nested model:\n"
            string += f"nug {nug}"
            for i in range(len(vmodels)):
                string += f"\n{vmodels[i]} model, s:{sills[i]}, r:{ranges[i]}"
    return string


""" ############# TEST #############"""
save = True
function = "Nested Spherical"
vcol = "Aperture"

params1 =    [
            [0],
            ["Spherical", 1, 2.5],
            ]
params2 =    [
            [0],
            ["Spherical", 1, 14],
            ]
params3 =    [
            [0],
            ["Spherical", 0.3, 2.5],
            ["Spherical", 0.7, 14]
            ]

params4 =    [
            [0],
            ["Spherical", 1, 8.5],
            ]

x, y1, p1 = semivar_v2(exp_var, params1, lagcol="lag", varcol="gamma")
x, y2, p2 = semivar_v2(exp_var, params2, lagcol="lag", varcol="gamma")
x, y3, p3 = semivar_v2(exp_var, params3, lagcol="lag", varcol="gamma")

x, y4, p4 = semivar_v2(exp_var, params4, lagcol="lag", varcol="gamma")
# =============================================================================
# function = "Spherical"
# nug = 0 #nugget
# r = 6 # Range
# isill = 1 #sill
# x, y = semivar(exp_var, function=function, nug=nug, r=r, isill=isill, lagcol="lag", varcol="gamma")
# plt.plot(x[:], y[:], label = f"{function} model fit to experimental variogram")
# =============================================================================

# # Works with ax.plt
# font = FontProperties()
# # font.set_family('serif')
# font.set_name('Arial')
# # font.set_style('italic')

# Works with plt
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
})

plt.plot(x[:], y1[:], label = paramslabelstrconstr(p1) )
plt.plot(x[:], y2[:], label = paramslabelstrconstr(p2) )
plt.plot(x[:], y3[:], label = paramslabelstrconstr(p3) )
plt.plot(exp_var['lag'], exp_var['gamma'], '.', markersize = 4, color = 'black', label = f"Experimental variogram")
plt.title(f'Freiberg gneiss Normalised Aperture Lag1 135deg Variogram')
plt.xticks(range(0, int(max(exp_var["lag"])), 10))
plt.grid(which="major", color='#CCCCCC', linestyle='--')
plt.grid(which='minor', color='#CCCCCC', linestyle=':')
plt.minorticks_on()
#plt.annotate(f"Nugget: {nug}\nRange: {r}\nSill: {isill}" , (exp_var["lag"].max()-15, 0.3))
plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
plt.ylabel(r'$\gamma \bf(h)$')
plt.legend(loc="lower right")
if save == True:
    plt.savefig(f"{wdsave}{inputfilename}{vcol}_{function}FIT_Arial{extension_png}", dpi=1000, bbox_inches = "tight")
plt.show()



plt.plot(x[:], y4[:], label = paramslabelstrconstr(p4) )
plt.plot(exp_var['lag'], exp_var['gamma'], '.', markersize = 4, color = 'black', label = f"Experimental variogram")
plt.title(f'Freiberg gneiss Normalised Aperture Lag1 45deg Variogram')
plt.xticks(range(0, int(max(exp_var["lag"])), 10))
plt.grid(which="major", color='#CCCCCC', linestyle='--')
plt.grid(which='minor', color='#CCCCCC', linestyle=':')
plt.minorticks_on()
#plt.annotate(f"Nugget: {nug}\nRange: {r}\nSill: {isill}" , (exp_var["lag"].max()-15, 0.3))
plt.xlabel(r'Lag Distance $\bf(h)$, (mm)')
plt.ylabel(r'$\gamma \bf(h)$')
plt.legend(loc="lower right")
if save == True:
    plt.savefig(f"{wdsave}{inputfilename}{vcol}_{function}FIT_Arial{extension_png}", dpi=1000, bbox_inches = "tight")
plt.show()



""" ############# TEST #############"""

####################################

def VarFit(df, *args):
    '''
    x:
    *args: Accepts only a list containing series of the following:
        function: The variogram model you wish to use for a specific interval. Options are "spherical", "power", "exponential" and "gaussian".
        nug: The nugget of the variogram model relating to the interval.
        r: The range of the variogram model relating to the interval.
        isill: The sill of the variogram model relating to the interval.
        end: The lag-value related to the end of the interval of this particular variogram model.
    
    *** All five parameters need to be provided, and in the specific order, for each variogram model interval, i.e. the list must have a size which is a multiple of five ***
    
    Example: VarFit(Dataframe, ['exponential', 0, 90, 2, 35, 'gaussian', 0, 110, 2, 100, 'power', 1, 150, 2.5, 120])
    
    formulas from Oliver & Webster 2015
    '''
    global comp, x, y
    
    try:
        type(df) == pd.DataFrame()
    except:
        print("First variable is not a pandas.core.frame.DataFrame. A pandas.core.frame.DataFrame is needed for the x values to be plotted.")
    try:
        (len(args) % 5 == 0 or len(args) == 0)
    except:
        print("The number of arguments is not concordant with the function chosen.")
    
    if df.values[-1][0] == 0.0:                                                                             #if last row has lag 0 (duplicate)
        df.drop(df.tail(1).index,inplace=True)                                                              # drop last n rows
    x = np.array(pd.concat([df["lag"], pd.DataFrame(arange(0,df["lag"].max(),0.1))]).sort_values(by=[0]))   # creates a 1-column DF of the lags in the input DF and merges it with a range from 0 to input-DF maximum ever 0.1 then sorts it.
    #x = x[np.argsort(x)]                                                                                   # Sorting the lag values (x-values)
    y = np.zeros((x.shape[0]), dtype=x.dtype).transpose()                                                   # adds a column with zeros

    args = args[0]                                                  # Extracts the list of arguments from the tuple (python creates a tuple with one single list of all the arguments, which needs to be extracted), i.e. args = (list,) which needs args[0] to extract list from the tuple.
    comp = len(args) // 5

    
    print(f'*args has size {len(args)}. arguments: {args}')
    print("Iterating over the optional arguments")
    for i in args:
        print(i)
    
    start = 0
    for w in range(comp):       
        print(f'\n#####\n{w+1}st model\n#####\n') if w == 0 else ( print(f'\n#####\n{w+1}nd model\n#####\n') if w == 1 else ( print(f'\n#####\n{w+1}rd model\n#####\n') if w == 2 else print(f'\n#####\n{i+w}th model\n#####\n')))
        
        for j in range(5): # Every 5 parameters, i.e.function, nugget, range, sill and limit
            #print(j+i*5)
            if j == 0: function = args[j+w*5]; print(f'function: {function}')
            if j == 1: nug = args[j+w*5]; print(f'nugget: {nug}')
            if j == 2: r = args[j+w*5] ; print(f'range: {r}')
            if j == 3: isill = args[j+w*5] ; print(f'sill: {isill}')
            if j == 4: limit = args[j+w*5]; print(f'limit: {limit}')
            
        print(f'function: {function}, nugget: {nug}, range: {r}, sill: {isill}, start: {start}, limit: {limit}')    
        
        for i in range(len(x)):
            if start < x[i] <= limit:
                
                if function == "spherical":
                    if x[i] == 0:
                        y[i] = 0
                    elif x[i] > 0 and x[i] <= r:
                        y[i] = nug + (isill - nug) * ((3*x[i])/(2*r) - 0.5*(x[i]/r)**3) # geostatspy uses y[i] = nug + isill  * (1.5*(x[i]/r) - 0.5*(x[i]/r)**3) which seems to give the same result
                    elif x[i] > r:
                        y[i] = isill
                elif function == "power":
                    if x[i] <= r:
                        y[i] = nug + g*x[i]**beta
                    else:
                        y[i] = isill
                elif function == "exponential":
                    if x[i] <= r:
                        y[i] =  nug + (isill + nug) * (1 - math.e**(-x[i]/r))
                    else:
                        y[i] = isill
                elif function == "gaussian":
                    if x[i] <= r:
                        y[i] =  nug + (isill + nug) * (1 - (math.e)**(-(x[i]**2)/r**2))
                    else:
                        y[i] = isill

        start = limit

    plt.plot(x[:], y[:], label = 'Azimuth', **pltfont)
    plt.plot(df['lag'], df['gamma'], 'x', markersize = 4, color = 'black')
    # plt.plot([0,2000], [df[vr].var(), df[vr].var()], color = 'black')
    plt.plot([0,2000], [isill, isill], color = 'black')                                 # Is this what was meant in the line above? plt.plot([0,2000], [df[vr].var(), df[vr].var()], color = 'black') ; vr is not defined anywhere else
    plt.xlim([0,df['lag'].max()+df['lag'].max()*0.1])
    plt.ylim([0,df['gamma'].max()+df['gamma'].max()*0.1])
    plt.title(f'Azimuth model: {function}', **pltfont)
    plt.show()


VarFit(exp_var, ['spherical', nug, r, isill, exp_var['lag'].max()])

"""  """
# import sys
# sys.path.append("c:\\users\\s2132627\\appdata\\local\\programs\\python\\python39\\lib\\site-packages")
# sys.path



""" For each of the experimental variograms """
# =============================================================================
# if function == "spherical":
#     for i in range(len(x)):
#         if x[i] == 0:
#             y[i] = 0
#         elif x[i] > 0 and x[i] <= r:
#             y[i] = nug + (isill - nug) * ((3*x[i])/(2*r) - 0.5*(x[i]/r)**3)
#         elif x[i] > r:
#             isill
# elif function == "power":
#     for i in range(len(x)):
#         y[i] = nug + g*x[i]**beta
# elif function == "exponential":
#     for i in range(len(x)):
#             y[i] =  nug + (isill + nug) * (1 - math.e**(-x[i]/r))
# elif function == "gaussian":
#     for i in range(len(x)):
#         y[i] =  nug + (isill + nug) * (1 - (math.e)**(-(x[i]**2)/r**2))
# 
# plt.plot(x[:], y[:], label = 'Azimuth ' + str(azi_mat[iazi]))
# plt.plot(lag[iazi,:],gamma[iazi,:], 'x', markersize = 4, color = 'black', label = 'Azimuth ' + str(azi_mat[iazi]))
# =============================================================================



# =============================================================================
# function = "gaussian"
# 
# x=data['lag']; x = x[np.argsort(x)] # Sorting the lag values (x-values)
# y = np.zeros((x.shape[0]), dtype=x.dtype).transpose() # adds a column with zeros
# 
# nug = 0 #nugget
# rang = 15 # Range
# isill = 1 #sill
# if function == "spherical":
#     for i in range(len(x)):
#         if x[i] == 0:
#             y[i] = 0
#         elif x[i] > 0 and x[i] <= rang:
#             y[i] = nug + (isill - nug) * ((3*x[i])/(2*rang) - 0.5*(x[i]/rang)**3)
#         elif x[i] > rang:
#             isill
# elif function == "power":
#     for i in range(len(x)):
#         y[i] = nug + g*x[i]**beta
# elif function == "exponential":
#     for i in range(len(x)):
#             y[i] =  nug + (isill + nug) * (1 - math.e**(-x[i]/rang))
# elif function == "gaussian":
#     for i in range(len(x)):
#         y[i] =  nug + (isill + nug) * (1 - (math.e)**(-(x[i]**2)/rang**2))
# 
# 
# plt.plot(data['lag'],data['gamma'], 'x', markersize = 4, color = 'black', label = 'Azimuth ')
# plt.plot(x[:], y[:])
# plt.show()
# =============================================================================

