# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:06:34 2024

@author: s2132627
"""

"""
This script maps initial point (x,y,z, aperture) data to mesh data (nodes, elements) so that we can create apertures and hydraulic conductivities (currently hardcoded for the cubic law) text files for the mesh elements (as opposed to spatial points) that OGS can read.

In this version, I use the elements' central points in the .msh file and krig the aperture from the .csv file. It requires therefore the spatial continuity parameters for the kriging algorithm. Alternatively, as with the last version, we can use the element's bounds x,y,z to filter all the .csv points that fall within that element and average them (works only for quad elements at the moment because I have been having issues with the winding number's algorithm). This way I don't have to faff about with KDtrees and such. Much simpler, probably more accurate.




Below are an example of the .csv file and the (GINA's MeshFormat 2.2) .msh file for reference.
***  .csv  ***
x	y	aperture	avg	var	std
-22.5	-32.5	0.49977896	0.515901725	0.005059804	0.071132302
-22.5	-31.5	0.48803121	0.508111098	0.004592228	0.067765978         
...     ...     ...         ...         ...         ...

***  GINA's .msh  ***
$NODES
  8186 
0 -22.5 -32.5 1
1 -22.5 -31.5 1
...
 $ELEMENTS
  8184 
0 0 tri 0 3977 4220
1 0 tri 4220 3977 4465

"""
import os
import pathlib
from tkinter import *               #Import everything from Tkinter
from tkinter import ttk
import tkinter.messagebox           #Import messagebox method
from tkinter import filedialog
import subprocess
from time import time, ctime
import time
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
from ast import literal_eval
import math
from scipy.spatial import KDTree
import skgstat as skg
import matplotlib.pyplot as plt
from pykrige.ok import OrdinaryKriging
import sys
import threading
from math import floor
'''User Defined Global Variables'''
output_extension = ".txt"
eamthreshold = 1E-10
prediction_method =  "Element Centre Kriging"
header = 0
readcols = [1,2,3]
skiprows = [""]
fringe_aperture = "max"
memory_mgmt = True # memory management. if Yes, it loops through each point, trims the df to a square around the point with side length rangeM
'''Global Variables'''
FractureMaterialGroup = 0
rangeM = 5
rangem = 2
angle = 157.5
variogram_model = "spherical"
sill = 1
nugget = 0
dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
ratio = rangem/rangeM
param = [sill, rangeM, nugget] #sill, range, nugget
# params3 =    [
#             [0],
#             ["Spherical", 0.3, 2.5],
#             ["Spherical", 0.7, 14]
#             ]
"""System Variables"""
cmap = plt.cm.plasma                    # color map
conductivity_method = "Cubic Law"
varmdl_method = "Spherical"
extension_csv = ".csv"
extension_txt = ".txt"


'''Methods'''
def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s
def eq(aperture):
    if conductivity_method == "Cubic Law":
        k = (aperture**2)/12; 
        conductivity_method_name = "Cubic Law";
    return k
def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s
def normal_round(n):
    if n - math.floor(n) < 0.5:
        return math.floor(n)
    return math.ceil(n)
def prec(df):
    precision = lambda x: (len(x)-1, len(x.split('.')[1]))
    precision = df.astype(str).applymap(precision).max().max()[1]
    #print(precision)
    return precision
def nsmall(a, n):
    return np.partition(a, n)[n]                                                #finds the n-th smallest value
class Spinner:
    busy = False
    delay = 0.1
    @staticmethod
    def spinning_cursor():
        while 1: 
            for cursor in '|/─\\': yield cursor
    def __init__(self, delay=None):
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay): self.delay = delay
    def spinner_task(self):
        while self.busy:
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()
    def __enter__(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()
    def __exit__(self, exception, value, tb):
        self.busy = False
        time.sleep(self.delay)
        if exception is not None:
            return False
def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f' days: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
# =============================================================================
#         # ############################# For testing ##############################################################################################################################################################################
# 
#         lc = 1
#         header = 0
#         readcols = [0,1,2]
#         skiprows = ""
#         CSVinputfilepath = pathlib.PureWindowsPath(r"C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Greywacke\GWQ4_Vfringes - ST - 1atmBP - TEST\z_ORIGvsORIG-rot(-0.02, 0.31, 0)_ApertureMap.csv").as_posix()
#         CSVinputfilename = CSVinputfilepath[CSVinputfilepath.rfind("/")+1:CSVinputfilepath.rfind(".")]
#         pathCSV = CSVinputfilepath[:CSVinputfilepath.rfind("/")+1]
#         CSVinput_extension = CSVinputfilepath[CSVinputfilepath.rfind("."):]
#         input_CSVfile = open(CSVinputfilepath, 'r')
#         # read the content of the file line by line 
#         CSVdata_input = input_CSVfile.readlines()
# 
# 
#         MSHinputfilepath = pathlib.PureWindowsPath(r"C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Greywacke\GWQ4_Vfringes - ST - 1atmBP - TEST\GW_Q4_90degOrientation.msh").as_posix()
#         MSHinputfilename = MSHinputfilepath[MSHinputfilepath.rfind("/")+1:MSHinputfilepath.rfind(".")]
#         pathMSH = MSHinputfilepath[:MSHinputfilepath.rfind("/")+1]
#         MSHinput_extension = MSHinputfilepath[MSHinputfilepath.rfind("."):]
#         input_MSHfile = open(MSHinputfilepath, 'r')
#         # read the content of the file line by line 
#         MSHdata_input = input_MSHfile.readlines()
# 
#         
#         header = int(0)
#         readcols = [1,2,3]
#         skiprows = ""
#         if "" in skiprows:
#             df = pd.read_csv(pathCSV + CSVinputfilename + extension_csv , usecols = readcols, header = header)
#         else:
#             skiprows = [int(i) for i in skiprows]
#             df = pd.read_csv(pathCSV + CSVinputfilename , usecols = readcols, header = header, skiprows = skiprows)
#         
#         dfMSH_OF = pd.read_csv(pathMSH_OF + MSH_OFinputfilename)
#         
#         with open(pathMSH_OF + MSH_OFinputfilename) as f:
#             contents = f.read()
#         f.close()
#         # Create a KDTree. We will need it for finding the closest points from df to the fringe points
#         nd = len(df)                                                                                                #  number of total points  
#         tree = KDTree(df.to_numpy())  # Create the K-D Tree used for a (very!) quick computation without the need for nested loops (how does this even work so fast??)
#         # ###########################################################################################################################################################################################################
# =============================================================================

def My_Custom_Function(params, dists):
    """
    Parameters
    ----------
    dataframe : Pandas dataframe
        Contains the experimental variogram with the lags and the corresponding semi-variances columns.
    params : List of Lists
        Provide a list containing, in this order, a list with the nugget (only!) followed a list containing the variogram type (string), variogram partial sill (float) and variogram range (float).
        For nested variograms (more than one structure and/or contribution), add any number of lists containing the same information ["type", partial sill, range] for each contribution of the nested structure.
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
                        if r2 != 0:      # Divide by zero exception catch all
                            variance += sills[n] * ( 1.5*(lag/ranges[n]) - 0.5*(lag/ranges[n])**3 )
                    elif vmodel[i] == "Power":
                        beta=2
                        g=0.01
                        variance = nug + g * (lag**beta)
                    elif vmodels[i] == "Exponential":
                        alpha = 1
                        FP = 1.1   # Oliver & Webster fitting parameter
                        variance =  (nug + sills[i]) * (alpha - math.e**(-lag/ (FP * r2) )) # Oliver & Webster 2015           cc
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

"""########
    Same as modelvariogrampoint but formated to gstools as a user defined covariance model
########## """
class UserDefinedVariogram(gs.CovModel):
    def cor(self, h, params):
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
                        print(f"Variogram model parameters are incomplete. Only {len(params)} parameters provided: 3 minimum are required (nugget, sill and range) in list format for a simple Spherical or Exponential Model. Follow the format described         below for other simp    le models\n \
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
                            if len(v) == 5 and (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[2]).__name__ == "int" or type(v[1]).__name__ == "float") or (type(v[1]).__name__ == "int" or type(v[1]).__name__ == "float")     and (type    (v[1]).__name__ == "int" or type(v[1]).__name__ == "float"):
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
            ###### Unpack params and format ######
            if 'nug' in locals():
                params_f = params[1:]
            vmodels, sills, ranges, *vargs = list(map(list, zip(*params_f)))
            sills.insert(0, nug)
            ranges.insert(0, 0)
            if nug != 0.0:
                ranges.insert(0, 0)
            print(vmodels, sills, ranges, *vargs)
            
            variance = 0
            for i, (r1, r2) in enumerate(zip(ranges[:len(ranges)-1], ranges[1:])):
                if h > r1 and lag <= r2:
                    for n in range(i+1):
                        variance += sills[n]
                    # ####### Specific to each variogram model #######
                    for n in range(i+1, len(ranges)):
                        if vmodels[i] == "Spherical":
                            # print(ranges, sills)
                            # print(n, r)
                            if r != 0:      # Divide by zero exception catch all
                                variance += sills[n] * ( 1.5*(h/ranges[n]) - 0.5*(h/ranges[n])**3 )
                        elif vmodel[i] == "Power":
                            beta=2
                            g=0.01
                            variance = nug + g * (h**beta)
                        elif vmodels[i] == "Exponential":
                            alpha = 1
                            FP = 1.1   # Oliver & Webster fitting parameter
                            variance =  (nug + sills[i]) * (alpha - math.e**(-h/ (FP * r2) )) # Oliver & Webster 2015           
                        elif vmodels[i] == "Gaussian":
                            variance =  nug + (sills[i] + nug) * (1 - (math.e)**(-(h**2)/r2**2))
                        elif vmodels[i] == "Cauchy":
                            beta = 0.5 
                            variance =  1 + ( (r2**2)/(h**2) )**(-beta/2)
                        else:
                            print("Variogram model not recognised.")
    
                # ####### Common for all variogram models #######        
                elif h <= 0:                            # if lag below or equal to 0
                    variance = 0
                elif h > r2 and i == len(ranges)-2:     # if lag beyond the last range
                    variance = sum(sills)
                    # print(f"i: {i},    r1: {r1},    r2: {r2},    lag: {lag},    vmodels: {vmodels[i]},   isills: {sills[i]}")
            # print(variance)
        return variance


def msh_of_SplitFunction_v4(path, input_file_name):
    """
    Initially built to create a pandas Data Frame from the .msh_of file from GINA's OGS pre-processor (hence the name), it now is able to distinguish between that and a Gmsh .msh file to create the same pandas df.
    
    Parameters
    ----------
    path : String
        system path address to the file.
    input_file_name : TYPE
        system file name provided without the extension. In the case of GINA's .msh_of, please provide the .msh file and the code will look for the .msh_of file itself.

    Returns
    -------
    dfN : Pandas Data Frame
        Pandas Data Frame with the nodes and their xy(z) locations.
    dfE : Pandas Data Frame
        Pandas Data Frame with the elements and their corresponding nodes.
    ENC_df : Pandas Data Frame
        Pandas Data Frame with the elements and their values.
    MGs : Pandas Data Frame
        Pandas Data Frame with the material groups of each element.
    Etype : String
        The type of element.

    """
    # #############################  Get the Nodes and Elements from the .msh (GINA's) file without having to export and read back .txt files #############################
    
    with open(path + input_file_name + ".msh", 'r') as file:
        
        # Check the file format
        head = file.readlines()[0:50] # Check the first 50 lines
        print(head)
        Format = "GINA"
        for l in head:    
            if "$MeshFormat" in l:
                Format = "Gmsh"
        if Format == "GINA":
            words = ["$NODES", "$ELEMENTS"]
        elif Format == "Gmsh":
            words = ["$Nodes", "$Elements"]
        print(Format)
        
    with open(path + input_file_name + ".msh", 'r') as file:                            # The below code block doesn't run it I put it in the same with open block above. I don't know why
        # This Loop is not working in tkinter for whatever reason so I hardcoded it
        # for word in words: 
        #     exec(f"{word[word.rfind('$')+1:]} = False") # Creates a False flag for each word in words
        NODES = False
        ELEMENTS = False
        for line in file:
            # print(line)
            casematch = next((word for word in words if word in line), False)
            if casematch:
                # for word in words: # Creates a False flag for each word in words  
                #     exec(f"{word[word.rfind('$')+1:]} = False")
                NODES = False
                ELEMENTS = False
                # exec(f"{match[match.rfind('$')+1:]} = True") # Sets the match word flag to True
                print(casematch)
                print(NODES)
                print(ELEMENTS)
                # if casematch.find("NODES") != -1 or casematch.find("Nodes") != -1: # If line contains "NODES"
                if "nodes" in casematch.lower():
                    NODES = True
                    Nnodes = int(next(file))
                    print(Nnodes)
                    Ns = []
                # elif casematch.find("ELEMENTS") != -1 or casematch.find("Elements") != -1: # If line contains "ELEMENTS"
                elif "elements" in casematch.lower():
                    ELEMENTS = True
                    Nelements = int(next(file))
                    print(Nelements)
                    Es = []
            else:
                if Format == "GINA":
                    if line.find("#STOP") != -1:
                        break
                    elif NODES == True:
                        Ns.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                    elif ELEMENTS == True:
                        Es.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                elif Format == "Gmsh":
                    if "$End" in line: # if line contains
                        NODES = False
                        ELEMENTS = False
                    if NODES == True:
                        Ns.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                    elif ELEMENTS == True:
                        Es.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
        
    dfN = pd.DataFrame(Ns, columns = ["NodeTag","x","y","z"]) # reads only the nodes
    # dfN = dfN.round(prec(df)) # This is necessary because of precision changes, later on when I am trying to find the aperture corresponding to the x,y,z values of a node, those x,y,z won't match.
    if Format == "GINA":
        Etype = Es[0][2] # Get E type (tri, quad)
        EsMaxLen = max([len(i) for i in Es])
        names = ["ElementNumber", "MaterialGroup", "ElementType"] 
        [names.append (f"Node{i+1}") for i in range(EsMaxLen - len(names))] # Add entries to the names list corresponding to the unknown number of nodes
        dfE = pd.DataFrame(Es, columns = names); #dfE = dfE.round(prec(df))
        EmaxN = max([(simplest_type(n[n.rfind("Node")-1:])) for n in list(dfE.columns) if n.find("Node") != -1]) # Count the number of columns with "Node"
        print(dfN)
        print(dfE)
        if os.path.isfile(path + input_file_name + ".msh_of"):
            with open(path + input_file_name + ".msh_of") as f:
                contents = f.read()
            f.close()
            Esplit = contents.split("\n\n")
            Esplit = [ s.split("\n") for s in Esplit] # Split contents into elements
            ENC = [[simplest_type(l2.strip()) for l2 in l1] for l1 in Esplit] # Elements' strings cleaned and converted into numbers    # Instead of for i, l1 in enumerate(f2split):
            EsMaxLen = max([len(i) for i in ENC])
            EsC_columns = ["MaterialGroup", "ElementNumber"] 
            [EsC_columns.append (f"Node{i+1}") for i in range(EsMaxLen - len(EsC_columns) - 5)] # Add entries to the names list corresponding to the unknown number of nodes. To get the unknown number of nodes, loop through a range to N where N is the number of total values read from the file minus the number of values that are known to make up the file.
            EsC_columns.extend (["ElementCentreX-Coordinate", "ElementCentreY-Coordinate", "ElementCentreZ-Coordinate", "Dunno", "Dunno2" ]) # Elements' Centres
            if "qua" in Etype.lower():
                ENC_df = pd.DataFrame(ENC[1:], columns = EsC_columns) # Start at row1 ignoring row 0 because the latter is not an Element
                # ENC_df_ElementCentrexCoordinate_round = round(ENC_df["Element Centre x-Coordinate"])
                ENC_df["GridX-coordinate"]  = (round(ENC_df["ElementCentreX-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))                                        # Round because there might be a precision issue. This rounding does not change the initial dataframe. -1 because the ranking function starts at 1
                ENC_df["GridX-coordinate"]  = ENC_df["GridX-coordinate"] / (np.sort(ENC_df["GridX-coordinate"].unique())[1] - np.sort(ENC_df["GridX-coordinate"].unique())[0])      # rank function assigns the same rank to duplicates but jumps that amount of duplicates when the next value is different, which we don't want. Thus we are calculating the difference between sequential values and dividing..
                ENC_df["GridY-coordinate"]  = (round(ENC_df["ElementCentreY-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))
                ENC_df["GridY-coordinate"]  = ENC_df["GridY-coordinate"] / (np.sort(ENC_df["GridY-coordinate"].unique())[1] - np.sort(ENC_df["GridY-coordinate"].unique())[0])
        else:      
            ENC_df = dfE
            ENC_df["ElementCentreX-Coordinate"] = 0
            ENC_df["ElementCentreY-Coordinate"] = 0
            ENC_df["ElementCentreZ-Coordinate"] = 0
            ENC_df["GridX-coordinate"] = 0
            ENC_df["GridY-coordinate"] = 0
            for row in ENC_df.itertuples():
# =============================================================================
#                 # Deprecated because below I don't need to hardcode
#                 if "tri" in Etype.lower():
#                     Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3") ]
#                 if "qua" in Etype.lower():
#                     Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3"), getattr(row, "Node4") ]
# =============================================================================
                Ns = [getattr(row, f"Node{i+1}") for i in range(EmaxN)]
                Ns = dfN[dfN["NodeTag"].isin(Ns)]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreX-Coordinate']  = Ns.x.describe()["mean"]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreY-Coordinate']  = Ns.y.describe()["mean"]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreZ-Coordinate']  = Ns.z.describe()["mean"]
            ENC_df["GridX-coordinate"]  = np.argsort(ENC_df["ElementCentreX-Coordinate"])
            ENC_df["GridY-coordinate"]  = np.argsort(ENC_df["ElementCentreY-Coordinate"])
            ENC_df["GridZ-coordinate"]  = np.argsort(ENC_df["ElementCentreZ-Coordinate"])
        MGs = np.sort(ENC_df["MaterialGroup"].unique() )# np.array containing all used Material Groups' numbers
    elif Format == "Gmsh":
        EsMaxLen = max([len(i) for i in Es])
        names = ["ElementNumber", "ElementTypeNumber", "ElementUnknown1", "ElementUnknown2", "MaterialGroup"] # Element number, Element type, Element unknown characteristic 1, Element unknown characteristic 2, Element Material Group
        [names.append (f"Node{i+1}") for i in range(EsMaxLen - len(names))] # Add entries to the names list corresponding to the unknown number of nodes
        dfE = pd.DataFrame(Es, columns = names)
        
        Gmshtypes_unique = dfE["ElementTypeNumber"].unique().tolist() # The second column has the element identifier. See https://gitlab.onelab.info/gmsh/gmsh/blob/master/src/common/GmshDefines.h as of 12.02.2024
        Gmshtypes = [
                    ['MSH_LIN_2', 1],         ['MSH_TRI_3', 2],         ['MSH_QUA_4', 3],         ['MSH_TET_4', 4],         ['MSH_HEX_8', 5],         ['MSH_PRI_6', 6],         ['MSH_PYR_5', 7],         ['MSH_LIN_3', 8],
                    ['MSH_TRI_6', 9],         ['MSH_QUA_9', 10],         ['MSH_TET_10', 11],         ['MSH_HEX_27', 12],         ['MSH_PRI_18', 13],         ['MSH_PYR_14', 14],         ['MSH_PNT', 15],         ['MSH_QUA_8', 16],
                    ['MSH_HEX_20', 17],         ['MSH_PRI_15', 18],         ['MSH_PYR_13', 19],         ['MSH_TRI_9', 20],         ['MSH_TRI_10', 21],         ['MSH_TRI_12', 22],         ['MSH_TRI_15', 23],         ['MSH_TRI_15I', 24],
                    ['MSH_TRI_21', 25],         ['MSH_LIN_4', 26],         ['MSH_LIN_5', 27],         ['MSH_LIN_6', 28],         ['MSH_TET_20', 29],         ['MSH_TET_35', 30],         ['MSH_TET_56', 31],         ['MSH_TET_22', 32],
                    ['MSH_TET_28', 33],         ['MSH_POLYG_', 34],         ['MSH_POLYH_', 35],         ['MSH_QUA_16', 36],         ['MSH_QUA_25', 37],         ['MSH_QUA_36', 38],         ['MSH_QUA_12', 39],         ['MSH_QUA_16I', 40],
                    ['MSH_QUA_20', 41],         ['MSH_TRI_28', 42],         ['MSH_TRI_36', 43],         ['MSH_TRI_45', 44],         ['MSH_TRI_55', 45],         ['MSH_TRI_66', 46],         ['MSH_QUA_49', 47],         ['MSH_QUA_64', 48],
                    ['MSH_QUA_81', 49],         ['MSH_QUA_100', 50],         ['MSH_QUA_121', 51],         ['MSH_TRI_18', 52],         ['MSH_TRI_21I', 53],         ['MSH_TRI_24', 54],         ['MSH_TRI_27', 55],         ['MSH_TRI_30', 56],
                    ['MSH_QUA_24', 57],         ['MSH_QUA_28', 58],         ['MSH_QUA_32', 59],         ['MSH_QUA_36I', 60],         ['MSH_QUA_40', 61],         ['MSH_LIN_7', 62],         ['MSH_LIN_8', 63],         ['MSH_LIN_9', 64],
                    ['MSH_LIN_10', 65],         ['MSH_LIN_11', 66],         ['MSH_LIN_B', 67],         ['MSH_TRI_B', 68],         ['MSH_POLYG_B', 69],         ['MSH_LIN_C', 70],         ['MSH_TET_84', 71],         ['MSH_TET_120', 72],
                    ['MSH_TET_165', 73],         ['MSH_TET_220', 74],         ['MSH_TET_286', 75],         ['MSH_TET_34', 79],         ['MSH_TET_40', 80],         ['MSH_TET_46', 81],         ['MSH_TET_52', 82],         ['MSH_TET_58', 83],
                    ['MSH_LIN_1', 84],         ['MSH_TRI_1', 85],         ['MSH_QUA_1', 86],         ['MSH_TET_1', 87],         ['MSH_HEX_1', 88],         ['MSH_PRI_1', 89],         ['MSH_PRI_40', 90],         ['MSH_PRI_75', 91],
                    ['MSH_HEX_64', 92],         ['MSH_HEX_125', 93],         ['MSH_HEX_216', 94],         ['MSH_HEX_343', 95],         ['MSH_HEX_512', 96],         ['MSH_HEX_729', 97],         ['MSH_HEX_1000', 98],         ['MSH_HEX_32', 99],
                    ['MSH_HEX_44', 100],         ['MSH_HEX_56', 101],         ['MSH_HEX_68', 102],         ['MSH_HEX_80', 103],         ['MSH_HEX_92', 104],         ['MSH_HEX_104', 105],         ['MSH_PRI_126', 106],         ['MSH_PRI_196', 107],
                    ['MSH_PRI_288', 108],         ['MSH_PRI_405', 109],         ['MSH_PRI_550', 110],         ['MSH_PRI_24', 111],         ['MSH_PRI_33', 112],         ['MSH_PRI_42', 113],         ['MSH_PRI_51', 114],         ['MSH_PRI_60', 115],
                    ['MSH_PRI_69', 116],         ['MSH_PRI_78', 117],         ['MSH_PYR_30', 118],         ['MSH_PYR_55', 119],         ['MSH_PYR_91', 120],         ['MSH_PYR_140', 121],         ['MSH_PYR_204', 122],         ['MSH_PYR_285', 123],
                    ['MSH_PYR_385', 124],         ['MSH_PYR_21', 125],         ['MSH_PYR_29', 126],         ['MSH_PYR_37', 127],         ['MSH_PYR_45', 128],         ['MSH_PYR_53', 129],         ['MSH_PYR_61', 130],         ['MSH_PYR_69', 131],
                    ['MSH_PYR_1', 132],         ['MSH_PNT_SUB', 133],         ['MSH_LIN_SUB', 134],         ['MSH_TRI_SUB', 135],         ['MSH_TET_SUB', 136],         ['MSH_TET_16', 137],         ['MSH_TRI_MINI', 138],         ['MSH_TET_MINI', 139],
                    ['MSH_TRIH_4', 140]
                    ] # Added the corresponding number of nodes in  some of the mesh element types.
        Gmshtypes_df = pd.DataFrame(Gmshtypes, columns = ["ElementType" , "ElementTypeNumber"])
        # Gmshtypes_df['No_Nodes'] = Gmshtypes_df['No_Nodes'].fillna(0).astype(int) # Not absolutely necessary. The objective is to convert theh number of nodes to int datatype but because most types don't have a number of nodes associated, pandas converts the whole column to float and those rows that don't have a node value become NaNs.
        for i in Gmshtypes_unique:
            exec ( f"Es{i} = dfE[dfE['ElementTypeNumber'] == {i}]" )  # Separate Gmsh elements (all points, lines, any time of element) into different pd dfs
        dfE = dfE[dfE["ElementTypeNumber"].isin([15, 1]) == False]    # Subset of the DataFrame excluding all mesh points (dfE["EtypeNo"] = 15) and mesh lines(dfE["EtypeNo"] = 1)
        ENC_df = dfE
        ENC_df["ElementCentreX-Coordinate"] = 0
        ENC_df["ElementCentreY-Coordinate"] = 0
        ENC_df["ElementCentreZ-Coordinate"] = 0
        ENC_df["GridX-coordinate"] = 0
        ENC_df["GridY-coordinate"] = 0
        for row in ENC_df.itertuples():
            Ns = [ getattr(row, f"Node{n+1}") for n in range( len([n for n in ENC_df.columns if n.find("Node") != -1]) )]       # Create list with the nodes' numbers of each element by looping through all the column names that contain "Node". This way we don't need to know a priori the amount of nodes and the type of element
            Ns = dfN[dfN["NodeTag"].isin(Ns)]                                                                                   # Creates a pd df with only the nodes from Ns list above
            ENC_df.at[getattr(row, "Index"), 'ElementCentreX-Coordinate']  = Ns.x.describe()["mean"]                          #   Calculates the mean to find the centre
            ENC_df.at[getattr(row, "Index"), 'ElementCentreY-Coordinate']  = Ns.y.describe()["mean"]
            ENC_df.at[getattr(row, "Index"), 'ElementCentreZ-Coordinate']  = Ns.z.describe()["mean"]
        ENC_df["GridX-coordinate"]  = np.argsort(ENC_df["ElementCentreX-Coordinate"])                                        # Sorts the elements by grid coordinates
        ENC_df["GridY-coordinate"]  = np.argsort(ENC_df["ElementCentreY-Coordinate"])
        MGs = np.sort(ENC_df["MaterialGroup"].unique() )# np.array containing all used Material Groups' numbers                # Creates a np array with all the Material Groups and shorts it
        Etype = Gmshtypes_df.loc[Gmshtypes_df['ElementTypeNumber'] == dfE["ElementTypeNumber"].unique()[0], 'ElementType'].iloc[0]
        
    return dfN, dfE, ENC_df, MGs, Etype
    

""" System variables treatment """
angOK=simplest_type(360-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html




if "" in skiprows:
    df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header)
else:
    skiprows = [simplest_type(i) for i in skiprows]
    df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header, skiprows = skiprows)           #  

dfcolnames = df.columns
xcol, ycol, vcol = dfcolnames    
# df = df.astype({f"{xcol}": np.float32, f"{ycol}": np.float32, f"{vcol}": np.float32})

        
print(dfcolnames)

# =============================================================================
#         # ############################# For testing ##############################################################################################################################################################################
#         lc = 1
#         fringe_aperture = 0.005
#         header = 0
#         readcols = [0,1,2]
#         skiprows = ""
#         CSVinputfilepath = pathlib.PureWindowsPath(r'M:\Documents\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture\Extracted M CNL_Aperture.csv').as_posix()
#         CSVinputfilename = CSVinputfilepath[CSVinputfilepath.rfind("/")+1:]
#         pathCSV = CSVinputfilepath[:CSVinputfilepath.rfind("/")+1]
# 
#         MSHinputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail-4x4\Hfringes\wMIDPOINTS\FB.msh').as_posix()
#         MSHinputfilename = MSHinputfilepath[MSHinputfilepath.rfind("/")+1:]
#         pathMSH = MSHinputfilepath[:MSHinputfilepath.rfind("/")+1]
# 
#         header = int(0)
#         readcols = [1,2,3]
#         skiprows = ""
#         if "" in skiprows:
#             df = pd.read_csv(pathCSV + CSVinputfilename , usecols = readcols, header = header)
#         else:
#             skiprows = [int(i) for i in skiprows]
#             df = pd.read_csv(pathCSV + CSVinputfilename , usecols = readcols, header = header, skiprows = skiprows)
#         # Create a KDTree. We will need it for finding the closest points from df to the fringe points
#         nd = len(df)                                                                                                #  number of total points  
#         tree = KDTree(df.to_numpy())  # Create the K-D Tree used for a (very!) quick computation without the need for nested loops (how does this even work so fast??)
#         # ###########################################################################################################################################################################################################
# =============================================================================

dfN, dfE, ENC_df, MGs, Etype = msh_of_SplitFunction_v4(pathMSH, MSHinputfilename)

if prediction_method == "Element Centre Kriging":
    print("# =====================================================================================================================================================\n \
    Kriging...\n \
    # =====================================================================================================================================================")
    with Spinner():
        if memory_mgmt == False:
            df2 = df.astype({"x": np.float16, "y": np.float16, vcol: np.float16})
            data = np.array(
                [
                    df2['x'],
                    df2['y'],
                    df2[vcol]
                ]
                ).transpose()
            
            print(df2)
            print("\n\n\n")
            print(data)
        
            if varmdl_method == "Custom":
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
            elif varmdl_method == "Spherical":
                OK = OrdinaryKriging(
                    data[:, 0],
                    data[:, 1],
                    data[:, 2],
                    variogram_model="spherical",
                    anisotropy_angle = angOK,
                    variogram_parameters = dic,
                    anisotropy_scaling = ratio,
                    exact_values  = True,
                    verbose=True,
                    enable_plotting=False,
                ) # scaling is the ratio between the major and minor directions' ranges
        
            ###############################################################################
            # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
            # grid of points, on a masked rectangular grid of points, or with arbitrary points.
            # (See OrdinaryKriging.__doc__ for more information.)
            # style (str) – Specifies how to treat input kriging points. Specifying ‘grid’ treats xpoints and ypoints as two arrays of x and y coordinates that define a rectangular grid.
            # Specifying ‘points’ treats xpoints and ypoints as two arrays that provide coordinate pairs at which to solve the kriging system.
            # Specifying ‘masked’ treats xpoints and ypoints as two arrays of x and y coordinates
        elif memory_mgmt == True:
            ENC_df[vcol] = 0
            # If memory is an issue, we can krig the centre of each element at a time. I think the traditional kriging system utilises the previous predictions to inform the rest of the predictions, which won't happen in this case but at least we can get some output......   
            for i, row in ENC_df.iterrows():
                if i%50 == 0:
                    print(f"#####\nKriging Element {i}/{len(ENC_df)}\n#####")
                df_filt = df[(df["x"] > row["ElementCentreX-Coordinate"] - rangeM) & (df["x"] < row["ElementCentreX-Coordinate"] + rangeM) & (df["y"] > row["ElementCentreY-Coordinate"] - rangeM) & (df["y"] < row["ElementCentreY-Coordinate"] + rangeM)]
                data = np.array(
                    [
                        df_filt['x'],
                        df_filt['y'],
                        df_filt[vcol]
                    ]
                    ).transpose()
                
                if varmdl_method == "Custom":
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
                elif varmdl_method == "Spherical":
                    OK = OrdinaryKriging(
                        data[:, 0],
                        data[:, 1],
                        data[:, 2],
                        variogram_model="spherical",
                        anisotropy_angle = angOK,
                        variogram_parameters = dic,
                        anisotropy_scaling = ratio,
                        exact_values  = True,
                        verbose=True,
                        enable_plotting=False,
                    ) # scaling is the ratio between the major and minor directions' ranges
    
                
                z, ss = OK.execute("points", row['ElementCentreX-Coordinate'], row['ElementCentreY-Coordinate'], backend = "loop")
                ENC_df.loc[ENC_df.ElementNumber == row["ElementNumber"], vcol] = z
    ENC_df = ENC_df.sort_values(by=['MaterialGroup', 'Element Number'])   
    
    # #############################       Write apertures and conductivities .txt files      ############################# 
    # Write apertures.txt
    outputfilepath = pathMSH + MSHinputfilename + "_Apertures_v4" + output_extension
    fringe_aperture = ENC_df[vcol].max() if fringe_aperture == "max" else fringe_aperture
    ENC_df.loc[ENC_df['MaterialGroup'] != FractureMaterialGroup, vcol] = fringe_aperture
    with open(outputfilepath, 'w') as output_file:        
        output_file.write("#MEDIUM_PROPERTIES_DISTRIBUTED\n")	
        output_file.write("$MSH_TYPE\n")
        output_file.write(" LIQUID_FLOW\n")
        output_file.write("$MMP_TYPE\n")
        output_file.write(" GEOMETRY_AREA\n")		
        output_file.write("$DIS_TYPE\n")		
        output_file.write(" ELEMENT\n")		
        output_file.write("$DATA\n")
        for i, element in ENC_df.iterrows():
            if element["MaterialGroup"] == FractureMaterialGroup:
                output_file.write(f"{element['ElementNumber']}\t{element[vcol]}\n")
            else:
                output_file.write(f"{element['ElementNumber']}\t{fringe_aperture}\n")
        output_file.write("$STOP")
    # Write conductivities.txt
    outputfilepath = pathMSH + MSHinputfilename + "_Conductivities_v4" + output_extension
    with open(outputfilepath, 'w') as output_file:        
        output_file.write("#MEDIUM_PROPERTIES_DISTRIBUTED\n")	
        output_file.write("$MSH_TYPE\n")
        output_file.write(" LIQUID_FLOW\n")
        output_file.write("$MMP_TYPE\n")
        output_file.write(" GEOMETRY_AREA\n")		
        output_file.write("$DIS_TYPE\n")		
        output_file.write(" ELEMENT\n")		
        output_file.write("$DATA\n")
        for i, element in ENC_df.iterrows():
            if element["MaterialGroup"] == FractureMaterialGroup:
                output_file.write(f"{element['ElementNumber']}\t{eq(element[vcol])}\n")
            else:
                output_file.write(f"{element['ElementNumber']}\t{eq(fringe_aperture)}\n")
        output_file.write("$STOP")
    
elif prediction_method == "Element Averaging":
    # #############################       Find the nodes for each element and calculate its aperture mean      ############################# 
    dfE = dfE.sort_values(by=['MaterialGroup', 'ElementNumber'])
    dfE['eam'] = 0 # Creates a zeros column
    dfE['ElementCentreX'] = 0
    dfE['ElementCentreY'] = 0
    if Etype == "tri":
# =============================================================================
#             for i, element in dfE.iterrows():
#                 e = [element['Node1'], element['Node2'], element['Node3']]; 
#                 ea = [] # Element Apertures
#                 for n in e:
#                     try:
#                         Nx = dfN._get_value(dfN[dfN['NodeTag'] == n].index[0],'x') # I am using dfN[dfN['NodeTag'] == e[0]].index[0] instead of e[0] in case the index of dfN does not match the 'NodeTag' value, in which case we would be sourcing different x,y,z values and later on different aprtures that would not correspond to the element.
#                         Ny = dfN._get_value(dfN[dfN['NodeTag'] == n].index[0],'y')
#                         Nz = dfN._get_value(dfN[dfN['NodeTag'] == n].index[0],'z')
#                         ea.append( df._get_value ( df[   (df['x']==Nx)   &   (df['y']==Ny) ].index[0], 'aperture' ))
#                     except:
#                         d, j = tree.query((dfN.iloc[n]['x'], dfN.iloc[n]['y'], 1), k=1) # using z == 1 because kdtree now requires a z-value and I believe GINA assigns z=1 when a z-value is not provided in the .geo file
#                         ea.append( df.iloc[j]['aperture'] )
#                 dfE.loc[i, 'eam'] =  sum(ea) / len(ea)
#                 print(dfE.loc[i])
# =============================================================================
        print("Doesn't work for tri elements just yet until I can figure out why the winding number functions returns bad results")
    elif Etype == "quad":
        for i, element in dfE.iterrows():
            e = [element['Node1'], element['Node2'], element['Node3'], element['Node4']]; 
            dfN_filt = dfN[dfN.NodeTag. isin(e)]
            df_filt = df[df.x.between(  dfN_filt['x'].min(), dfN_filt['x'].max() )]
            df_filt = df_filt[df_filt.y.between(  dfN_filt['y'].min(), dfN_filt['y'].max() )]
            dfE.loc[i, 'ElementCentreX'] =  df_filt['x'].mean()
            dfE.loc[i, 'ElementCentreY'] =  df_filt['y'].mean()
            if df_filt.size==0:                 # Control in case there are no samples in within the element. We were having issues when using upscaled data because when upscaling the xy df datapoints shrink.
                eam = eamthreshold
            else:
                eam = df_filt[vcol].mean() # Element Apertures
                # if eam < eamthreshold:
                #     eam = eamthreshold
            dfE.loc[i, 'eam'] =  eam
            print(dfE.loc[i])
    dfE.loc[dfE["eam"]<=eamthreshold, "eam"] = nsmall(dfE['eam'].unique(), 2) # Replace eam column <= eamthreshold with the minimum aperture of dfE after all points are average into their respective elements. This is to manage the spread in the apertures and conductivities txt files for OGS not to complain.
    print(dfE)
    fringe_aperture = dfE["eam"].max() if fringe_aperture == "max" else fringe_aperture
    
    dfE.loc[dfE['MaterialGroup'] != FractureMaterialGroup, "eam"] = fringe_aperture
    
# =============================================================================
#         # Doesn't seem to run without the "root." outside the msh_of_SplitFunction()
#         MSH_OFinputfilepath = root.MSHinputfilepath
#         MSH_OFinputfilename = MSH_OFinputfilepath[MSH_OFinputfilepath.rfind("/")+1:MSH_OFinputfilepath.rfind(".")] + ".msh_of"
#         pathMSH_OF = MSH_OFinputfilepath[:MSH_OFinputfilepath.rfind("/")+1]
#         ENC_df, MGs = msh_of_SplitFunction(pathMSH_OF, MSH_OFinputfilename, Etype)
# =============================================================================
    # #############################       Write apertures and conductivities .txt files      ############################# 
    # Write apertures.txt
    outputfilepath = pathMSH + MSHinputfilename + "_Apertures_v4" + output_extension
    with open(outputfilepath, 'w') as output_file:        
        output_file.write("#MEDIUM_PROPERTIES_DISTRIBUTED\n")	
        output_file.write("$MSH_TYPE\n")
        output_file.write(" LIQUID_FLOW\n")
        output_file.write("$MMP_TYPE\n")
        output_file.write(" GEOMETRY_AREA\n")		
        output_file.write("$DIS_TYPE\n")		
        output_file.write(" ELEMENT\n")		
        output_file.write("$DATA\n")
        for mg in np.sort(dfE["MaterialGroup"].unique() ): # For each Material Group
            if mg == FractureMaterialGroup:
                for i, element in dfE.iterrows():
                    if element["MaterialGroup"] == FractureMaterialGroup:
                        output_file.write(f"{element['ElementNumber']}\t{element['eam']}\n")
                    else:
                        output_file.write(f"{element['ElementNumber']}\t{fringe_aperture}\n")
            else:
                for i, element in dfE.iterrows():
                    if element["MaterialGroup"] != FractureMaterialGroup:
                        output_file.write(f"{element['ElementNumber']}\t{fringe_aperture}\n")
        
        output_file.write("$STOP")
    
    # Write conductivities.txt
    outputfilepath = pathMSH + MSHinputfilename + "_Conductivities_v4" + output_extension
    with open(outputfilepath, 'w') as output_file:        
        output_file.write("#MEDIUM_PROPERTIES_DISTRIBUTED\n")	
        output_file.write("$MSH_TYPE\n")
        output_file.write(" LIQUID_FLOW\n")
        output_file.write("$MMP_TYPE\n")
        output_file.write(" GEOMETRY_AREA\n")		
        output_file.write("$DIS_TYPE\n")		
        output_file.write(" ELEMENT\n")		
        output_file.write("$DATA\n")
        for mg in np.sort(dfE["MaterialGroup"].unique() ): # For each Material Group
            if mg == FractureMaterialGroup:
                for i, element in dfE.iterrows():
                    if element["MaterialGroup"] == FractureMaterialGroup:
                        output_file.write(f"{element['ElementNumber']}\t{eq(element['eam'])}\n")
                    else:
                        output_file.write(f"{element['ElementNumber']}\t{eq(fringe_aperture)}\n")
            else:
                for i, element in dfE.iterrows():
                    if element["MaterialGroup"] != FractureMaterialGroup:
                        output_file.write(f"{element['ElementNumber']}\t{eq(fringe_aperture)}\n")                   
        output_file.write("$STOP")
        output_file.write(f"// Averaging method: {prediction_method}")
        output_file.write(f"// rangeM: {rangeM}, rangem: {rangem}, angle: {angle}, variogram_model: {variogram_model}, sill: {sill}, nugget{nugget}")
        
dfE.to_csv( pathMSH + MSHinputfilename + "_dfE.csv")
dfN.to_csv( pathMSH + MSHinputfilename + "_dfN.csv")
ENC_df.to_csv( pathMSH + MSHinputfilename + "_ENC_df.csv")
# =============================================================================
# """
# Plot for QC
# """
# # Plot Nodes, Element centres and base points
# plt.scatter(ENC_df['Element Centre x-Coordinate'], ENC_df['Element Centre y-Coordinate'], c = "blue", label = "Element centres")
# plt.scatter(df['x'],df['y'], c = "orange", label = "Aperture data points")
# plt.scatter(dfN['x'],dfN['y'], c = "red", s = 5, label = "Nodes")
# plt.xlim(0,5)
# plt.ylim(-20, -15)
# plt.legend()
# plt.xlabel('X (mm)'); plt.ylabel('Y (mm)')
# plt.tight_layout()
# plt.show()
# # Plot kriged aperture at element centres
# sc = plt.scatter(ENC_df['Element Centre x-Coordinate'], ENC_df['Element Centre y-Coordinate'], c = ENC_df['aperture'], cmap=cmap, label = "Element centres")
# plt.colorbar(sc)
# plt.legend()
# plt.xlabel('X (mm)'); plt.ylabel('Y (mm)')
# plt.tight_layout()
# plt.show()
# =============================================================================

# =============================================================================
#          # Write README.txt
#         treeN = KDTree(dfN[["x","y", "z"]].to_numpy()) 
#         POLYLINESpath = root.MSHinputfilepath[:root.MSHinputfilepath.rfind(".")] + "_POLYLINES" + ".csv"
#         PLdf = pd.read_csv(POLYLINESpath, header = 0)
#         readmepath = root.MSHinputfilepath[:root.MSHinputfilepath.rfind(".")] + "_README" + ".csv"
#         PLdf["Node1"] = np.NaN
#         PLdf["Node2"] = np.NaN
#         PLdf["MidNode"] = np.NaN
#         with open(readmepath, 'w') as readme_file:  
#             readme_file.write(f"{','.join([i for i in PLdf.columns])}\n")
#             for i, pl in PLdf.iterrows():
#                 try:
#                     PLdf.iloc[i, "Node1"] = int(dfN[(dfN["x"] == PLdf.iloc[i]["Point1x"]) & (dfN["y"] == PLdf.iloc[i]["Point1y"])]["NodeTag"])
#                     PLdf.iloc[i, "Node2"] = int(dfN[(dfN["x"] == PLdf.iloc[i]["Point2x"]) & (dfN["y"] == PLdf.iloc[i]["Point2y"])]["NodeTag"])
#                     PLdf.iloc[i, "MidNode"] = int(dfN[(dfN["x"] == PLdf.iloc[i]["MidPointx"]) & (dfN["y"] == PLdf.iloc[i]["MidPointy"])]["NodeTag"])
#                 except:
#                     d1, j1 = treeN.query((PLdf.iloc[i]["Point1x"], PLdf.iloc[i]["Point1y"], 0), k=1)
#                     d2, j2 = treeN.query((PLdf.iloc[i]["Point2x"], PLdf.iloc[i]["Point2y"], 0), k=1)
#                     dm, jm = treeN.query((PLdf.iloc[i]["MidPointx"], PLdf.iloc[i]["MidPointy"], 0), k=1)
#                     PLdf.loc[i, "Node1"] = int(dfN.iloc[j1]["NodeTag"])
#                     PLdf.loc[i, "Node2"] = int(dfN.iloc[j2]["NodeTag"])
#                     PLdf.loc[i, "MidNode"] = int(dfN.iloc[jm]["NodeTag"])
#                 readme_file.write(f"{', '.join([str(s) for s in PLdf.iloc[i]])}\n")
#             print(PLdf)
#                 
#                 
#                 
#             readme_file.write(f"Aperture Minimum: {dfE['eam'].min()}\n")
#             readme_file.write(f"Aperture Maximum: {dfE['eam'].max()}\n")
#             readme_file.write(f"Aperture Mean: {dfE['eam'].mean()}\n")
# =============================================================================
        # plt.hist(df['aperture'], bins = 1000)
# =============================================================================
#         # ############################# For testing #############################
#         lc = 1
#         header = 0
#         readcols = [0,1,2]
#         skiprows = ""
#         inputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture\upscale\Extracted M CNL_Aperture_LR_Upscaled1_aperture_FractureMap.csv').as_posix()
#         outputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Freiberg Gneiss\MidSquare\aperture\upscale\Extracted M CNL_Aperture_LR_Upscaled1_aperture_FractureMap_APERTURES.txt').as_posix()
#         # Import necessary data
#         df = pd.read_csv(inputfilepath , usecols = readcols, header = header); dfcolnames = list(df.columns)
#         dfN = pd.read_csv(pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture_Upscaled1\Extracted M CNL_Aperture_LR_Upscaled1_aperture_FractureMap Gmsh MSHNodes.txt').as_posix(), sep=' ', skiprows = 2, skipfooter = 1, usecols = [0,1,2,3], names = ["NodeTag","x","y","aperture"]) # reads only the nodes
#         # Read Elements differently because it changes the number of columns in the middle of the file and confusespd.read_csv
#         f = open(pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture_Upscaled1\Extracted M CNL_Aperture_LR_Upscaled1_aperture_FractureMap.msh').as_posix(), 'r');
#         data_input = f.readlines();
#         f.close()
#         lE = []
#         for iline, line in enumerate(data_input):
#             if any(word in line for word in words):
#                 continue
#             elif str(line).isnumeric():
#                 l =  int( line[:line.rfind("\n")].split(" ")[0] )
#             else:
#                 tline = line[:line.rfind("\n")].split(" ") #trimmed line
#                 e = [int(tline[0]), int(tline[-3]), int(tline[-2]), int(tline[-1])]
#                 lE.append(e)
#         dfE = pd.DataFrame( lE, columns = ["ElementTag", "Node1", "Node2", "Node3"] )
#         
#         
#         input_file = open(inputfilepath, 'r')
#         output_file  = open(outputfilepath, 'w') 
#         data_input = input_file.readlines()["MidPointx"]
#         # #############################
# =============================================================================
        
        

# =============================================================================
#         # ============================================================================= Section for files with $MeshFormat 4.1 0 8 (version, filetype and datasize) =============================================================================
#         # Read Nodes .txt file and create pd.df
#         i=-1
#         while i < len(lines):
#             i+=1
#             if lines[i].find("$Nodes") != -1:
#                 Nodeslineno = i; print("yes")
#                 
#             elif i == Nodeslineno + 1:
#                 line = lines[i][:lines[i].rfind("\n")] # Gets rid of the \n at the end 
#                 s = [int(n) for n in line.split(" ")]; numEntityBlocks = s[0]; numNodes = s[1]; minNodeTag = s[2]; maxNodeTag = s[3]; print(s);
#                 
#             elif lines[i].find("$EndNodes") = -1:
#                 line = lines[i][:lines[i].rfind("\n")] # Gets rid of the \n at the end 
#                 s = [int(n) for n in line.split(" ")];  entityDim = s[0]; entityTag = s[1]; parametric = s[2]; numNodesInBlock = s[3];
#                 for j in range(1, numNodesInBlock * 2 + 1): # *2 because it lists each node tag and then lists each node location. +1 to accommodate for python's 0 start
#                     ListNodesInBlock = []
#                     ListNodesCoordsInBlock = []
#                     if j <= numNodesInBlock:
#                         line= lines[i+j][:lines[i+j].rfind("\n")]
#                         nodeTag = [int(n) for n in line.split(" ")][0]; ListNodesInBlock.append(nodeTag);
#                     else:            
#                     line= lines[i+j][:lines[i+j].rfind("\n")]
#                     nodeCoords = [int(n) for n in line.split(" ")];  nodeX = nodeCoords[0]; nodeY = nodeCoords[1]; nodeZ = nodeCoords[2]; ListNodesCoordsInBlock.append(nodeCoords);
#                 
#                 i+=numNodesInBlock * 2 + 1
#         # =============================================================================
# =============================================================================

# =============================================================================
#         #Submit input and output files' contents (.csv and .geo, respectively) onto the respective text boxes for review by the user
#         CSV_file_Text.delete(0.0, END)
#         CSV_file_Text.insert(0.0, data_input)
#         MSH_file_Text.delete(0.0, END)
#         MSH_file_Text.insert(0.0, data_output)
# 
#         #disable text boxes to edit (As these edits won't be replicated in the files)
#         CSV_file_Text.config(state=DISABLED)
#         MSH_file_Text.config(state=DISABLED) 
# =============================================================================
