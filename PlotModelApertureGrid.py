# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:09:39 2023

@author: s2132627

Single Fracture OGS Model plotting

Requires .geo file to get the grid configuration (x-elements x y-elements) to be used in the matplotlib.imshow function

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
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import math
from ast import literal_eval

import time
from time import ctime
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
extension_txt = ".txt"
extension_geo = ".geo"

from bisect import bisect
import numpy as np
cmap = plt.cm.plasma                    # color map


#set working directory and filename
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Manuscript\CodeTest\Tkinter_GmshGINA_UpscaleMapping\Testing').as_posix() #flips the backslashes to plug in open()
inputfilename = pathlib.PureWindowsPath(r'AttributePointCloud').as_posix()
wdmsh = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Manuscript\CodeTest\Tkinter_GmshGINA_UpscaleMapping\Testing').as_posix() #flips the backslashes to plug in open()
inputfilename_msh = pathlib.PureWindowsPath(r'OGSGINAmsh').as_posix()

bs="\\"; wd=wd+bs                               # Use this instead in Windows
wdmsh=wdmsh+bs
wdsave=wd 

save = False
plot = True

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
def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s
def msh_of_SplitFunction_v2(path, input_file_name):
    # #############################  Get the Nodes and Elements from the .msh (GINA's) file without having to export and read back .txt files #############################
    words = ["$NODES", "$ELEMENTS"]
    with open(path + input_file_name + ".msh", 'r') as file:
        # This Loop is not working in tkinter for whatever fucking reason so I hardcoded it
        # for word in words: 
        #     exec(f"{word[word.rfind('$')+1:]} = False") # Creates a False flag for each word in words
        NODES = False
        ELEMENTS = False
        for line in file:
            match = next((word for word in words if word in line), False)
            if match:
                print(line)
                # for word in words: # Creates a False flag for each word in words  
                #     exec(f"{word[word.rfind('$')+1:]} = False")
                NODES = False
                ELEMENTS = False
                # exec(f"{match[match.rfind('$')+1:]} = True") # Sets the match word flag to True
                print(match)
                print(NODES)
                print(ELEMENTS)
                if match.find("NODES") != -1: # If line contains "NODES"
                    NODES = True
                    Nnodes = int(next(file))
                    print(Nnodes)
                    Ns = []
                elif match.find("ELEMENTS") != -1: # If line contains "ELEMENTS"
                    ELEMENTS = True
                    Nelements = int(next(file))
                    print(Nelements)
                    Es = []
            else:
                if line.find("#STOP") != -1:
                    break
                elif NODES == True:
                    Ns.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                elif ELEMENTS == True:
                    Es.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
    # Get E type (tri, quad)
    Etype = Es[0][2]
    dfN = pd.DataFrame(Ns, columns = ["NodeTag","x","y","z"]) # reads only the nodes
    # dfN = dfN.round(prec(df)) # This is necessary because of precision changes, later on when I am trying to find the aperture corresponding to the x,y,z values of a node, those x,y,z won't match.
    if Etype == "tri":
        dfE = pd.DataFrame(Es, columns = ['Element Number', 'Material Group', 'ElementType', 'Node1', 'Node2', 'Node3']); #dfE = dfE.round(prec(df))
    elif Etype == "quad":
        dfE = pd.DataFrame(Es, columns = ['Element Number', 'Material Group', 'ElementType', 'Node1', 'Node2', 'Node3', 'Node4']); #dfE = dfE.round(prec(df))
    EmaxN = max([(simplest_type(n[n.rfind("Node")-1:])) for n in list(dfE.columns) if n.find("Node") != -1])
    print(dfN)
    print(dfE)
    
    
    
    if os.path.isfile(path + input_file_name + ".msh_of"):
        with open(path + input_file_name + ".msh_of") as f:
            contents = f.read()
        f.close()
        
        Esplit = contents.split("\n\n")
        Esplit = [ s.split("\n") for s in Esplit] # Split contents into elements
        ENC = [[simplest_type(l2.strip()) for l2 in l1] for l1 in Esplit] # Elements' strings cleaned and converted into numbers    # Instead of for i, l1 in enumerate(f2split):
                                                                                                                                    #for j, l2 in enumerate(l1):
                                                                                                                                    #    f2split[i][j] = simplest_type(l2.strip())
        if Etype == "tri":
            EsC_columns = ["Material Group", "Element Number", "Node1", "Node2", "Node3", "Element Centre x-Coordinate", "Element Centre y-Coordinate", "Element Centre z-Coordinate", "Dunno", "Dunno2" ] # Elements' Centres
            ENC_df = pd.DataFrame(ENC[1:], columns = EsC_columns) # Start at row1 ignoring row 0 because the latter is not an Element
        elif Etype == "quad":
            EsC_columns = ["Material Group", "Element Number", "Node1", "Node2", "Node3","Node4", "Element Centre x-Coordinate", "Element Centre y-Coordinate", "Element Centre z-Coordinate", "Dunno", "Dunno2" ] # Elements' Centres
            ENC_df = pd.DataFrame(ENC[1:], columns = EsC_columns) # Start at row1 ignoring row 0 because the latter is not an Element
            # ENC_df_ElementCentrexCoordinate_round = round(ENC_df["Element Centre x-Coordinate"])
            ENC_df["grid x-coordinate"]  = floor(round(ENC_df["Element Centre x-Coordinate"]).rank(method = "min") - 1)                                                             # Round because there might be a precision issue. This rounding does not change the initial dataframe. -1 because the ranking function starts at 1
            ENC_df["grid x-coordinate"]  = ENC_df["grid x-coordinate"] / (np.sort(ENC_df["grid x-coordinate"].unique())[1] - np.sort(ENC_df["grid x-coordinate"].unique())[0])      # rank function assigns the same rank to duplicates but jumps that amount of duplicates when the next value is different, which we don't want. Thus we are calculating the difference between sequential values and dividing..
            ENC_df["grid y-coordinate"]  = floor(round(ENC_df["Element Centre y-Coordinate"]).rank(method = "min") - 1)
            ENC_df["grid y-coordinate"]  = ENC_df["grid y-coordinate"] / (np.sort(ENC_df["grid y-coordinate"].unique())[1] - np.sort(ENC_df["grid y-coordinate"].unique())[0])
            # ENC_df["grid x-coordinate"]
            # ENC_df["grid y-coordinate"]
    # =============================================================================
    #     # Get N lists, one for each Material Group, containing all Elements' numbers of that Material Group
    #     for mg in MGs:
    #         exec(f"MG{mg}Es =  ENC_df['Element Number'][ ENC_df['Material Group'] == {mg} ].to_list()   ")
    # =============================================================================

    else:      
        ENC_df = dfE
        ENC_df["Element Centre x-Coordinate"] = 0
        ENC_df["Element Centre y-Coordinate"] = 0
        ENC_df["Element Centre z-Coordinate"] = 0
        ENC_df["grid x-coordinate"] = 0
        ENC_df["grid y-coordinate"] = 0
        
        for row in ENC_df.itertuples():
            if Etype == "tri":
                Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3") ]
            if Etype == "quad":
                Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3"), getattr(row, "Node4") ]
            Ns = dfN[dfN["NodeTag"].isin(Ns)]
            
            ENC_df.at[getattr(row, "Index"), 'Element Centre x-Coordinate']  = Ns.x.describe()["mean"]
            ENC_df.at[getattr(row, "Index"), 'Element Centre y-Coordinate']  = Ns.y.describe()["mean"]
            ENC_df.at[getattr(row, "Index"), 'Element Centre z-Coordinate']  = Ns.z.describe()["mean"]
            
        ENC_df["grid x-coordinate"]  = np.argsort(ENC_df["Element Centre x-Coordinate"])
        ENC_df["grid y-coordinate"]  = np.argsort(ENC_df["Element Centre y-Coordinate"])
    MGs = np.sort(ENC_df["Material Group"].unique() )# np.array containing all used Material Groups' numbers
        
    return dfN, dfE, ENC_df, MGs, Etype

def Get_dotgeo_grid_info(path, input_file_name):
    keywords = ["Point(","Point (", "Line(", "Line (", "Transfinite Line(", "Transfinite Line ("]
    keywordsvars = [keyword.replace(" ","").replace("(","") for keyword in keywords]
    # for i, keywordvar in enumerate(keywordsvars):         # Hardcoded because it's not working
    #     exec(f"{keywordvar} = []")
    #     exec(f'{keywordvar}_df = []')
    #     exec(f'{keywordvar}_df')
    P=-1
    L=-1
    TL=-1
    Point_df = []
    Line_df = []
    TransfiniteLine_df = []
    
    with open(path + input_file_name + ".geo", 'r') as file:
        for line in file:
           #print(line)
           line[0:line.find("//")]  # Ignore comments
           #for i, keyword in enumerate(keywords):
           if keywords[0] in line or keywords[1] in line:
               P+=1
               #exec(f"{keywordsvars[0]}.append(line)")
               pointn = simplest_type(line[line.find("(")+1:line.find(")")])
               x = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[0].replace(" ",""))
               y = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[1].replace(" ",""))
               z = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[2].replace(" ",""))
               lc = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[3].replace(" ",""))
               Point_df.append([pointn, x, y, z, lc])
               #print(linen, p1, p2)
               #break
           elif keywords[4] in line or keywords[5] in line:                         # 2 comes first because "Transfinite Line" string contains "Line" string
               TL+=1
               #exec(f"{keywordsvars[2]}.append(line)")
               tlinen = int(line[line.find("(")+1:line.find(")")].replace(" ",""))
               nom = simplest_type(line[line.find("=")+1:line.find("Using Progression")].replace(" ",""))
               denom = simplest_type(line[line.find("Using Progression")+len("Using Progression"):line.find(";")].replace(" ",""))
               TransfiniteLine_df.append([tlinen, nom/denom, ""])
               #break
           elif line.startswith(keywords[2] or keywords[3]):
               L+=1
               #exec(f"{keywordsvars[1]}.append(line)")
               linen = simplest_type(line[line.find("(")+1:line.find(")")].replace(" ",""))
               p1 = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[0].replace(" ",""))
               p2 = simplest_type(line[line.find("{")+1:line.find("}")].split(",")[1].replace(" ",""))
               Line_df.append([linen, p1, p2])
               #break
           
    Line_df = pd.DataFrame(Line_df, columns = ["linen", "point1", "point2"])
    Point_df = pd.DataFrame(Point_df, columns = ["pointn", "x", "y", "z", "lc"])
    TransfiniteLine_df = pd.DataFrame(TransfiniteLine_df, columns = ["tlinen", "elementn", "alongaxis"])
    
    for i, row in TransfiniteLine_df.iterrows():
        row
        p1 = Line_df[Line_df["linen"] == row[0]]["point1"].values[0]
        p2 = Line_df[Line_df["linen"] == row[0]]["point2"].values[0]
        p1x = Point_df[Point_df["pointn"] == p1]["x"].values[0]
        p1y = Point_df[Point_df["pointn"] == p1]["y"].values[0]
        #p1z = Point_df[Point_df["pointn"] == p1]["z"].values[0]
        p2x = Point_df[Point_df["pointn"] == p2]["x"].values[0]
        p2y = Point_df[Point_df["pointn"] == p2]["y"].values[0]
        #p2z = Point_df[Point_df["pointn"] == p2]["z"].values[0]
        alongaxis = ["x" if p1y==p2y else "y"][0]
        TransfiniteLine_df.loc[i, "alongaxis"] = alongaxis
    
    M = TransfiniteLine_df["elementn"].max()
    m = TransfiniteLine_df["elementn"].min()
    mid = TransfiniteLine_df.describe().transpose()["50%"]["elementn"]
    fringeaxis = TransfiniteLine_df[TransfiniteLine_df["elementn"]==m]["alongaxis"].values[0]
    longaxis = TransfiniteLine_df[TransfiniteLine_df["elementn"]==M]["alongaxis"].values[0]
    if fringeaxis == "y":
        if longaxis == "y":
            leny = (M-1)+2*(m-1)
            lenx = mid-1
        elif longaxis == "x":
            leny = (mid-1)+2*(m-1)
            lenx = M-1
    elif fringeaxis == "x":
        if longaxis == "x":
            lenx = (M-1)+2*(m-1)
            leny = mid-1
        elif longaxis == "y":
            lenx = (mid-1)+2*(m-1)
            leny = M-1

    return lenx, leny,  fringeaxis, longaxis

#================================================================================================================================================================================ 
'''###########################################################################################################################################################
################################################################# IMPORT DATAFRAME AND PLOT ##################################################################
###########################################################################################################################################################'''
"""
Import and calculate necessary data
"""
# Apertures_df = pd.read_csv(wd + inputfilename + extension_txt, names = ["ElementNumber", "Aperture"], skiprows=8, sep="\t", skipfooter=1, engine='python'); Apertures_df.reset_index(drop=False)  # For reading the Apertures.txt file     #This new line is because reading the .csv file with a header and an index column will create some issue
Apertures_df = pd.read_csv(wd + inputfilename + extension_csv, names = ["ElementNumber", "Aperture"], skiprows=8, sep="\t", skipfooter=1, engine='python'); Apertures_df.reset_index(drop=False)  # For reading the initial .csv point cloud
Apertures_df = pd.read_csv(wd + inputfilename + extension_csv, header = 0, index_col = 0);
dfN, dfE, ENC_df, MGs, Etype = msh_of_SplitFunction_v4(wdmsh, inputfilename_msh)
lenx, leny,  fringeaxis, longaxis = Get_dotgeo_grid_info(wdmsh, "Gmsh")
ENC_df = pd.merge(ENC_df, Apertures_df, on='ElementNumber')

x_coords = ENC_df['GridX-coordinate']
y_coords = ENC_df['GridY-coordinate']

Z = np.zeros((int(leny),int(lenx)))
for i, (MG, En, Node1, Node2, Node3, Node4, ECx, ECy, ECz, dunno1, dunno2, Gx, Gy, Aperture) in ENC_df.iterrows(): # ECx Element x-coordinate, Gx, Element Grid x-coordinate
    Z[int(Gy)][int(Gx)] = Aperture




"""
Plot Upscaled attribute and mesh
"""
im = plt.imshow(Z, interpolation = None, aspect = 1,
            origin="lower", extent=[0, lenx, 0, leny],
            cmap = cmap) #
fmt = mpl.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))
cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(Z.min(), Z.max(), 10),  format=fmt)       ## contour map, colour bar and colourbar label??? ##
cbar.set_label(f'Aperture (mm)', rotation=90, labelpad=20)
plt.title("Greywacke OpenGeoSys Model")
ax1 = plt.gca()
ax1.set_aspect(1)
ax1.set_xticks( np.arange(0, lenx, 1), minor=True)
ax1.set_yticks(np.arange(0, leny, 1), minor=True)
# ax1.set_xticklabels(np.arange(0, max(x_coords)+1, 20))
# ax1.set_yticklabels(np.arange(0, max(y_coords)+1, 20))
plt.grid(True, which='both', color='black', linewidth=0.25)
plt.gcf().set_dpi(1000)
if save == True:
    plt.savefig(f"{wdsave}ModelAttribute{extension_png}", dpi=1000, bbox_inches = "tight")
if plot == True:
    plt.show()








"""
Plot attribute point cloud and mesh
"""

def norm(val, vmin, vmax):
    nval = (val - vmin)/(vmax-vmin)
    return nval

edited = array([-43., -42., -41., -40., -39., -38., -37., -36., -35., -34.,
       -33., -32., -31., -30., -29., -28., -27., -26., -25., -24., -23.,
       -22., -21., -20., -19., -18., -17., -16., -15., -14., -13., -12.,
       -11., -10.,  -9.,  -8.,  -7.,  -6.,  -5.,  -4.,  -3.,  -2.,  -1.,
        -0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
        11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,  20.,  21.,
        22.,  23.,  24.,  25.,  26.,  27.,  28.,  29.,  30.,  31.,  32.,
        33.,  34.,  35.,  36.,  37.,  38.,  39.,  40.,  41.,  42.,  43.,
        44.,  45.,  46.,  47.,  48.,  49.,  50.,  51.,  52.,  53.,  54.,
        55.,  56.,  57.,  58.,  59.,  60.,  61.,  62.,  63.,  64.,  65.,
        66.,  67., 68.])



vcol = "Aperture"
plt.scatter(Apertures_df['x'], Apertures_df['y'], c=Apertures_df[vcol], cmap=cmap, s=2, marker = ".", alpha=1 ) 
for h in round(dfN['y']).unique():
    plt.axhline(h, 0.01, 0.98, color = "black", linestyle = "-", linewidth = 0.5) 
for v in edited:
    plt.axvline(v, 0.04, 0.97, color = "black", linestyle = "-", linewidth = 0.5) 
# for v in list(round(dfN['x']).unique()):
#     plt.axvline(v, dfN['y'].min(), dfN['y'].max(), color = "black", linestyle = "-", linewidth = 0.5)    
# for h in list(round(dfN['y']).unique()):
#     plt.axhline(h, round(dfN['x'].min()), round(dfN['x'].max()), color = "black", linestyle = "-", linewidth = 0.5)  
plt.colorbar(mappable = None, label = f'Attribute', orientation="vertical", ticks=np.linspace(amin(Apertures_df[vcol]), amax(Apertures_df[vcol]), 10))
# plt.scatter(dfN['x'], dfN['y'], c="black",  s=1, marker = "+", alpha=1 ) 
plt.xlabel(r'$X$', fontsize=15)
plt.ylabel(r'$Y$', fontsize=15)
plt.title(f'Gmsh mesh and Attribute distribution')
if save == True:
    plt.savefig(f"{wdsave}GmshAttributeDist{extension_png}", dpi=1000) #, bbox_inches = "tight")
if plot == True:
    plt.show()
# gridx = np.arange(0, int(lenx), 1) # Based on original points' coordinates
# gridy = np.arange(int(leny), 0, -1)
# xx, yy = np.meshgrid(gridx, gridy)
# data = np.array([  x_coords, y_coords, Z  ]).transpose()


# im = plt.contourf(xx, yy, data, cmap = cmap, vmin = Apertures_df["Aperture"].min(), vmax = Apertures_df["Aperture"].max(), levels = np.linspace(Apertures_df["Aperture"].min(), Apertures_df["Aperture"].max(), 100),)                                               ####################################################
# cbar = plt.colorbar(im, orientation="vertical", ticks=np.linspace(z.min(), z.max(), 10))                                                                              ## contour map, colour bar and colourbar label??? ##
# cbar.set_label(f'{vcol} (mm)', rotation=90, labelpad=20)                                                                                                              ####################################################
# plt.imshow(Apertures_df["Aperture"], interpolation = None, vmin = Apertures_df["Aperture"].min(), vmax = Apertures_df["Aperture"].max(), cmap = cmap) # extent = [dfu['x'].min(), dfu['x'].max(), dfu['y'].min(), dfu['y'].max()]
# plt.title(f"title"); plt.xlabel('X (mm)'); plt.ylabel('Y (mm)')
# plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.2, wspace=0.2, hspace=0.2)
# plt.tight_layout()
# if save == True:
#     plt.savefig(f"{wdsave}ModelAperture{extension_png}", dpi=300, bbox_inches = "tight")
# if plot == True:
#     plt.show()














