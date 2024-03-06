# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:28:20 2023
Created on Mon May  8 10:32:54 2023
Created on Wed Apr 26 11:46:02 2023
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
from msh_of_SplitFunction_v3 import msh_of_SplitFunction_v3
import matplotlib.pyplot as plt
from pykrige.ok import OrdinaryKriging
import sys
import threading
'''User Defined Global Variables'''
output_extension = ".txt"
eamthreshold = 1E-10
interval = 1000 # Use 1 if no interval is wanted

"""System Variables"""
cmap = plt.cm.plasma                    # color map
'''Tkinter Variables'''
root = Tk()
BROWSE_frame = LabelFrame(root, text = "Browse Frame")
BROWSE_frame.grid(row = 1, column = 1, columnspan = 3, sticky = EW, pady = 10, padx = 10)
CSV_frame = LabelFrame(root, highlightbackground="black", highlightthickness=1)
CSV_frame.grid(row = 3, column = 1, pady = 20, padx = 10)
MSH_frame = Frame(root, highlightbackground="black", highlightthickness=1)
MSH_frame.grid(row = 3, column = 3, pady = 20, padx = 10)
PRED_frame = LabelFrame(root, text = "Prediction Frame")
PRED_frame.grid(row = 2, column = 1, pady = 20, padx = 10)
VARMODEL_frame = LabelFrame(root, text = "Variogram Model Frame (ignore if Element Averaging is selected in the Prediction Frame)")
VARMODEL_frame.grid(row = 2, column = 2, pady = 20, padx = 10)
prediction_method = StringVar(value = "Element Centre Kriging")
conductivity_method = StringVar(value = "Cubic Law")
varmdl_method = StringVar(value = "Spherical")

'''Methods'''
def eq(aperture):
    if conductivity_method.get() == "Cubic Law":
        k = (aperture**2)/12; conductivity_method_name = "Cubic Law";
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

'''Tkinter Methods'''
def prediction_method_selection():
    prediction_method_text = f"**{str(prediction_method.get())}** Prediction method selected"
    prediction_method_label.config(text = prediction_method_text)
def conductivity_method_selection():
    conductivity_method_text = f"**{str(conductivity_method.get())}** Conductivity calculation function selected"
    conductivity_method_label.config(text = conductivity_method_text)
def varmdl_method_selection():
    varmdl_method_text = f"**{str(varmdl_method.get())}** Variogram Model selected"
    varmdl_method_label.config(text = varmdl_method_text)
def description():
    tkinter.messagebox.showinfo(title = "Description", message = " Simply use the 'Browse' button to browse to the output .CSV file from Gmsh.\n\n Then click 'Format', or alternatively 'File' > 'Format' and a new file will be created in the same directory formatted and ready to input in GINA.")
def browseCSV():
    root.CSVinputfilepath = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".csv files","*.csv"), (".txt files", "*.txt"), ("All files", "*.*") ))
    CSV_Text.delete(0.0, END)
    CSV_Text.insert(0.0, root.CSVinputfilepath)
    CSVinputfilename = root.CSVinputfilepath[root.CSVinputfilepath.rfind("/")+1:]
    pathCSV = root.CSVinputfilepath[:root.CSVinputfilepath.rfind("/")+1]
    CSVinput_extension = root.CSVinputfilepath[root.CSVinputfilepath.rfind("."):]
    try:
        f = open(root.CSVinputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_CSVfile = open(root.CSVinputfilepath, 'r')
        # read the content of the file line by line 
        CSVdata_input = input_CSVfile.readlines()
        #Submit input file's contents(.CSV) onto the respective text box for review by the user
        CSV_file_Text.delete(0.0, END)
        CSV_file_Text.insert(0.0, CSVdata_input)
def browseMSH():
    root.MSHinputfilepath = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".msh files","*.msh"), (".txt files", "*.txt"), (".csv files","*.csv"), ("All files", "*.*") ))
    MSH_Text.delete(0.0, END)
    MSH_Text.insert(0.0, root.MSHinputfilepath)
    MSHinputfilename = root.MSHinputfilepath[root.MSHinputfilepath.rfind("/")+1:root.MSHinputfilepath.rfind(".")]
    pathMSH = root.MSHinputfilepath[:root.MSHinputfilepath.rfind("/")+1]
    MSHinput_extension = root.MSHinputfilepath[root.MSHinputfilepath.rfind("."):]
    try:
        f = open(root.MSHinputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_MSHfile = open(root.MSHinputfilepath, 'r')
        # read the content of the file line by line 
        MSHdata_input = input_MSHfile.readlines()
        #Submit input file's contents(.CSV) onto the respective text box for review by the user
        MSH_file_Text.delete(0.0, END)
        MSH_file_Text.insert(0.0, MSHdata_input)
def createApertnPermtxt():
    
    CSVinputfilename = root.CSVinputfilepath[root.CSVinputfilepath.rfind("/")+1:root.CSVinputfilepath.rfind(".")]
    pathCSV = root.CSVinputfilepath[:root.CSVinputfilepath.rfind("/")+1]
    CSVinput_extension = root.CSVinputfilepath[root.CSVinputfilepath.rfind("."):]
    MSHinputfilename = root.MSHinputfilepath[root.MSHinputfilepath.rfind("/")+1:root.MSHinputfilepath.rfind(".")]
    pathMSH = root.MSHinputfilepath[:root.MSHinputfilepath.rfind("/")+1]
    MSHinput_extension = root.MSHinputfilepath[root.MSHinputfilepath.rfind("."):]
    try:
        # path = root.inputfilepath[:root.inputfilepath.rfind("/")+1]
        header = int(header_Entry.get())
        readcols = [int(i) for i in readcols_Entry.get().split(",")]
        skiprows = [i for i in skiprows_Entry.get().split(",")]
        fringe_aperture = simplest_type(fringeaperture_Entry.get())
        '''Global Variables'''
        FractureMaterialGroup = simplest_type(FractureMaterialGroup_Entry.get())
        rangeM = simplest_type(rangeM_Entry.get())
        rangem = simplest_type(rangem_Entry.get())
        angle = simplest_type(angle_Entry.get())
        variogram_model = simplest_type(varmdl_method.get())
        sill = simplest_type(sill_Entry.get())
        nugget = simplest_type(nugget_Entry.get())
        dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
        angOK=simplest_type(-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html
        angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
        ratio = rangem/rangeM
        param = [sill, rangeM, nugget] #sill, range, nugget
        if "" in skiprows:
            df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header)
        else:
            skiprows = [simplest_type(i) for i in skiprows]
            df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header, skiprows = skiprows)           #  
    except:
        tkinter.messagebox.showerror(title = "Error", message = "There was an error reading the input file.\nThe most likely reason is that the file does not conform with the required format.\n\nMake sure the file has only 3 columns (x, y and aperture, respectively) and a header on row 1.")
    finally:
        if "" in skiprows:
            df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header)
        else:
            skiprows = [int(i) for i in skiprows]
            df = pd.read_csv(pathCSV + CSVinputfilename + CSVinput_extension, usecols = readcols, header = header, skiprows = skiprows)           #    
        
        dfcolnames = df.columns
        xcol, ycol, vcol = dfcolnames
        df = df.astype({xcol: np.float32, ycol: np.float32, vcol: np.float32})
        FractureMaterialGroup = simplest_type(FractureMaterialGroup_Entry.get())
        rangeM = simplest_type(rangeM_Entry.get())
        rangem = simplest_type(rangem_Entry.get())
        angle = simplest_type(angle_Entry.get())
        partition_interval = simplest_type(partition_Entry.get())
        variogram_model = simplest_type(varmdl_method.get())
        sill = simplest_type(sill_Entry.get())
        nugget = simplest_type(nugget_Entry.get())
        dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
        angOK=simplest_type(-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html
        angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
        ratio = rangem/rangeM
        param = [sill, rangeM, nugget] #sill, range, nugget
        
        
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
#         MSHinputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail-4x4\Hfringes\wMIDPOINTS\FB').as_posix()
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
        
        
        dfN, dfE, ENC_df, MGs, Etype = msh_of_SplitFunction_v3(pathMSH, MSHinputfilename)
        
        
        if prediction_method.get() == "Element Centre Kriging":
            
            # =====================================================================================================================================================
            print("Kriging...")
            # =====================================================================================================================================================
            df2 = df.astype({"x": np.float16, "y": np.float16, vcol: np.float16})
            data = np.array(
                [
                    df2['x'],
                    df2['y'],
                    df2[vcol]
                ]
                ).transpose()
            if varmdl_method.get() == "Custom":
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
                    enable_plotting=False,
                ) # scaling is the ratio between the major and minor directions' ranges
            elif varmdl_method.get() == "Spherical":
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
    
            with Spinner():
                
                if partition_interval == 1:     # Krig the elements' centres in chunks if decimation is != 1
                    z, ss = OK.execute("points", ENC_df['Element Centre x-Coordinate'], ENC_df['Element Centre y-Coordinate'], backend = "loop")
                else:
                    print(f"Kriging in intervals of {interval}...")
                    z = np.ndarray((0))
                    ss = np.ndarray((0))
                    intervals = [[list(ENC_df.index)[i:i + partition_interval][0], list(ENC_df.index)[i:i + partition_interval][-1]] for i in range(0, len(ENC_df), partition_interval)] # Calculate a list of lists with the intervals
                    for i, (n, m) in enumerate(intervals):
                        print(f"Kriging interval [{n}/{m}]...")
                        z_concat, ss_concat = OK.execute("points", ENC_df['Element Centre x-Coordinate'][n:m], ENC_df['Element Centre y-Coordinate'][n:m], backend = "loop")
                        z = np.concatenate ((z, z_concat))        
                # z is the value, ss should stand for sigma squared which is the variance. Backend options are vectorized or loop: vectorized faster but memory intensive whereas loop is slower but also less memory-intensive.
                ###############################################################################
                ENC_df[vcol] = z
                ENC_df = ENC_df.sort_values(by=['Material Group', 'Element Number'])
                # #############################       Write apertures and conductivities .txt files      ############################# 
                # Write apertures.txt
                outputfilepath = pathMSH + MSHinputfilename + "_Apertures_v4" + output_extension
                
                fringe_aperture = ENC_df[vcol].max() if fringe_aperture == "max" else fringe_aperture
                ENC_df.loc[ENC_df['Material Group'] != FractureMaterialGroup, vcol] = fringe_aperture
                
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
                        if element["Material Group"] == FractureMaterialGroup:
                            output_file.write(f"{element['Element Number']}\t{element[vcol]}\n")
                        else:
                            output_file.write(f"{element['Element Number']}\t{fringe_aperture}\n")
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
                        if element["Material Group"] == FractureMaterialGroup:
                            output_file.write(f"{element['Element Number']}\t{eq(element[vcol])}\n")
                        else:
                            output_file.write(f"{element['Element Number']}\t{eq(fringe_aperture)}\n")
                    output_file.write("$STOP")
            
        elif prediction_method.get() == "Element Averaging":
            # #############################       Find the nodes for each element and calculate its aperture mean      ############################# 
            dfE = dfE.sort_values(by=['Material Group', 'Element Number'])
            dfE['eam'] = 0 # Creates a zeros column
            dfE['Element Centre x'] = 0
            dfE['Element Centre y'] = 0
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
                    dfE.loc[i, 'Element Centre x'] =  df_filt['x'].mean()
                    dfE.loc[i, 'Element Centre y'] =  df_filt['y'].mean()
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
            
            dfE.loc[dfE['Material Group'] != FractureMaterialGroup, "eam"] = fringe_aperture
            
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
                for mg in np.sort(dfE["Material Group"].unique() ): # For each Material Group
                    if mg == FractureMaterialGroup:
                        for i, element in dfE.iterrows():
                            if element["Material Group"] == FractureMaterialGroup:
                                output_file.write(f"{element['Element Number']}\t{element['eam']}\n")
                            else:
                                output_file.write(f"{element['Element Number']}\t{fringe_aperture}\n")
                    else:
                        for i, element in dfE.iterrows():
                            if element["Material Group"] != FractureMaterialGroup:
                                output_file.write(f"{element['Element Number']}\t{fringe_aperture}\n")
                
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
                for mg in np.sort(dfE["Material Group"].unique() ): # For each Material Group
                    if mg == FractureMaterialGroup:
                        for i, element in dfE.iterrows():
                            if element["Material Group"] == FractureMaterialGroup:
                                output_file.write(f"{element['Element Number']}\t{eq(element['eam'])}\n")
                            else:
                                output_file.write(f"{element['Element Number']}\t{eq(fringe_aperture)}\n")
                    else:
                        for i, element in dfE.iterrows():
                            if element["Material Group"] != FractureMaterialGroup:
                                output_file.write(f"{element['Element Number']}\t{eq(fringe_aperture)}\n")                   
                output_file.write("$STOP")
                output_file.write(f"// Averaging method: {prediction_method}")
                output_file.write(f"// rangeM: {rangeM}, rangem: {rangem}, angle: {angle}, variogram_model: {variogram_model}, sill: {sill}, nugget{nugget}")
                
        dfE.to_csv( pathMSH + MSHinputfilename + "_dfE.csv")
        dfN.to_csv( pathMSH + MSHinputfilename + "_dfN.csv")
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

'''Labels'''
#Create labels
CSV_Label = Label(BROWSE_frame, text = ".CSV file", pady = 5)
MSH_Label = Label(BROWSE_frame, text = ".MSH file", pady = 5)
header_Label = Label(BROWSE_frame, text = "Header row:", pady = 5)
readcols_Label = Label(BROWSE_frame, text = "Columns to read (separete by comma):", pady = 5)
skiprows_Label = Label(BROWSE_frame, text = "Rows to skip (separete by comma):\nLeave blank if no rows to skip.", pady = 5)
csv_file_label = Label(CSV_frame, text = ".CSV file content", pady = 5)
msh_file_label = Label(MSH_frame, text = ".MSH file content", pady = 5)

fringeaperture_Label = Label(PRED_frame, text = "Fringe Aperture:\nLeave blank if using element average.", pady = 5)
FractureMaterialGroup_label = Label(PRED_frame, text = "Fracture Material Group", pady = 5)
rangeM_label = Label(VARMODEL_frame, text = "Major continuity direction range", pady = 5)
rangem_label = Label(VARMODEL_frame, text = "Minor continuity direction range", pady = 5)
angle_label = Label(VARMODEL_frame, text = "Major contiuity direction angle (CCW from North)", pady = 5)
variogram_model_label = Label(VARMODEL_frame, text = "Variogram_model (spherical)", pady = 5)
partition_label = Label(VARMODEL_frame, text = "Partition (if no partition required, use 1)", pady = 5)
sill_label = Label(VARMODEL_frame, text = "Sill", pady = 5)
nugget_label = Label(VARMODEL_frame, text = "Nugget", pady = 5)

prediction_method_label = Label(PRED_frame, pady = 5, text = "Prediction Method")
conductivity_method_label = Label(PRED_frame, pady = 5, text = "Conductivity Function")
varmdl_method_label = Label(VARMODEL_frame, pady = 5, text = "Variogram Model")

'''Text boxes and scrolls'''
#Create scrollbars
CSVscroll = Scrollbar(CSV_frame)
MSHscroll = Scrollbar(MSH_frame)
#Create text boxes
CSV_Text = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
MSH_Text = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
header_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); header_Entry.insert(-1, "0")
readcols_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); readcols_Entry.insert(-1, "1,2,3")
skiprows_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); skiprows_Entry.insert(-1, "")
CSV_file_Text = Text(CSV_frame, height=20, state=NORMAL, yscrollcommand = CSVscroll.set)
MSH_file_Text = Text(MSH_frame, height=20, state=NORMAL, yscrollcommand = MSHscroll.set)

fringeaperture_Entry = Entry(PRED_frame, width = 10, state=NORMAL); fringeaperture_Entry.insert(-1, "max")
FractureMaterialGroup_Entry = Entry(PRED_frame, width = 10, state=NORMAL); FractureMaterialGroup_Entry.insert(-1, "0")
rangeM_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangeM_Entry.insert(-1, "5") #rangeM_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangeM_Entry.insert(-1, "12")
rangem_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangem_Entry.insert(-1, "2") #rangem_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangem_Entry.insert(-1, "8")
angle_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); angle_Entry.insert(-1, "157.5") #angle_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); angle_Entry.insert(-1, "135")
variogram_model_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); variogram_model_Entry.insert(-1, "spherical")
partition_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); partition_Entry.insert(-1, 1)
sill_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); sill_Entry.insert(-1, "1")
nugget_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); nugget_Entry.insert(-1, "0")
#Configure scrollbar
CSVscroll.config(command=CSV_file_Text.yview)
MSHscroll.config(command=MSH_file_Text.yview)

'''Buttons'''
BrowseCSVButt= Button(BROWSE_frame, text="Browse .CSV file", fg="black", font=("Ariel", 9, "bold"), command=browseCSV)
BrowseMSHButt= Button(BROWSE_frame, text="Browse .MSH file", fg="black", font=("Ariel", 9, "bold"), command=browseMSH)
FormatButt= Button(BROWSE_frame, text="Create Aperture & Permeabilty .txt files", fg="black", font=("Ariel", 9, "bold"), command=createApertnPermtxt)
"Radio Buttons"
avg_RButt = Radiobutton(PRED_frame, text = "Element Averaging", variable = prediction_method, value = "Element Averaging", command = prediction_method_selection)
krig_RButt = Radiobutton(PRED_frame, text = "Element Centre Kriging", variable = prediction_method, value = "Element Centre Kriging", command = prediction_method_selection)
conductivity_RButt = Radiobutton(PRED_frame, text = "Cubic Law", variable = conductivity_method, value = "Cubic Law", command = conductivity_method_selection)
sph_varmdl_RButt = Radiobutton(VARMODEL_frame, text = "Spherical", variable = varmdl_method, value = "Spherical", command = varmdl_method_selection)


'''Allocate widgets'''
# BROWSE_frame
BrowseCSVButt.grid(row = 0, column = 2, sticky = W)
BrowseMSHButt.grid(row = 0, column = 5, sticky = W)
header_Label.grid(row = 1, column = 0, sticky = W)
header_Entry.grid(row = 1, column = 1, sticky = W)
readcols_Label.grid(row = 2, column = 0, sticky = W)
readcols_Entry.grid(row = 2, column = 1, sticky = W)
skiprows_Label.grid(row = 3, column = 0, sticky = W)
skiprows_Entry.grid(row = 3, column = 1, sticky = W)
FormatButt.grid(row = 0, column = 6, sticky = W)
# MSH_frame
MSH_Label.grid(row = 0, column = 3, sticky = W)
MSH_Text.grid(row = 0, column = 4, sticky = W)
msh_file_label.grid(row = 1, column = 1, sticky = W)
MSH_file_Text.grid(row = 2, column = 1)
MSHscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
#CSV_frame
CSV_Label.grid(row = 0, column = 0, sticky = W)
CSV_Text.grid(row = 0, column = 1, sticky = W)
csv_file_label.grid(row = 1, column = 1, sticky = W)
CSV_file_Text.grid(row = 2, column = 1)
CSVscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
# PRED_frame
FractureMaterialGroup_label.grid(row = 0, column = 0, sticky = W)
FractureMaterialGroup_Entry.grid(row = 0, column = 1, sticky = W)
prediction_method_label.grid(row = 1, column = 0, sticky = W)
avg_RButt.grid(row = 1, column = 1, sticky = W)
krig_RButt.grid(row = 2, column = 1, sticky = W)
conductivity_method_label.grid(row = 3, column = 0, sticky = W)
conductivity_RButt.grid(row = 3, column = 1, sticky = W)
fringeaperture_Label.grid(row = 4, column = 0, sticky = W)
fringeaperture_Entry.grid(row = 4, column = 1, sticky = W)
#VARMODEL_frame
partition_label.grid(row = 0, column = 0, sticky = W)
varmdl_method_label.grid(row = 1, column = 0, sticky = W)
rangeM_label.grid(row = 2, column = 0, sticky = W)
rangem_label.grid(row = 3, column = 0, sticky = W)
angle_label.grid(row = 4, column = 0, sticky = W)
sill_label.grid(row = 5, column = 0, sticky = W)
nugget_label.grid(row = 6, column = 0, sticky = W)
partition_Entry.grid(row = 0, column = 1, sticky = W)
sph_varmdl_RButt.grid(row = 1, column = 1, sticky = W)
rangeM_Entry.grid(row = 2, column = 1, sticky = W)
rangem_Entry.grid(row = 3, column = 1, sticky = W)
angle_Entry.grid(row = 4, column = 1, sticky = W)
sill_Entry.grid(row = 5, column = 1, sticky = W)
nugget_Entry.grid(row = 6, column = 1, sticky = W)

'''Menu bar'''
menubar = Menu(root)
root.config(menu = menubar)

file_menu = Menu(menubar)
menubar.add_cascade(label = "File", menu = file_menu)
file_menu.add_command(label = "Create Aperture & Permeability .txt files", command = createApertnPermtxt)
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

options_menu = Menu(menubar)
menubar.add_cascade(label = "Options", menu = options_menu)
options_menu.add_command(label = "Description", command = description)

'''Prompt on open'''
#Uncomment below if you wish the description() function (info box explaining how the app works) to be prompted at app opening
#description()

'''Add Window title, CSVmetry and create Window's main loop'''
root.title("Tkinter - CSV and GINA .msh Aperture Mapping v4")
root.geometry("2000x800")
root.mainloop()