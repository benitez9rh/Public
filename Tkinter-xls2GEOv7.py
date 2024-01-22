# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:12:03 2023

@author: s2132627
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 09:19:08 2022

@author: s2132627
"""

import os
import pathlib
from tkinter import *               #Import everything from Tkinter
from tkinter import ttk
import tkinter.messagebox           #Import messagebox method
from tkinter import filedialog
import subprocess
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
from time import time, ctime
import time
from scipy.spatial import KDTree
import math
'''Parameteres'''
fringes = "WE"
use_aperture_as_z = False
# square = True
# Dimensions = 2
tri_Elements = True #Use triangular or square elements
fringe_aperture = 0.005 # type value in m (1mm=0.001m). If None is used, it will grab the aperture from the nearest point
Transfinite = True
h_threshold = 0.1
'''Create main window, call it root'''
root = Tk()
BROWSE_frame = Frame (root)
BROWSE_frame.grid(row = 1, column = 1, columnspan = 3, sticky = EW, pady = 10, padx = 10)
CSV_frame = Frame(root)
CSV_frame.grid(row = 3, column = 1, pady = 20, padx = 10)
GEO_frame = Frame(root)
GEO_frame.grid(row = 3, column = 3, pady = 20, padx = 10)

'''Methods'''
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
        d, j = tree.query(pt, k=[2])                                                                            # d = distance, i = index. k = 2 to get the first and second nearest points because the first nearest point will be the point itself, hence dist = 0.
        dist_df.loc[i,'dist_to_closest_pt'] = d                                                                     # Assign the value of the distance between each point in the loop and it's closest point (apart from itself)
    avg_min_dist = dist_df['dist_to_closest_pt'].mean()
    min_dist = dist_df['dist_to_closest_pt'].min()
    max_dist = dist_df['dist_to_closest_pt'].max()
    print( f"The average minimum distance between closest points is {avg_min_dist}, the minimum distance between closest points is {min_dist} and the maximum distance between closest points  is {max_dist}.")  
    return avg_min_dist, min_dist, max_dist, dist_df
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
        https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def ang3p(p0,p1,p2):    
    """
    Computes the angle for the p0p1p2 corner
    Inputs:
        p0,p1,p2 - points in the form of [x,y]
    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python

    Parameters
    ----------
    p0 : List, tuple
        DESCRIPTION.
    p1 : List, tuple
        DESCRIPTION.
    p2 : List, tuple
        DESCRIPTION.

    Returns
    -------
    radians : Float
        DESCRIPTION.
    degrees : Float
        DESCRIPTION.

    """
    #p0 = [-22.5, -32.5]
    #p1 = [df2['x'].mean(), df2['y'].mean()]
    #p2 = [df2['x'].mean() + 100, df2['y'].mean()]
    
    try:
        len(p0) == len(p1) == len(p2)
        type(p0) == type(p1) == type(p2)
    except:
        print("Error")
    finally:
        if len(p0) < 3:
            v0 = np.array(p0) - np.array(p1)
            v1 = np.array(p2) - np.array(p1)
            
            radians = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
            degrees = np.degrees(radians)
            print(f"Radians: {radians}, Degrees: {degrees}")
        else:
            v0 = np.array(p0) - np.array(p1)
            v1 = np.array(p2) - np.array(p1)
            radians = angle_between(v0, v1)
            degrees = np.degrees(radians)
            print(f"Radians: {radians}, Degrees: {degrees}")
    return radians, degrees

def angle_trunc(a):
    while a < 0.0:
        a += math.pi * 2
    return a

def xaxis_vs_lineseg_angle(p1, p2):
    ''' https://stackoverflow.com/questions/7586063/how-to-calculate-the-angle-between-a-line-and-the-horizontal-axis '''
    deltaY = p2[1] - p1[1]
    deltaX = p2[0] - p1[0]
    rads = angle_trunc(math.atan2(deltaY, deltaX))
    degs = math.degrees(rads)
    return rads, degs
def splitlist(listtosplit, n): # https://stackoverflow.com/questions/1624883/alternative-way-to-split-a-list-into-groups-of-n
    split_list = [listtosplit[i:i+n] for i in range(0, len(listtosplit), n)]
    return split_list

#description() function (info box explaining how the app works)
def description():
    tkinter.messagebox.showinfo(title = "Description", message = " Simply use the 'Browse' button to browse to the output .CSV file from Gmsh.\n\n Then click 'Format', or alternatively 'File' > 'Format' and a new file will be created in the same directory formatted and ready to input in GINA.")

#browse() function to grab the .CSV file and assign respective file path
def browse():
    root.inputfilepath = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".csv files","*.csv"), (".txt files", "*.txt"), ("All files", "*.*") ))
    CSV_Text.delete(0.0, END)
    CSV_Text.insert(0.0, root.inputfilepath)
    try:
        f = open(root.inputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_file = open(root.inputfilepath, 'r')
        # read the content of the file line by line 
        data_input = input_file.readlines()
        #Submit input file's contents(.CSV) onto the respective text box for review by the user
        CSV_file_Text.delete(0.0, END)
        CSV_file_Text.insert(0.0, data_input)
        

def createGEO():
    try:
        f = open(root.inputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        #from input file path, reads the path to the folder where the file is located, file extension, file name and assigns respective variables
        inputfilename = root.inputfilepath[root.inputfilepath.rfind("/")+1:root.inputfilepath.rfind(".")]
        path = root.inputfilepath[:root.inputfilepath.rfind("/")+1]
        input_extension = root.inputfilepath[root.inputfilepath.rfind("."):]
        output_extension = ".geo"
        outputfilepath = root.inputfilepath[:root.inputfilepath.rfind(".")] + output_extension
        polylinesfilepath = root.inputfilepath[:root.inputfilepath.rfind("/")+1]  + inputfilename + "_POLYLINES" + ".csv"
        
        input_file = open(root.inputfilepath, 'r')  # open file in read mode 
        output_file  = open(outputfilepath, 'w')    # open other file in write mode 
        polylines_file = open(polylinesfilepath, 'w')
        data_input = input_file.readlines()     # read the content of the file line by line

        try:
            lc = float(LC_Entry.get())
            header = int(header_Entry.get())
            readcols = [int(i) for i in readcols_Entry.get().split(",")]
            skiprows = [i for i in skiprows_Entry.get().split(",")]
            if "" in skiprows:
                df = pd.read_csv(path + inputfilename + input_extension , usecols = readcols, header = header)
            else:
                skiprows = [int(i) for i in skiprows]
                df = pd.read_csv(path + inputfilename + input_extension , usecols = readcols, header = header, skiprows = skiprows) #   
        except:
            tkinter.messagebox.showerror(title = "Error", message = "There was an error reading the input file.\nThe most likely reason is that the file does not conform with the required format.\n\nMake sure the file has only 3 columns (x, y and aperture, respectively) and a header on row 1.")
        finally:
            lc = float(LC_Entry.get())
            header = int(header_Entry.get())
            readcols = [int(i) for i in readcols_Entry.get().split(",")]
            skiprows = [i for i in skiprows_Entry.get().split(",")]
            if "" in skiprows:
                df = pd.read_csv(path + inputfilename + input_extension, usecols = readcols, header = header)
                print(df)
            else:
                skiprows = [int(i) for i in skiprows]
                df = pd.read_csv(path + inputfilename + input_extension, usecols = readcols, header = header, skiprows = skiprows) # 
            colnames = list(df.columns)
            
            #Check dimensions
            if ("x" in colnames and "y" in colnames and "z" in colnames):
                Dimensions = 3
            elif ("x" in colnames and "y" in colnames):
                Dimensions = 2
            else:
                Dimensions = None
            

# =============================================================================
#             # ############################# For testing
#             lc = 1
#             header = 0
#             readcols = [1,2,3]
#             skiprows = ""
#             inputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail-2x2\Kriged_aperture_Based_on_Extracted_M_CNL_Aperture_Dens2x2.csv').as_posix()
#             outputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail-2x2\Kriged_aperture_Based_on_Extracted_M_CNL_Aperture_Dens2x2.geo').as_posix()
#             readmefilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail-2x2\Kriged_aperture_Based_on_Extracted_M_CNL_Aperture_Dens2x2README.csv').as_posix()
#             vcol = "aperture"
#             if "" in skiprows:
#                 df = pd.read_csv(inputfilepath , usecols = readcols, header = header)
#             else:
#                 skiprows = [int(i) for i in skiprows]
#                 df = pd.read_csv(inputfilepath , usecols = readcols, header = header, skiprows = skiprows) # 
#             df[vcol] = df[vcol]*0.001 # make it meters instead of mm
#             
#             input_file = open(inputfilepath, 'r')
#             output_file  = open(outputfilepath, 'w') 
#             readme_file = open(readmefilepath, 'w') 
#             data_input = input_file.readlines()
#             # #############################
# =============================================================================

            # Calculate mins, maxs, deltas
            ptx = (df['x'].max() - df['x'].min())
            pty = (df['y'].max() - df['y'].min())
            dx = ptx / (len(df['x'].unique())-1)    # Because the points area equidistant, calculate the x-spacing
            dy = pty / (len(df['y'].unique())-1)    # Because the points area equidistant, calculate the y-spacing)
            # Get bound mins and maxs
            bxmin = df['x'].min()
            bxmax = df['x'].max()
            bymin = df['y'].min()
            bymax = df['y'].max()
            surfpts_ctrmass = [df['x'].mean(), df['y'].mean()]

            
            # Calculate convex hull (boundary around points)     
            points =  np.array( list(map(list, zip( df['x'], df['y'] ) )) ) # colnames[0] and colnames[1] should be 'x' and 'y' respectively
            hull = ConvexHull(points, qhull_options = None)
            vindex = tuple(hull.vertices); print(vindex)     # list of df indexes in order for the convex hull ; 
            bindex = list(vindex); bindex.append(bindex[0]); print(bindex) 
            blines = [list(a) for a in zip(bindex[:-1], bindex[1:])]; print(blines) # list comprehension converts tuples from zip into lists #print(line);  
            # Create points on convex hull's vertices and expand them outwards (to create and outer boundary) by an amount
            df2 = df[df.index.isin(vindex)] # Create a df with only the points of vindex but ***** NOT ***** in order

            avg_min_dist, min_dist, max_dist, dist_df = grid_fit(df)
            tree = KDTree(df[["x","y"]].to_numpy())
            
            df["name"] = "" # add an empty string column to pd df
            colnames = list(df.columns)
            blinesN = len(blines)# number of boundary lines
            square = True if blinesN == 4 else False
            surfptsN = len(df) # number of points in original surface
            ip = surfptsN -1# current point index
            fpts = [] # Fringe points
            bmpts = [] # Boundary mid points
            # fptsdf = pd.DataFrame( columns = colnames) # Fringe points dataframe empty
            LsN = blinesN # Number of lines
            Ls = [] # Lines
            Ls_index = [] # Lines index list (for constructing the DataFrame in correct order so that the lines corresponding to the boundary come in first but we have only one main loop)
            LLsN = 0
            LLs = [] # Line Loops
            LLs_index = []
            
            # Create DataFrames containing all the points, lines, lineloops ,etc
            LLsN +=1
            LLs_index.append(LLsN-1)
            LLs.append( [L for L in range(len(blines)+1)] ) # Add boundary LineLoop first and foremost
            for b, bline in enumerate (blines):
                Ls_index.append(b+1)
                Ls.append([bline[0], bline[1]]) # Create boundary line first and foremost

                # if sqrt((df[bline[0]]['x'] + bline[1]]['x']) ** 2 + (df[bline[0]]['y'] + bline[1]]['y']) ) > h_threshold: # Control if we want to create fringe when distance between two points is very short < threshold
                # Create fringe points after controls checks
                theta_rad, theta_deg = xaxis_vs_lineseg_angle([df.loc[bline[0]]['x'], df.loc[bline[0]]['y']], [df.loc[bline[1]]['x'], df.loc[bline[1]]['y']] )
                exec(f"pt{ip+1} = [  df.loc[{bline[0]}]['x'] + np.cos({theta_rad-math.radians(90)}) * {avg_min_dist},   df.loc[{bline[0]}]['y'] + np.sin({theta_rad-math.radians(90)}) * {avg_min_dist}, {fringe_aperture}, 'L{len(blines)+b+1}_LEFT']  ")
                exec(f"pt{ip+2} = [  df.loc[{bline[1]}]['x'] + np.cos({theta_rad-math.radians(90)}) * {avg_min_dist},   df.loc[{bline[1]}]['y'] + np.sin({theta_rad-math.radians(90)}) * {avg_min_dist}, {fringe_aperture}, 'L{len(blines)+b+1}_RIGHT']  ")
                exec(f"pt{ip+3} = [  (pt{ip+1}[0] + pt{ip+2}[0])/2,   (pt{ip+1}[1] + pt{ip+2}[1])/2, {fringe_aperture}, 'L{b+1}_MID']  ") # fringe mid point
                exec(f"fpts.append( pt{ip+1} )  ")
                exec(f"fpts.append( pt{ip+2} )  ")
                exec(f"fpts.append( pt{ip+3} )  ")
                LsN += 1
                Ls_index.append(LsN-1)
                exec(f"Ls.append(     [{bline[1]},  {ip+2}] )  ")
                LsN += 1
                Ls_index.append(LsN-1)
                exec(f"Ls.append(     [{ip+2},      {ip+1}] )  ")
                LsN += 1
                Ls_index.append(LsN-1)
                exec(f"Ls.append(     [{ip+1},      {bline[0]}] )  ")
                LLsN +=1
                LLs_index.append(LLsN-1)
                exec(f"LLs.append( [{b+1}, {LsN-2}, {LsN-1}, {LsN} ])  ")
                
                
                ip += 3
            FLs = [l for l in Ls if l not in blines] # Remove all blines from Ls
            
            
            FLLs = list( map( list, zip(splitlist(FLs, int(len(FLs)/len(blines)))) ) )
            LLPoints = [ np.unique(linelooppoints).tolist() for linelooppoints in FLLs] # https://stackoverflow.com/questions/716477/join-list-of-lists-in-python
            
            fptsdf = pd.DataFrame(fpts, columns = colnames) # Fringe points dataframe 
            df = pd.concat([df, fptsdf], ignore_index = True)
            if use_aperture_as_z == False:
                df['aperture'] = 0
            
            
            if fringes == False:

                #   O      O
                #    ######
                #    ######
                #    ######
                #   O      O

                # Write Line Loop for the boundary
                Nlineloops = 1
                output_file.write("\n\n")
                output_file.write( "Line Loop(1)={" + ",".join( [str(i) for i in range(len(bindex)-1)]) + "};" )
                
                # Write Plane Surface
                Nplanesurfaces = 1
                output_file.write("\n\n")
                output_file.write( "Plane Surface(1)={1};")
                
                # Add Points to Surface
                output_file.write("\n\n")
                output_file.write( "Point{" + ",".join( [str(i) for i in range(Nsurfpoints)]) + "} In Surface{1};" )
                
                # Write Physical Surface
                Nphysicalsurfaces = 1
                output_file.write("\n\n")
                output_file.write( "Physical Surface(1)={1};")
                
                # Add Points to Surface
                output_file.write("\n\n")
                output_file.write( "Point{" + ",".join( [str(i) for i in range(Nsurfpoints)]) + "} In Surface{1};" )
                
                # Write Physical Surface
                Nphysicalsurfaces = 1
                output_file.write("\n\n")
                output_file.write( "Physical Surface(1)={1};")
                
            else:
                if square == True:
                    if Dimensions == 2:

# =============================================================================
#                         
#                         # Create fringes points
#                         #     O    O
#                         #    O######O
#                         #     ######
#                         #    O######O
#                         #     O    O
#                         p1df2i = df.loc[(df['x']==df['x'].min()) & (df['y']==df['y'].min())].index[0] # Get df2 index of (xmin, ymin)
#                         p2df2i = df.loc[(df['x']==df['x'].max()) & (df['y']==df['y'].min())].index[0] # Get df2 index of (xmax, ymin)
#                         p3df2i = df.loc[(df['x']==df['x'].max()) & (df['y']==df['y'].max())].index[0] # Get df2 index of (xmax, ymax)
#                         p4df2i = df.loc[(df['x']==df['x'].min()) & (df['y']==df['y'].max())].index[0] # Get df2 index of (xmin, ymax)
#                         # Get tuples of all necessary boundary points for the fringes
#                         p1 = (df2.loc[p1df2i]['x'], df2.loc[p1df2i]['y'], df2.loc[p1df2i]['aperture'])
#                         p2 = (df2.loc[p2df2i]['x'], df2.loc[p2df2i]['y'], df2.loc[p2df2i]['aperture'])
#                         p3 = (df2.loc[p3df2i]['x'], df2.loc[p3df2i]['y'], df2.loc[p3df2i]['aperture'])
#                         p4 = (df2.loc[p4df2i]['x'], df2.loc[p4df2i]['y'], df2.loc[p4df2i]['aperture'])
#                         
#                         if fringe_aperture == None:
#                             p5 = [bxmax, bymin - dy, df2.loc[p1df2i]['aperture']] # Copies the aperture from the closest points hardcoded
#                             p6 = [bxmin, bymin - dy, df2.loc[p1df2i]['aperture']]
#                             p7 = [bxmax + dx, bymax, df2.loc[p2df2i]['aperture']]
#                             p8 = [bxmax + dx, bymin, df2.loc[p2df2i]['aperture']]
#                             p9 = [bxmin, bymax + dy, df2.loc[p3df2i]['aperture']] 
#                             p10 = [bxmax, bymax + dy, df2.loc[p3df2i]['aperture']]
#                             p11 = [bxmin - dx, bymin, df2.loc[p4df2i]['aperture']] 
#                             p12 = [bxmin - dx, bymax, df2.loc[p4df2i]['aperture']]
#                         else:
#                             p5 = [bxmax, bymin - dy, fringe_aperture] # Uses the aperture user defined above
#                             p6 = [bxmin, bymin - dy, fringe_aperture]
#                             p7 = [bxmax + dx, bymax, fringe_aperture]
#                             p8 = [bxmax + dx, bymin, fringe_aperture]
#                             p9 = [bxmin, bymax + dy, fringe_aperture] 
#                             p10 = [bxmax, bymax + dy, fringe_aperture]
#                             p11 = [bxmin - dx, bymin, fringe_aperture] 
#                             p12 = [bxmin - dx, bymax, fringe_aperture]
#                         
#                         # Boundary points' indexes. Left and Right are from the perspective of the viewer looking to the polyline from the outside
#                         Left_Boundary_SOUTHi = df.loc[((df['x']==df['x'].min()) & (df['y']==df['y'].min()))].index[0]
#                         Left_Boundary_EASTi = df.loc[((df['x']==df['x'].max()) & (df['y']==df['y'].min()))].index[0]
#                         Left_Boudary_NORTHi = df.loc[((df['x']==df['x'].max()) & (df['y']==df['y'].max()))].index[0]
#                         Left_Boundary_WESTi =   df.loc[((df['x']==df['x'].min()) & (df['y']==df['y'].max()))].index[0]
#                         Right_Boundary_SOUTHi = df.loc[((df['x']==df['x'].max()) & (df['y']==df['y'].min()))].index[0]
#                         Right_Boundary_EASTi = df.loc[((df['x']==df['x'].max()) & (df['y']==df['y'].max()))].index[0]
#                         Right_Boudary_NORTHi = df.loc[((df['x']==df['x'].min()) & (df['y']==df['y'].max()))].index[0]
#                         Right_Boundary_WESTi = df.loc[((df['x']==df['x'].min()) & (df['y']==df['y'].min()))].index[0]
#                         
#                         #Mid points' lists
#                         # Boundary Mid points
#                         Mid_Boundary_SOUTH =  [ (p1[0]+p2[0])/2,        p1[1],              (p1[2]+p2[2])/2,  "Mid_Boundary_SOUTH"]   # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_EAST =   [ p2[0],                  (p2[1]+p3[1])/2,    (p2[2]+p3[2])/2,  "Mid_Boundary_EAST"]  # EAST Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_NORTH =  [ (p3[0]+p4[0])/2,        p3[1],              (p3[2]+p4[2])/2, "Mid_Boundary_NORTH"]  # NORTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_WEST =   [ p4[0],                  (p4[1]+p1[1])/2,    (p4[2]+p1[2])/2,  "Mid_Boundary_WEST"]
#                         # Fringe Mid points
#                         Mid_Fringe_SOUTH = [(p5[0]+p6[0])/2, p5[1], (p5[2]+p6[2])/2,  "Mid_Fringe_SOUTH"]    # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Fringe_EAST = [p7[0], (p7[1]+p8[1])/2 , (p7[2]+p8[2])/2,  "Mid_Fringe_EAST"]   # EAST Boundary Mid point (between p2 and p3)
#                         Mid_Fringe_NORTH = [(p9[0]+p10[0])/2, p9[1], (p9[2]+p10[2])/2,  "Mid_Fringe_NORTH"]   # NORTH Boundary Mid point (between p3 and p4)
#                         Mid_Fringe_WEST = [p11[0], (p11[1]+p12[1])/2 , (p11[2]+p12[2])/2,  "Mid_Fringe_WEST"]   # WEST Boundary Mid point (between p4 and p1)
#                         
#                         # Boundary Mid points
#                         Mid_Boundary_SOUTHdf = pd.DataFrame(  [ Mid_Boundary_SOUTH] , columns = colnames)   # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_EASTdf = pd.DataFrame(   [ Mid_Boundary_EAST ] , columns = colnames)   # EAST Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_NORTHdf = pd.DataFrame(  [ Mid_Boundary_NORTH ] , columns = colnames)  # NORTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_WESTdf = pd.DataFrame(   [ Mid_Boundary_WEST ] , columns = colnames)   
#                         # Fringe Mid points
#                         Mid_Fringe_SOUTHdf = pd.DataFrame(  [ Mid_Fringe_SOUTH] , columns = colnames)   # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Fringe_EASTdf = pd.DataFrame(   [ Mid_Fringe_EAST ] , columns = colnames)   # EAST Boundary Mid point (between p1 and p2)
#                         Mid_Fringe_NORTHdf = pd.DataFrame(  [ Mid_Fringe_NORTH ] , columns = colnames)  # NORTH Boundary Mid point (between p1 and p2)
#                         Mid_Fringe_WESTdf = pd.DataFrame(   [ Mid_Fringe_WEST ] , columns = colnames)   
# =============================================================================
                        
                        #Concatenate
                        # df = pd.concat([df, Mid_Boundary_SOUTHdf, Mid_Boundary_EASTdf, Mid_Boundary_NORTHdf, Mid_Boundary_WESTdf, Mid_Fringe_SOUTHdf, Mid_Fringe_EASTdf, Mid_Fringe_NORTHdf, Mid_Fringe_WESTdf], ignore_index = True)         
                        
# =============================================================================
#                         # Boundary Mid points' indexes
#                         Mid_Boundary_SOUTHi = df.loc[(df['x']==(p1[0]+p2[0])/2) & (df['y']==p1[1])].index[0]   # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_EASTi  = df.loc[(df['x']==p2[0])           & (df['y']==(p2[1]+p3[1])/2)].index[0]  # EAST Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_NORTHi = df.loc[(df['x']==(p3[0]+p4[0])/2) & (df['y']==p3[1])].index[0]  # NORTH Boundary Mid point (between p1 and p2)
#                         Mid_Boundary_WESTi  = df.loc[(df['x']==p4[0])           & (df['y']==(p4[1]+p1[1])/2)].index[0]
#                         # Fringe Mid points' indexes
#                         Mid_Fringe_SOUTHi   =   df.loc[(df['x']==(p5[0]+p6[0])/2)   & (df['y']== p5[1])              ].index[0]    # SOUTH Boundary Mid point (between p1 and p2)
#                         Mid_Fringe_EASTi    =   df.loc[(df['x']==p7[0]  )           & (df['y']==(p7[1]+p8[1])/2)     ].index[0] # EAST Boundary Mid point (between p2 and p3)
#                         Mid_Fringe_NORTHi   =   df.loc[(df['x']==(p9[0]+p10[0])/2)  & (df['y']==p9[1]      )         ].index[0]  # NORTH Boundary Mid point (between p3 and p4)
#                         Mid_Fringe_WESTi    =   df.loc[(df['x']==p11[0] )           & (df['y']==(p11[1]+p12[1])/2)   ].index[0]  # WEST Boundary Mid point (between p4 and p1)
# =============================================================================
                        
                        
                        
                        """ Don't need this anymore?"""
# =============================================================================
#                         # Add names to df's boundary points                        
#                         mid_bps = {Mid_Boundary_SOUTHi, Mid_Boundary_EASTi, Mid_Boudary_NORTHi, Mid_Boundary_WESTi} # These are indexes, not (x,y,z)s
#                         mid_bps_names = ['Mid_Boundary_SOUTH', 'Mid_Boundary_EAST', 'Mid_Boudary_NORTH', 'Mid_Boundary_WEST']
#                         corner_bps = {Left_Boundary_SOUTHi, Right_Boundary_SOUTHi, Left_Boudary_NORTHi, Right_Boudary_NORTHi} # These are indexes, not (x,y,z)s
#                         corner_bps_names = ['Bottom_Left_Boundary', 'Bottom_Right_Boundary', 'Top_Left_Boundary', 'Top_Right_Boundary']
#                         for i, n in list(map(list,zip( mid_bps, mid_bps_names))):
#                             df.loc[i, 'name'] = n
#                         for i, n in list(map(list,zip( corner_bps, corner_bps_names))):
#                             df.loc[i, 'name'] = n
# =============================================================================
                        
# =============================================================================
#                         # If you want to actually create the points but they already exist and you will have duplicates, thus why we are getting their corresponding index above so we can name them in the .geo file
#                         p13 = [(p1[0]+p2[0])/2, p1[1], (p1[2]+p2[2])/2]    # SOUTH Boundary Mid point (between p1 and p2)
#                         p14 = [p2[0], (p2[1]+p3[1])/2 , (p2[2]+p3[2])/2]   # EAST Boundary Mid point (between p1 and p2)
#                         p15 = [(p3[0]+p4[0])/2, p3[1], (p3[2]+p4[2])/2]   # NORTH Boundary Mid point (between p1 and p2)
#                         p16 = [p4[0], (p4[1]+p1[1])/2 , (p4[2]+p1[2])/2]   # WEST Boundary Mid point (between p1 and p2)
#                         # Add points to df
#                         dataframe = [
#                                     [p5[0], p5[1], p5[2]],
#                                     [p6[0], p6[1], p6[2]],
#                                     [p7[0], p7[1], p7[2]],
#                                     [p8[0], p8[1], p8[2]],
#                                     [p9[0], p9[1], p9[2]],
#                                     [p10[0], p10[1], p10[2]],
#                                     [p11[0], p11[1], p11[2]],
#                                     [p12[0], p12[1], p12[2]],
#                                     [p13[0], p13[1], p13[2]],
#                                     [p14[0], p14[1], p14[2]],
#                                     [p15[0], p15[1], p15[2]],
#                                     [p16[0], p16[1], p16[2]],
#                                     [p17[0], p17[1], p17[2]],
#                                     [p18[0], p18[1], p18[2]],
#                                     [p19[0], p19[1], p19[2]],
#                                     [p20[0], p20[1], p20[2]]
#                                     ]
# =============================================================================
                        
# =============================================================================
#                         # Add points to df
#                         Fringe_corner_pts_lst = [
#                                     [p5[0], p5[1], p5[2], "Right_Fringe_SOUTH"],
#                                     [p6[0], p6[1], p6[2], "Left_Fringe_SOUTH"],
#                                     [p7[0], p7[1], p7[2], "Right_Fringe_EAST"],
#                                     [p8[0], p8[1], p8[2], "Left_Fringe_EAST"],
#                                     [p9[0], p9[1], p9[2], "Right_Fringe_NORTH"],
#                                     [p10[0], p10[1], p10[2], "Left_Fringe_NORTH"],
#                                     [p11[0], p11[1], p11[2], "Right_Fringe_WEST"],
#                                     [p12[0], p12[1], p12[2], "Left_Fringe_WEST"],
#                                     ]
#                         Fringe_mid_pts_lst = [
#                                     [Mid_Fringe_SOUTH[0], Mid_Fringe_SOUTH[1], Mid_Fringe_SOUTH[2], "Mid_Fringe_SOUTH"],
#                                     [Mid_Fringe_EAST[0], Mid_Fringe_EAST[1], Mid_Fringe_EAST[2], "Mid_Fringe_EAST"],
#                                     [Mid_Fringe_NORTH[0], Mid_Fringe_NORTH[1], Mid_Fringe_NORTH[2], "Mid_Fringe_NORTH"],
#                                     [Mid_Fringe_WEST[0], Mid_Fringe_WEST[1], Mid_Fringe_WEST[2], "Mid_Fringe_WEST"]
#                                     ]
#                         fringes_long = []
#                         if fringes.rfind("S") != -1:
#                             fringes_long.append("SOUTH")
#                         if fringes.rfind("E") != -1:
#                             fringes_long.append("EAST") 
#                         if fringes.rfind("N") != -1:
#                             fringes_long.append("NORTH")
#                         if fringes.rfind("W") != -1:
#                             fringes_long.append("WEST")
#                         print(fringes, fringes_long)
#                         Fringe_corner_pts_df = pd.DataFrame(Fringe_corner_pts_lst, columns=colnames)
#                         Fringe_corner_pts_df = Fringe_corner_pts_df[Fringe_corner_pts_df['name'].str.contains("|".join([str(s) for s in fringes_long])) == True]
#                         
#                         findex = list(range(len(df), len(df) + len(Fringe_corner_pts_df))) # index of the points making up the fringes continuing from df index
#                         Nsurfpoints = df.shape[0] # number of points before the fringes
#                         ip = Nsurfpoints # Keep track of the tag number of points because we are creating more
#                         Fringe_mid_pts_df = pd.DataFrame(Fringe_mid_pts_lst, columns=colnames)
#                         Fringe_mid_pts_df = Fringe_mid_pts_df[Fringe_mid_pts_df['name'].str.contains("|".join([str(s) for s in fringes_long])) == True]
#                         Fringe_pts_df = pd.concat([Fringe_corner_pts_df, Fringe_mid_pts_df], ignore_index=True)
#                         df = pd.concat([df, Fringe_pts_df], ignore_index=True) # add Fringe_pts_df to df continuing with existing df's indexing
# =============================================================================
                        
                        
                    # # Under construction
                    # elif Dimensions == 3: 
                    #     ptz = (df['z'].max() - df['z'].min())
                    
                    
            output_file.write("Mesh.MshFileVersion = 2.2;")
            output_file.write(f"lc={lc};//Characteristic length")
            output_file.write("\n\n")
            # Write Points whilst $NAMEing the mid and edge points
            for row in df.itertuples():
                if df.iloc[row[0]]['name'] != "" :  # The below way of concatenating the string accommodates for 2+ dimensional points
                    output_file.write( f"Point({row[0]})=" + "{" + ", ".join([str(s) for s in row[1:-1]]) + ", lc};" + f"\t//$NAME {row[-1]}\n")
                else:
                    output_file.write( f"Point({row[0]})=" + "{" + ", ".join([str(s) for s in row[1:-1]]) + ", lc};\n") # Not using z-values, i.e. rendering a flat mesh
                

                    
# =============================================================================
#             # Work in progress.... Decided to hardcode below for now
#             c = [i for i in range(0, len(df2)+1, 3)] # hardcoded 3 because each line is 4 corners
#             d = [[i,j] for i in c[:-1] for j in c[1:]][::len(c)]
#             lines = []
#             for i,j in d:
#                 l=[]
#                 for t in range(i, j+1):
#                     l.append(t)
#                 l.append(i); print(l)
#                 line = list(zip(l[1:], l[:-1])); print(line);
# =============================================================================
            Nlines = 0
            Ntransline = 0
            # polylines_file.write("#POLYLINES\n")
            polylines_file.write("Polyline,Point1,Point2,Point1x,Point1y,Point2x,Point2y,MidPoint,MidPointx,MidPointy,Length\n")
            for line in blines: # I could put the two loops in one but the main boundary lines would not be 1-4 anymore so i have left it.
                Nlines += 1
                Ntransline += 1
                output_file.write( f"Line({Nlines})=" + "{" + f"{line[0]}," + f"{line[1]}" + "};\n")
                d, j = tree.query(  (   (df.iloc[line[0]]['x'] + df.iloc[line[1]]['x'])/2    , (df.iloc[line[0]]['y'] + df.iloc[line[1]]['y'])/2), k=1)
                Point1,Point2,Point1x,Point1y,Point2x,Point2y,MidPoint,MidPointx,MidPointy,Length = \
                line[0], \
                line[1], \
                df.iloc[line[0]]['x'], \
                df.iloc[line[0]]['y'], \
                df.iloc[line[1]]['x'], \
                df.iloc[line[1]]['y'], \
                j, \
                df.iloc[j]['x'], \
                df.iloc[j]['y'], \
                math.sqrt( (df.iloc[line[0]]['x'] - df.iloc[line[1]]['x'])**2 + (df.iloc[line[0]]['y'] - df.iloc[line[1]]['y'])**2 )
                
                polylines_file.write( f"{Nlines}, {Point1},{Point2},{Point1x},{Point1y},{Point2x},{Point2y},{MidPoint},{MidPointx},{MidPointy},{Length}\n" )
                output_file.write( f"Transfinite Line ({Ntransline}) = {max(abs(Point1x - Point2x), abs(Point1y - Point2y))+1} Using Progression 1.0 ;\n" )

            for i, fline in enumerate(FLs): # for each boundary line of the initial data, start at the last point of that line and work your way around the fringe clockwise
                Nlines += 1
                Ntransline += 1
                output_file.write( f"Line({Nlines})=" + "{" + f"{fline[0]}, {fline[1]}" + "};\n")
                d, j = tree.query(  (   (df.iloc[fline[0]]['x'] + df.iloc[fline[1]]['x'])/2    , (df.iloc[fline[0]]['y'] + df.iloc[fline[1]]['y'])/2), k=1)
                Point1,Point2,Point1x,Point1y,Point2x,Point2y,MidPoint,MidPointx,MidPointy,Length = \
                fline[0], \
                fline[1], \
                df.iloc[fline[0]]['x'], \
                df.iloc[fline[0]]['y'], \
                df.iloc[fline[1]]['x'], \
                df.iloc[fline[1]]['y'], \
                j, \
                df.iloc[j]['x'], \
                df.iloc[j]['y'], \
                math.sqrt( (df.iloc[fline[0]]['x'] - df.iloc[fline[1]]['x'])**2 + (df.iloc[fline[0]]['y'] - df.iloc[fline[1]]['y'])**2 )
                
                polylines_file.write( f"{Nlines}, {Point1},{Point2},{Point1x},{Point1y},{Point2x},{Point2y},{MidPoint},{MidPointx},{MidPointy},{Length}\n" )
                output_file.write( f"Transfinite Line ({Ntransline}) = {max(abs(Point1x - Point2x), abs(Point1y - Point2y))+1} Using Progression 1.0;\n" )

                                
            # Write Line Loop for the boundary
            output_file.write("\n\n")
            for i, LL in enumerate(LLs):
                Nlineloops = i+1
                output_file.write( f"Line Loop({Nlineloops})=" + "{" f"{','.join(str(l) for l in LL)}" + "};\n" )
            # output_file.write( "Line Loop(2)={1,5,6,7};\n" )
            # output_file.write( "Line Loop(3)={2,8,9,10};\n" )
            # output_file.write( "Line Loop(4)={3,11,12,13};\n" )
            # output_file.write( "Line Loop(5)={4,14,15,16};\n" )
            
            # Write Plane Surface
            Nplanesurfaces = Nlineloops
            output_file.write("\n\n")
            for i in range(1,Nplanesurfaces+1):
                output_file.write(f"Plane Surface({i})="+"{"+f"{i}"+"};\n")
            
            # Write Transfinite Surfaces
            # Write Line Loop for the boundary
            Ntranssurf = Nplanesurfaces
            output_file.write("\n\n")
            output_file.write(f"Transfinite Surface(1)=" + "{" + ", ".join(str(s) for s in np.unique(blines).tolist()) +"};\n") # Boundary surface
            for TS, LL in enumerate(LLPoints):
                output_file.write(f"Transfinite Surface({TS+2})=" + "{" + ", ".join(str(s) for s in np.unique(LL).tolist()) +"};\n")
            
            if tri_Elements == False:
                for i in range(1, Ntranssurf+1):
                    output_file.write("\n")
                    output_file.write("Recombine Surface{" + f"{i}" + "}; //Makes it square")
            else:
                for i in range(1, Ntranssurf+1):
                    output_file.write("\n")
                    output_file.write("//Recombine Surface{" + f"{i}" + "}; //Makes it square")
            #Close output file to submit changes and open again to read
            output_file.close()
            polylines_file.close()
            output_file  = open(outputfilepath, 'r')
            data_output = output_file.read()
            #Submit input and output files' contents (.csv and .geo, respectively) onto the respective text boxes for review by the user
            CSV_file_Text.delete(0.0, END)
            CSV_file_Text.insert(0.0, data_input)
            GEO_file_Text.delete(0.0, END)
            GEO_file_Text.insert(0.0, data_output)
    
            # close all files 
            input_file.close() 
            output_file.close()
            #disable text boxes to edit (As these edits won't be replicated in the files)
            CSV_file_Text.config(state=DISABLED)
            GEO_file_Text.config(state=DISABLED) 







'''Labels'''
#Create labels
CSV_Label = Label(BROWSE_frame, text = ".CSV file", pady = 5)
LC_Label = Label(BROWSE_frame, text = "Characteristic length:", pady = 5)
header_Label = Label(BROWSE_frame, text = "Header row:", pady = 5)
readcols_Label = Label(BROWSE_frame, text = "Columns to read (separete by comma):", pady = 5)
skiprows_Label = Label(BROWSE_frame, text = "Rows to skip (separete by comma):", pady = 5)
csv_file_label = Label(CSV_frame, text = ".CSV file content", pady = 5)
geo_file_label = Label(GEO_frame, text = ".GEO file content", pady = 5)

'''Text boxes and scrolls'''
#Create scrollbars
CSVscroll = Scrollbar(CSV_frame)
GEOscroll = Scrollbar(GEO_frame)
#Create text boxes
CSV_Text = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
LC_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); LC_Entry.insert(-1, 1)
header_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); header_Entry.insert(-1, "0")
readcols_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); readcols_Entry.insert(-1, "1,2,3")
skiprows_Entry = Entry(BROWSE_frame, width = 10, state=NORMAL); skiprows_Entry.insert(-1, "")
CSV_file_Text = Text(CSV_frame, height=20, state=NORMAL, yscrollcommand = CSVscroll.set)
GEO_file_Text = Text(GEO_frame, height=20, state=NORMAL, yscrollcommand = GEOscroll.set)
#Configure scrollbar
CSVscroll.config(command=CSV_file_Text.yview)
GEOscroll.config(command=GEO_file_Text.yview)

'''Buttons'''
BrowseButt= Button(BROWSE_frame, text="Browse .CSV file", fg="black", font=("Ariel", 9, "bold"), command=browse)
FormatButt= Button(BROWSE_frame, text="CreateGEO", fg="black", font=("Ariel", 9, "bold"), command=createGEO)

'''Allocate widgets'''
CSV_Label.pack(side = LEFT)
CSV_Text.pack(side = LEFT)
BrowseButt.pack(side = LEFT, fill = Y)
LC_Label.pack(side = LEFT)
LC_Entry.pack(side = LEFT)

header_Label.pack(side = LEFT)
header_Entry.pack(side = LEFT)
readcols_Label.pack(side = LEFT)
readcols_Entry.pack(side = LEFT)
skiprows_Label.pack(side = LEFT)
skiprows_Entry.pack(side = LEFT)




FormatButt.pack(side = BOTTOM, fill = X)

csv_file_label.grid(row = 1, column = 1, sticky = W)
geo_file_label.grid(row = 1, column = 1, sticky = W)
CSV_file_Text.grid(row = 2, column = 1)
GEO_file_Text.grid(row = 2, column = 1)
CSVscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
GEOscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)

'''Menu bar'''
menubar = Menu(root)
root.config(menu = menubar)

file_menu = Menu(menubar)
menubar.add_cascade(label = "File", menu = file_menu)
file_menu.add_command(label = "CreateGEO", command = createGEO)
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

options_menu = Menu(menubar)
menubar.add_cascade(label = "Options", menu = options_menu)
options_menu.add_command(label = "Description", command = description)

'''Prompt on open'''
#Uncomment below if you wish the description() function (info box explaining how the app works) to be prompted at app opening
#description()

'''Add Window title, CSVmetry and create Window's main loop'''
root.title("Python .csv file Gmsh to Gina Formater")
root.geometry("1600x500")
root.mainloop()