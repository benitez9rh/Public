# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:33:29 2023

@author: s2132627
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import pathlib
from ast import literal_eval
import pandas as pd
import numpy as np
from math import floor
'''Methods'''
def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s

# =============================================================================
#         # ############################# For testing ##############################################################################################################################################################################
#         import pathlib
#         lc = 1
#         header = 0
#         readcols = [0,1,2]
#         skiprows = ""
#         CSVinputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\M CNL_Aperture-KrigFromFullDetail\Kriged_aperture_Based_on_Extracted_M_CNL_Aperture.csv').as_posix()
#         CSVinputfilename = CSVinputfilepath[CSVinputfilepath.rfind("/")+1:]
#         pathCSV = CSVinputfilepath[:CSVinputfilepath.rfind("/")+1]
# 
#         MSHinputfilepath = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\Main Model\FB.msh').as_posix()
#         MSHinputfilename = MSHinputfilepath[MSHinputfilepath.rfind("/")+1:]
#         pathMSH = MSHinputfilepath[:MSHinputfilepath.rfind("/")+1]
#         
#         MSH_OFinputfilepath = MSHinputfilepath
#         MSH_OFinputfilename = MSH_OFinputfilepath[MSH_OFinputfilepath.rfind("/")+1:] + "_of"
#         pathMSH_OF = MSH_OFinputfilepath[:MSH_OFinputfilepath.rfind("/")+1]
# 
#         
#         header = int(0)
#         readcols = [1,2,3]
#         skiprows = ""
#         if "" in skiprows:
#             df = pd.read_csv(pathCSV + CSVinputfilename , usecols = readcols, header = header)
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
        
def msh_of_SplitFunction_v3(path, input_file_name):
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
            ENC_df["grid x-coordinate"]  = (round(ENC_df["Element Centre x-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))                                                             # Round because there might be a precision issue. This rounding does not change the initial dataframe. -1 because the ranking function starts at 1
            ENC_df["grid x-coordinate"]  = ENC_df["grid x-coordinate"] / (np.sort(ENC_df["grid x-coordinate"].unique())[1] - np.sort(ENC_df["grid x-coordinate"].unique())[0])      # rank function assigns the same rank to duplicates but jumps that amount of duplicates when the next value is different, which we don't want. Thus we are calculating the difference between sequential values and dividing..
            ENC_df["grid y-coordinate"]  = (round(ENC_df["Element Centre y-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))
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