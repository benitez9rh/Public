# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:13:05 2024

@author: s2132627
"""

import os
import pathlib
import pandas as pd
import numpy as np
from ast import literal_eval



extension_tec = ".tec"

wd = pathlib.PureWindowsPath(r"C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 2 - HM Shear\Freiberg Model\Models comparison").as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux

inputfilenames = [
    pathlib.PureWindowsPath(r'Trans33_avg_FB_POLYLINE_L4_TIM_LIQUID_FLOW').as_posix(),
    pathlib.PureWindowsPath(r'Trans65_avg_FB_POLYLINE_L4_TIM_LIQUID_FLOW').as_posix(),
    pathlib.PureWindowsPath(r'Trans33_krig_FB_POLYLINE_L4_TIM_LIQUID_FLOW').as_posix(),
    pathlib.PureWindowsPath(r'Trans65_krig_FB_POLYLINE_L4_TIM_LIQUID_FLOW').as_posix()
    ]



os.chdir(wd)                                   # set the working directory
wdsave = wd
bs="//"; wdsave=wdsave+bs  



def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s



for i, inputfilename in enumerate(inputfilenames):
    with open(wd + inputfilename + extension_tec, 'r') as file:
        
        # Check the file format
        # head = file.readlines() # Check the first 50 lines
        words = ["VARIABLES", "ZONE"]
        
        Nzones = 0
        ZONE = False
        VARIABLES = False
        nlines = 0
        for line in file:
            # print(line)
            casematch = next((word for word in words if word in line), False)
            if casematch:
                # for word in words: # Creates a False flag for each word in words  
                #     exec(f"{word[word.rfind('$')+1:]} = False")
                ZONE = False
                VARIABLES = False
                print(casematch)
                print(VARIABLES)
                print(ZONE)
                # if casematch.find("NODES") != -1 or casematch.find("Nodes") != -1: # If line contains "NODES"
                if "zone" in casematch.lower():
                    ZONE = True
                    print(ZONE)
                    Nzones += 1
                if "variables" in casematch.lower():
                    VARIABLES = True
                    print(VARIABLES)
                    line = line.replace('"', '')
                    line = line.replace('VARIABLES = ', '')
                    line = line.replace('\n', '')
                    columns = line.split(",")
                    # columns.append("Zone")
                    try: #  Check if df has been defined
                        print(df.head())
                    except AttributeError: # catch when df1 is None
                        pass
                    except NameError: # catch when it hasn't even been defined
                        df = pd.DataFrame( columns = columns)
            else:
                line = [float(x) for x in line.split()]
                line.append(Nzones)
                res_df = {columns[j]: line[j] for j in range(len(columns))}
                df = df.append(res_df, ignore_index = True)
            
        nlines +=1
        print(nlines)                
                # Nnodes = int(next(file))
                
        # k = [f"Zone {key}"for key in range(1, Nzones + 1)] # Check multi index or keys in pandas
        # result = pd.concat(frames, keys=k)
        exec(f"df{i} = df.set_index('Time')")
        
    
    if i == 0:
        result_df = df0
    else:
        # df.drop('Time', axis=1)
        # result_df = pd.concat([result_df, df], axis=1, join='inner', ignore_index=True)
        result_df = pd.concat([result_df, df], axis = "columns")
    
    
# df = pd.concat(exec(f"{[df{i} for i in {range(len(inputfilenames))}]}"), axis=0, join='inner')

# for i in range(len(inputfilenames)):
    
#     exec( f"df = df0.merge(df{i}, on='A', how='inner') ")
        
        
            