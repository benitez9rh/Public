# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:00:58 2022

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
import geostatspynscore
from time import time, ctime
import time
from scipy.spatial import KDTree
from matplotlib import cm
from sklearn.metrics import r2_score

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
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals\Upscale').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs                              # Use this instead in linux
inputfilename = pathlib.PureWindowsPath(r'Fracs').as_posix()
os.chdir(wd)                                   # set the working directory
# =============================================================================
mainpd = pd.read_csv(wd + "Fracs" + extension_csv, index_col = 0)



x = "var"; y = "std"; #avg, std,var

deg=2
plot="False"
scalekw = "Scale"

if deg == 1:
    for scale in mainpd[scalekw].unique():
        a = mainpd[mainpd[scalekw] == scale]
        m, b = np.polyfit(a[x],a[y], deg)
        r_squared_value = r2_score(a[y], [x*m+b for x in a[x]])
        plt.scatter(a[x],a[y], s = 1)
        plt.axline((0, b), slope=m, color='red', linestyle='dashed', linewidth=1, label = f"R^2 = {r_squared_value}")
        plt.xlabel(f'{x}', fontsize=15)
        plt.ylabel(f'{y}', fontsize=15)
        plt.title(f'{inputfilename} {x} vs {y} at Scale: {scale}')
        plt.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(f'{wd}{inputfilename}_corel_{x}_vs_{y}_scale_{int(scale)}{extension_png}')
        plt.show()
    
elif deg == 2:
    
    for scale in mainpd[scalekw].unique():
        a = mainpd[mainpd[scalekw] == scale]; a = a.sort_values(by=[x])
        p2,p1,p0 = np.polyfit(a[x],a[y], deg) # p2 is the value to be multiplied by the x raised to the power with the same index (2), p1 to x, etc.
        print(p0,p1,p2)
        pred_y = [p2*x**2 + p1*x + p0 for x in a[x]]
        r_squared_value = r2_score(a[y], pred_y)
        plt.scatter(a[x], a[y], s = 2)
        plt.yscale('log', base=10)
        plt.plot(np.linspace(a[x].min(), a[x].max(),300), [p2*x**2 + p1*x + p0 for x in np.linspace(a[x].min(), a[x].max(),300)], color='red', linestyle='dashed', linewidth=1, label = f"R^2 = {r_squared_value}")
        plt.xlabel(f'{x}', fontsize=15)
        plt.ylabel(f'{y}', fontsize=15)
        plt.title(f'{inputfilename} {x} vs log({y}) at Scale: {scale}')
        plt.legend(loc="upper left")
        plt.tight_layout()
        if plot=="True":
            plt.savefig(f'{wd}{inputfilename}_correl_{x}_vs_log({y})_scale_{int(scale)}{extension_png}')
        else:
            plt.savefig(f'{wd}{inputfilename}_correl_{x}_vs_{y}_scale_{int(scale)}{extension_png}')
        plt.show()    