# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:13:26 2022

@author: s2132627
"""

import os
import pathlib
import pandas as pd
from math import sin, cos, e, sqrt, floor, ceil, radians
import matplotlib.pyplot as plt
from sympy import limit, Symbol, oo
import numpy as np
import math

""" """
extension_csv = ".csv"
extension_png = ".png"
LSAC = 'LSAC'; zero = 'zero' # Last Scale Average Centre, i.e. centre the harmonic motion and envelopes around the last scale average
""" ####################################################################################################################################### """


""" ####################################################### User-defined Variables ###################################################### """
vcol = 'dipdir' # 'dipdir', 'aperture'
m = 5                 # mass (kg)
k =  50                  # spring stiffness
c = 0.5                # damping coefficient
v0 = 0                  
save = True
centre = None # Options are None, LSAC and zero
radians_ticks = True

wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\DipDirPlot').as_posix() #flips the backslashes to plug in open()
bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
inputfilename = pathlib.PureWindowsPath(r'GW1Q4_z_ORIGvsORIG-tr(-0.02, 0.31, 0)_ApertureMap_aperture_DipDirSpread').as_posix()
""" ####################################################################################################################################### """

os.chdir(wd)                                   # set the working directory
# df = pd.read_csv(wd + inputfilename + extension_csv)                  #This new line is because reading the .csv file with a header and an index column will create some issue
df = pd.read_csv(wd + inputfilename + extension_csv, usecols=[1,2])   # For readingdip and dipdir spread csv files

# Change working directory to save
wd = pathlib.PureWindowsPath(wd+r'\PolarPlot').as_posix() #flips the backslashes to plug in open()
bs="//"; wd=wd+bs 


""" ######################################################## Define Functions ##########################################################  """
def cdc(dforlist, vcol):    # Converts circular (cyclic) data (). # Inspired from https://en.wikipedia.org/wiki/Circular_mean#Example and https://ianlondon.github.io/blog/encoding-cyclical-features-24hour-time/
    if type(dforlist) == "list":
        # https://datascience.stackexchange.com/questions/5990/what-is-a-good-way-to-transform-cyclic-ordinal-attributes
        sins = ([math.sin(i) for i in dforlist]) # sines
        coss = ([math.cos(i) for i in dforlist]) # cosines
    elif isinstance(df, pd.DataFrame):
        
        sins = list(dforlist[vcol].apply(math.sin)) # sines
        coss = list(dforlist[vcol].apply(math.cos)) # cosines

    atans = []
    for s, c in zip(sins, coss):
        if (s > 0 and c > 0): # Check if quadrant 1
            atans.append(math.atan(s/c)); #print("Quadrant 1"); 
        elif c < 0: # Check if quadrant 3 or 4
            atans.append(math.atan(s/c)+ 180) ; #print("Quadrant 3 or 4"); 
        elif (s < 0 and c > 0): # Check if quadrant 2
            atans.append(math.atan(s/c) + 360) ; #print("Quadrant 2"); 
        
    return atans

def circmean(dforlist, degrees = True):    # Calculates the mean of circular (cyclic) data (). # https://en.wikipedia.org/wiki/Circular_mean#Example https://datascience.stackexchange.com/questions/5990/what-is-a-good-way-to-transform-cyclic-ordinal-attributes
    if type(dforlist) == "list":
        if degrees == True:
            sins = sum([math.sin(i) for i in dforlist])/len(dforlist) # Average of the sines
            coss = sum([math.cos(i) for i in dforlist])/len(dforlist) # Average of the cosines
            arctand = math.degrees(math.atan(sins/coss))
            if (sins > 0 and coss > 0):
                circmean = arctand; print("sins > 0 and coss > 0"); print("Circular Mean: ", circmean)
            elif coss < 0:
                circmean = arctand + 180; print("coss < 0"); print("Circular Mean: ", circmean)
            elif (sins < 0 and coss > 0):
                circmean = arctand + 360; print("sins < 0 and coss > 0"); print("Circular Mean: ", circmean)
        else:
            sins = sum([math.sin(math.radians(i)) for i in dforlist])/len(dforlist) # Average of the sines
            coss = sum([math.cos(math.radians(i)) for i in dforlist])/len(dforlist) # Average of the cosines
            arctanr = math.atan(sins/coss)
            if (sins > 0 and coss > 0):
                circmean = arctanr; print("sins > 0 and coss > 0"); print("Circular Mean: ", circmean)
            elif coss < 0:
                circmean = arctanr + math.pi; print("coss < 0"); print("Circular Mean: ", circmean)
            elif (sins < 0 and coss > 0):
                circmean = arctanr + 2*math.pi; print("sins < 0 and coss > 0"); print("Circular Mean: ", circmean)
""" ####################################################################################################################################### """


    
""" ################################################## Centre plot around y-value ############################################################ """
if centre == 'LSAC': # centre the y-axis at the Last Scale Average
    circavg = pd.DataFrame(list(zip(df['scale'].unique(), df.groupby(['scale']).apply(circmean))), columns = ['scale', vcol+" circmean"]) # Calculate Circular average per scale
    convavg = circavg; convavg.loc[convavg[vcol+" circmean"] > 180, vcol+" circmean"] = convavg[vcol+" circmean"] - 360
    dfc = df; dfc.loc[(dfc[vcol] > 180), vcol] = dfc[vcol] - 360 # https://thispointer.com/get-last-value-of-a-column-in-pandas-dataframe/
    convlastavg = convavg.iloc[-1, convavg.columns.get_loc(vcol+" circmean")]
elif centre == 'zero': # centre at 0deg  i.e. zero equals north. Degrees count clockwise are positive up to 180deg and counterclockwise are negative up to -180
    circavg = pd.DataFrame(list(zip(df['scale'].unique(), df.groupby(['scale']).apply(circmean))), columns = ['scale', vcol+" circmean"]) # Circular average
    convavg = circavg; convavg.loc[convavg[vcol+" circmean"] > 180, vcol+" circmean"] = convavg[vcol+" circmean"] - 360
    dfc = df; dfc.loc[(dfc[vcol] > 180), vcol] = dfc[vcol] - 360 # https://thispointer.com/get-last-value-of-a-column-in-pandas-dataframe/
    clastavg = convavg.iloc[-1, convavg.columns.get_loc(vcol+" circmean")]
else:
    dfc = df
    circavg = pd.DataFrame(list(zip(dfc['scale'].unique(), dfc.groupby(['scale']).apply(circmean))), columns = ['scale', vcol+" circmean"]) # Circular average
    convavg = circavg
""" ####################################################################################################################################### """

""" Calculated Variables """
diff = (abs(df[vcol].max()) - abs(df[vcol].min()))/2
tmax = df['scale'].max()              # how long to run for. TO DO: Use 0 to calculate the limit
y0 = max(abs(df[vcol].max()),abs(df[vcol].min()))                  # initial releasing position (initial condition)
dr = c/(2*sqrt(m*k))    # damping ratio (ζ)
                        # undamped ( ζ=0 , blue)
                        # underdamped ( 0<ζ<1 , orange)
                        # critically damped ( ζ=1 , green)
                        # overdamped ( ζ>1 , red)
alpha = c/(2*m); alphas = alpha**2
w0 = sqrt(k/m)          # natural (undamped) angular frequency
w0s = w0**2             # 
w = sqrt(w0s - alphas)  # angular frequency
if dr == 0:
    print(f'Undamped system: dr = {dr}')
elif (dr>0 and dr<1):
    print(f'Underdamped system: dr = {dr}')
elif dr == 1:
    print(f'Critically damped system: dr = {dr}')
elif dr > 1:
    print(f'Overdamped system: dr = {dr}')

# Was trying to get the limit of the function so we don't have to set a tmax
if tmax == 0:
    x = Symbol('x')
    y=y0*(e**(-alpha*x))*cos(w*x+v0)
    limit = limit(y, x, oo)
else:
    t = list(np.linspace(0, tmax, 500))


if centre == 'LSAC':
    # Harmonic motion graph
    y = [y0*(e**(-alpha*t))*cos(w*t+v0) + diff + convlastavg for t in t]
    # Envelopes
    envplus = [y0*e**(-alpha*t) + diff + convlastavg  for t in t]
    envminus = [-y0*e**(-alpha*t) + diff  + convlastavg for t in t]
else:
    # If statements try to adjust the harmonic motion and envelopes to the data. If no adjustment is required, used the commented code below without "diff" 
    # Harmonic motion graph
    if abs(df[vcol].max()) > abs(df[vcol].min()):
        y = [y0*(e**(-alpha*t))*cos(w*t+v0)+diff for t in t]
    else:
        y = [y0*(e**(-alpha*t))*cos(w*t+v0)-diff for t in t]

    # Envelopes
    if abs(df[vcol].max()) > abs(df[vcol].min()):
        envplus = [y0*e**(-alpha*t)+diff for t in t]
        envminus = [-y0*e**(-alpha*t)+diff for t in t]
    else:
        envplus = [y0*e**(-alpha*t)-diff for t in t]
        envminus = [-y0*e**(-alpha*t)-diff for t in t]


model = np.poly1d(np.polyfit(t, y, 2))
z = [model(np.linspace(0, tmax, 1)) for t in t]
# z = [y0*e**(-alpha*t) + diff + convlastavg  for t in t]
# =============================================================================
#     # Harmonic motion graph
#     y = [y0*(e**(-alpha*t))*cos(w*t+v0) - diff  for t in t]
#     # Envelopes
#     envplus = [y0*e**(-alpha*t) + diff   for t in t]
#     envminus = [-y0*e**(-alpha*t) + diff   for t in t]
# =============================================================================

# Plot
# harmonic = vcol #  options are values, std, avg, var
plt.scatter(dfc['scale'], dfc[vcol], s=5, c= '#80afe8', label = f'{vcol}')
plt.scatter(convavg['scale'], convavg[vcol+' circmean'], s=15, c = '#043d82')
plt.plot(t,envplus,color='red', linestyle='dashed', label = 'Envelope +')
plt.plot(t,envminus,color='red', linestyle='dashed', label = 'Envelope -')
plt.plot(t,y,color = 'green', linewidth=0.5, linestyle='solid', label = f'DHM m={m}, k={k}, c={c}')
plt.plot(t, z)
plt.xlabel(f'Scale')
if (vcol == "std" or vcol == "var"):
    plt.ylabel(f'{vcol}')
if (vcol == "dip" or vcol == "dipdir"):
    plt.ylabel(f'{vcol} (°)')
else:
    plt.ylabel(f'{vcol} (mm)')

plt.legend(loc = 1, fontsize='small') # Location Code:'best'0, 'upper right'1, 'upper left'2, 'lower left'3, 'lower right'4, 'right'5, 'center left'6, 'center right'7, 'lower center'8, 'upper center'9, 'center'10. fontsizeint or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
if save == True:
    if centre == 'LSAC':
        plt.title(f'Damped Harmonic Motion Last Scale Average Centred\n{inputfilename}', wrap = True)
        plt.savefig(f'{wd}DHM-LSAC_{inputfilename}{extension_png}', bbox_inches = 'tight')
    elif centre == 'zero':
        plt.title(f'Damped Harmonic Motion 180deg centred\n{inputfilename}', wrap = True)
        plt.savefig(f'{wd}DHM-180C_{inputfilename}{extension_png}', bbox_inches = 'tight')
    else:
        plt.title(f'Damped Harmonic Motion\n{inputfilename}', wrap = True)
        plt.savefig(f'{wd}DHM_{inputfilename}{extension_png}', bbox_inches = 'tight')
plt.show()

# =============================================================================
# # Polar Plot
# df = df.sort_values(by=['scale', vcol])
# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# # Convert ticks to radians https://stackoverflow.com/questions/21172228/python-matplotlib-polar-plots-with-angular-labels-in-radians
# if radians_ticks == True:
#     xT=plt.xticks()[0]
#     xL=['0', 
#         r'$\frac{\pi}{4}$', 
#         r'$\frac{\pi}{2}$',
#         r'$\frac{3\pi}{4}$',
#         r'$\pi$',
#         r'$\frac{5\pi}{4}$',
#         r'$\frac{3\pi}{2}$',
#         r'$\frac{7\pi}{4}$']
#     plt.xticks(xT, xL)
# for scale in df['scale'].unique():
#     f = df[df['scale'] == scale]
#     f = f.sort_values(by=['scale', vcol])
#     f[vcol] = f[vcol].apply(radians) # Convert to radians so it plots in the correct position in the polar plot.
#     ax.scatter(f[vcol], f['scale'], s = 15, label = scale)
#     # DHM line
# rads = [radians(i) for i in y]
# rads.sort(); t.sort()
# ax.plot(rads, t, 'red', lw=1)
# ax.set_rmax(2)
# ax.set_rticks(df['scale'])
# ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
# ax.grid(True)
# ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
# ax.set_theta_zero_location("N")  # theta=0 at the top
# ax.set_theta_direction(-1)  # theta increasing clockwise
# ax.set_title(f'Damped Harmonic Motion-Polar\n{inputfilename}', va='bottom')
# 
# if save == True:
#     plt.savefig(f'{wd}DHM-Polar_{inputfilename}{extension_png}', bbox_inches = 'tight')
# plt.show()
# 
# =============================================================================
