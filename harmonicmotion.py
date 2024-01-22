# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 17:51:34 2022

@author: s2132627
"""
import os
import pathlib
import pandas as pd
from math import cos, e, sqrt, floor, ceil
import matplotlib.pyplot as plt
from sympy import limit, Symbol, oo
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit

extension_png = ".png"
extension_csv = ".csv"
wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals\Upscale').as_posix() #flips the backslashes to plug in open()
bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
inputfilename = pathlib.PureWindowsPath(r'Fracs').as_posix()
os.chdir(wd)                                   # set the working directory

wdsave = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals\Upscale\DHM').as_posix() #flips the backslashes to plug in open()
bs="\\"; wdsave=wdsave+bs  

df = pd.read_csv(wd + inputfilename + extension_csv)                  #This new line is because reading the .csv file with a header and an index column will create some issue


# =============================================================================
# # If df does not contain scale 0
# wd = pathlib.PureWindowsPath(r'C:\Users\s2132627\OneDrive - University of Edinburgh\The University of Edinburgh\PhD\Step 1 - Variograms\Greywacke scans\Q4Aperture\Residuals').as_posix() #flips the backslashes to plug in open()
# bs="\\"; wd=wd+bs                               # Use this instead in Windows                             
# inputfilename = pathlib.PureWindowsPath(r'z_ORIGvsORIG-rot(-0.02, 0.31, 0)_ApertureMap').as_posix()
# df0 = pd.read_csv(wd + inputfilename + extension_csv, header=0, index_col=0)                  #This new line is because reading the .csv file with a header and an index column will create some issue
# df0['ApertureR'] = df0['ApertureR']*1000 # raw data is in m whereas the rest is in mm  
# df0['Scale'] = 0
# df = pd.concat([df, df0], axis=0)
# =============================================================================


""" User-defined Variables """
vcol = "ApertureR"
m = 1                # mass (kg)
k =  10                  # spring stiffness
c = 3*df[df['Scale']==0][vcol].std()                  # damping coefficient # Use mean/std ratio for the initial raw data
v0 = 0              # initial velocity? (initial condition)
y0 = 3*df[df["Scale"]==0][vcol].std()            # initial releasing position (initial condition)
ctr = df.groupby(["Scale"])[vcol].mean().iloc[-1] # The centre DHM of the axis and corresponding envelopes. Default is 0. Only active when auto_fit_centre_to_data == False
save = False
auto_fit_centre_to_data = False # Automatically fit the DHM axis to the centre of the data. Otherwise, the centre is controlled by ctr
auto_y0 = False
"""
################################################################################
"""
""" Calculated Variables """
tmax = df['Scale'].max()              # how long to run for. TO DO: Use 0 to calculate the limit
if auto_y0 == True:
    y0 = max(abs(df[vcol].max()),abs(df[vcol].min()))                  # initial releasing position (initial condition)
dr = c/(2*sqrt(m*k))    # damping ratio (ζ)
                        # undamped ( ζ=0 , blue)
                        # underdamped ( 0<ζ<1 , orange)
                        # critically damped ( ζ=1 , green)
                        # overdamped ( ζ>1 , red)
alpha = c/(2*m)
alphas = alpha**2
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
    t = list(np.linspace(0, tmax, 1000))
if auto_fit_centre_to_data == True:
    ctr = (abs(df[vcol].max()) - abs(df[vcol].min()))/2
    # ctr = df.groupby(["scale"])[vcol].mean().iloc[-1] # If you wish to get the first or last scale average of vcol
    # If statements try to adjust the harmonic motion and envelopes to the data. If no adjustment is required, used the commented code below without "diff" 
    # Harmonic motion graph
    if abs(df[vcol].max()) > abs(df[vcol].min()):
        y = [y0*(e**(-alpha*t))*cos(w*t+v0)+ctr for t in t]
    else:
        y = [y0*(e**(-alpha*t))*cos(w*t+v0)-ctr for t in t]
    # Envelopes
    if abs(df[vcol].max()) > abs(df[vcol].min()):
        envplus = [y0*e**(-alpha*t)+ctr for t in t]
        envminus = [-y0*e**(-alpha*t)+ctr for t in t]
    else:
        envplus = [y0*e**(-alpha*t)-ctr for t in t]
        envminus = [-y0*e**(-alpha*t)-ctr for t in t]
else:
    # Harmonic motion graph
    y = [y0*(e**(-alpha*t))*cos(w*t+v0)+ctr for t in t]
    # Envelopes
    envplus = [y0*e**(-alpha*t)+ctr for t in t]
    envminus = [-y0*e**(-alpha*t)+ctr for t in t]
# Calculate scales' averages, stdevs, variances
scaleASV = pd.DataFrame(list(zip(df.groupby(['Scale'])['Scale'].min(), df.groupby(['Scale'])[vcol].mean(), df.groupby(['Scale'])[vcol].std(), df.groupby(['Scale'])[vcol].var(), df.groupby(['Scale'])[vcol].mean() / df.groupby(['Scale'])[vcol].std(), df.groupby(['Scale'])[vcol].count() )), columns = ['Scale', vcol+" Mean", vcol+" StandardDeviation", vcol+" Variance", "ratiomstd", "Count"]) # Calculate Circular average per scale
stddf = scaleASV[["Scale", f"{vcol} Mean", f"{vcol} StandardDeviation"]]
# Plot Damped Harmonic Motion (DHM)
harmonic = vcol #  options are values, std, avg, var
plt.scatter(df['Scale'], df[vcol], s=2, c= '#80afe8', label = f'{vcol}')
plt.scatter(scaleASV['Scale'], scaleASV[vcol+" Mean"], s=15, c = '#043d82', label = f'{vcol} mean')
for i in range(1,4):    # Include -3 to +3 standard deviations
    stddf[f"m plus {i} StD"] = stddf[f"{vcol} Mean"] + i * stddf[f"{vcol} StandardDeviation"]
    stddf[f"m minus {i} StD"] = stddf[f"{vcol} Mean"] - i * stddf[f"{vcol} StandardDeviation"]
    plt.scatter(stddf["Scale"], stddf[f"m plus {i} StD"], s=5, c = 'yellow')
    plt.scatter(stddf["Scale"], stddf[f"m minus {i} StD"], s=5, c = 'yellow')
plt.plot(t, envplus,color='red', linestyle='dashed', label = 'Envelope +')
plt.plot(t, envminus,color='red', linestyle='dashed', label = 'Envelope -')
plt.plot(t, y, color = 'green', linewidth=0.5, linestyle='solid', label = f'DHM m={m}, k={k}, c={c:.2f}, y0={y0:.2f}')
plt.xticks(df['Scale'].unique())
plt.xlabel(f'Scale')
if (vcol == "std" or vcol == "var"):
    plt.ylabel(f'{vcol}')
else:
    plt.ylabel(f'{vcol} (mm)')
plt.title(f'Damped Harmonic Motion\n{inputfilename}')
plt.legend(loc = 1, fontsize='small') # Location Code:'best'0, 'upper right'1, 'upper left'2, 'lower left'3, 'lower right'4, 'right'5, 'center left'6, 'center right'7, 'lower center'8, 'upper center'9, 'center'10. fontsizeint or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}
if save == True:
    plt.savefig(f"{wdsave}{inputfilename}{vcol}DHM{extension_png}", dpi=300, bbox_inches = "tight")
plt.show()


#  Plot Linear best fit with reducing datapoints until r-value is reached (requires sorted data) # https://stackoverflow.com/questions/71056937/fit-a-straight-line-ignoring-data-points-having-ceiling-effect
function = "ExponentialDecay" # Linear or ExponentialDecay
prop = "Mean" # property variance, standard deviation, ratiomstd 
xdata = list(scaleASV['Scale'])
ydata = list(scaleASV[f'{vcol} {prop}'])
xc = np.linspace(df['Scale'].min() , df['Scale'].max(), 1000) # x-continuous
if function == "Linear":
    # It needs to be increasing or else the rvalue is negative and it never reaches the limit set in the while loop
    xdata = list(scaleASV['Scale'])
    ydata = list(scaleASV[f'{vcol} {prop}'])
    fit = linregress(xdata, ydata)
    r_squaredi = fit.rvalue
    # Using an r.value threshold
# =============================================================================
#     i=0
#     while fit.rvalue < 0.995:
#         i += 1
#         print(i)
#         fit = linregress(xdata[:-i], ydata[:-i])
#         print(fit.rvalue)
# =============================================================================
    
    # Using best rvalue (except when only having 2-points which would obviously have an r-value = 1)
    try:
        for i in range(1, len(ydata)-2): # Data needs at least 2 data points to make a model
            #print(i)
            #print(fit.rvalue)
            fit = linregress(xdata[:-i], ydata[:-i])
            if abs(fit.rvalue) > abs(r_squaredi): # absolute because sometimes the r-value is negative because the slope is negative
                iopt = i
                r_squaredf = fit.rvalue
                r_squaredi = fit.rvalue
                print(iopt)
                print(r_squaredf)
    except Exception: # curve_fit sometimes has a runtime error because it can't find a solution after a certain number of iterations and it's stopping the code, hence the exception catcher.
        pass
    finally:
        x2data = xdata[:-iopt]
        y2data = ydata[:-iopt]
        fit = linregress(x2data, y2data)
        yc = [fit.slope*xi+fit.intercept for xi in xc]
        m, b, r_squared, _, _ = fit
elif function == "ExponentialDecay": # https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    gamma = 2 # first guess at exponent
    c = y[0]
    b = y[-1]
    def f(x, c, gamma, b): # c is the initial condition, i.e. y(x=0). b is the assymptote y tends to as lim x->infinity
        return c * np.exp(-gamma * x) + b
        #return e**(-x/gamma)
    popt, pcov = curve_fit(f, np.array(xdata), np.array(ydata) )
    residuals = ydata - f(np.array(xdata), *popt)
    ss_tot = np.sum((ydata - np.mean(ydata))**2)
    ss_res = np.sum(residuals**2)
    r_squaredi = 1 - (ss_res / ss_tot); r_squared = r_squaredi;
    print(f"initial r^2: {r_squaredi}")
    try:
        for i in range(1, len(ydata)-1): # Data needs at least 2 data points to make a model
            popt, pcov = curve_fit(f, np.array(xdata[:-i]), np.array(ydata[:-i]) )
            residuals = ydata[:-i] - f( np.array(xdata[:-i]), *popt)
            ss_tot = np.sum( (ydata[:-i] - np.mean(ydata[:-i]))**2)  
            ss_res = np.sum(residuals**2)
            r_squared = 1 - (ss_res / ss_tot)
            print("i: " + str(i))
            print("r^2: " + str(r_squared))
            if r_squared > r_squaredi:
                iopt = i
                r_squaredf = r_squared
                r_squaredi = r_squared
    except Exception: # curve_fit sometimes has a runtime error because it can't find a solution after a certain number of iterations and it's stopping the code, hence the exception catcher.
        pass
    finally:
        popt, pcov = curve_fit(f, np.array(xdata[:-iopt]), np.array(ydata[:-iopt]) )
        residuals = ydata[:-iopt] - f( np.array(xdata[:-iopt]), *popt)
        ss_tot = np.sum( (ydata[:-iopt] - np.mean(ydata[:-iopt]))**2)  
        ss_res = np.sum(residuals**2)
        r_squaredf = 1 - (ss_res / ss_tot)
        x2data = xdata[:-iopt]
        y2data = ydata[:-iopt]
        c, gamma, b = popt
        yc = [f(xi, *popt) for xi in xc]

plt.plot(xdata, ydata, 'o', label = 'Complete Dataset')
plt.plot(x2data, y2data, '.', label = 'Filtered dataset to fit highest r-value')
plt.plot(xc, yc, '-', label = 'Fit to filtered dataset')
plt.xticks(df['Scale'].unique())
if function == "Linear":
    plt.ylabel(f'Mean/{prop} Ratio');
    if r_squaredf<0:
        plt.annotate("y = {:.2f}x + {:.2f}".format(m, b), (min(xc), min(yc)+4))
        plt.annotate("r-squared = {:.3f}".format(r_squaredf), ((min(xc), min(yc)+2))  )
    else:
        plt.annotate("y = {:.2f}x + {:.2f}".format(m, b), (min(xc), max(yc)-2))
        plt.annotate("r-squared = {:.3f}".format(r_squaredf), ((min(xc), max(yc)-4))  )
elif function == "ExponentialDecay":
    plt.ylabel(f'{prop}')
    plt.annotate("y =  {:.2f} * exp (-{:.2f}*x)+{:.2f}".format(c, gamma, b), (xc.mean(), 0.7*sqrt(max(yc)**2+min(yc)**2)) )
    plt.annotate("r-squared = {:.3f}".format(r_squaredf), (xc.mean(), 0.6*sqrt(max(yc)**2+min(yc)**2))  )
plt.xlabel(f'Scale')
plt.title(f"GW1_Q4 Original Aperture")    #plt.title(f"{inputfilename}")
plt.legend(loc = 1, fontsize='small')
if save == True:
    plt.savefig(f"{wdsave}{inputfilename}{vcol}DHM{extension_png}", dpi=600, bbox_inches = "tight")

# =============================================================================
# #Hist plot
# dfplot = df[df['scale']==0]
# plt.hist(dfplot[vcol], facecolor='red',bins=np.linspace(dfplot[vcol].min(),dfplot[vcol].max(),1000), histtype="stepfilled",alpha=0.1,density=True,cumulative=False,edgecolor='black',label='Original')
# # plt.xlim([0.05,0.25
# # plt.ylim([0,5.0])
# plt.xlabel(f'{vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} {vcol}')
# plt.legend(loc='upper left')
# plt.grid(True)
# 
# scale = 1
# plt.hist(df[df['scale']==scale][vcol], facecolor='blue',bins=np.linspace(df[df['scale']==scale][vcol].min(),df[df['scale']==scale][vcol].max(),1000), histtype="stepfilled",alpha=0.2,density=True,cumulative=False,edgecolor='black',label=f'Scale {scale}')
# # plt.xlim([0.05,0.25
# # plt.ylim([0,5.0])
# plt.xlabel(f'{vcol}'); plt.ylabel('Frequency'); plt.title(f'{inputfilename} {vcol}')
# plt.legend(loc='upper left')
# plt.grid(True)
# =============================================================================

#  Plot Basic Statistics as function of scale
vcolstr = "Residual Aperture" #  options are values, std, avg, var
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax3 = ax1.twinx()
dist = ax1.scatter(df['Scale'], df[vcol], s=2, c= '#80afe8', label = f'{vcolstr}')
mean = ax1.scatter(scaleASV['Scale'], scaleASV[vcol+" Mean"], marker="_", s=50, c = '#043d82', label = f'Mean')
var = ax2.scatter(scaleASV['Scale'], scaleASV[vcol+" Variance"], s=15, c = '#008000', label = f'Variance')
std = ax2.scatter(scaleASV['Scale'], scaleASV[vcol+" StandardDeviation"], s=15, c = '#580F41', label = f'StandardDeviation')
LOGstd = ax3.scatter(scaleASV['Scale'], scaleASV[vcol+" StandardDeviation"], marker="s", s=15, c = '#800080', label = f'Log Standard Deviation')
LOGvar = ax3.scatter(scaleASV['Scale'], scaleASV[vcol+" Variance"], marker="s", s=15, c = '#15B01A', label = f'Log Variance')
for i in range(1,4):    # Include -3 to +3 standard deviations
    stddf[f"m plus {i} StD"] = stddf[f"{vcol} Mean"] + i * stddf[f"{vcol} StandardDeviation"]
    stddf[f"m minus {i} StD"] = stddf[f"{vcol} Mean"] - i * stddf[f"{vcol} StandardDeviation"]
    ax1.scatter(stddf["Scale"], stddf[f"m plus {i} StD"], s=5, c = 'yellow')
    ax1.scatter(stddf["Scale"], stddf[f"m minus {i} StD"], s=5, c = 'yellow')
plt.plot(t, envplus,color='red', linestyle='dashed', label = 'Envelope +')
plt.plot(t, envminus,color='red', linestyle='dashed', label = 'Envelope -')
plt.plot(t, y, color = 'green', linewidth=0.5, linestyle='solid', label = f'DHM m={m}, k={k}, c={c:.2f}, y0={y0:.2f}')
#ax3.set_yscale('log')
#ax3.spines['right'].set_position(('outward', 60)) # Move the 3rd axis away from the 2nd so it's legible
#ax3.set_ylabel(f'Log Variance/StandardDeviation')
ax1.set_xticks(df['Scale'].unique())
ax1.set_xlabel(f'Scale')
if (vcol == "std" or vcol == "var"):
    ax1.set_ylabel(f'{vcolstr}')
else:
    ax1.set_ylabel(f'{vcolstr} (mm)')

ax2.set_ylabel(f'Variance/Standard Deviation')
ax1.set_title(f'GW1Q4 {vcolstr} Upscaling Effects on Basic Statistics')

# Merge all legends from all axes into one
handles, labels = ax1.get_legend_handles_labels() # https://gis.stackexchange.com/questions/379706/python-matplotlib-combine-legend-from-histogram-and-lines
handles.append(var)
handles.append(std)
handles.append(LOGstd)
handles.append(LOGvar)
ax1.legend(handles = handles, loc = "best", fontsize='x-small') # Location Code:'best'0, 'upper right'1, 'upper left'2, 'lower left'3, 'lower right'4, 'right'5, 'center left'6, 'center right'7, 'lower center'8, 'upper center'9, 'center'10. fontsizeint or {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}

if save == True:
    plt.savefig(f"{wdsave}{inputfilename}{vcol}DHM{extension_png}", dpi=1000, bbox_inches = "tight")
plt.show()