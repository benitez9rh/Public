# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 23:15:23 2022

@author: s2132627
"""
# ODE Pendulum https://www.youtube.com/watch?v=p_di4Zn4wz4&t=59s
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from time import ctime

fast = False
save = False
g = 9.8             # acceleration of gravity m/s^2
L = 2               # Length of pendulum (m)
# W = 2 # width
# H = 2 # height

mu = 0.1            # air resistance constant
delta_t = 0.01      # time interval in seconds
theta_0 = 60        # Initial angle position. type value in degrees
omega_0 = 0         # initial angular velocity (m/s), angle theta first derivative. type value in degrees
time_steps = 10     # number of time steps if using linspace
tfinal = 100


theta_0 = math.radians(theta_0)
omega_0 = math.radians(omega_0)
theta = theta_0
omega = omega_0

def duration():
    finish = time.time()
    days = math.floor( (finish-stop0)/86400)
    hours = math.floor( (finish-stop0)/3600 - days*24 )
    minutes = math.floor( (finish-stop0)/60 - (days*24+hours)*60)
    seconds = math.floor( (finish-stop0) - ((days*24+hours)*60+minutes)*60 )
    print(f'\n**Duration:**\ndays: {days} \nhours: {hours} \nminutes: {minutes} \nseconds: {seconds}')
    
def pendulum_acceleration(theta, omega):        # second derivative of theta, theta_double_dot
    pa = -mu * omega - (g/L) * np.sin(theta)
    return pa

def theta_calcs(t):
    global theta, omega, g, L, mu
    pa = pendulum_acceleration(theta, omega)
    theta += omega * delta_t
    omega += pa * delta_t
    return theta, omega

def theta_calcf(t):
    global theta, omega, g, L, mu
    pa = pendulum_acceleration(theta, omega)
    theta = omega * t
    omega = pa * t
    return theta, omega

def arc(theta, L): # x is the arc the pendulum makes with the position at rest
    arcx = L*theta
    return arcx

def acceleration(theta, g): # acceleration = x second derivative, x double dot
    global a
    a = -g*np.sin(theta)
    return a

stop0 = time.time()

t = np.arange(0, tfinal, delta_t)
# t = np.linspace(0, delta_t, time_steps)

fast = True

""" Fast way """
if fast == True:
    th_c = np.vectorize(theta_calcf); th, o = th_c(t) # thc becomes a function to apply the function fed to np.vectorize and applies it to t https://thispointer.com/apply-a-function-to-every-element-in-numpy-array/
    arc_c = np.vectorize(arc); arcdist = arc_c(th, L)
    acceleraction_c = np.vectorize(acceleration); ac = acceleraction_c(th, g)
    f = np.vstack(( (t, o, th, ac, arcdist) ))
    speed = "Fast"
    """ Slow way """
else:
    speed = "Slow"
    th = []
    ac = []
    arcdist = []
    o = []
    
    for time in t:
        theta, omega = theta_calcs(time)
        th.append(theta)
        o.append(omega)
        pa = pendulum_acceleration(theta, omega)
        ac.append(pa)
        arcdist.append(arc(theta, L))


plt.plot( t, o, color = 'blue',label = 'Omega (m/s)' )
plt.plot( t, th, color = 'red',label = 'Theta (radians)' )
plt.plot( t, ac, color = 'green',label = r'$Acceleration (m/s^2)$' ) # https://matplotlib.org/2.0.2/users/mathtext.html
plt.plot( t, arcdist, color = 'black',label = 'Arc Distance (m)' )


plt.xlabel(r'$\bf time (s)$')
# plt.ylabel(r'$\gamma \bf(h)$')
plt.title(f'PendulumODE {speed}')
plt.xlim([0, 10 ])
# plt.ylim([0, np.amax(gamf[gamf!=tmax])+0.25 ])
plt.legend(loc='upper right')
plt.grid(True)
#plt.subplots_adjust(left=0.0, bottom=0.0, right=2.0, top=4.2, wspace=0.2, hspace=0.3)
if save == True:
    plt.savefig(f'{wdsave}{inputfilename}_Res_Semi-variograms_lag{dxlag}_{int(azi_mat[iazi])}deg{extension_png}',dpi=300, bbox_inches = "tight")
plt.show()
    
# duration()
    
    
    
    





    