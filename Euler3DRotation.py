# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:20:59 2022

@author: s2132627
"""

import numpy as np
import math as m
  
def Rx(theta):
  return np.matrix([[ 1, 0           , 0           ],
                   [ 0, m.cos(theta),-m.sin(theta)],
                   [ 0, m.sin(theta), m.cos(theta)]])
  
def Ry(theta):
  return np.matrix([[ m.cos(theta), 0, m.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-m.sin(theta), 0, m.cos(theta)]])
  
def Rz(theta):
  return np.matrix([[ m.cos(theta), -m.sin(theta), 0 ],
                   [ m.sin(theta), m.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])


# =============================================================================
# phi = 0 # phi is the rotation on the x axis
# theta = 5 # theta is the rotation on the y axis
# psi = 0 # psi is the rotation on the z axis
# =============================================================================
def rotate_df(dataframe, phi, theta, psi):
    phi =  phi * m.pi / 180.0       # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    theta = theta * m.pi / 180.0   # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    psi = psi * m.pi / 180.0       # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    print("phi = ", phi)
    print("theta = ", theta)
    print("psi = ", psi)
      
    R = Rz(psi) * Ry(theta) * Rx(phi)
    print(np.round(R, decimals=2))

    array = dataframe.to_numpy()
    array = array * R
    dataframe = pd.DataFrame(array, columns = ['x', 'y', 'z'])
    return dataframe

def rotate_np(array, phi, theta, psi):
    phi =  phi * m.pi / 180.0       # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    theta = theta * m.pi / 180.0   # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    psi = psi * m.pi / 180.0       # azimuth converted to conventional (i.e. not mathematical) and converted to radians
    print("phi = ", phi)
    print("theta  = ", theta)
    print("psi = ", psi)
      
    R = Rz(psi) * Ry(theta) * Rx(phi)
    print("R = ")
    print(np.round(R, decimals=2))
    
    array2 = np.zeros(array.shape)
    for i,l in enumerate(array):
        array2[i] = l * R.T

    return array2





# =============================================================================
# """ Reverse Operation """
# ###################################################
# import sys
# tol = sys.float_info.epsilon * 10
#   
# if abs(R.item(0,0))< tol and abs(R.item(1,0)) < tol:
#    eul1 = 0
#    eul2 = m.atan2(-R.item(2,0), R.item(0,0))
#    eul3 = m.atan2(-R.item(1,2), R.item(1,1))
# else:   
#    eul1 = m.atan2(R.item(1,0),R.item(0,0))
#    sp = m.sin(eul1)
#    cp = m.cos(eul1)
#    eul2 = m.atan2(-R.item(2,0),cp*R.item(0,0)+sp*R.item(1,0))
#    eul3 = m.atan2(sp*R.item(0,2)-cp*R.item(1,2),cp*R.item(1,1)-sp*R.item(0,1))
#   
# print("phi =", eul1)
# print("theta  =", eul2)
# print("psi =", eul3)
# ####################################################
# =============================================================================
