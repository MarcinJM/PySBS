#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Tis file is part of PySBS.
    PySBS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    PySBS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    
    If you use it in research please cite an appropriate paper from README file.
    
    If you have any questions email: marcinjmalinowski@gmail.com
    
    author: Marcin Malinowski
"""
import numpy as np
from math import cos, sin, pi
import itertools


Voigt_index_maping = {0:(0,0), 1:(1,1) , 2:(2,2), 3:(1,2), 4: (0,2), 5:(0,1)}



def RotationMatrixYaxis(theta):
    """ define the 3D rotation matrix around y axis """
    R = np.array([[cos(theta), 0.0, sin(theta)],
              [0.0, 1.0, 0.0],
              [-sin(theta), 0.0, cos(theta)]])
    return R

Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

def full_3x3_to_Voigt_6_index(i, j):
    if i == j:
        return i
    return 6-i-j


def Voigt_6x6_to_full_3x3x3x3(C):
    
    C = np.asarray(C)
    C_out = np.zeros((3,3,3,3), dtype=float)
    for i, j, k, l in itertools.product(range(3), range(3), range(3), range(3)):
        Voigt_i = full_3x3_to_Voigt_6_index(i, j)
        Voigt_j = full_3x3_to_Voigt_6_index(k, l)
        C_out[i, j, k, l] = C[Voigt_i, Voigt_j]
        
    return C_out
    
def full_3x3x3x3_to_Voigt_6x6(C):

    C = np.asarray(C)
    Voigt = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            Voigt[i,j] = C[k,l,m,n]

    return Voigt



def Rotate_Voigt_tensor(T, angle):
    TT =  Voigt_6x6_to_full_3x3x3x3(T)
    R = RotationMatrixYaxis(angle)
    PP = np.einsum('ai,bj,ck,dl,abcd->ijkl', R, R, R, R, TT)
    P = full_3x3x3x3_to_Voigt_6x6(PP)
    return P





if __name__ == "__main__":
    from material import CubicStiffness
    
    c44 = 79.6
    c12 = 63.9
    c11 = 165.7    

    C = np.array([[194.5, 64.1, 35.7, 0, 0, 0],
                       [64.1, 165.7, 64.1, 0, 0, 0],
                       [35.7, 64.1, 194.5, 0, 0, 0],
                       [0  , 0  , 0  , 79.6, 0  ,   0],
                       [0  , 0  , 0  , 0  , 50.9,   0],
                       [0  , 0  , 0  , 0  ,   0, 79.6]])



                   
    C100 = CubicStiffness(c11, c12,c44)
    C110 = Rotate_Voigt_tensor(C100, pi/4.0)
    print(C110)
    print(C)                          
                                         
                                         
                                         
