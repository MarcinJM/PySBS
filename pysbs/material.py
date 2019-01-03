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
from dolfin import as_matrix
import numpy as np

# stiffness tensor in Voigt notation
def IsotropicStiffness(M,G):
    C = np.array([[M, M-2*G, M-2*G, 0, 0, 0],
                   [M-2*G, M, M-2*G, 0, 0, 0],
                   [M-2*G, M-2*G, M, 0, 0, 0],
                   [0  , 0  , 0  , G, 0  ,   0],
                   [0  , 0  , 0  , 0  , G,   0],
                   [0  , 0  , 0  , 0  ,   0, G]])
    return C

# stiffness tensor in Voigt notation
def CubicStiffness(c11, c12, c44):
    C = np.array([[c11, c12, c12, 0, 0, 0],
                   [c12, c11, c12, 0, 0, 0],
                   [c12, c12, c11, 0, 0, 0],
                   [0  , 0  , 0  , c44, 0  ,   0],
                   [0  , 0  , 0  , 0  , c44,   0],
                   [0  , 0  , 0  , 0  ,   0, c44]])
    return C



# photoelastic tensor in Voigt notation
def IsotropicPhotoelasticity(p11, p12):
    C = np.array([[p11, p12, p12, 0, 0, 0],
                   [p12, p11, p12, 0, 0, 0],
                   [p12, p12, p11, 0, 0, 0],
                   [0  , 0  , 0  , 0.5*(p11-p12)  , 0,  0],
                   [0  , 0  , 0  , 0  ,  0.5*(p11-p12) ,   0],
                   [0  , 0  , 0  , 0  ,   0,  0.5*(p11-p12) ]])
    return C

# photoelastic tensor in Voigt notation
def CubicPhotoelasticity(p11, p12, p44):
    C = np.array([[p11, p12, p12, 0, 0, 0],
                   [p12, p11, p12, 0, 0, 0],
                   [p12, p12, p11, 0, 0, 0],
                   [0  , 0  , 0  , p44, 0  ,   0],
                   [0  , 0  , 0  , 0  , p44,   0],
                   [0  , 0  , 0  , 0  ,   0, p44]])
    return C



class ElasticProperties():
    """ container for elastic properties of a material """

    def __init__(self, rho, stiffness):
        self.rho = rho
        self.C = stiffness
    
class ElectromagneticProperties():
    """ containter for electromagnetic properties """
    
    def __init__(self, u_r, e_r, p):
        self.u_r = u_r
        self.e_r = e_r
        self.p = p

class Material():
    """ class for storing all material properties 
        and associated domain """

    def __init__(self, domain):
        self.em = ElectromagneticProperties(1.0, 1.0, CubicPhotoelasticity(1.0, 1.0, 1.0))
        self.el = ElasticProperties(1.0, IsotropicStiffness(3.0, 1.0))
        self.domain = domain
        



