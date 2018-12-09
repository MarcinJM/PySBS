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
from scipy.constants import c as c0
from math import pi
from dolfin import (FunctionSpace, project, as_vector, interpolate, Function)




def project_Efield_on_submesh(em_solver, eigenvalue_index, submesh):
    """ project E field onto a submesh 
        multiplies Ez from solver to gamma Ez (true field) 
        assume Ez is imaginary as in the em solver"""
    r, c, rx, cx = em_solver.esolver.get_eigenpair(eigenvalue_index)
    e_raw = Function(em_solver._combined_space, rx)
    Esubspace = FunctionSpace(submesh, em_solver._finite_element)
    E_raw_sub = interpolate(e_raw, Esubspace)
    gamma = (abs(r))**0.5
    E_sub = as_vector([E_raw_sub[0], E_raw_sub[1], gamma*E_raw_sub[2]])
    E_sub = project(E_sub)
    return E_sub