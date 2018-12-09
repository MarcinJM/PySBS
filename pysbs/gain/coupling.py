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
from dolfin import (FunctionSpace, project, VectorElement, as_vector, dx,
                    interpolate, Function, triangle, Constant, dot, assemble)

def displacement_at_boundary(boundary, u):
    """ calculate the displacement at the boundary, return numpy vector """
    (ur, ui) = u.split()
    ur_b = np.zeros((boundary.Nfacets, 3))
    ui_b = np.zeros((boundary.Nfacets, 3))
    for idx, point in enumerate(boundary.midpoints):
        ur_b[idx] = ur(point)
        ui_b[idx] = ui(point)

    return (ur_b, ui_b)

def boundary_force_coupling(fr, fi, ur, ui, boundary):
    """ caclulate f*u* """
    coupling = np.multiply((ur - 1j*ui), (fr + 1j*fi))
    coupling = np.sum(coupling, axis = 1) #sum vector components
    coupling = np.dot(coupling, boundary.facet_lengths)
    #coupling = np.square((np.absolute(coupling)))
    return coupling


def bulk_force_coupling(fr_bulk, fi_bulk, u):
    """ caclulate |f*u*|^2 for bulk electrostriction"""
    (ur, ui) = u.split()
    coup_r = dot(ur,fr_bulk)*dx + dot(ui,fi_bulk)*dx
    coup_i = -dot(ui,fr_bulk)*dx + dot(ur,fi_bulk)*dx
    int_coup_r = assemble(coup_r)
    int_coup_i = assemble(coup_i)
    #coupling_bulk = int_coup_r**2 + int_coup_i**2
    coupling = int_coup_r + 1j*int_coup_i
    return coupling


def bulk_force_gain(Q, power_opt, power_mech, omega_opt, fr_bulk, fi_bulk, u):
    """ calculate gain due to bulk electrostriction 
        returns gain in units of  W-1 m-1 """
    coupling_bulk = bulk_force_coupling(fr_bulk, fi_bulk, u)
    gain_bulk = Q*omega_opt/(4.0*power_opt**2*power_mech)*10**21*np.absolute(coupling_bulk)**2
    return gain_bulk


def boundary_force_gain(Q, power_opt, power_mech, omega_opt, fr_bdr, fi_bdr, ur_bdr, ui_bdr, boundary):
    """ calculate gain due to boundary forces
         returns gain in units of   W-1 m-1 """
    gain_bdr =  (Q*omega_opt/(4.0*power_opt**2*power_mech)*10**21*
            np.absolute(boundary_force_coupling(fr_bdr, fi_bdr, ur_bdr, ui_bdr, boundary))**2)
    return gain_bdr






