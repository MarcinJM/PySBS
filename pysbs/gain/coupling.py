#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:31:44 2018
coupling calculation 
need the E field and the U field from both solvers
# run soi in air first
@author: marcin

"""

import numpy as np
from scipy.constants import c as c0
from dolfin import assemble, dot
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

def bulk_force_coupling(fr_bulk, fi_bulk, u):
    """ caclulate |f*u*|^2 for bulk electrostriction"""
    (ur, ui) = u.split()
    coup_r = dot(ur,fr_bulk)*dx + dot(ui,fi_bulk)*dx
    coup_i = -dot(ui,fr_bulk)*dx + dot(ur,fi_bulk)*dx
    int_coup_r = assemble(coup_r)
    int_coup_i = assemble(coup_i)
    #coupling_bulk = int_coup_r**2 + int_coup_i**2
    coupling_bulk = int_coup_r + 1j*int_coup_i
    return coupling_bulk
