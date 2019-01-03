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
from pysbs.gain.coupling import displacement_at_boundary, boundary_force_gain
from pysbs.em.forces.electrostriction import (bulk_electrostriction, boundary_electrostriction)
from pysbs.em.forces.stress import boundary_stress
from pysbs.gain.plot import plot_bulk_electrostriction, plot_boundary_force


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


class Forces():
    """ a class for calculating forces acting on a boundary 
       input units same as in the solvers lam = [um], freq_mech = [GHz]
       TODO implement multiple boundaries, right now only one
       is possible for boundary force calculation """
    
    def __init__(self, E, u, q_b, Q, power_opt, power_mech, direction, freq_mech, lam, boundary):
        self.E = E
        self.u = u
        self.boundary = boundary
        self.direction = direction
        self.power_opt = power_opt
        self.power_mech = power_mech
        self.q_b = q_b
        self.Q = Q
        self.omega_mech = freq_mech*2.0*pi*1e9
        self.omega_opt = 2.0*pi*c0/(lam*1e-6)
        self.u_bdr = None
        self.f_bdr_electst = None
        self.f_bdr_stress = None
        self.u_bdr = None
        self.interior_material = None
        self.exterior_material = None
        self.gain_bdr_elcst = None
        self.gain_bdr_stress = None
        self.displacement_at_boundary()
        self.assign_materials()
        
    def displacement_at_boundary(self):
        self.u_bdr = displacement_at_boundary(self.boundary, self.u)
        
    def assign_materials(self, materials):
        interior = self.boundary.interior
        exterior = self.boundary.exterior
        for material in materials:
            if material.domain.subdomain_id() == interior:
                self.interior_material = material
            elif material.domain.subdomain_id() == exterior:
                self.exterior_material = material
            else:
                pass
            
    def calculate_boundary_electrostriction(self):
        self.f_bdr_electst = boundary_electrostriction(self.E, 
                                                       self.exterior_material.em.p,
                                                       self.interior_material.em.p, 
                                                       self.exterior_material.em.e_r, 
                                                       self.interior_material.em.e_r, 
                                                       self.direction, 
                                                       self.boundary, 
                                                       offset = 1e-12)
        
    def calculate_boundary_electrostriction_gain(self):
        self.gain_bdr_elcst =  boundary_force_gain(self.Q,
                                                   self.power_opt, 
                                                   self.power_mech, 
                                                   self.omega_opt, 
                                                   self.f_bdr_elctst[0], 
                                                   self.f_bdr_elctst[0], 
                                                   self.u_bdr[0], 
                                                   self.u_bdr[1],
                                                   self.boundary)
        
        
    def calculate_boundary_stress(self):
        self.f_bdr_stress = boundary_stress(self.E,
                                  self.exterior_material.em.e_r,
                                  self.interior_material.em.e_r,
                                  self.direction, 
                                  self.boundary, 
                                  offset = 1e-12)
    
    def calculate_boundary_stress_gain(self):
        self.gain_bdr_stress =  boundary_force_gain(self.Q,
                                                    self.power_opt,
                                                    self.power_mech,
                                                    self.omega_opt,
                                                    self.f_bdr_stress[0], 
                                                    self.f_bdr_stress[1],
                                                    self.u_bdr[0], 
                                                    self.u_bdr[1], 
                                                    self.boundary)
        
    def calculate_bulk_electrostriction(self):
        pass
    
    def calculate_bulk_gain(self):
        pass
    
    
    def calculate_forces(self):
        self.calculate_boundary_electrostriction()
        self.calculate_boundary_stress()
        self.calculate_bulk_electrostriction()
    
    
    def plot_forces():
        plot_boundary_force(self.f_bdr_stress[0], self.boundary, self.power_opt)
        plot_boundary_force(self.f_bdr_stress[0], self.boundary, self.power_opt)
    















