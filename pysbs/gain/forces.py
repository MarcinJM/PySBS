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


from pysbs.gain.coupling import displacement_at_boundary, boundary_force_gain
from pysbs.gain.electrostriction import (bulk_electrostriction, boundary_electrostriction)
from pysbs.gain.stress import boundary_stress
from pysbs.gain.plot import plot_bulk_electrostriction, plot_boundary_force
from pysbs.gain.coupling import (displacement_at_boundary, boundary_force_gain, 
                                 bulk_force_gain, bulk_force_coupling)
from math import pi
from scipy.constants import c as c0


class Forces():
    """ a class for calculating forces acting on a waveguide 
       Input units same as in the solvers lam = [um], freq_mech = [GHz]
       The class initialization requires boundary class from which it finds
       relevant materials.
       For waveguides surrounded by air set exclude_exterior = True,
       in that case also project the E field on the same mesh as u using
       project_Efield_on_submesh().
       
       TODO implement multiple boundaries, right now only one
       is possible for boundary force calculation """
    
    def __init__(self, E, u, q_b, Q, power_opt, power_mech, direction, 
                 freq_mech, lam, boundary, materials, exclude_exterior = True):
        
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
        self.f_bulk_elcst= None
        self.f_bdr_elcst = None
        self.f_bdr_stress = None
        self.u_bdr = None
        self.interior_material = None
        self.exterior_material = None
        self.gain_bdr_elcst = None
        self.gain_bdr_stress = None
        self.gain_bulk_elcst = None
        self.exclude_exterior = exclude_exterior
        self.displacement_at_boundary()
        self.assign_materials(materials)
        
    def displacement_at_boundary(self):
        self.u_bdr = displacement_at_boundary(self.boundary, self.u)
        
    def project_Efield(self):
        pass
        
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
        self.f_bdr_elcst = boundary_electrostriction(self.E, 
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
                                                   self.f_bdr_elcst[0], 
                                                   self.f_bdr_elcst[0], 
                                                   self.u_bdr[0], 
                                                   self.u_bdr[1],
                                                   self.boundary)
        return self.gain_bdr_elcst
        
        
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
        return self.gain_bdr_stress
        
    def calculate_bulk_electrostriction(self):
        if self.exclude_exterior == True:
            self.f_bulk_electst = bulk_electrostriction(self.E,
                                                      self.interior_material.em.e_r,
                                                      self.interior_material.em.p, 
                                                      self.direction, 
                                                      self.q_b)
        else: 
            pass #TODO implement!
    
    def calculate_bulk_gain(self):
        if self.exclude_exterior == True:
            self.gain_bulk_elcst = bulk_force_gain(self.Q, 
                                                  self.power_opt, 
                                                  self.power_mech, 
                                                  self.omega_opt, 
                                                  self.f_bulk_elcst[0], 
                                                  self.f_bulk_elect[1], 
                                                  self.u, 
                                                  dx)
            return self.gain_b
        else:
            pass
    
    def calculate_total_gain(self):
        pass
    
    def calculate_forces(self):
        self.calculate_boundary_electrostriction()
        self.calculate_boundary_stress()
        self.calculate_bulk_electrostriction()
    
    
    def plot_forces():
        plot_boundary_force(self.f_bdr_elcst[0], self.boundary, self.power_opt)
        plot_boundary_force(self.f_bdr_stress[0], self.boundary, self.power_opt)
        if self.exclude_exterior == True:
            pass
        else:
            pass















