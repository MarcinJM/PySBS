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
# fully tensorial elastic mode solver

from dolfin import (FunctionSpace, TrialFunction, TestFunction, VectorElement, 
assemble, Function, split, SLEPcEigenSolver, dot, as_matrix)
from pysbs.el.weak_form import RHS, LHS
from pysbs.el.clamped_walls import build_clamped_walls
from pysbs.mode_solver import ModeSolver
from pysbs.material import Material
from pysbs.el.clamped_walls import build_clamped_walls
import collections
from pysbs.el.plot import plot_displacement
from math import pi

# log level
#set_log_level(10) # for debugging


class ELSolver(ModeSolver):
    """ fully tensorial elastic mode solver  """

    def __init__(self,mesh,materials):
        super(ELSolver, self).__init__(mesh, materials)
        self.q_b = 0.0
        self.lagrange_order = 2
        self._VComplex = None

    def assemble_matrices(self):
        V = VectorElement("Lagrange", self.mesh.ufl_cell(),
                          self.lagrange_order, dim = 3)
        self._VComplex = FunctionSpace(self.mesh, V*V)
        u = TrialFunction(self._VComplex)
        (ur, ui) = split(u)
        v = TestFunction(self._VComplex)
        (vr, vi) = split(v)

        A_ij = None
        B_ij = None

        # construct matrices for each domain
        if not isinstance(self.materials, collections.Iterable):
            material = self.materials
            C = as_matrix(material.el.C)
            rho = material.el.rho
            DX  = material.domain
            A_ij = LHS(self.q_b, C, ur, ui, vr, vi, DX)
            B_ij = RHS(rho, ur, ui, vr, vi, DX)


        else:
            counter = 0 # a work around for summing forms
            for material in self.materials:
                C = as_matrix(material.el.C)
                rho = material.el.rho
                DX  = material.domain
                
                if counter == 0:
                    A_ij = LHS(self.q_b, C, ur, ui, vr, vi, DX)
                    B_ij = RHS(rho, ur, ui, vr, vi, DX)
                    counter += 1
                else:
                    a_ij = LHS(self.q_b, C, ur, ui, vr, vi, DX)
                    b_ij = RHS(rho, ur, ui, vr, vi, DX)

                    A_ij += a_ij
                    B_ij += b_ij


        # assemble the matrices
        assemble(A_ij, tensor=self._A)
        assemble(B_ij, tensor=self._B)

    def setup_solver(self):
        shift = (2.0*pi*self.eigenmode_guess)**2
        self.esolver = SLEPcEigenSolver(self._A, self._B) 
        self.esolver.parameters["solver"] = "krylov-schur"
        self.esolver.parameters["tolerance"] = 1e-12
        self.esolver.parameters["spectral_transform"] = "shift-and-invert"
        self.esolver.parameters["spectral_shift"] = shift
        self.esolver.parameters["spectrum"] = "target magnitude"


    def set_clamped_walls(self):
        bcs = build_clamped_walls(self._VComplex)
        bcs.apply(self._A)
        bcs.apply(self._B)

    def set_free_walls(self):
        # this is the natural BC when no changes are made
        # to the assembled matrices
        pass
    
    
    def calculate_power(self, eigenvalue_index):
        
        r, c, rx, cx = self.esolver.get_eigenpair(eigenvalue_index)
        u = Function(self._VComplex, rx)
        omega_sqrd = r*(1e9)**2

         # construct matrices for each domain
        if not isinstance(self.materials, collections.Iterable):
            material = self.materials
            rho = material.el.rho
            DX  = material.domain
            norm = assemble(dot(u, rho*u)*DX)
            p = 0.5*omega_sqrd*norm
            return p
           

        else:
            p = 0.0
            for material in self.materials:
                rho = material.el.rho
                DX  = material.domain
                norm = assemble(dot(u, rho*u)*DX)
                p +=  0.5*omega_sqrd*norm
            return p


    def compute_eigenvalues(self):
        self.esolver.solve(self.n_modes)
        if self.esolver.get_number_converged()==0:
            print( 'Eigensolver did not calculate eigenmodes. Increase log level')
        for i in range(0,self.esolver.get_number_converged(),1):
            r, c = self.esolver.get_eigenvalue(i)
            try:
                f_b = (r)**0.5/(2*pi)
                print('Frequency {} GHz'.format(f_b))

                if self.plot_eigenmodes == True:
                    # TODO save plots to a folder
                    r, c, xr, xi = self.esolver.get_eigenpair(i)
                    f = Function(self._VComplex, xr)
                    plot_displacement(self.mesh, f)
            except:
                print("Found negative frequency")

    def extract_field(self, eigenvalue_index):
         r, c, xr, xi = self.esolver.get_eigenpair(eigenvalue_index)
         f = Function(self._VComplex, xr)
         f_b = (r)**0.5/(2*pi)
         return (f, f_b)


