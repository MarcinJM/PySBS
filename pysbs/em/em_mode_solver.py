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
# based on FEniCS documentation by Evan Lezar

from dolfin import (FunctionSpace, TrialFunctions, TestFunctions, FiniteElement, 
assemble, Function, split, SLEPcEigenSolver, set_log_level, as_vector, dot )
from pysbs.em.weak_form import em_weak_form, curl_t
from pysbs.em.electric_walls import build_electric_walls
from pysbs.mode_solver import ModeSolver
from pysbs.geometry import EmbeddedWaveguide
from pysbs.material import Material
from pysbs.em.plot import plot_transverse_field
from scipy.constants import epsilon_0 as eps0
from scipy.constants import c as c0
from math import pi
import collections


# log level
#set_log_level(10) # for debugging

class EMSolver(ModeSolver):
    """ electromagnetic mode solver """

    def __init__(self,mesh,materials,  wavelength):
        super(EMSolver, self).__init__(mesh, materials)   
        self.wavelength = wavelength
        self.mesh_unit = 1e-6
        self._combined_space = None
        self._finite_element = None
        self.nedelec_order = 1
        self.lagrange_order = 1
        self._k_o_squared = (2.0*pi/self.wavelength)**2

    def assemble_matrices(self):
        nedelec = FiniteElement( "Nedelec 1st kind H(curl)",
                                 self.mesh.ufl_cell(), self.nedelec_order)
        lagrange = FiniteElement( "Lagrange", self.mesh.ufl_cell(),
                                  self.lagrange_order)
        self._finite_element = nedelec*lagrange
        self._combined_space = FunctionSpace(self.mesh, nedelec*lagrange)
        # define the test and trial functions from the combined space

        (N_i, L_i) = TestFunctions(self._combined_space)
        (N_j, L_j) = TrialFunctions(self._combined_space)

        # construct matrices for each domain
        if not isinstance(self.materials, collections.Iterable):
            material = self.materials
            u_r = material.em.u_r
            e_r = material.em.e_r
            DX  = material.domain
            (A_ij, B_ij) = em_weak_form(N_i, N_j, L_i, L_j, u_r, e_r,
                                            self._k_o_squared,DX)

        else:
            counter = 0 # a work around for summing forms
            for material in self.materials:
                u_r = material.em.u_r
                e_r = material.em.e_r
                DX  = material.domain
                if counter == 0:
                    (A_ij, B_ij) = em_weak_form(N_i, N_j, L_i, L_j, u_r, e_r,
                                            self._k_o_squared, DX)
                    counter += 1
                else:
                    (a_ij, b_ij) = em_weak_form(N_i, N_j, L_i, L_j, u_r, e_r,
                                            self._k_o_squared, DX)
                    A_ij += a_ij
                    B_ij += b_ij

        # assemble the matrices
        assemble(A_ij, tensor=self._A)
        assemble(B_ij, tensor=self._B)


    def setup_solver(self):
        self.esolver= SLEPcEigenSolver(self._A, self._B) 
        shift = -(2.0*pi*self.eigenmode_guess/self.wavelength)**2
        self.esolver.parameters["solver"] = "krylov-schur"
        self.esolver.parameters["tolerance"] = 1e-12
        self.esolver.parameters["spectral_transform"] = "shift-and-invert"        
        self.esolver.parameters["spectrum"] = "target magnitude"
        self.esolver.parameters["spectral_shift"] = shift  


    def set_electric_walls(self):
        bcs = build_electric_walls(self._combined_space)
        bcs.apply(self._A)
        bcs.apply(self._B)

    def extract_normalized_field(self, eigenvalue_index):
        """ extract and rescale the E field from the solver  
            note that Ex, Ey, are real, Ez is imaginary """
        r, c, rx, cx = self.esolver.get_eigenpair(eigenvalue_index)
        e_raw = Function(self._combined_space, rx)
        n_eff = ((-r)**0.5)*self.wavelength/(2*pi)
        gamma = (abs(r))**0.5
        E = as_vector([e_raw[0], e_raw[1], gamma*e_raw[2]])
        return (E, n_eff)


    def calculate_power(self, eigenvalue_index, ng):
        """ 1. note the normalization of Ez in the equations
            2. assume u_r = 1 everywhere  """

        r, c, rx, cx = self.esolver.get_eigenpair(eigenvalue_index)
        f = Function(self._combined_space, rx)
        gamma_sqrd = abs(r)
        (Et, Ez) = split(f)        

         # construct matrices for each domain
        if not isinstance(self.materials, collections.Iterable):
            material = self.materials
            e_r = material.em.e_r
            DX  = material.domain

            int_t = assemble(dot(Et,Et)*DX)
            int_z = assemble(Ez*Ez*DX)*gamma_sqrd 
            power = 0.5*(c0/ng*eps0*e_r*(int_t+ int_z))
            return power

        else:
            power = 0.0
            for material in self.materials:
                e_r = material.em.e_r
                DX  = material.domain
                int_t = assemble(dot(Et,Et)*DX)
                int_z = assemble(Ez*Ez*DX)*gamma_sqrd 
                power += 0.5*(c0/ng*eps0*e_r*(int_t + int_z))
            return power


    def compute_eigenvalues(self):
        self.esolver.solve(self.n_modes)
        if self.esolver.get_number_converged()==0:
            print('Eigensolver did not converge')
        for i in range(0,self.esolver.get_number_converged(),1):
            r, c = self.esolver.get_eigenvalue(i)
            try:
                f_b = ((-r)**0.5)*self.wavelength/(2*pi)
                print('Effective index: {} '.format(f_b))
                if self.plot_eigenmodes == True:
                    # TODO save plots to a folder
                    r, c, xr, xi = self.esolver.get_eigenpair(i)
                    e_raw = Function(self._combined_space, xr)
                    plot_transverse_field(e_raw)
            except:
                print("Found negative frequency")


if __name__ == "__main__":
    from dolfin import RectangleMesh, Point, dx, plot
    from plot import plot_vector_in_xy_plane
    width = 2.0
    height = 2.0
    n_modes = 10
    lam = 1.55
    mesh = RectangleMesh(Point(0, 0), Point(width, height), 10, 10)
    core = Material(dx)
    core.em.e_r = 1.0   
    mode_solver = EMSolver(mesh, core, lam)    
    mode_solver.eigenmode_guess = 1.1   
    mode_solver.plot_eigenmodes = True
    mode_solver.assemble_matrices()    
    mode_solver.set_electric_walls()       
    mode_solver.setup_solver()   
    mode_solver.compute_eigenvalues()
    (E, neff) = mode_solver.extract_normalized_field(0)    
    plot_vector_in_xy_plane(E)
    #plot(dot(E,E))
    
 