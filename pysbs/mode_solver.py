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
from dolfin import PETScMatrix, SLEPcEigenSolver


class ModeSolver(object):
    """ this is the base class for the two mode solvers """

    def __init__(self, mesh, materials):
        """ init class always require mesh and materials to
        build the weak form """
        self.mesh = mesh
        self.materials = materials #FIXME throw error if not all domains are initilized
        self.eigenmode_guess = 1.0
        self.n_modes = 10 # number of modes to solve for      
        self.plot_eigenmodes = False
        self._A = PETScMatrix()
        self._B = PETScMatrix()	
        self.esolver = SLEPcEigenSolver(self._A, self._B)
   

    def build_weak_form(self):
        pass

    def setup_solver(self):
        pass
    
    def set_boundary_conditions(self):
        pass
    
    def compute_eigenvalues(self):
        pass

    
