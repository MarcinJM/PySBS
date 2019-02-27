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

from pysbs.el.el_mode_solver import ELSolver
from pysbs.el.plot import plot_displacement, plot_projection
from pysbs.geometry import coupled_ridge_waveguides
from pysbs.material import Material, CubicStiffness,  CubicPhotoelasticity
from dolfin import Function, SubMesh, dx, as_matrix, plot, Point
from mshr import Polygon, generate_mesh
import matplotlib.pyplot as plt
from math import pi



plt.close('all')
""" find the fundamental optical mode """


direction = 'forward'
# sizes in microns
w_wg1 = 1.0
w_wg2 = 1.0
w_gap = 1.0
w_wing1 = 1.0
w_wing2 = 1.0
h_total = 0.22
h_mem = 0.12
q_b = 0.0

# numerical parameters
res = 100


# stiffness tensor in Voigt notation
def Si_stiffness_110():
    C = as_matrix(((194.5, 64.1, 35.7, 0, 0, 0),
                   (64.1, 165.7, 64.1, 0, 0, 0),
                   (35.7, 64.1, 194.5, 0, 0, 0),
                   (0  , 0  , 0  , 79.6, 0  ,   0),
                   (0  , 0  , 0  , 0  , 50.9,   0),
                   (0  , 0  , 0  , 0  ,   0, 79.6)))
    return C

domain = coupled_ridge_waveguides(w_wg1,w_wg2,w_gap,w_wing1, w_wing2, h_total,h_mem)
mesh = generate_mesh(domain,res)
#plot(mesh)



""" find the fundamental elastic mode """
core = Material(dx)
core.el.C = Si_stiffness_110()
core.el.rho = 2.328
el_solver = ELSolver(mesh, core)
el_solver.n_modes = 10
el_solver.q_b = q_b
el_solver.plot_eigenmodes = False
el_solver.assemble_matrices()
el_solver.eigenmode_guess =  4.26
el_solver.setup_solver()
el_solver.compute_eigenvalues()
(u, freq_mech) = el_solver.extract_field(0)
plot_projection(mesh, u, 'X')
plot_displacement(mesh, u, scale = 0.1)








