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
from pysbs.el.plot import plot_displacement
from math import pi
from pysbs.geometry import EmbeddedWaveguide
from pysbs.material import IsotropicStiffness, Material
from mshr import *

"""
Note that there are plenty of modes, so if you have trouble locating the 
funamental it might be better to run the scalar solver first to find the 
good guess for the eigenvalue solver

"""
w_sim = 2.0
h_sim = 2.0
w_wg = 1.0
h_wg = 1.0    
res = 40
q_b = 2.0*2.0*pi*2.24/1.544


# waveguide
wg = EmbeddedWaveguide(w_sim, h_sim, w_wg, h_wg, res)

# plot domains to see how to assign materials
# wg.plot_domains()
materials = []
core = Material(wg.dx(1))
# As2S3
core.el.C = IsotropicStiffness(22.19, 6.20) # in GPa
core.el.rho = 3.2 # g/cm3
materials.append(core)
# thermal oxide
cladding = Material(wg.dx(0))
cladding.el.C = IsotropicStiffness(78.0, 31.6) # in GPa
cladding.rho = 2.2 # g/cm3
materials.append(cladding)
mode_solver = ELSolver(wg.mesh, materials)
mode_solver.q_b = q_b
mode_solver.n_modes = 10
mode_solver.plot_eigenmodes = False
mode_solver.assemble_matrices()
mode_solver.set_clamped_walls()
mode_solver.eigenmode_guess = 7.7
mode_solver.setup_solver()
mode_solver.compute_eigenvalues()
(u, freq_mech) = mode_solver.extract_field(0)
plot_displacement(wg.mesh, u)
