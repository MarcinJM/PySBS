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

from pysbs.em.em_mode_solver import EMSolver
from pysbs.el.el_mode_solver import ELSolver
from pysbs.el.plot import plot_displacement, plot_projection
from pysbs.em.plot import plot_transverse_field
from pysbs.geometry import EmbeddedWaveguide
from pysbs.material import Material, CubicStiffness,  CubicPhotoelasticity
from pysbs.tensor_rotation import Rotate_Voigt_tensor
from dolfin import Function, SubMesh, dx
import matplotlib.pyplot as plt
from math import pi


plt.close('all')
""" find the fundamental optical mode """

# lengths in um
w_sim = 1.5
h_sim = 1.5
w_wg = 0.45
h_wg = 0.23
lam = 1.55
res = 50
n_si = 3.48
n_air = 1.0

direction = 'forward'
#direction = 'backward'

""" geometery """
wg = EmbeddedWaveguide(w_sim, h_sim, w_wg, h_wg, res)


""" material parameters """

stiffness0 = CubicStiffness(164.0, 64.0,  79.0)
photoelasticity0 = CubicPhotoelasticity(-0.09, 0.017, -0.051) 
stiffness = Rotate_Voigt_tensor(stiffness0, 0)
photoelasticity = Rotate_Voigt_tensor(photoelasticity0, 0)

# silicon
materials = []
core = Material(wg.dx(2))
core.em.e_r = n_si**2
core.el.C = stiffness
core.em.p = photoelasticity
core.el.rho = 2.328
materials.append(core)
# air
cladding = Material(wg.dx(1))
cladding.em.er = n_air**2
materials.append(cladding)

""" find the fundamental electromagnetic mode """

# simulation is run twice to calculate group refractive index
dlam = 0.001
em_solver = EMSolver(wg.mesh, materials, lam+dlam)
em_solver.n_modes = 2
em_solver.eigenmode_guess = core.em.e_r 
em_solver.compute_eigenvalues()

(E2, n_eff2) = em_solver.extract_normalized_field(0)
em_solver = EMSolver(wg.mesh, materials, lam)
em_solver.n_modes = 2
em_solver.eigenmode_guess = core.em.e_r 
em_solver.compute_eigenvalues()
(E, n_eff1) = em_solver.extract_normalized_field(0)
# calculate ng
ng = n_eff1 - lam*(n_eff2-n_eff1)/dlam
print("Group refractive index is %f " % (ng))
power_opt = em_solver.calculate_power(0, ng)
plot_transverse_field(E)

    """find the fundamental elastic mode """


if direction == 'forward':

    q_b = 0.0
    # need to extract submesh for silicon layer only
    submesh = SubMesh(wg.mesh, wg.domains, 2)
    core.domain = dx # change domain to whole submesh
    el_solver = ELSolver(submesh, core)
    el_solver.n_modes = 10
    el_solver.q_b = q_b
    el_solver.eigenmode_guess =  9.3
    el_solver.compute_eigenvalues()
    (u, freq_mech) = el_solver.extract_field(0)
    plot_displacement(submesh, u)
    power_mech = el_solver.calculate_power(0)
    
else:

    q_b = 2.0*2.0*pi/lam*n_eff1
    # need to extract submesh for silicon layer only
    submesh = SubMesh(wg.mesh, wg.domains, 2)
    core.domain = dx # change domain to whole submesh
    el_solver = ELSolver(submesh, core)
    el_solver.n_modes = 10
    el_solver.q_b = q_b
    el_solver.eigenmode_guess =  24.6
    el_solver.compute_eigenvalues()
    (u, freq_mech) = el_solver.extract_field(0)
    plot_displacement(submesh, u)
    power_mech = el_solver.calculate_power(0)



""" calculate forces and plot forces | calculate gain  """

from pysbs.gain.internal_boundary import InternalBoundary
from pysbs.gain.forces import Forces
from pysbs.misc.projection import project_Efield_on_submesh
# elastic simulation is done on a subdomain hence we need to project onto submesh
Esub = project_Efield_on_submesh(em_solver, 0, submesh)
core.domain = wg.dx(2) # change to original
Q = 249
# build internal boundary with its normal
boundary = InternalBoundary(wg.domains,2,1)
forces =  Forces(E, u, q_b, Q, power_opt, power_mech, direction, freq_mech, lam, boundary, materials, Esub)
gain_bdr_electrostriction = forces.calculate_boundary_electrostriction_gain()
gain_stress = forces.calculate_boundary_stress_gain()
gain_bulk = forces.calculate_bulk_electrostriction_gain()
gain_total = forces.calculate_total_gain()
forces.plot_boundary_forces()
# plot xy electrostriction component
from pysbs.gain.plot import plot_bulk_electrostriction
core.domain = dx
plot_bulk_electrostriction(Esub, direction, q_b, core, submesh, power_opt)

print("Radiation pressure gain is %d W-1m-1" % (gain_stress))
print("Bulk electrostriction gain is %d W-1m-1" % (gain_bulk))
print("Total gain is %d W-1m-1" % (gain_total))

