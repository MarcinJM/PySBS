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
from dolfin import Function, SubMesh, dx
import matplotlib.pyplot as plt
from math import pi
from pysbs.material import IsotropicStiffness, Material, IsotropicPhotoelasticity


"""
Note that there are plenty of modes, so if you have trouble locating the 
funamental it might be better to run the scalar solver first to find the 
good guess for the eigenvalue solver

"""
w_sim = 2.0
h_sim = 2.0
w_wg = 1.0
h_wg = 0.5
res = 40
lam = 1.544
q_b = 2.0*2.0*pi*2.24/lam
n_chg = 2.37
n_ox = 1.44
direction = 'backward'


# waveguide
wg = EmbeddedWaveguide(w_sim, h_sim, w_wg, h_wg, res)

# plot domains to see how to assign materials
# wg.plot_domains()
materials = []
core = Material(wg.dx(1))
# As2S3
core.el.C = IsotropicStiffness(22.19, 6.20) # in GPa
core.em.p = IsotropicPhotoelasticity(0.25, 0.24)
core.em.e_r = n_chg**2
core.el.rho = 3.2 # g/cm3
materials.append(core)
# thermal oxide
cladding = Material(wg.dx(0))
cladding.el.C = IsotropicStiffness(78.0, 31.6) # in GPa
cladding.em.p = IsotropicPhotoelasticity(0.12, 0.27)
cladding.em.e_r = n_ox**2
cladding.rho = 2.2 # g/cm3
materials.append(cladding)




""" find the fundamental electromagnetic mode """

# simulation is run twice to calculate group refractive index
dlam = 0.001
em_solver = EMSolver(wg.mesh, materials, lam+dlam)
em_solver.n_modes = 2
em_solver.eigenmode_guess = core.em.e_r 
em_solver.plot_eigenmodes = False
em_solver.assemble_matrices()
em_solver.set_electric_walls()
em_solver.setup_solver()
em_solver.compute_eigenvalues()

(E2, n_eff2) = em_solver.extract_normalized_field(0)
em_solver = EMSolver(wg.mesh, materials, lam)
em_solver.n_modes = 2
em_solver.eigenmode_guess = core.em.e_r 
em_solver.plot_eigenmodes = False
em_solver.assemble_matrices()
em_solver.set_electric_walls()
em_solver.setup_solver()
em_solver.compute_eigenvalues()
(E, n_eff1) = em_solver.extract_normalized_field(0)
# calculate ng
ng = n_eff1 - lam*(n_eff2-n_eff1)/dlam
print("Group refractive index is %f " % (ng))
power_opt = em_solver.calculate_power(0, ng)
plot_transverse_field(E)



""" find the fundamental elastic mode """

el_solver = ELSolver(wg.mesh, materials)
el_solver.n_modes = 10
el_solver.q_b = q_b
el_solver.plot_eigenmodes = False
el_solver.assemble_matrices()
el_solver.eigenmode_guess = 7.79
el_solver.setup_solver()
el_solver.compute_eigenvalues()
(u, freq_mech) = el_solver.extract_field(0)
power_mech = el_solver.calculate_power(0)
plot_projection(wg.mesh, u, 'Z')


""" calculate forces and plot forces | calculate gain  """

from pysbs.gain.internal_boundary import InternalBoundary
from pysbs.gain.forces import Forces

# Q from  "On-chip stimulated Brillouin scattering"
Q = 7700.0/34
# build internal boundary with its normal
boundary = InternalBoundary(wg.domains, 1, 0)
forces =  Forces(E, u, q_b, Q, power_opt, power_mech, direction, freq_mech, lam, boundary, materials)
forces.calculate_boundary_electrostriction()
gain_bdr_electrostriction = forces.calculate_boundary_electrostriction_gain()
forces.calculate_boundary_stress()
gain_stress = forces.calculate_boundary_stress_gain()
forces.calculate_bulk_electrostriction()
gain_bulk = forces.calculate_bulk_electrostriction_gain()
gain_total = forces.calculate_total_gain()

print("Radiation pressure gain is %d W-1m-1" % (gain_stress))
print("Bulk electrostriction gain is %d W-1m-1" % (gain_bulk))
print("Total gain is %d W-1m-1" % (gain_total))


""" analytical formula for fibers """
from scipy.constants import c as c0
p12 = 0.24
Aeff = w_wg*1e-6*h_wg*1e-6
df = freq_mech*1e9/Q
rho = 3200
overlap = 0.95
gain_estimate = overlap*4*pi*n_chg**8*p12**2/(rho*c0*lam**3*1e-18*df*freq_mech*1e9*Aeff)
print("Bulk electrostriction gain estimate is %d W-1m-1" % (gain_estimate))

#





