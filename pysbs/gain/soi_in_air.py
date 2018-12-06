#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:03:58 2018

@author: marcin
"""

from src.em.em_mode_solver import EMSolver
from src.el.el_mode_solver import ELSolver
from src.el.plot import plot_displacement, plot_projection
from src.em.plot import plot_transverse_field
from src.geometry import EmbeddedWaveguide
from src.material import Material, CubicStiffness,  CubicPhotoelasticity
from dolfin import Function, SubMesh, dx
import matplotlib.pyplot as plt
from math import pi


plt.close('all')
""" find the fundamental optical mode """

# lengths in um
w_sim = 1.5
h_sim = 1.5
w_wg = 0.45
h_wg = 0.22
lam = 1.55
res = 100
n_si = 3.48
n_air = 1.0

direction = 'forward'

wg = EmbeddedWaveguide(w_sim, h_sim, w_wg, h_wg, res)

# silicon
materials = []
core = Material(wg.dx(1))
core.em.e_r = n_si**2
core.el.C = CubicStiffness(165.7, 63.9,  79.6) #100 direction GPa
core.em.p = CubicPhotoelasticity(-0.09, 0.017, -0.051) #100 direction
core.el.rho = 2.328
materials.append(core)
# thermal oxide
cladding = Material(wg.dx(0))
cladding.em.er = n_air**2
cladding.em.p =  CubicPhotoelasticity(0.0, 0.0, 0.0) # set to zero for air
materials.append(cladding)

""" find the fundamental electromagnetic mode """

dlam = 0.001
em_solver = EMSolver(wg.mesh, materials, lam+dlam)
em_solver.n_modes = 2
em_solver.eigenmode_guess = core.em.e_r 
em_solver.plot_eigenmodes = False
em_solver.assemble_matrices()
em_solver.set_electric_walls()
em_solver.setup_solver()
em_solver.compute_eigenvalues()

# make sure all values correspond to the first solution
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
print(ng)

plot_transverse_field(E)

if direction == 'forward':
    """ find the fundamental elastic mode """
    q_b = 0.0
    # need to extract submesh for silicon layer only
    submesh = SubMesh(wg.mesh, wg.domains, 1)
    core_el_sim = core
    core_el_sim.domain = dx
    el_solver = ELSolver(submesh, core_el_sim)
    el_solver.n_modes = 6
    el_solver.q_b = q_b
    el_solver.plot_eigenmodes = False
    el_solver.assemble_matrices()
    el_solver.eigenmode_guess = 9
    el_solver.setup_solver()
    el_solver.compute_eigenvalues()
    (u, freq_mech) = el_solver.extract_field(0)
    plot_displacement(submesh, u)
    power_mech = el_solver.calculate_power(0)
    
elif direction == 'backward':
    
    """ find the fundamental elastic mode """
    q_b = 2.0*2.0*pi/lam*n_eff1
    # need to extract submesh for silicon layer only
    submesh = SubMesh(wg.mesh, wg.domains, 1)
    core_el_sim = core
    core_el_sim.domain = dx
    el_solver = ELSolver(submesh, core_el_sim)
    el_solver.n_modes = 10
    el_solver.q_b = q_b
    el_solver.plot_eigenmodes = False
    el_solver.assemble_matrices()
    el_solver.eigenmode_guess = 25
    el_solver.setup_solver()
    el_solver.compute_eigenvalues()
    (u, freq_mech) = el_solver.extract_field(0)
#    plot_displacement(submesh, u)
    plot_projection(submesh, u, 'Z')
    power_mech = el_solver.calculate_power(0)