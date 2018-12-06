#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:31:44 2018
coupling calculation 
need the E field and the U field from both solvers
# run soi in air first
@author: marcin

#"""
##
#from src.gain.coupling import displacement_at_boundary, boundary_force_coupling, project_Efield_on_submesh
#from src.gain.internal_boundary import InternalBoundary
#from src.em.forces.electrostriction import (bulk_electrostriction,
#                    boundary_electrostriction, plot_bulk_electrostriction)
#from src.em.forces.stress import boundary_stress
import numpy as np
from scipy.constants import c as c0
from dolfin import assemble, dot
from math import pi

# project field onto the domain of computation
power_opt = em_solver.calculate_power(0, ng)
E_sub = project_Efield_on_submesh(em_solver,0, submesh)

# build internal boundary with its normal
boundary = InternalBoundary(wg.domains, 1,0)
omega_mech = freq_mech*2.0*pi*1e9
omega_opt = 2.0*pi*c0/(lam*1e-6)
Q = 1000

#E = project(E)

boundary = InternalBoundary(wg.domains, 1,0)
omega_opt = 2.0*pi*c0/(lam*1e-6)    
scaling = 1.0/(power_opt*1e3)*1e12 ##pN/(um^2mW)
Q = 1000
(ur_bdr, ui_bdr) = displacement_at_boundary(boundary, u)

# boundary electrostriction
(fr_bdr_elctst, fi_bdr_elctst) = boundary_electrostriction(E, cladding.em.p, core.em.p,
                             cladding.em.e_r, core.em.e_r, direction, boundary, offset = 1e-12)


gain_bdr_elcst =  (Q*omega_opt/(4.0*power_opt**2*power_mech)*10**21*
            np.absolute(boundary_force_coupling(fr_bdr_elctst, fi_bdr_elctst, ur_bdr, ui_bdr, boundary))**2)
print(gain_bdr_elcst)
# radiation pressure
(fr_bdr_stress, fi_bdr_stress) = boundary_stress(E, cladding.em.e_r, core.em.e_r,
                                        direction, boundary, offset = 1e-12)

gain_bdr_stress =  (Q*omega_opt/(4.0*power_opt**2*power_mech)*10**21*
                    np.absolute(boundary_force_coupling(fr_bdr_stress, fi_bdr_stress, ur_bdr, ui_bdr, boundary))**2)
print(gain_bdr_stress)



########################## plots #############################################
green = '#0e6600'
red = '#c81515'
blue = '#0e2fb4'

quiver_scale = 100
scaling = 1.0/(power_opt*1e3)*1e12 ##pN/(um^2mW)
fr_bdr_elctst = fr_bdr_elctst*scaling
fr_bdr_stress = fr_bdr_stress*scaling
plt.figure(0)
Q1 = plt.quiver(boundary.midpoints[:,0],boundary.midpoints[:,1],
           fr_bdr_elctst[:,0],fr_bdr_elctst[:,1], color=red, scale=quiver_scale, scale_units='inches')
plt.quiverkey(Q1, 0.61, 0.85, 50.0, r' $50 \frac{pN}{\mu m^2 mW}, ES$', labelpos='E',
                   coordinates='figure', fontproperties={'size': 24})
arrow1 = max(fr_bdr_elctst[:,0])
Q2 = plt.quiver(boundary.midpoints[:,0],boundary.midpoints[:,1],
           fr_bdr_stress[:,0],fr_bdr_stress[:,1], color=blue,  scale=quiver_scale, scale_units='inches' )
Q2 = plt.quiver(boundary.midpoints[:,0],boundary.midpoints[:,1],
           fr_bdr_stress[:,0],fr_bdr_stress[:,1], color=blue,  scale=quiver_scale, scale_units='inches' )
arrow2 = max(fr_bdr_stress[:,0])
plt.quiverkey(Q2, 0.19, 0.85, 50.0, r'$50 \frac{pN}{\mu m^2 mW}, RP$', labelpos='E',
                   coordinates='figure', fontproperties={'size':24})
ax_scale = 2.5
plt.axis('equal')
plt.xlim((-w_wg/2.0*ax_scale, w_wg/2.0*ax_scale))
plt.ylim((-h_wg/2.0*ax_scale, h_wg/2.0*ax_scale*1.3))
plt.xlabel('x [$\mu$m]', fontsize=24)
plt.ylabel('y [$\mu$m]', fontsize=24)
ax = plt.gca()
ax.yaxis.set_label_coords(0.08, 0.65)
ax.xaxis.set_label_coords(0.85, 0.10)
plt.tight_layout(pad=1.5)
plt.tick_params(labelsize=24)

plt.show()

