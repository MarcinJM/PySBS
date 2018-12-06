#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 11:37:28 2018
coupling calculation 
need the E field and the U field from both solvers
# run soi in air first


@author: marcin
"""
from dolfin import *
from scipy.constants import epsilon_0
from scipy.constants import c as c0
import numpy as np
import matplotlib.tri as tri
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from src.gain.coupling import bulk_force_coupling
#from src.em.forces.electrostriction import bulk_electrostriction


def mesh2triang(mesh):
    xy = mesh.coordinates() 
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())


def plot_component(mesh, u, component):
    # ploting the projection of the real part of the mode
    try:
        w0 = u.compute_vertex_values(mesh)
    except AttributeError:
        uu = project(u)
        w0 = uu.compute_vertex_values(mesh)
    
    n = mesh.num_vertices()
    start = component*n
    stop = (component+1)*n
    tt = w0[start:stop]

    top = max(tt)
    bottom = min(tt)
    fig = plt.figure()    
    cax = plt.tripcolor(mesh2triang(mesh), tt, cmap = cm.coolwarm)
    plt.xlabel('x [$\mu$m]', fontsize=18)
    plt.ylabel('y [$\mu$m]', fontsize=18)
    plt.colorbar()
    plt.tight_layout()
    plt.tick_params(labelsize=16)
    plt.show()
    return fig


######################extract E field ##########################

power_opt = em_solver.calculate_power(0, ng)
E_sub = project_Efield_on_submesh(em_solver,0, submesh)
scaling = 1.0/(power_opt*1e3)*1e12 ##pN/(um^2mW)
(f_r, f_i) = bulk_electrostriction(E_sub, core_el_sim.em.e_r, core_el_sim.em.p, direction, q_b)
Q = 1000
omega_opt = 2.0*pi*c0/(lam*1e-6)
coupling_bulk = bulk_force_coupling(f_r, f_i, u)
gain_bulk = Q*omega_opt/(4.0*power_opt**2*power_mech)*10**21*np.absolute(coupling_bulk)**2
print(gain_bulk)


f_r = f_r*scaling
f_i = f_i*scaling

##################### pretty plot ##############################

plt.close('all')
quiver_scale = 4000
scaling = 1.0/(power_opt*1e3)*1e12
ff = project(f_r)
mesh = ff.function_space().mesh()
w0 = ff.compute_vertex_values(mesh)

nv = mesh.num_vertices()
w0 = [w0[i * nv: (i + 1) * nv] for i in range(3)]
U = w0[0]
V = w0[1]
XY = mesh.coordinates()
X = XY[:,0]
Y = XY[:,1]




if direction == 'backward':
    uu = project(f_i[2])
    w0 = uu.compute_vertex_values(mesh)    
    
    #Z = np.zeros(nv)    
    fig = plt.figure()
    ax = fig.gca()
    cax = plt.tripcolor(mesh2triang(mesh), w0,cmap = cm.Blues_r)
    cbar = plt.colorbar(orientation='horizontal',  pad=0.2)
    tick_locator = ticker.MaxNLocator(nbins=3)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label(r'$Im(f_z) [\frac{pN}{\mu m^3 mW}]$',x = 0.8, y = 1.08, fontsize=24,  color ='#08306B')
    cbar.ax.tick_params(labelsize=24)
    Q1 = ax.quiver(X,Y, U,V, scale=quiver_scale, scale_units='inches')
    plt.quiverkey(Q1, 0.05, 0.1, 1000.0, r'$  1000 \frac{pN}{\mu m^3 mW} Re(f_{x,y})$', labelpos='E',
                       coordinates='figure', fontproperties={'size': 24})
    ax_scale = 1.2
    plt.axis('equal')
    plt.xlim((-w_wg/2.0*ax_scale, w_wg/2.0*ax_scale))
    plt.ylim((-h_wg/2.0*ax_scale, h_wg/2.0*ax_scale))
    plt.xlabel('x [$\mu$m]', fontsize=24, rotation = 0)
    plt.ylabel('y [$\mu$m]', fontsize=24)
    #ax = plt.gca()
    ax.yaxis.set_label_coords(0.08, 0.55)
    ax.xaxis.set_label_coords(0.65, -0.10)
    plt.tight_layout(pad=2.0)
    plt.tick_params(labelsize=24)
    plt.show()
    
elif direction == 'forward':
    
    #Z = np.zeros(nv)    
    fig = plt.figure()
    ax = fig.gca()

    Q1 = ax.quiver(X,Y, U,V, scale=quiver_scale, scale_units='inches')
    plt.quiverkey(Q1, 0.35, 0.85, 1000.0, r'$  1000 \frac{pN}{\mu m^3 mW} Re(f_{x,y})$', labelpos='E',
                       coordinates='figure', fontproperties={'size': 24 })
    ax_scale = 1.5
    plt.axis('equal')
    plt.xlim((-w_wg/2.0*ax_scale, w_wg/2.0*ax_scale))
    plt.ylim((-h_wg/2.0*ax_scale, h_wg/2.0*ax_scale))
    plt.xlabel('x [$\mu$m]', fontsize=24, rotation = 0)
    plt.ylabel('y [$\mu$m]', fontsize=24)
    #ax = plt.gca()
    ax.yaxis.set_label_coords(0.08, 0.55)
    ax.xaxis.set_label_coords(0.65, -0.10)
    plt.tight_layout(pad=2.0)
    plt.tick_params(labelsize=24)
    plt.show()

