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
from dolfin import *
from scipy.constants import epsilon_0
from scipy.constants import c as c0
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import ticker, cm
from mpl_toolkits.mplot3d import Axes3D
import collections
from pysbs.em.forces.electrostriction import bulk_electrostriction


def mesh2triang(mesh):
    xy = mesh.coordinates() 
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())


def plot_component(mesh, u, component):
    """ plots a single component of a vector field as a heat maps """
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


def plot_bulk_electrostriction(E, direction, q_b, materials, mesh, power_opt):
    """ plots bulk electrostriction in 2D for the whole simulation space
        the solution is first projected into Lagrange basis """
        
    from dolfin import (VectorElement, FunctionSpace, TrialFunction, Function,
                        TestFunction, split, dot, inner, lhs, rhs, solve, plot)
    import matplotlib.pyplot as plt
    
    # bulk electrostriction
    V = VectorElement("Lagrange", mesh.ufl_cell(), 1, dim = 3)
    VComplex = FunctionSpace( mesh, V*V)
    u = TrialFunction(VComplex)
    (ur, ui) = split(u)
    v = TestFunction(VComplex)
    (vr, vi) = split(v)
    f_bulk = Function(VComplex)
    
    if not isinstance(materials, collections.Iterable):
        (fr_bulk, fi_bulk) = bulk_electrostriction(E, materials.em.e_r, materials.em.p, direction, q_b)       
        F = dot(vr,fr_bulk)*materials.domain + dot(vi,fi_bulk)*materials.domain
        F -= inner(vr,ur)*materials.domain + inner(vi,ui)*materials.domain
    
    
    else:
        for idx, material in enumerate(materials):
            if idx == 0:
                (fr_bulk, fi_bulk) = bulk_electrostriction(E, material.em.e_r, material.em.p, direction, q_b)       
                F = dot(vr,fr_bulk)*material.domain + dot(vi,fi_bulk)*material.domain
                F -= inner(vr,ur)*material.domain + inner(vi,ui)*material.domain
            else:
                (fr_bulk, fi_bulk) = bulk_electrostriction(E, material.em.e_r, material.em.p, direction, q_b)
                F += dot(vr,fr_bulk)*material.domain + dot(vi,fi_bulk)*material.domain
                F -= inner(vr,ur)*material.domain + inner(vi,ui)*material.domain
                


    scaling = 1.0/(power_opt*1e3)*1e12 ##pN/(um^2mW)
    a = lhs(F)
    L = rhs(F)    
    solve(a==L, f_bulk)
    
    
    
    w0 = f_bulk.compute_vertex_values(mesh)
    nv = mesh.num_vertices()
    w0 = [w0[i * nv: (i + 1) * nv] for i in range(3)]
    U = w0[0]*scaling
    V = w0[1]*scaling
    #W = w0[2]
    XY = mesh.coordinates()
    X = XY[:,0]
    Y = XY[:,1]
    #Z = np.zeros(nv)    
    
    # make a pretty plot
    
    fig = plt.figure()
    ax = fig.gca()
    Q1 = ax.quiver(X,Y, U,V, scale=4000, scale_units='inches')
    plt.quiverkey(Q1, 0.4, 0.9, 1000.0, r'$  1000 \frac{pN}{\mu m^3 mW} Re(f_{x,y})$', labelpos='E',
                       coordinates='figure', fontproperties={'size': 24})
    plt.xlabel('x [$\mu$m]', fontsize=24, rotation = 0)
    plt.ylabel('y [$\mu$m]', fontsize=24)
    plt.tick_params(labelsize=24)
    plt.show()
    
    return fig

def plot_boundary_force(fr_bdr, boundary, power_opt):
    quiver_scale = 100
    scaling = 1.0/(power_opt*1e3)*1e12 ##pN/(um^2mW)
    fr_bdr = fr_bdr*scaling
    fig = plt.figure()
    Q1 = plt.quiver(boundary.midpoints[:,0],boundary.midpoints[:,1],
               fr_bdr[:,0],fr_bdr[:,1], scale=quiver_scale, scale_units='inches')
    plt.quiverkey(Q1, 0.5, 0.85, 50.0, r' $50 \frac{pN}{\mu m^2 mW}$', labelpos='E',
                       coordinates='figure', fontproperties={'size': 24})
    plt.axis('equal')
    plt.xlabel('x [$\mu$m]', fontsize=24)
    plt.ylabel('y [$\mu$m]', fontsize=24)
    ax = plt.gca()
    plt.tick_params(labelsize=24)    
    plt.show()
    return fig

