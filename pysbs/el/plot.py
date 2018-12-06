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


import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def mesh2triang(mesh):
    xy = mesh.coordinates() 
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())


def plot_projection(mesh, u, component):
    # ploting the projection of the real part of the mode
    w0 = u.compute_vertex_values(mesh)
    n = mesh.num_vertices()
    tt =0
    if component == 'X':
        tt = w0[0:n]
    elif component == 'Y':
        tt = w0[n:2*n]
    elif component == 'Z':
        tt = w0[2*n:3*n]
    else:
        print('Specify component as X, Y or Z')

    top = max(tt)*0.9
    bottom = min(tt)*0.9
    fig = plt.figure()    
    plt.tripcolor(mesh2triang(mesh), tt, cmap = cm.coolwarm)
    plt.xlabel('x [$\mu$m]', fontsize=18)
    plt.ylabel('y [$\mu$m]', fontsize=18)
    cbar = plt.colorbar(ticks=[bottom, 0, top])
    label = component + ' displacement A.U.'
    cbar.set_label(label, fontsize = 18)
    cbar.ax.set_yticklabels(['-1', '0', '1'], fontsize = 16) 
    plt.tight_layout()
    plt.axis('equal')
    plt.tick_params(labelsize=16)
    plt.show()
    return fig
   


def plot_displacement(mesh, u):
    # extract displacement values
    w0 = u.compute_vertex_values(mesh)
    n = mesh.num_vertices()
    wx = w0[0:n]
    wy = w0[n:2*n]
    wz = w0[2*n:3*n]
    rr = [0] * n
    for ii in range(0,n):
        rr[ii] = (wx[ii]**2 + wy[ii]**2 + wz[ii]**2)**0.5

    mm = max(rr)
    for ii in range(0,n):
        rr[ii] = rr[ii]/mm


    # scale wx, wy for ploting
    # convert displacement to true xy position  
    mm = max( max(wx), max(wy))
    xy = mesh.coordinates() 
    scale = max(max(xy[:,0]), max(xy[:,1]))
    for ii in range(0,n):
        wx[ii] = scale/mm*wx[ii]*0.1
        wy[ii] = scale/mm*wy[ii]*0.1
        wx[ii] = wx[ii] + xy[ii,0]
        wy[ii] = wy[ii] + xy[ii,1]

    # define new displaced mesh
    displaced_mesh = tri.Triangulation(wx, wy, mesh.cells())    
          
     
    fig = plt.figure(figsize=(8,6))   
    plt.tripcolor(displaced_mesh, rr, cmap = cm.coolwarm)    
    plt.xlabel('x [$\mu$m]', fontsize=24)
    plt.ylabel('y [$\mu$m]', fontsize=24)    
    cbar = plt.colorbar()    
    label = '|displacement|  A.U.'
    cbar.set_label(label, fontsize = 24)
    cbar.ax.tick_params(labelsize= 22)    
    cbar.set_ticks([0.1, 0.95]) # rr is normalized
    cbar.ax.set_yticklabels([ '0', '1'], fontsize = 24)    
    plt.tight_layout()
    plt.tick_params(labelsize=22)
    plt.axis('equal')    
    plt.show()
    return fig


