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
from dolfin import plot, split, dot, project, facets, as_vector
from dolfin.cpp.mesh import MeshFunctionBool
import matplotlib.pyplot as plt
import matplotlib.tri as tri



class MyDict(dict):
    def get(self, key):
        return dict.get(self, sorted(key))


def mesh2triang(mesh):    
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def plot_transverse_field(f):
    transverseE = as_vector([f[0], f[1]])
    normE = dot(transverseE, transverseE)
    plt.figure()
    plot(normE, backend = "matplotlib")
    plt.show()

def plot_vector_in_xy_plane(f):
    try:
        mesh = f.function_space().mesh()
        w0 = f.compute_vertex_values(mesh)
    except AttributeError:
        ff = project(f)
        mesh = ff.function_space().mesh()
        w0 = ff.compute_vertex_values(mesh)

    nv = mesh.num_vertices()
    w0 = [w0[i * nv: (i + 1) * nv] for i in range(3)]
    U = w0[0]
    V = w0[1]
    #W = w0[2]
    XY = mesh.coordinates()
    X = XY[:,0]
    Y = XY[:,1]
    #Z = np.zeros(nv)    
    fig = plt.figure()
    ax = fig.gca()
    ax.quiver(X,Y, U,V)
    plt.show()
    

def plot_marked_edges(f):    
    assert isinstance(f,MeshFunctionBool)
    assert (f.dim() == 1)    
    
    fm_2_vm = MyDict((facet.index(), tuple(facet.entities(0)))
             for facet in facets(f.mesh()))    
    ### plot
    plt.triplot(mesh2triang(f.mesh()), color = 'b')
    for idx,edge in enumerate(f.array()):
        if edge == True:
            (v1, v2) = fm_2_vm[idx]
            (x1, y1) = f.mesh().coordinates()[v1]
            (x2, y2) = f.mesh().coordinates()[v2]
            plt.plot([x1, x2], [y1, y2], color='r', linewidth = 2)        
    
    plt.show()
        
    
    
    