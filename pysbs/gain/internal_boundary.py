#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 15:10:57 2018

@author: marcin

"""

from dolfin import MeshFunction, facets, Vertex
import numpy as np


def facet_length(facet):
    """ calculate FEniCS facet length """
    v_idx = facet.entities(0)
    v0 = Vertex(facet.mesh(), v_idx[0])
    v1 = Vertex(facet.mesh(), v_idx[1])
    x0 = v0.point().x()
    x1 = v1.point().x()
    y0 = v0.point().y()
    y1 = v1.point().y()
    return ((x1-x0)**2 + (y1-y0)**2)**0.5
    


class InternalBoundary():
    """ class for calculating the normal of an internal boundary
        The normal is stored as numpy array because operating
        withing the fenics framework is cumbersome i.e. direction of internal
        normals is not defined"""
        
    def __init__(self, domains, interior, exterior):
        self.interior = interior
        self.exterior = exterior
        self.facet_lengths = None
        (self.edges, self.Nfacets) = self.create_boundary(domains, interior, exterior)
        (self.normal, self.midpoints) = self.calculate_normal(domains, interior, exterior)
        
    def create_boundary(self, domains, interior, exterior):
        """ mark commmon facets between two domains with 1 and a return FacetFunction"""
    
        edges = MeshFunction("size_t", domains.mesh(), domains.mesh().topology().dim() - 1)
        edges.set_all(0)
        
        # count number of boundary facets
        Nfacets = 0
        for facet in facets(domains.mesh()):
            cells = facet.entities(2)
            if facet.exterior() == False:
                
                cell1 = cells[0]
                cell2 = cells[1]
                domain1 = domains.array()[cells[0]]
                domain2 = domains.array()[cells[1]]
                if (sorted([domain1, domain2]) == 
                    sorted([interior, exterior])):
                    idx = facet.index()
                    edges.array()[idx] = True
                    Nfacets += 1                    

        return (edges, Nfacets)

    def calculate_normal(self, domains, interior, exterior):
        """calculate normal pointing towards the exterior domain at the 
        middle of each facet. Return the normal and its position """
       
        self.normal = np.zeros((self.Nfacets, 2))
        self.midpoints = np.zeros((self.Nfacets, 2))
        self.facet_lengths = np.zeros(self.Nfacets)
        counter = 0

        for facet in facets(domains.mesh()):
            cells = facet.entities(2)
            if facet.exterior() == False:
                domain1 = domains.array()[cells[0]]
                domain2 = domains.array()[cells[1]]
                if (sorted([domain1, domain2]) == 
                    sorted([interior, exterior])):
                    self.facet_lengths[counter] = facet_length(facet)
                    if domain1 ==interior:
                        self.normal[counter,0] =  facet.normal().x()
                        self.normal[counter,1] =  facet.normal().y()
                        self.midpoints[counter,0] =  facet.midpoint().x()
                        self.midpoints[counter,1] =  facet.midpoint().y()

                    elif domain1 == exterior:
                        self.normal[counter,0] =  -facet.normal().x()
                        self.normal[counter,1] =  -facet.normal().y()
                        self.midpoints[counter,0] =  facet.midpoint().x()
                        self.midpoints[counter,1] =  facet.midpoint().y()
                        
                    counter += 1
                    
        return (self.normal, self.midpoints)

                    
            
if __name__ == "__main__":


    from dolfin import *
    from scipy.constants import epsilon_0
    import matplotlib.pyplot as plt
    from mshr import *
    import numpy as np
  
    w_sim = 2.0
    h_sim = 2.0
    w_wg = 1.0
    h_wg = 1.0
    res = 20



    class Waveguide(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[0], (-w_wg/2,w_wg/2)) and between(x[1], (-h_wg/2, h_wg/2 )))
            
        
    domain = Rectangle(Point(-w_sim/2,-h_sim/2),
               Point(w_sim/2,h_sim/2))
    domain.set_subdomain(1, Rectangle(Point(-w_sim/2,-h_sim/2),
                              Point(w_sim/2,h_sim/2)))
    domain.set_subdomain(2, Rectangle(Point(-w_wg/2,-h_wg/2),
                              Point(w_wg/2,h_wg/2)))
    mesh = generate_mesh(domain, res)
    mesh.init()

    # Initialize mesh function for interior domains
    waveguide = Waveguide()
    domains = MeshFunction("size_t", mesh, 2)
    domains.set_all(0)
    waveguide.mark(domains, 1)
    
    boundary = InternalBoundary(domains, 1,0)
    normal = boundary.normal
    midpoints = boundary.midpoints
    
            
    plt.quiver(midpoints[:,0],midpoints[:,1],normal[:,0],normal[:,1], color='r')
    plt.show()

        
    
        
        