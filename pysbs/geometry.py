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
from dolfin import SubDomain, between, plot, MeshFunction, Point, Measure 
from mshr import Rectangle, generate_mesh, Polygon
import matplotlib.pyplot as plt
from dolfin import Point


def coupled_ridge_waveguides(w_wg1,w_wg2,w_gap,w_wing1, w_wing2,h_total,h_mem):
    # manually define all points
    domain_vertices = []
    x = 0.0
    y= 0.0
    domain_vertices.append(Point(x,y))
    x = 0.0
    y= h_mem
    domain_vertices.append(Point(x,y))
    x += w_wing1
    y= h_mem
    domain_vertices.append(Point(x,y))
    x = x
    y= h_total
    domain_vertices.append(Point(x,y))
    x += w_wg1
    y= h_total
    domain_vertices.append(Point(x,y))
    x = x
    y= h_mem
    domain_vertices.append(Point(x,y))
    x += w_gap
    y= h_mem
    domain_vertices.append(Point(x,y))
    x = x
    y= h_total
    domain_vertices.append(Point(x,y))
    x += w_wg2
    y = h_total
    domain_vertices.append(Point(x,y))
    x = x
    y = h_mem
    domain_vertices.append(Point(x,y))
    x += w_wing2
    y = h_mem
    domain_vertices.append(Point(x,y))
    x = x
    y = 0.0
    domain_vertices.append(Point(x,y))
    domain_vertices.reverse() # need counter clockwise order
    
    
    domain = Polygon(domain_vertices)
    return domain

        
        
class RidgeWaveguide():
    """ a ridge waveguide with cladding class 
        builds the mesh and marks domains """
    
    def __init__(self, w_membrane, w_ridge, h_ridge, h_memebrane, w_sim, h_sim, res):
        self.res = res
        domain_vertices = [Point(w_membrane/2, 0.0),
                           Point(w_membrane/2, h_membrane),
                           Point(w_ridge/2.0, h_membrane),
                           Point(w_ridge/2.0 , h_membrane + h_ridge),
                           Point(-w_ridge/2.0, h_membrane + h_ridge),
                           Point(-w_ridge/2.0, h_membrane),
                           Point(-w_membrane/2.0, h_membrane),
                           Point(-w_membrane/2.0, 0.0)]
        waveguide = Polygon(domain_vertices)
        domain = Rectangle(Point(-w_sim/2,-h_sim/2),
                           Point(w_sim/2,h_sim/2))
        domain.set_subdomain(1, Rectangle(Point(-w_sim/2,-h_sim/2),
                                          Point(w_sim/2,h_sim/2)))
        domain.set_subdomain(2, waveguide)
        self.mesh = generate_mesh(domain, res)
        self.mesh.init()
        self.domains = MeshFunction('size_t', self.mesh, 2, self.mesh.domains())
        self.dx = Measure("dx")(subdomain_data=self.domains)
        
        



class EmbeddedWaveguide():
    """ a rectangular waveguide with cladding class 
        builds the mesh and marks domains """

    def __init__(self, w_sim, h_sim, w_wg, h_wg, res):
        # Create mesh with two domains, waveguide + cladding
        self.res = res
        domain = Rectangle(Point(-w_sim/2,-h_sim/2),
                   Point(w_sim/2,h_sim/2))
        domain.set_subdomain(1, Rectangle(Point(-w_sim/2,-h_sim/2),
                                  Point(w_sim/2,h_sim/2)))
        domain.set_subdomain(2, Rectangle(Point(-w_wg/2,-h_wg/2),
                                  Point(w_wg/2,h_wg/2)))
        self.mesh = generate_mesh(domain, self.res)
        self.mesh.init()
        self.domains = MeshFunction('size_t', self.mesh, 2, self.mesh.domains())
        self.dx = Measure("dx")(subdomain_data=self.domains)
