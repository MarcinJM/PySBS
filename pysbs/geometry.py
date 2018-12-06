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
from mshr import Rectangle, generate_mesh
import matplotlib.pyplot as plt

class Waveguide(SubDomain):
    """ returns True if point is inside the waveguide """
    def __init__(self, w_wg, h_wg):
        super(Waveguide, self).__init__()
        self.w_wg = w_wg
        self.h_wg = h_wg
    def inside(self, x, on_boundary):
        return (between(x[0], (-self.w_wg/2, self.w_wg/2)) and
                    between(x[1], (-self.h_wg/2, self.h_wg/2 )))

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
        # Initialize mesh function for interior domains
        self.waveguide = Waveguide(w_wg, h_wg)
        self.domains = MeshFunction("size_t", self.mesh, 2)
        self.domains.set_all(0)
        self.waveguide.mark(self.domains, 1)
        # Define new measures associated with the interior domains 
        self.dx = Measure("dx")(subdomain_data=self.domains)
