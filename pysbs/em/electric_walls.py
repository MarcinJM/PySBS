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
from dolfin import SubDomain, Constant, DirichletBC

class ElectricWalls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


def build_electric_walls(combined_space):
    c0vec = Constant((0.0, 0.0, 0.0))              
    bcs = DirichletBC(combined_space, c0vec, ElectricWalls())
    return bcs
