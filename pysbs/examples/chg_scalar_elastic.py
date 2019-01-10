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

from __future__ import print_function
from dolfin import *
from mshr import *
from pylab import show,triplot
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from math import pi
import numpy as np
from pysbs.geometry import EmbeddedWaveguide


# physical parameters
rho_chg = 3.2 #g/cm3
rho_ox = 2.2
G_chg = Constant(6.2) # shear modulus GPa
M_chg = Constant(22.19) # P-wave modulus GPa
G_ox = Constant(31.6)
M_ox = Constant(78.0)
q_b = 2.0*2.0*pi*2.24/1.544  #sizes [um]
q_b_sqrd = Constant(q_b*q_b)

# sizes [um]
w_sim = 2.0
h_sim = 2.0
w_wg = 1.0
h_wg = 0.5




# numerical parameters
n_modes = 3
res = 50


# waveguide
wg = EmbeddedWaveguide(w_sim, h_sim, w_wg, h_wg, res)
# define function space
V = FunctionSpace(wg.mesh, "Lagrange", 1)
dx = wg.dx


# define basis and bilinear form
u = TrialFunction(V)
v = TestFunction(V)
a = G_ox*inner(grad(u), grad(v))*dx(0) + q_b_sqrd*M_ox*inner(u,v)*dx(0) \
+ G_chg*inner(grad(u), grad(v))*dx(1) + q_b_sqrd*M_chg*inner(u,v)*dx(1)

b = rho_ox*inner(u,v)*dx(0) + rho_chg*inner(u,v)*dx(1)

# Define Dirichlet boundary conditions

def boundary(x, on_boundary):
    return on_boundary
              
bcs = DirichletBC(V, 0.0,boundary)


# Assemble stiffness form
A = PETScMatrix()
B = PETScMatrix()
assemble(a, tensor=A)
assemble(b, tensor=B)
bcs.apply(A)
bcs.apply(B)



# Create eigensolver
eigensolver = SLEPcEigenSolver(A,B)
eigensolver.parameters["solver"] = "krylov-schur"
eigensolver.parameters["tolerance"] = 1e-12
eigensolver.parameters["spectral_shift"] = (7.0*2*pi)**2
eigensolver.parameters["spectrum"] = "target magnitude"
eigensolver.parameters["spectral_transform"] = "shift-and-invert"

# Compute all eigenvalues of A x = \lambda x
print("Computing eigenvalues. This can take a minute.")
eigensolver.solve(n_modes)

# Extract largest (first) eigenpair
for i in range(0,n_modes):
    r, c, rx, cx = eigensolver.get_eigenpair(i)
    if r > 0.0:
        f_b = (r)**0.5/(2*pi)
        print("Frequency ", f_b, "GHz")
    else:
        print("Negative Frequency")


r, c, rx, cx = eigensolver.get_eigenpair(0)
# Initialize function and assign eigenvector
u = Function(V)
u.vector()[:] = rx
plot(u)



