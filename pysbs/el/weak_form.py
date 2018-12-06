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
from ufl import indices


def eps_T_Vgt(u):    
    # included factor of 2 for e_ij where i!=j
    epsilon = as_vector([Dx(u[0],0),
                        Dx(u[1],1),
                        0.0,
                        Dx(u[2],1),
                        Dx(u[2],0),
                        Dx(u[0],1) + Dx(u[1],0)])               
    return epsilon


def eps_Z_Vgt(u):
    epsilon = as_vector([0.0,0.0,u[2],u[1],u[0],0.0])               
    return epsilon


def LHS(q_b, C, ur, ui, vr, vi, DX):
# define the weak form with norm u v*
    i, j = indices(2) # define to avoid name collisions
        
    TT = C[i,j]*eps_T_Vgt(ur)[j]*eps_T_Vgt(vr)[i]*DX + \
         C[i,j]*eps_T_Vgt(ui)[j]*eps_T_Vgt(vi)[i]*DX + \
         C[i,j]*eps_T_Vgt(ui)[j]*eps_T_Vgt(vr)[i]*DX - \
         C[i,j]*eps_T_Vgt(ur)[j]*eps_T_Vgt(vi)[i]*DX

    TZ = q_b*(C[i,j]*eps_T_Vgt(ur)[j]*eps_Z_Vgt(vr)[i]*DX + \
          C[i,j]*eps_T_Vgt(ui)[j]*eps_Z_Vgt(vi)[i]*DX - \
          C[i,j]*eps_T_Vgt(ui)[j]*eps_Z_Vgt(vr)[i]*DX + \
          C[i,j]*eps_T_Vgt(ur)[j]*eps_Z_Vgt(vi)[i]*DX)

    ZT = q_b*(C[i,j]*eps_Z_Vgt(ur)[j]*eps_T_Vgt(vr)[i]*DX + \
          C[i,j]*eps_Z_Vgt(ui)[j]*eps_T_Vgt(vi)[i]*DX - \
          C[i,j]*eps_Z_Vgt(ui)[j]*eps_T_Vgt(vr)[i]*DX + \
          C[i,j]*eps_Z_Vgt(ur)[j]*eps_T_Vgt(vi)[i]*DX)

    ZZ = q_b*q_b*(C[i,j]*eps_Z_Vgt(ur)[j]*eps_Z_Vgt(vr)[i]*DX + \
              C[i,j]*eps_Z_Vgt(ui)[j]*eps_Z_Vgt(vi)[i]*DX + \
              C[i,j]*eps_Z_Vgt(ui)[j]*eps_Z_Vgt(vr)[i]*DX - \
              C[i,j]*eps_Z_Vgt(ur)[j]*eps_Z_Vgt(vi)[i]*DX)
    a = - TT + TZ + ZT - ZZ
    return a


def RHS(rho, ur, ui, vr, vi, DX):
    b = -rho*(inner(ur,vr)*DX + inner(ui,vi)*DX + \
              inner(ui,vr)*DX - inner(ur,vi)*DX)
    return b
