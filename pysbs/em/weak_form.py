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
# this is a modified mode solver from Elan Lezar Ph.D. thesis

from dolfin import Dx, dot, grad

def curl_t(w):
    return Dx(w[1], 0) - Dx(w[0], 1)

def em_weak_form(N_i, N_j, L_i, L_j, u_r, e_r, k_o_squared, DX):
   
    s_tt_ij = 1.0 / u_r * dot ( curl_t ( N_i ) , curl_t ( N_j ) )
    t_tt_ij = e_r * dot ( N_i , N_j )
    s_zz_ij = 1.0 / u_r * dot ( grad ( L_i ) , grad ( L_j ) )
    t_zz_ij = e_r * L_i * L_j
   
    s_ij = s_tt_ij*DX + s_zz_ij*DX
    t_ij = t_tt_ij*DX + t_zz_ij*DX    

    b_tt_ij = 1.0/u_r * dot ( N_i , N_j )
    b_tz_ij = 1.0/u_r * dot ( N_i , grad ( L_j ) )
    b_zt_ij = 1.0/u_r * dot ( grad ( L_i ) , N_j )
   
    a_ij = ( s_tt_ij - k_o_squared * t_tt_ij ) * DX
    b_ij = ( s_zz_ij - k_o_squared * t_zz_ij + b_tt_ij + b_tz_ij + b_zt_ij ) * DX

    return (a_ij, b_ij)
