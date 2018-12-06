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
    You should have received a copy of the GNU General Public License
    along with pyNLO.  If not, see <http://www.gnu.org/licenses/>.
    
    If you use it in research please cite an appropriate paper from README file.
    
    If you have any questions email: marcinjmalinowki@gmail.com
    
    author: Marcin Malinowski
"""
# calculate the Maxwells stress tensor

from scipy.constants import epsilon_0
import numpy as np
import matplotlib.pyplot as plt


def MaxwellsStressTensor(E, e_r, direction):
    """ calculate Maxwells stress tensor for the domain
        assumes e_r = constant across domains thus only
        the electric component is calculated i.e. H_i
        terms are neglected in
        T_ij = eps0 epsR E_i E_j* - 0.5 d_ij E_i E_j*
               + mu0 mur H_i H_j* - 0.5 d_ij H_i H_j*
        
        Only the phase matched components are retained
        exp iWt = epx i(w_p-w_s)t
        
        The E_x,y is assumed to be real while 
        E_z complex as taken from mode solver
        projected and scaled by gamma
        check em solver for details
        """

    # for forward scattering Ep = Es, E_i = E_p, E_j = E_s
    if direction == 'forward':
        norm = 0.5*(E[0]*E[0] + E[1]*E[1] + E[2]*E[2])
        t_r = (e_r*epsilon_0*np.array(
            [[E[0]*E[0]- norm, E[0]*E[1], 0.0],
             [E[1]*E[0], E[1]*E[1] - norm, 0.0],
             [0.0, 0.0,   E[2]*E[2] - norm]]))

        t_i = 0.5*(e_r*epsilon_0*np.array(
            [[0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0]]))

    # for backward scattering Es = Ep*, E_i = E_p, E_j = E_s
    elif direction == 'backward':
        norm = 0.5*(E[0]*E[0] + E[1]*E[1] - E[2]*E[2])
        t_r = (e_r*epsilon_0*np.array(
            [[E[0]*E[0]- norm, E[0]*E[1], 0.0],
             [E[1]*E[0], E[1]*E[1] - norm, 0.0],
             [0.0, 0.0, E[2]*E[2] - norm]]))

        t_i = 0.5*(e_r*epsilon_0*np.array(
            [[0.0, 0.0, E[0]*E[2]],
             [0.0, 0.0, E[1]*E[2]],
             [E[2]*E[0], E[2]*E[1], 0.0]]))
    else:
        raise ValueError('Specify scattering direction as forward or backward')

    return (t_r, t_i)


def boundary_stress(E, oe_r, ie_r, direction, boundary , offset = 1e-12):

    fr = np.zeros((boundary.Nfacets, 3))
    fi = np.zeros((boundary.Nfacets,3))

    for idx, point in enumerate( boundary.midpoints):
        n = boundary.normal[idx]
        ipoint = point - offset*n
        opoint = point + offset*n
        iElocal = np.array([E[0](ipoint), E[1](ipoint), E[2](ipoint)])
        oElocal = np.array([E[0](opoint), E[1](opoint), E[2](opoint)])
        (itr, iti) = MaxwellsStressTensor(iElocal, ie_r, direction)
        (otr, oti) = MaxwellsStressTensor(oElocal, oe_r, direction)
        tr = otr-itr 
        ti = oti-iti
        
        f_r = np.array([tr[0,0]*n[0] + tr[0,1]*n[1],
                        tr[1,0]*n[0] + tr[1,1]*n[1],
                        tr[2,0]*n[0] + tr[2,1]*n[1]])
    
        f_i = np.array([ti[0,0]*n[0] + ti[0,1]*n[1],
                        ti[1,0]*n[0] + ti[1,1]*n[1],
                        ti[2,0]*n[0] + ti[2,1]*n[1]])
        fr[idx,:] = f_r 
        fi[idx,:] = f_i

    return (fr, fi)


