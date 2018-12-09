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
from dolfin import as_vector, Dx, Constant, triangle, as_matrix
from scipy.constants import epsilon_0
from pysbs.gain.internal_boundary import InternalBoundary
import collections
import numpy as np


def sigma_boundary(E,p, e_r, direction):
    """ returns the stress tensor as a numpy matrix
        to be used in boundary electrostriction calculation """
        
    assert isinstance(p,np.ndarray)
        
            # for forward scattering Ep = Es
    if direction == 'forward':
        EEr = np.array([E[0]*E[0], E[1]*E[1], E[2]*E[2],
                        0.0, 0.0, 2.0*E[0]*E[1]])
        #EEi = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sigma_r =  -0.5*epsilon_0*e_r**2*p.dot(EEr)
        sigma_i = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # for backward scattering Es = Ep *
    elif direction == 'backward':
        EEr = np.array([E[0]*E[0], E[1]*E[1],-E[2]*E[2],
                        0.0, 0.0, 2.0*E[0]*E[1]])
        EEi = np.array([0.0, 0.0, 0.0, 2.0*E[1]*E[2], 2.0*E[0]*E[2], 0.0])
        sigma_r =  -0.5*epsilon_0*e_r**2*p.dot(EEr)
        sigma_i =  -0.5*epsilon_0*e_r**2*p.dot(EEi)
    else:
        raise ValueError('Specify scattering direction as forward or backward')
    
    return (sigma_r, sigma_i)
        

def facet_electrostriction_force(E, p, e_r, direction, normal):
    """ calcuate the boundary electrostricion for one material
        return a tuple of real and imaginary force 
         The E_x,y is assumed to be real while 
        E_z complex as taken from mode solver
        projected and scaled by gamma """
    
    (sigma_r, sigma_i) = sigma_boundary(E, p, e_r, direction)
    fr = np.array([sigma_r[0]*normal[0]+sigma_r[5]*normal[1],
                   sigma_r[5]*normal[0]+sigma_r[1]*normal[1],
                   sigma_r[4]*normal[0]+sigma_r[3]*normal[1]])
    fi = np.array([sigma_i[0]*normal[0]+sigma_i[5]*normal[1],
                   sigma_i[5]*normal[0]+sigma_i[1]*normal[1],
                   sigma_i[4]*normal[0]+sigma_i[3]*normal[1]])

    return (fr, fi)



def boundary_electrostriction(E, po, pi, oe_r, ie_r, direction, boundary , offset = 1e-12):
    """ calculate the boundary electrostriction at the interface of two materials """
    fr = np.zeros((boundary.Nfacets, 3))
    fi = np.zeros((boundary.Nfacets, 3))
 
    for idx, point in enumerate( boundary.midpoints):
        n = boundary.normal[idx]
        ipoint = point- offset*n
        opoint = point + offset*n
        iElocal = np.array([E[0](ipoint), E[1](ipoint), E[2](ipoint)])
        oElocal = np.array([E[0](opoint), E[1](opoint), E[2](opoint)])
        (ifr, ifi) = facet_electrostriction_force(iElocal, pi, ie_r, direction, n)
        (ofr, ofi) = facet_electrostriction_force(oElocal, po, oe_r, direction, n)
        fr[idx,:] = (ifr-ofr)
        fi[idx,:] = (ifi-ofi)
    return (fr, fi)



def f_elst_r(sigma_r, sigma_i, q_b):
    # assumes sigma is a 6 vector and sigma exp(i*q_b*z)
    f = as_vector([-Dx(sigma_r[0], 0) - Dx(sigma_r[5], 1) + q_b*sigma_i[4],
                   -Dx(sigma_r[5], 0) - Dx(sigma_r[1], 1) + q_b*sigma_i[3],
                   -Dx(sigma_r[4], 0) - Dx(sigma_r[3], 1) + q_b*sigma_i[2]])
    return f



def f_elst_i(sigma_r, sigma_i, q_b):
    # assumes sigma is a 6 vector and sigma exp(i*q_b*z)
    f = as_vector([-Dx(sigma_i[0], 0) - Dx(sigma_i[5], 1) - q_b*sigma_r[4],
                   -Dx(sigma_i[5], 0) - Dx(sigma_i[1], 1) - q_b*sigma_r[3],
                   -Dx(sigma_i[4], 0) - Dx(sigma_i[3], 1) - q_b*sigma_r[2]])
    return f


def bulk_electrostriction(E, e_r,p, direction, q_b):
        """ calculate the bulk electrostriction force """
        pp = as_matrix(p) # photoelestic tensor is stored as as numpy array
        if direction == 'backward':
            EE_6vec_r = as_vector([E[0]*E[0], E[1]*E[1], -E[2]*E[2],
                        0.0, 0.0, 2.0*E[0]*E[1]])
            EE_6vec_i = as_vector([0.0, 0.0, 0.0, 2.0*E[1]*E[2],
                                   2.0*E[0]*E[2], 0.0])
            sigma_r = -0.5*epsilon_0*e_r**2*pp*EE_6vec_r
            sigma_i = -0.5*epsilon_0*e_r**2*pp*EE_6vec_i
            f_r = f_elst_r(sigma_r, sigma_i, q_b)
            f_i = f_elst_i(sigma_r, sigma_i, q_b)
            
        elif direction == 'forward':
            EE_6vec_r = as_vector([E[0]*E[0], E[1]*E[1], E[2]*E[2],
                        0.0, 0.0, 2.0*E[0]*E[1]])
            # EE_6vec_i, sigma_i, q_b is zero
            sigma_r = -0.5*epsilon_0*e_r**2*pp*EE_6vec_r
            # no need to multiply zeros ...
            f_r = as_vector([-Dx(sigma_r[0], 0) - Dx(sigma_r[5], 1),
                             -Dx(sigma_r[5], 0) - Dx(sigma_r[1], 1),
                             -Dx(sigma_r[4], 0) - Dx(sigma_r[3], 1)])
            f_i = Constant((0.0, 0.0, 0.0), cell = triangle)
            #f_i = as_vector([ - q_b*sigma_r[4],- q_b*sigma_r[3],- q_b*sigma_r[2]])

        else:
            raise ValueError('Specify scattering direction as forward or backward')
            
        return (f_r, f_i)




if __name__ == "__main__":

    """ Test electrostriction calculation. For uniform field the force
        is the same as in fibers f = i0.5 q p12 n^4 E^2 """
    from dolfin import *
    from scipy.constants import epsilon_0
    from pysbs.material import IsotropicPhotoelasticity
    

    
    direction = 'backward'
    e_r = 2.2
    q_b = 1.1
    mesh = UnitSquareMesh(32, 32)
    # define some E - field
    V = VectorElement("Lagrange", mesh.ufl_cell(), 2, dim = 3)
    em_space  = FunctionSpace(mesh, V)    
    E_analytical = Expression(('0.0', '1.0', '0.0'),
                               element = em_space.ufl_element())
    E = interpolate(E_analytical, em_space)
    
    p = IsotropicPhotoelasticity(3.0, 2.0)
    (f_r,f_i) = bulk_electrostriction(E, e_r,p, direction, q_b)
    
    m1 = f_i[2]*dx
    v1 = assemble(m1)
    # analytical force per unit length
    f = 0.5*epsilon_0*e_r**2*q_b*2.0*1**2 # f = i0.5 q p12 n^4 E^2
    error = abs((f-v1)/f)
    print(error)
    
    
