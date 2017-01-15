#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import math
import unittest
import numpy as np_M

import dolfin as df_M

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

class TestChargeDensity(unittest.TestCase):
    """Test the particle-to-field functions in class Field_C"""
    
    def setUp(self):
        # initializations for each test go here...

        return

    def test_1_interpolate_particleCharge_to_1Dmesh(self):
        """Compute the mesh charge-density due to a planar charge.
        """

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'
        
        # Create a 1D mesh
        mesh = df_M.IntervalMesh(2, -0.5, 0.5)

        mesh1DCI = Mesh_C(Mesh=mesh, computeDictionaries=True, computeTree=True, plotFlag=False)

# Should I make a Mesh_C?        

        # Put 3 DT particles on it

        x0 = -0.25; y0 = 0.0; z0 = 0.0
        ux0 = 0.0; uy0 = 0.0; uz0 = 0.0
        weight0 = 2.0e10 # number of electrons per macroparticle
        bitflags0 = 0b0
        cell_index0 = 0

        p0 = (x0,y0,z0, ux0,uy0,uz0, weight0,bitflags0,cell_index0)


        # 2nd particle
        x1 = 0.25; y1 = 0.0; z1 = 1.0
        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 3.0e10
        bitflags1 = 0b0
        cell_index1 = 1

        p1 = (x1,y1,z1, ux1,uy1,uz1, weight1,bitflags1,cell_index1)

        # 3nd particle
        x2 = 0.0; y2 = 0.0; z2 = 1.0
        ux2 = uy2 = 0.0; uz2 = -uz0
        weight2 = 4.0e10
        bitflags2 = 0b0
        cell_index2 = 0 # Particle lies on boundary between 0 and 1 

        p2 = (x2,y2,z2, ux2,uy2,uz2, weight2,bitflags2,cell_index2)

        # Create the DT particle record type
        pvars = ['x', 'y', 'z', 'ux', 'uy', 'uz', 'weight', 'bitflags', 'cell_index']
        pvartypes = [np_M.float64]*7
        pvartypes.append(np_M.int32) # bitflags
        pvartypes.append(np_M.int32) # cell_index

        p_dtype = {'names' : pvars, 'formats': pvartypes}

        # Put the particles into an ndarray with the above type
        nparticles = 3
        particles = np_M.empty(nparticles, dtype=p_dtype)
        particles[0] = p0
        particles[1] = p1
        particles[2] = p2

        # Allocate storage for the number-density values.  Then
        # number-density array stores the integral of the physical
        # charge-distribution times the element basis functions.

        number_density_element_type = 'Lagrange'
        number_density_element_degree = 1
        number_density_field_type = 'scalar'
        number_density = Field_C(meshCI=mesh1DCI,
                                 element_type=number_density_element_type,
                                 element_degree=number_density_element_degree,
                                 field_type=number_density_field_type)

#        print "size of number_density:", number_density.function.vector().size()

        # The expected number-density values
        n_expected = np_M.empty(number_density.function.vector().size(), dtype=np_M.float64)

        n_expected[0] = 1.5e10
        n_expected[1] = 6.5e10
        n_expected[2] = 1.0e10

        for p in particles:
            number_density.integrate_delta_function(p)

        n_calc = number_density.function.vector().array()
#        print n_calc

        # Compare results
        function_space = number_density.function_space
        gdim = number_density.mesh_gdim
        # Reshape the coordinates to get (x), or (x,y), or (x,y,z) tuples.
        if df_M.DOLFIN_VERSION_STRING > "1.5.0":
            dofcoords = function_space.tabulate_dof_coordinates().reshape((-1, gdim))
        else:
            print '\n!!!WARNING!!!: ', fncname, ": DOLFIN too old.  Skipping rest of test"
            return


#        print "dofcoords=", dofcoords

        for idof in range(len(dofcoords)):
#            print "n_calc, n_expected =", n_calc[idof], n_expected[idof]
            self.assertAlmostEqual(n_calc[idof], n_expected[idof], places=3, msg="Wrong value of number_density")
        return

if __name__ == '__main__':
    unittest.main()
