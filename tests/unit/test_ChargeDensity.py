#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest
import numpy as np_M

import dolfin as df_M

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

# Here's the user's mesh definition for the 2D test
from UserMesh_FE_XYZ_Module import *


class TestChargeDensity(unittest.TestCase):
    """Test the particle-to-field functions in class Field_C

       See test_ChargeDensity.ods
"""
    
    def setUp(self):
        # initializations for each test go here...

        return

#class TestChargeDensity(unittest.TestCase):
    def test_1_interpolate_particleDensity_to_1Dmesh(self):
        """Compute the number-density generated by planar particles.

           Macroparticles are created within a 1D meshed region and
           are weighted to nodal points on the mesh.
           
           No species data is defined for the particles.  No segmented-array
           particle storage is used for the particles, just one numpy array is
           used, with a particle dtype.

           See test_ChargeDensity.ods:Test1 for calculated values of the density.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'
        
        # Create a 1D mesh from -0.5 to 0.5
        mesh = df_M.IntervalMesh(2, -0.5, 0.5)

        mesh1DCI = Mesh_C(Mesh=mesh, compute_dictionaries=True, compute_tree=True, plot_flag=False)

        # Put 3 DT particles in the meshed region.

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
        # particle weight-distribution times the element basis
        # functions.

        number_density_element_type = 'Lagrange'
        number_density_element_degree = 1
        number_density_field_type = 'scalar'
        number_density = Field_C(meshCI=mesh1DCI,
                                 element_type=number_density_element_type,
                                 element_degree=number_density_element_degree,
                                 field_type=number_density_field_type)

#        print "size of number_density:", number_density.function.vector().size()

        # The expected number-density values are put into a numpy array
        n_expected = np_M.empty(number_density.function.vector().size(), dtype=np_M.float64)

        # Note that these are in the vertex numbering order, not DoF order
        n_expected[0] = 1.0e10
        n_expected[1] = 6.5e10
        n_expected[2] = 1.5e10

        for p in particles:
            number_density.integrate_delta_function(p)
#            number_density.interpolate_delta_function_to_dofs(p)

        # array() returns a numpy array that has a copy of the values in vector().
        n_calc = number_density.function.vector().array()
#        print n_calc

        # Compare results
        function_space = number_density.function_space
        gdim = number_density.mesh_gdim
        # Reshape the coordinates to get (x), or (x,y), or (x,y,z) tuples. The
        # index "-1" is short-hand for the last element of the array.
        if df_M.DOLFIN_VERSION_STRING > "1.5.0":
            dofcoords = function_space.tabulate_dof_coordinates().reshape((-1, gdim))
        else:
            print '\n!!!WARNING!!!: ', fncName, ": DOLFIN too old.  Skipping rest of test"
            return


#        print "dofcoords=", dofcoords

        # Convert vertex indices to DoF indices.  For CG1 elements, there are
        # the same number of vertices as DoFs.
        v2d=df_M.vertex_to_dof_map(function_space)

        for ivert in range(len(n_expected)):
            idof = v2d[ivert]
#            print "vertex", ivert, "is DoF index", idof
#            print "n_calc, n_expected =", n_calc[idof], n_expected[ivert]
            self.assertAlmostEqual(n_calc[idof], n_expected[ivert], places=3, msg="Wrong value of number_density")

        return
#    def test_1_interpolate_particleCharge_to_1Dmesh(self):ENDDEF

#class TestChargeDensity(unittest.TestCase):
    def test_2_interpolate_particleCharge_to_2Dmesh(self):
        """Compute the number-density generated by line particles on a 2D mesh

           No species data is defined for the particles.  No segmented-array
           particle storage is used for the particles, just one numpy array is
           used, with a particle dtype.

           See test_ChargeDensity.ods:Test2 for calculated values of the
           density.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'
        
        # Describe a 2D mesh from (-10,-10) to (10,10) with 2 cells on a side.
        mi2DCI = UserMeshInput_C()
        mi2DCI.pmin = df_M.Point(-10.0, -10.0)
        mi2DCI.pmax = df_M.Point(10.0, 10.0)
        mi2DCI.cells_on_side = (2, 2)
        mi2DCI.diagonal = 'left'

        # UserMesh_FE_XYZ_Module can make the mesh from the above input.
        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"

        pmesh2DCI = UserMesh_C(mi2DCI, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle)
#        pmesh2DCI.compute_cell_vertex_dict()
#        pmesh2DCI.compute_cell_dict()

        # Put 3 particles inside the meshed region

        x0 = -5.0; y0 = -5.0; z0 = 0.0
        ux0 = 0.0; uy0 = 0.0; uz0 = 0.0
        weight0 = 2.0e10 # number of electrons per macroparticle
        bitflags0 = 0b0
        cell_index0 = 1

        p0 = (x0,y0,z0, ux0,uy0,uz0, weight0,bitflags0,cell_index0)

        # 2nd particle
        x1 = 1.0; y1 = 1.0; z1 = 1.0
        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 3.0e10
        bitflags1 = 0b0
        cell_index1 = 6

        p1 = (x1,y1,z1, ux1,uy1,uz1, weight1,bitflags1,cell_index1)

        # 3nd particle
        x2 = -9.0; y2 = 1.0; z2 = 1.0
        ux2 = uy2 = 0.0; uz2 = -uz0
        weight2 = 4.0e10
        bitflags2 = 0b0
        cell_index2 = 4 # Particle lies on boundary between 0 and 1 

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
        number_density = Field_C(meshCI=pmesh2DCI,
                                 element_type=number_density_element_type,
                                 element_degree=number_density_element_degree,
                                 field_type=number_density_field_type)

#        print "size of number_density:", number_density.function.vector().size()

        # The expected number-density values from test_ChargeDensity.ods:Test2
        n_expected = np_M.empty(number_density.function.vector().size(), dtype=np_M.float64)

        # These are in vectex, not DOF order
        n_expected[0] = 0.0
        n_expected[1] = 1.0e10
        n_expected[2] = 0.0
        n_expected[3] = 4.2e10
        n_expected[4] = 2.8e10
        n_expected[5] = 3.0e9
        n_expected[6] = 4.0e9
        n_expected[7] = 3.0e9
        n_expected[8] = 0.0

        for p in particles:
#            number_density.integrate_delta_function(p)
            number_density.interpolate_delta_function_to_dofs(p)

        n_calc = number_density.function.vector().array()
#        print fncName, n_calc

        # Compare results
        function_space = number_density.function_space
        gdim = number_density.mesh_gdim
        # Reshape the coordinates to get (x), or (x,y), or (x,y,z) tuples.
        if df_M.DOLFIN_VERSION_STRING > "1.5.0":
            dofcoords = function_space.tabulate_dof_coordinates().reshape((-1, gdim))
        else:
            print '\n!!!WARNING!!!: ', fncName, ": DOLFIN too old.  Skipping rest of test"
            return

#        print "dofcoords=", dofcoords

        # Plot the result

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": number density"
        df_M.plot(number_density.function, title=plotTitle)
        df_M.interactive()

        ## Check the values vs. those computed in test_ChargeDensity.ods:Test2

        # Convert vertex indices to DoF indices.  For CG1 elements, there are
        # the same number of vertices as DoFs.
        v2d=df_M.vertex_to_dof_map(function_space)

        for ivert in range(len(n_expected)):
            idof = v2d[ivert]
#            print "vertex", ivert, "is DoF index", idof
#            print "n_calc, n_expected =", n_calc[idof], n_expected[ivert]
            # "places" means numbers after the decimal point:
            self.assertAlmostEqual(n_calc[idof], n_expected[ivert], places=4, msg="Wrong value of number_density")
        return
#    def test_2_interpolate_particleCharge_to_2Dmesh(self):ENDDEF

#class TestChargeDensity(unittest.TestCase):ENDCLASS


if __name__ == '__main__':
    unittest.main()
