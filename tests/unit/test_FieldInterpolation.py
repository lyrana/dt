#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import math
import unittest
import numpy as np_m

import dolfin as df_m

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

# Use the C++ functions in the dolfin_functions_solib.so library
import dolfin_functions_solib as df_so

class TestFieldInterpolation(unittest.TestCase):
    """Test the field-to-particle interpolation functions in Field_C class"""
    
    def setUp(self):
        # initializations for each test go here...

        # Create mesh from a file
        coordinateSystem = 'Cartesian'
        mesh2D_M = Mesh_C(mesh_file="mesh_quarter_circle_crossed.xml", coordinate_system=coordinateSystem, compute_dictionaries=True, compute_cpp_arrays=False, compute_tree=True, plot_flag=False)

#        df_m.plot(self.mesh, title='cylindrical mesh', axes=True)
#        df_m.interactive()

        phi_element_type = 'Lagrange'
        phi_element_degree = 1
        self.phi = Field_C(mesh2D_M,
                           element_type=phi_element_type,
                           element_degree=phi_element_degree,
                           field_type='scalar')

        # Read the potential from a file
        file = df_m.File("phi_test_2_2D.xml")
        file >> self.phi.function

        # Plot phi

#        df_m.plot(self.phi)
#        df_m.interactive()

        if phi_element_degree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electric_field_element_type = "DG"
        else:
            electric_field_element_type = "Lagrange"

        self.neg_electric_field = Field_C(mesh2D_M,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        # Read the electric field from a file
        file = df_m.File("negE_test_2_2D.xml")
        file >> self.neg_electric_field.function

        # Create the points for interpolation of the E-field
        
        # 1st point
        x0 = 1.0; y0 = 0.0; z0 = 0.0
        ux0 = 3000.0; uy0 = 2000.0; uz0 = 1000.0
        weight0 = 1.0 # number of electrons per macroparticle

#        p0 = (x0,y0,z0, ux0,uy0,uz0, weight0)

# Can these be named? Yes, with a different dtype.
        p0 = np_m.array([x0,y0,z0, ux0,uy0,uz0, weight0], dtype=float)

        # 2nd point
        x1 = 5.0; y1 = 0.0; z1 = 1.0
        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 2.0

#        p1 = (x1,y1,z1, ux1,uy1,uz1, weight1)
        p1 = np_m.array([x1,y1,z1, ux1,uy1,uz1, weight1], dtype=float)

        # 3nd point: same cell as 2nd point, so same E for DG0 a element.
        x2 = 4.5; y2 = 0.0; z2 = 1.0
        ux2 = uy2 = 0.0; uz2 = -uz0
        weight2 = 2.0

#        p1 = (x1,y1,z1, ux1,uy1,uz1, weight1)
        p2 = np_m.array([x2,y2,z2, ux2,uy2,uz2, weight2], dtype=float)

        self.points = np_m.array([p0, p1, p2])

        # dtype of E field at the particles
        self.Ecomps = ['x', 'y',]
        Eexpected = np_m.empty([len(self.points), len(self.Ecomps)], dtype=float)
        Eexpected[0] = [-0.84919658, -0.03336502]
        Eexpected[1] = [-0.19694748, -0.00773809]
        Eexpected[2] = [-0.19694748, -0.00773809]
        self.E_expected = Eexpected
        
        flpoint = np_m.float64
        self.E_points = np_m.empty(self.points.shape[0], dtype={'names': self.Ecomps, 'formats': (flpoint, flpoint)})
        
        return

    def test_1_interpolate_vectorField_to_points(self):
        """Provide 3 points and compute the values of a vector field at these points.

           The point data has spatial coordinates and other data, such as velocities.
           The number of spatial coordinates can be greater than the spatial
           dimension of the vector field. E.g., a point can have (x, y, z) spatial
           coordinates, while the vector may be defined only on an (x, y) space. The
           extra particle coordinates are ignored in evaluating the vector field.

        """
        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')
        
        self.neg_electric_field.interpolate_field_to_points(self.points, self.E_points)

        # Check the interpolated electric field against the expected values
        for ip in range(len(self.points)):
            for ic in range(len(self.Ecomps)):
#                print "Ecalc, Eexpect =", self.E_points[ip][ic], self.E_expected[ip][ic]
                self.assertAlmostEqual(self.E_points[ip][ic], self.E_expected[ip][ic], places=3, msg="Wrong value of E")

        return
#    def test_1_interpolate_vectorField_to_points(self): ENDDEF

    def test_2_cpp_interpolate_vectorField_to_points(self):
        """Provide 3 points and compute the values of a vector field at these points.

           The point data has spatial coordinates and other data, such as velocities.
           The number of spatial coordinates can be greater than the spatial
           dimension of the vector field. E.g., a point can have (x, y, z) spatial
           coordinates, while the vector may be defined only on an (x, y) space. The
           extra particle coordinates are ignored in evaluating the vector field.

        """
        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

  # template <Ptype PT, typename Ftype>
  # void interpolate_field_to_points(dolfin::Function& field,
  #                                  py::array_t<Pstruct<PT>, 0> points,
  #                                  py::ssize_t npoints,
  #                                  Ftype& field_at_points)
        # Create the name of the specialized MeshEntityArrays class with the right
        # number of cell-facets
        tDim = mesh_df.topology().dim()
        nFacets = tDim + 1
        meaClass = "MeshEntityArrays_" + str(nFacets) + "_facets"
        meaCtor = getattr(mea_so, meaClass)

        
        self.neg_electric_field.interpolate_field_to_points(self.points, self.E_points)

        # Check the interpolated electric field against the expected values
        for ip in range(len(self.points)):
            for ic in range(len(self.Ecomps)):
#                print "Ecalc, Eexpect =", self.E_points[ip][ic], self.E_expected[ip][ic]
                self.assertAlmostEqual(self.E_points[ip][ic], self.E_expected[ip][ic], places=3, msg="Wrong value of E")

        return
#    def test_2_cpp_interpolate_vectorField_to_points(self): ENDDEF

if __name__ == '__main__':
    unittest.main()
