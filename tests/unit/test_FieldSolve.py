#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import importlib as im_M

import numpy as np_M
import dolfin as df_M

from UserMesh_y_Fields_Spherical1D_Module import *

class TestPoissonSolve(unittest.TestCase):
    """Test the field solvers in UserMesh_y_Fields_FE2D_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print "\n", fncName, ": This is DOLFIN Version", df_M.DOLFIN_VERSION_STRING, '\n'

        return
#    def setUp(self):ENDDEF

#class TestPoissonSolve(unittest.TestCase):
    def test_1D_poisson_solver(self):
        """Test a 1D Laplace equation in spherical coordinates.

           This is a Laplace solve as there is no source term, only
           boundary-conditions.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ## Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        miCI = UserMeshInput_C()

        ## Make the mesh

        # radial
        miCI.rmin, miCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        miCI.nr = 10 # Number of divisions in r direction
        miCI.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        miCI.field_boundary_dict = fieldBoundaryDict

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        meshCI = UserMesh_C(meshInputCI=miCI, compute_tree=False, plot_flag=True, plot_title=plotTitle)

        ## Storage for the potential and electric field

        # For this calculation, we use the Dolfin framework for finite-element fields.
        field_infrastructure_module = 'Dolfin_Module'
        fI_M = im_M.import_module(field_infrastructure_module)

        phiElementType = 'Lagrange'
        phiElementDegree = 2
        phiFieldType = 'scalar'

        phi = fI_M.Field_C(meshCI=meshCI,
                           element_type=phiElementType,
                           element_degree=phiElementDegree,
                           field_type=phiFieldType)

        if phiElementDegree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electricFieldElementType = 'DG'
        else:
            electricFieldElementType = 'Lagrange'

        negElectricField = fI_M.Field_C(meshCI=meshCI,
                                      element_type=electricFieldElementType,
                                      element_degree=phiElementDegree-1,
                                      field_type='vector')

        ## The Poisson solver parameters

        linearSolver = 'lu'
        preconditioner = None

        # Dirichlet boundaries are marked with a boundary marker mesh function
        fieldBoundaryMarker = meshCI.field_boundary_marker
       
        # Boundary values of the potential
        phiVals = {'rmin':  0.0,
                    'rmax': -1.5,
                    }

        phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        # Compute the electrostatic field from phi
        # Have to pass "mesh" because a SpatialCoordinate is used in Poisson's equation.
        poissonsolveCI = UserPoissonSolve_C(phi,
                                            linearSolver, preconditioner,
#                                            boundary_marker, phi_rmin, phi_rmax,
                                            fieldBoundaryMarker, phiBCs,
                                            neg_electric_field=negElectricField)

        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolveCI.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_M.File('phi1D.pvd')
        # phi_name goes into the output file; phi_label doesn't
        phi.function.rename("phi1D", "phi_label")
        file << phi.function
        file = df_M.File('phi1D.xml')
        file << phi.function
        
        # Write -E to a file in VTK and XML formats
        # if (fsiCI.computeEflag == True):
        #     file = df_M.File('negE1D.pvd')
        #     file << fieldsolveCI.negE
        #     file = df_M.File('negE1D.xml')
        #     file << fieldsolveCI.negE
        # Doesn't work for a 1D vector field!

        # Check the 1/r^2 fall-off in E:

        Emin, Emax = negElectricField.function(miCI.rmin), negElectricField.function(miCI.rmax)

#        print "Emin, Emax =", Emin, Emax
        ratio = Emin*miCI.rmin**2/(Emax*miCI.rmax**2)
#        print "ratio:", ratio
        
        decimal_places = 1 # Note: can get 2 decimal places for degree-3 elements.
        self.assertAlmostEqual(ratio, 1.0, places=decimal_places, msg="Falloff in E is not 1/r^2 to %d decimals" % decimal_places)

        r_arr = meshCI.mesh.coordinates()[:,0]
        phi_arr = phi.function.vector().array()
        nnodes = len(phi_arr)
        nvertices = len(r_arr)
#        print "nnode =", nnodes
#        print "nvertex =", nvertex
#        print "r_arr =", r_arr
#        print "phi_arr =", phi_arr

        # Compute E from a simple derivative, for linear elements
        if nnodes == nvertices:
            # Check the derivative: (nvertices*dim array, so need the ,0)

            # Set up temp arrays
            ncells = len(r_arr)-1
            negEexp = np_M.zeros(ncells)
            dr = np_M.zeros(ncells)

            dr[0:ncells] = r_arr[1:ncells+1] - r_arr[0:ncells]

            # Have to reverse phi_arr to get it in same order as r_arr
            phi_arr = phi_arr[::-1]

            # Compute dphi/dr
            negEexp[0:ncells] = (phi_arr[1:ncells+1] - phi_arr[0:ncells])/dr[0:ncells]
     #       print "negEexp =", negEexp

            negEget = poissonsolveCI.negE.vector().array()
    #        print "negEget =", fieldsolveCI.negE.vector().array()

            # Check that the derivative equals the stored negE
            for ic in range(ncells):
                self.assertAlmostEqual(negEget[ic], negEexp[ic], msg="Expected and stored values of negE are not the same")

        return
#    def test_1D_poisson_solver(self):ENDDEF

#class TestPoissonSolve(unittest.TestCase):
    def test_2D_poisson_solver(self):
        """Test a 2D Laplace equation in cylindrical coordinates.

           This is a Laplace solve as there is no source term, only
           boundary-conditions.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        # Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        miCI = UserMeshInput_C()

        # Make the mesh
        # radial
        miCI.rmin, miCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        miCI.nr = 10 # Number of divisions in r direction
        miCI.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        miCI.field_boundary_dict = fieldBoundaryDict

        # theta, starts at 0
        miCI.tmax = math.pi/2 # quarter-circle
        miCI.nt = 20  # Number of divisions in theta direction

        # The diagonal that makes the triangular mesh
        # Options: 'left, 'right', 'left/right', 'crossed'
        miCI.diagonal = 'crossed'

        from UserMesh_y_Fields_FE2D_Module import UserMesh_C
        meshCI = UserMesh_C(meshInputCI=miCI, compute_tree=False, plot_flag=False)

        # For this calculation, we use the Dolfin framework for finite-element fields.
        field_infrastructure_module = 'Dolfin_Module'
        fI_M = im_M.import_module(field_infrastructure_module)

        # Storage for the potential and electric field
        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi = fI_M.Field_C(meshCI=meshCI,
                           element_type=phiElementType,
                           element_degree=phiElementDegree,
                           field_type=phiFieldType)

        if phiElementDegree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To get these values, we need Discontinuous Galerkin
            # elements.
            electricFieldElementType = "DG"
        else:
            electricFieldElementType = "Lagrange"

        negElectricField = fI_M.Field_C(meshCI=meshCI,
                                      element_type=electricFieldElementType,
                                      element_degree=phiElementDegree-1,
                                      field_type='vector')

        # The Poisson solver parameters

        linearSolver = 'lu'
        preconditioner = None
# Don't get exactly the same solutions with the following
#        linearSolver = 'cg'
#        preconditioner = 'ilu'

        # Dirichlet Boundaries are marked with a boundary marker mesh
        # function
        fieldBoundaryMarker = meshCI.field_boundary_marker
       
        # Boundary values of the potential.  These names have to be
        # the same as those assigned to the boundaries above.
        phiVals = {'rmin':  0.0,
                    'rmax': -1.5,
                    }

        phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        computeEflag = True

        from UserMesh_y_Fields_FE2D_Module import UserPoissonSolve_C

        poissonsolveCI = UserPoissonSolve_C(phi,
                                            linearSolver, preconditioner,
                                            fieldBoundaryMarker, phiBCs,
                                            neg_electric_field=negElectricField)

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolveCI.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_M.File("phi2D_crossed.pvd")
        phi.function.rename("phi2D", "phi_label")
        file << phi.function
        file = df_M.File("phi2D_crossed.xml")
        file << phi.function
        
        # Write -E to a file in VTK and XML formats
        if negElectricField is not None:
            negElectricField.function.rename("E2D", "E_label")
            file = df_M.File("negE2D_crossed.pvd")
            file << negElectricField.function
            file = df_M.File("negE2D_crossed.xml")
            file << negElectricField.function

        return
#    def test_2D_poisson_solver(self):ENDDEF

#class TestPoissonSolve(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
