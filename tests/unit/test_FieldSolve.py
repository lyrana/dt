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
import dolfin as df_m

from Dolfin_Module import Field_C
from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C


class TestFieldSolve(unittest.TestCase):
    """Test the field solvers in UserMesh_y_Fields_FE2D_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print "\n", fncName, ": This is DOLFIN Version", df_m.DOLFIN_VERSION_STRING, '\n'

        return
#    def setUp(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_1_1D_spherical_laplace_solver(self):
        """Test a 1D Laplace equation in spherical coordinates.

           This is a Laplace solve as there is no source term, only
           boundary-conditions.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ## Specialized modules for the mesh and field solver
        from UserMesh_y_Fields_Spherical1D_Module import UserMeshInput_C
        from UserMesh_y_Fields_Spherical1D_Module import UserMesh_C
        from UserMesh_y_Fields_Spherical1D_Module import UserPoissonSolve_C

        ## Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        # For this calculation, we use the Dolfin framework for finite-element fields.
#        field_infrastructure_module = 'Dolfin_Module'
#        fI_M = im_M.import_module(field_infrastructure_module)

        umiCI = UserMeshInput_C()

        ## Make the mesh

        # radial
        umiCI.rmin, umiCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umiCI.nr = 10 # Number of divisions in r direction
        umiCI.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        umiCI.field_boundary_dict = fieldBoundaryDict

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        meshCI = UserMesh_C(meshInputCI=umiCI, compute_tree=False, plot_flag=True, plot_title=plotTitle)

        ## Storage for the potential and electric field

        phiElementType = 'Lagrange'
        phiElementDegree = 2
        phiFieldType = 'scalar'

        phi = Field_C(meshCI=meshCI,
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

        negElectricField = Field_C(meshCI=meshCI,
                                      element_type=electricFieldElementType,
                                      element_degree=phiElementDegree-1,
                                      field_type='vector')

        ## The Laplace solver parameters
        # Gives PETSc LU solver, (null). (Direct solver).
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
        laplacesolveCI = UserPoissonSolve_C(phi,
                                            linearSolver, preconditioner,
#                                            boundary_marker, phi_rmin, phi_rmax,
                                            fieldBoundaryMarker, phiBCs,
                                            neg_electric_field=negElectricField)

        # Set the source term: it's zero for Laplace's equation.
        # This sets the b vector
        laplacesolveCI.assemble_source_expression(0.0)

        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        laplacesolveCI.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_m.File('phi1D.pvd')
        # phi_name goes into the output file; phi_label doesn't
        phi.function.rename("phi1D", "phi_label")
        file << phi.function
        file = df_m.File('phi1D.xml')
        file << phi.function
        
        # Write -E to a file in VTK and XML formats
        # if (fsiCI.computeEflag == True):
        #     file = df_m.File('negE1D.pvd')
        #     file << fieldsolveCI.negE
        #     file = df_m.File('negE1D.xml')
        #     file << fieldsolveCI.negE
        # Doesn't work for a 1D vector field!

        # Check the 1/r^2 fall-off in E:

        Emin, Emax = negElectricField.function(umiCI.rmin), negElectricField.function(umiCI.rmax)

#        print "Emin, Emax =", Emin, Emax
        ratio = Emin*umiCI.rmin**2/(Emax*umiCI.rmax**2)
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

            negEget = laplacesolveCI.negE.vector().array()
    #        print "negEget =", fieldsolveCI.negE.vector().array()

            # Check that the derivative equals the stored negE
            for ic in range(ncells):
                self.assertAlmostEqual(negEget[ic], negEexp[ic], msg="Expected and stored values of negE are not the same")

        return
#    def test_1_1D_spherical_laplace_solver(self):ENDDEF


#class TestFieldSolve(unittest.TestCase):
    def test_2_2D_cyl_laplace_solver(self):
        """Test a 2D Laplace equation in cylindrical coordinates.

           This is a Laplace solve as there is no source term, only
           boundary-conditions.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ## Specialized modules for the mesh and field solver
        from UserMesh_y_Fields_FE2D_Module import UserMesh_C
        from UserMesh_y_Fields_FE2D_Module import UserMeshInput_C
        from UserMesh_y_Fields_FE2D_Module import UserPoissonSolve_C

        # Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        # For this calculation, we use the Dolfin framework for finite-element fields.
#        field_infrastructure_module = 'Dolfin_Module'
#        fI_M = im_M.import_module(field_infrastructure_module)

        umiCI = UserMeshInput_C()

        # Make the mesh
        # radial
        umiCI.rmin, umiCI.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umiCI.nr = 10 # Number of divisions in r direction
        umiCI.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        umiCI.field_boundary_dict = fieldBoundaryDict

        # theta, starts at 0
        umiCI.tmax = math.pi/2 # quarter-circle
        umiCI.nt = 20  # Number of divisions in theta direction

        # The diagonal that makes the triangular mesh
        # Options: 'left, 'right', 'left/right', 'crossed'
        umiCI.diagonal = 'crossed'

        meshCI = UserMesh_C(meshInputCI=umiCI, compute_tree=False, plot_flag=False)

        # Storage for the potential and electric field
        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi = Field_C(meshCI=meshCI,
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

        negElectricField = Field_C(meshCI=meshCI,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        # The Laplace solver parameters

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
        
        phiBCs = dict((bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        computeEflag = True

        laplacesolveCI = UserPoissonSolve_C(phi,
                                            linearSolver, preconditioner,
                                            fieldBoundaryMarker, phiBCs,
                                            neg_electric_field=negElectricField)

        laplacesolveCI.assemble_source_expression(0.0)

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        laplacesolveCI.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_m.File("phi2D_crossed.pvd")
        phi.function.rename("phi2D", "phi_label")
        file << phi.function
        file = df_m.File("phi2D_crossed.xml")
        file << phi.function
        
        # Write -E to a file in VTK and XML formats
        if negElectricField is not None:
            negElectricField.function.rename("E2D", "E_label")
            file = df_m.File("negE2D_crossed.pvd")
            file << negElectricField.function
            file = df_m.File("negE2D_crossed.xml")
            file << negElectricField.function

        return
#    def test_2_2D_cyl_laplace_solver(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_3_1D_poisson_solver(self):
        """Test a 1D Poisson equation in cartesian coordinates.

           Note: the mesh, mesh BC markers, and the charge-density are read from
           external files written by test_4_compute_charge_density_on_1Dmesh() in
           test_ChargeDensity.py. Boundary values are specified here.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ### Specialized mesh and field-solver modules for this test ###

        from UserMesh_FE_XYZ_Module import UserMesh_C
        from UserMesh_FE_XYZ_Module import UserMeshInput_C
        from UserMesh_FE_XYZ_Module import UserPoissonSolve1D_C

        ### Set plot flag ###

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        ### Read the mesh and boundary-condition markers ###

        # Create and read in the mesh object
        mesh1d_M = Mesh_C(mesh_file='mesh-2cell-1D.xml', field_boundary_marker_file='mesh-2cell-1D_Fbcs.xml', compute_dictionaries=True, compute_tree=True, plot_flag=False)

        # Read the field boundary marker
        field_boundary_marker_file = df_m.File('mesh-2cell-1D_Fbcs.xml')
        field_boundary_marker_file >> mesh1d_M.field_boundary_marker

        ### Create the charge-density vector ###

        # Copy over from test_ChargeDensity.py:test_4
        # This associates the field marker integers with boundaries.
        # IS THIS NEEDED?
        xminIndx = 1
        xmaxIndx = 2
        fieldBoundaryDict = {'xmin': xminIndx,
                             'xmax': xmaxIndx,
                             }

        chargeDensityElementType = 'Lagrange'
        chargeDensityElementDegree = 1
        chargeDensityFieldType = 'scalar'

        chargeDensity_F = Field_C(meshCI=mesh1d_M,
                                  element_type=chargeDensityElementType,
                                  element_degree=chargeDensityElementDegree,
                                  field_type=chargeDensityFieldType)

        ### Set the charge density values ###

        # Read in the charge-density written by
        # test_ChargeDensity.py:test_4_compute_charge_density_on_1Dmesh
        file = df_m.File("charge-1D.xml")
        file >> chargeDensity_F.function

        ### Set the Poisson solver parameters ###

        linearSolver = 'lu'
        preconditioner = None

        ### Set the field boundary conditions ###

        # Dirichlet boundaries are applied using a mesh function that marks the outer
        # edges of the mesh.
        fieldBoundaryMarker = mesh1d_M.field_boundary_marker
       
        # Specify the boundary values of the potential
        phiVals = {'xmin':  -2.0,
                   'xmax': 1.0,
                   }

        phiBCs = dict((bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        ### Create vectors for for the potential and electric field ###

        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi_F = Field_C(meshCI=mesh1d_M,
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

        negElectricField = Field_C(meshCI=mesh1d_M,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        ### Create the Poisson-solver object ###

        poissonsolveCI = UserPoissonSolve1D_C(phi_F,
                                              linearSolver, preconditioner,
                                              fieldBoundaryMarker,
                                              phiBCs,
                                              assembled_charge=chargeDensity_F.function,
                                              neg_electric_field=negElectricField)

        # Create the source term from a given density
        # (In this test, the charge comes from kinetic point particles, rather
        # than from a density function, so the following isn't needed.)
#        poissonsolveCI.assemble_source_expression(charge_density)

#        self.b = df_m.assemble(self.L)


        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolveCI.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

        return
#    def test_3_1D_poisson_solver(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
