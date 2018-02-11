#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import math
import unittest

import importlib as im_m

import numpy as np_m
import dolfin as df_m

from Dolfin_Module import Field_C
from Dolfin_Module import Mesh_C
from Dolfin_Module import PoissonSolve_C


class TestFieldSolve(unittest.TestCase):
    """Test the field solvers in UserMesh_y_Fields_FE2D_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...

        fncName = '(' + __file__ + ') ' + sys._getframe().f_code.co_name + "():"
        print "\n", fncName, "This is DOLFIN Version", df_m.DOLFIN_VERSION_STRING

        return
#    def setUp(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_1_1D_spherical_laplace_solver(self):
        """Test a 1D Laplace equation in spherical coordinates.

           This is a Laplace solve as there is no source term, only
           boundary-conditions.

           NB: The potential is represented with quadratic elements are used instead
           of the usual linear ones.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ## Specialized modules for the mesh and field solver
        from UserMesh_y_Fields_Spherical1D_Module import UserMeshInput1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserMesh1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserPoissonSolve1DS_C

        ## Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        # For this calculation, we use the Dolfin framework for finite-element fields.
#        field_infrastructure_module = 'Dolfin_Module'
#        fI_M = im_m.import_module(field_infrastructure_module)

        umi = UserMeshInput1DS_C()

        ## Make the mesh

        # radial
        umi.rmin, umi.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umi.nr = 10 # Number of divisions in r direction
        umi.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        umi.field_boundary_dict = fieldBoundaryDict

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        mesh_M = UserMesh1DS_C(umi, compute_tree=False, plot_flag=True, plot_title=plotTitle)
#        yesno = raw_input("Just called UserMesh1DS_C() in test_1_1D_spherical_laplace_solver")

        ## Storage for the potential and electric field

        phiElementType = 'Lagrange'
        phiElementDegree = 2
        phiFieldType = 'scalar'

        phi_F = Field_C(mesh_M=mesh_M,
                           element_type=phiElementType,
                           element_degree=phiElementDegree,
                           field_type=phiFieldType)

        if phiElementDegree == 1:
            # For linear elements, grad(phi_F) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electricFieldElementType = 'DG'
        else:
            electricFieldElementType = 'Lagrange'

        negElectricField_F = Field_C(mesh_M=mesh_M,
                                      element_type=electricFieldElementType,
                                      element_degree=phiElementDegree-1,
                                      field_type='vector')

        ## The Laplace solver parameters
        # Gives PETSc LU solver, (null). (Direct solver).
        linearSolver = 'lu'
        preconditioner = None

        # Dirichlet boundaries are marked with a boundary marker mesh function
        fieldBoundaryMarker = mesh_M.field_boundary_marker
       
        # Boundary values of the potential
        phiVals = {'rmin':  0.0,
                    'rmax': -1.5,
                    }

        phiBCs = dict( (bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        # Compute the electrostatic field from phi
        laplacesolve = UserPoissonSolve1DS_C(phi_F,
                                               linearSolver, preconditioner,
                                               fieldBoundaryMarker, phiBCs,
                                               neg_electric_field=negElectricField_F)

        # Set the source term: it's zero for Laplace's equation.
        # This sets the b vector
        laplacesolve.assemble_source_expression(0.0)

        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        laplacesolve.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_m.File('phi_test_1_1D.pvd')
        # phi_name goes into the output file; phi_label doesn't
        phi_F.function.rename("phi1D", "phi_label")
        file << phi_F.function
        file = df_m.File('phi_test_1_1D.xml')
        file << phi_F.function
        
        # Write -E to a file in VTK and XML formats
        # !!! VTK cannot write a vector field on a 1D grid.
        # file = df_m.File('negE_test_1_1D.pvd')
        # file << fieldsolveCI.negE
        file = df_m.File('negE_test_1_1D.xml')
        file << negElectricField_F.function
        
        # Check the 1/r^2 fall-off in E:

        (Emin, Emax) = (negElectricField_F.function(umi.rmin), negElectricField_F.function(umi.rmax))
#        print "Emin, Emax =", Emin, Emax
        ratio = Emin*umi.rmin**2/(Emax*umi.rmax**2)
#        print "ratio:", ratio
        
        decimal_places = 1 # Note: can get 2 decimal places for degree-3 elements.
        self.assertAlmostEqual(ratio, 1.0, places=decimal_places, msg="Falloff in E is not 1/r^2 to %d decimals" % decimal_places)

        r_arr = mesh_M.mesh.coordinates()[:,0]
        phi_arr = phi_F.function.vector().get_local()
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
            negEexp = np_m.zeros(ncells)
            dr = np_m.zeros(ncells)

            dr[0:ncells] = r_arr[1:ncells+1] - r_arr[0:ncells]

            # Have to reverse phi_arr to get it in same order as r_arr
            phi_arr = phi_arr[::-1]

            # Compute dphi/dr
            negEexp[0:ncells] = (phi_arr[1:ncells+1] - phi_arr[0:ncells])/dr[0:ncells]
     #       print "negEexp =", negEexp

            negEget = laplacesolve.negE.vector().get_local()
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
        print '\ntest: ', fncName

        ## Specialized modules for the mesh and field solver
        from UserMesh_y_Fields_FE2D_Module import UserMesh2DCirc_C
        from UserMesh_y_Fields_FE2D_Module import UserMeshInput2DCirc_C
        from UserMesh_y_Fields_FE2D_Module import UserPoissonSolve2DCirc_C

        # Plotting
        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        # For this calculation, we use the Dolfin framework for finite-element fields.
#        field_infrastructure_module = 'Dolfin_Module'
#        fI_M = im_m.import_module(field_infrastructure_module)

        umi = UserMeshInput2DCirc_C()

        # Make the mesh
        # radial
        umi.rmin, umi.rmax = 1.0, 5.0 # Mesh goes from rmin to rmax in radius
        umi.nr = 10 # Number of divisions in r direction
        umi.stretch = 1.3 # Stretch parameter

        # Name the Dirichlet boundaries and assign integers to them.
        # These are the boundary-name -> int pairs used to mark mesh
        # facets:
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        umi.field_boundary_dict = fieldBoundaryDict

        # theta, starts at 0
        umi.tmax = math.pi/2 # quarter-circle
        umi.nt = 20  # Number of divisions in theta direction

        # The diagonal that makes the triangular mesh
        # Options: 'left, 'right', 'left/right', 'crossed'
        umi.diagonal = 'crossed'

        mesh_M = UserMesh2DCirc_C(umi, compute_tree=False, plot_flag=False)

        # Storage for the potential and electric field
        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi_F = Field_C(mesh_M=mesh_M,
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

        negElectricField_F = Field_C(mesh_M=mesh_M,
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
        fieldBoundaryMarker = mesh_M.field_boundary_marker
       
        # Boundary values of the potential.  These names have to be
        # the same as those assigned to the boundaries above.
        phiVals = {'rmin':  0.0,
                   'rmax': -1.5,
                   }
        
        phiBCs = dict((bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        computeEflag = True

        laplacesolve = UserPoissonSolve2DCirc_C(phi_F,
                                            linearSolver, preconditioner,
                                            fieldBoundaryMarker, phiBCs,
                                            neg_electric_field=negElectricField_F)

        laplacesolve.assemble_source_expression(0.0)

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        laplacesolve.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

#        yesno = raw_input("Looks OK [Y/n]?")
#        self.assertNotEqual(yesno, 'n', "Problem with mesh")

        # Write the potential to a file in VTK and XML formats
        file = df_m.File("phi_test_2_2D.pvd")
        phi_F.function.rename("phi2D", "phi_label")
        file << phi_F.function
        file = df_m.File("phi_test_2_2D.xml")
        file << phi_F.function
        
        # Write -E to a file in VTK and XML formats
        if negElectricField_F is not None:
            negElectricField_F.function.rename("E2D", "E_label")
            file = df_m.File("negE_test_2_2D.pvd")
            file << negElectricField_F.function
            file = df_m.File("negE_test_2_2D.xml")
            file << negElectricField_F.function

        return
#    def test_2_2D_cyl_laplace_solver(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_3_1D_poisson_solver(self):
        """Test a 1D Poisson equation in cartesian coordinates.

           Note: the mesh, mesh BC markers, and the charge-density are read from
           external files written by test_4_compute_charge_density_on_1Dmesh() in
           test_ChargeDensity.py. Boundary values are specified here.

           See Calc file test_FieldSolve.ods for calculation of electric field.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ### Specialized mesh and field-solver modules for this test ###

        from UserMesh_y_Fields_FE_XYZ_Module import UserMesh_C
        from UserMesh_y_Fields_FE_XYZ_Module import UserMeshInput_C
        from UserMesh_y_Fields_FE_XYZ_Module import UserPoissonSolve1D_C

        ### Set plot flag ###

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        ### Read the mesh and boundary-condition markers ###

        # Create the mesh object and read in the mesh and boundary markers
        mesh1d_M = Mesh_C(mesh_file='mesh_2cell_1D.xml', field_boundary_marker_file='mesh_2cell_1D_Fbcs.xml', compute_dictionaries=True, compute_tree=True, plot_flag=False)

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

        chargeDensity_F = Field_C(mesh_M=mesh1d_M,
                                  element_type=chargeDensityElementType,
                                  element_degree=chargeDensityElementDegree,
                                  field_type=chargeDensityFieldType)

        ### Set the charge density values ###

        # Read in the charge-density written by
        # test_ChargeDensity.py:test_4_compute_charge_density_on_1Dmesh
        file = df_m.File("charge_1D.xml")
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

        phi_F = Field_C(mesh_M=mesh1d_M,
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

        negElectricField_F = Field_C(mesh_M=mesh1d_M,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        ### Create the Poisson-solver object ###

        poissonsolve = UserPoissonSolve1D_C(phi_F,
                                              linearSolver, preconditioner,
                                              fieldBoundaryMarker,
                                              phiBCs,
                                              assembled_charge=chargeDensity_F.function,
                                              neg_electric_field=negElectricField_F)

        # Create the source term from a given density
        # (In this test, the charge comes from kinetic point particles, rather
        # than from a density function, so the following isn't needed.)
#        poissonsolve.assemble_source_expression(charge_density)

#        self.b = df_m.assemble(self.L)


        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolve.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

        return
#    def test_3_1D_poisson_solver(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_4_1D_poisson_solver(self):
        """Test a Poisson solve in 1D radial spherical coordinates starting at r = 1.

           The mesh and mesh BC markers used here read from files written by
           test_1D_spherical_mesh in test_UserMesh.py.  

           The charge-density is from a file written by
           test_5_compute_charge_density_on_1Dmesh() in test_ChargeDensity.py.

           Boundary values are specified here.  The potential is set to zero at the
           outer radius at 5 meters.  The potential at the inner radius (1 meter) is
           set to the analytic potential due to the charge.  If the solution is
           correct, then the potential inside the radial location of the particle
           should be flat.

           See Calc file test_FieldSolve.ods:test_4 for calculation of potential

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ### Specialized mesh and field-solver modules for this test ###

        from UserMesh_y_Fields_Spherical1D_Module import UserMesh1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserMeshInput1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserPoissonSolve1DS_C

        ### Set plot flag ###

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        ### Read the mesh and boundary-condition markers ###

        # Create a mesh object and read in the mesh.
        mesh1d_M = Mesh_C(mesh_file='mesh_1D_radial.xml', field_boundary_marker_file='mesh_1D_radial_Fbcs.xml', compute_dictionaries=True, compute_tree=True, plot_flag=False)

        ### Create the charge-density vector ###

        chargeDensityElementType = 'Lagrange'
        chargeDensityElementDegree = 1
        chargeDensityFieldType = 'scalar'

        chargeDensity_F = Field_C(mesh_M=mesh1d_M,
                                  element_type=chargeDensityElementType,
                                  element_degree=chargeDensityElementDegree,
                                  field_type=chargeDensityFieldType)

        ### Set the charge density values ###

        # Read in the charge-density written by
        # test_ChargeDensity.py:test_5_compute_charge_density_on_1Dmesh
        file = df_m.File("charge_1DS.xml")
        file >> chargeDensity_F.function

        ### Set the Poisson solver parameters ###

        linearSolver = 'lu'
        preconditioner = None

        ### Set the field boundary conditions ###

        # Copy over from test_ChargeDensity.py:test_5 for convenience in using the
        # integer markers to set boundary values.
        # 
        # This associates the field marker integers with boundary names.
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        # Dirichlet boundaries are applied using a mesh function that marks the outer
        # edges of the mesh.
        fieldBoundaryMarker = mesh1d_M.field_boundary_marker
       
        # Specify the boundary values of the potential
        phiVals = {'rmin': -2.696e9, # Analytic potential inside the particle radius
#        phiVals = {'rmin': -2.72e9, # Better
                   'rmax': 0.0,
                   }

        phiBCs = dict((bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        ### Create vectors for for the potential and electric field ###

        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi_F = Field_C(mesh_M=mesh1d_M,
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

        negElectricField_F = Field_C(mesh_M=mesh1d_M,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        ### Create the Poisson-solver object ###

        poissonsolve = UserPoissonSolve1DS_C(phi_F,
                                             linearSolver, preconditioner,
                                             fieldBoundaryMarker,
                                             phiBCs,
                                             assembled_charge=chargeDensity_F.function,
                                             neg_electric_field=negElectricField_F)

        # Create the source term from a given density
        # (In this test, the charge comes from kinetic point particles, rather
        # than from a density function, so the following isn't needed.)
#        poissonsolve.assemble_source_expression(charge_density)

#        self.b = df_m.assemble(self.L)


        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolve.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

        # Write the potential to a file in VTK and XML formats
        file = df_m.File('phi_test_4_1D.pvd')
        # phi_name goes into the output file; phi_label doesn't
        phi_F.function.rename("phi1D", "phi_label")
        file << phi_F.function
        file = df_m.File('phi_test_4_1D.xml')
        file << phi_F.function

        # Write -E to a file in VTK and XML formats
        if negElectricField_F is not None:
            negElectricField_F.function.rename("E1D", "E_label")
            # !!! VTK cannot write a vector field on a 1D grid.
#            file = df_m.File("negE_test_4_1D.pvd")
#            file << negElectricField_F.function
            file = df_m.File("negE_test_4_1D.xml")
            file << negElectricField_F.function

        return
#    def test_4_1D_poisson_solver(self):ENDDEF

#class TestFieldSolve(unittest.TestCase):
    def test_5_1D_poisson_solver(self):
        """Test a Poisson solve in 1D radial spherical coordinates starting at r = 0.

           The mesh and mesh BC markers used here read from files written by
           test_1D_spherical_mesh in test_UserMesh.py.  

           The charge-density is from a file written by
           test_6_compute_charge_density_on_1Dmesh() in test_ChargeDensity.py.

           Boundary values are specified here.  The potential is set to zero at the
           outer radius at 5 meters.  The potential at the inner radius (1 meter) is
           set to the analytic potential due to the charge.  If the solution is
           correct, then the potential inside the radial location of the particle
           should be flat.

           See Calc file test_FieldSolve.ods:test_5 for calculation of potential

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ### Specialized mesh and field-solver modules for this test ###

        from UserMesh_y_Fields_Spherical1D_Module import UserMesh1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserMeshInput1DS_C
        from UserMesh_y_Fields_Spherical1D_Module import UserPoissonSolve1DS_C

        ### Set plot flag ###

        if os.environ.get('DISPLAY') is None:
            plotFlag=False
        else:
            plotFlag=True

        ### Read the mesh and boundary-condition markers ###

        # Create a mesh object and read in the mesh.
        mesh1d_M = Mesh_C(mesh_file='mesh_1D_radial_r0.xml', field_boundary_marker_file='mesh_1D_radial_r0_Fbcs.xml', compute_dictionaries=True, compute_tree=True, plot_flag=False)

        ### Create the charge-density vector ###

        chargeDensityElementType = 'Lagrange'
        chargeDensityElementDegree = 1
        chargeDensityFieldType = 'scalar'

        chargeDensity_F = Field_C(mesh_M=mesh1d_M,
                                  element_type=chargeDensityElementType,
                                  element_degree=chargeDensityElementDegree,
                                  field_type=chargeDensityFieldType)

        ### Set the charge density values ###

        # Read in the charge-density written by
        # test_ChargeDensity.py:test_5_compute_charge_density_on_1Dmesh
        file = df_m.File("charge_1DS_r0.xml")
        file >> chargeDensity_F.function

        ### Set the Poisson solver parameters ###

        linearSolver = 'lu'
        preconditioner = None

        ### Set the field boundary conditions ###

        # Copy over from test_ChargeDensity.py:test_5 for convenience in using the
        # integer markers to set boundary values.
        # 
        # This associates the field marker integers with boundary names.
        rminIndx = 1
        rmaxIndx = 2
        fieldBoundaryDict = {'rmin': rminIndx,
                             'rmax': rmaxIndx,
                             }

        # Dirichlet boundaries are applied using a mesh function that marks the outer
        # edges of the mesh.
        fieldBoundaryMarker = mesh1d_M.field_boundary_marker
       
        # Specify the boundary values of the potential
        phiVals = {'rmin': 'unset', # Analytic potential inside the particle radius
                   'rmax': 0.0,
                   }

        phiBCs = dict((bnd, [fieldBoundaryDict[bnd], phiVals[bnd]]) for bnd in fieldBoundaryDict.keys())

        ### Create vectors for for the potential and electric field ###

        phiElementType = 'Lagrange'
        phiElementDegree = 1
        phiFieldType = 'scalar'

        phi_F = Field_C(mesh_M=mesh1d_M,
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

        negElectricField_F = Field_C(mesh_M=mesh1d_M,
                                   element_type=electricFieldElementType,
                                   element_degree=phiElementDegree-1,
                                   field_type='vector')

        ### Create the Poisson-solver object ###

        poissonsolve = UserPoissonSolve1DS_C(phi_F,
                                             linearSolver, preconditioner,
                                             fieldBoundaryMarker,
                                             phiBCs,
                                             assembled_charge=chargeDensity_F.function,
                                             neg_electric_field=negElectricField_F)

        # Create the source term from a given density
        # (In this test, the charge comes from kinetic point particles, rather
        # than from a density function, so the following isn't needed.)
#        poissonsolve.assemble_source_expression(charge_density)

#        self.b = df_m.assemble(self.L)


        # Solve for the potential
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name
        poissonsolve.solve_for_phi(plot_flag=plotFlag, plot_title=plotTitle)

        # Write the potential to a file in VTK and XML formats
        file = df_m.File('phi_test_5_1D.pvd')
        # phi_name goes into the output file; phi_label doesn't
        phi_F.function.rename("phi1D", "phi_label")
        file << phi_F.function
        file = df_m.File('phi_test_5_1D.xml')
        file << phi_F.function

        # Write -E to a file in VTK and XML formats
        if negElectricField_F is not None:
            negElectricField_F.function.rename("E1D", "E_label")
            # !!! VTK cannot write a vector field on a 1D grid.
#            file = df_m.File("negE_test_5_1D.pvd")
#            file << negElectricField_F.function
            file = df_m.File("negE_test_5_1D.xml")
            file << negElectricField_F.function

        return
#    def test_5_1D_poisson_solver(self):ENDDEF


#class TestFieldSolve(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
