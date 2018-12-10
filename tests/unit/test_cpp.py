#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2018 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import unittest

import dolfin as df_m
import matplotlib.pyplot as mplot_m

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C
from Particle_Module import *
from SegmentedArrayPair_Module import SegmentedArray_C

# Here's the user's mesh definition for test_2
from UserMesh_y_Fields_FE_XYZ_Module import *

from UserUnits_Module import MyPlasmaUnits_C

# Use the C++ functions in .so libraries
import pseg_cpp
import dolfin_cpp

class TestCPP(unittest.TestCase):
    """Test the C++ functions in pseg_cpp.so and dolfin_cpp.so"""
    
    def setUp(self):

        # initializations for each test go here...

        self.segment_length = 100
        self.delete_flag = 0b01 # the lowest bit is 1
        self.trajectory_flag = 0b10 # the second lowest bit is 1
        precision = np_m.float64

        # 1D
        pcoords1D = ['x', 'x0', 'ux',]

        pvars = pcoords1D
        pvars.append('weight')
        pvartypes = [precision for var in pvars]
        pvars.append('bitflags')
        pvartypes.append(np_m.int32)
        pvars.append('cell_index')
        pvartypes.append(np_m.int32) # The size determines how many local cells you can have.
        pvars.append('unique_ID')
        pvartypes.append(np_m.int32)
        pvars.append('crossings')
        pvartypes.append(np_m.int32)

        metadata = {'names' : pvars, 'formats': pvartypes}
        self.seg_array_obj1D = SegmentedArray_C(self.segment_length, metadata)

        self.pvars1D = pvars
        self.particle_dtype1D = {'names' : pvars, 'formats': pvartypes}

        # 3D
        pcoords3D = ['x','y','z', 'x0','y0','z0', 'ux','uy','uz']

        pvars = pcoords3D
        pvars.append('weight')
        pvartypes = [precision for var in pvars]
        pvars.append('bitflags')
        pvartypes.append(np_m.int32)
        pvars.append('cell_index')
        pvartypes.append(np_m.int32) # The size determines how many local cells you can have.
        pvars.append('unique_ID')
        pvartypes.append(np_m.int32)
        pvars.append('crossings')
        pvartypes.append(np_m.int32)

        metadata = {'names' : pvars, 'formats': pvartypes}
        self.seg_array_obj3D = SegmentedArray_C(self.segment_length, metadata)

        self.pvars3D = pvars
        self.particle_dtype3D = {'names' : pvars, 'formats': pvartypes}

        return

#class TestCPP(unittest.TestCase):
    def test_1_print_pseg(self):
        """Check the print_pseg function in pseg_cpp.so.

           Make two 1D particle tuples and put them into a SegmentedArray_C
           object. Pass the first segment to pseg_cpp and check that
           the printout from print_pseg1D() is correct.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        particle_dimension = 1
        x=1.5; x0=1.0; ux=3.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL

        unique_ID = 7
        crossings = 5

        # Put two particles into the SegmentedArray_C
        
        pseg_arr = self.seg_array_obj1D
        putparticle = (x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings)
        pseg_arr.put(putparticle)
        putparticle = (x+0.5, x0+0.5, ux+0.5, weight, bitflags, cell_index, unique_ID+1, crossings)
        pseg_arr.put(putparticle)

        # Get the first segment and print it
        (npSeg, pseg) = pseg_arr.init_out_loop()        
        #print("npSeg=", npSeg)
        if particle_dimension == 1:
            print("pseg_cpp returns", pseg_cpp.print_pseg1D(pseg))
        
        return
#    def test_1_print_pseg:ENDDEF

#class TestCPP(unittest.TestCase):
    def test_2_particle_densities_on_2D_mesh(self):
        """Compute the number-density generated by line particles on a 2D mesh,
           using a C++ function.

           Macroparticles are created within a 2D meshed region and are weighted
           to nodal points on the mesh.
           
           No species data is defined for the particles.  
           No segmented-array particle storage is used for the particles, just
           one np_m array is used, with a particle dtype.

           See test_ChargeDensity.ods:test_2 for calculated values of the
           density.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')
        
        # Describe a 2D mesh from (-10,-10) to (10,10) with 2 cells on a side.
        umi2d_I = UserMeshInput_C()
        umi2d_I.pmin = df_m.Point(-10.0, -10.0)
        umi2d_I.pmax = df_m.Point(10.0, 10.0)
        umi2d_I.cells_on_side = (2, 2)
        umi2d_I.diagonal = 'left'

        # UserMesh_FE_XYZ_Module can make the mesh from the above input.
        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"

        pmesh2d_M = UserMesh_C(umi2d_I, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle)
#        pmesh2d_M.compute_cell_vertex_dict()
#        pmesh2d_M.compute_cell_dict()

        # Put 3 particles inside the meshed region

        x0 = -5.0; y0 = -5.0; z0 = 0.0
        ux0 = 0.0; uy0 = 0.0; uz0 = 0.0
        weight0 = 2.0e10 # number of electrons per macroparticle
        bitflags0 = 0b0
        cell_index0 = 1
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p0 = (x0,y0,z0, ux0,uy0,uz0, weight0, bitflags0, cell_index0, unique_ID, crossings)

        # 2nd particle
        x1 = 1.0; y1 = 1.0; z1 = 1.0
        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 3.0e10
        bitflags1 = 0b0
        cell_index1 = 6
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p1 = (x1,y1,z1, ux1,uy1,uz1, weight1, bitflags1, cell_index1, unique_ID, crossings)

        # 3nd particle
        x2 = -9.0; y2 = 1.0; z2 = 1.0
        ux2 = uy2 = 0.0; uz2 = -uz0
        weight2 = 4.0e10
        bitflags2 = 0b0
        cell_index2 = 4 # Particle lies on boundary between 0 and 1
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p2 = (x2,y2,z2, ux2,uy2,uz2, weight2, bitflags2, cell_index2, unique_ID, crossings)

        # Create the DT particle record type
        pvars = ['x', 'y', 'z', 'ux', 'uy', 'uz', 'weight', 'bitflags', 'cell_index', 'unique_ID', 'crossings']
        pvartypes = [np_m.float64]*7
        pvartypes.append(np_m.int32) # bitflags
        pvartypes.append(np_m.int32) # cell_index
        pvartypes.append(np_m.int32) # unique_ID
        pvartypes.append(np_m.int32) # crossings

        p_dtype = {'names' : pvars, 'formats': pvartypes}

        # Put the particles into an ndarray with the above type
        nparticles = 3
        particles = np_m.empty(nparticles, dtype=p_dtype)
        particles[0] = p0
        particles[1] = p1
        particles[2] = p2

        ###### DoF density vector
        
        # Allocate storage for the number-density values.  The number-density array
        # stores the integral of the physical density-distribution times the element
        # basis functions.

        dofNumberDensityElementType = 'Lagrange'
        dofNumberDensityElementDegree = 1
        dofNumberDensityFieldType = 'scalar'
        dofNumberDensity_F = Field_C(pmesh2d_M,
                                  element_type=dofNumberDensityElementType,
                                  element_degree=dofNumberDensityElementDegree,
                                  field_type=dofNumberDensityFieldType)

#        print "size of dofNumberDensity:", dofNumberDensity_F.function.vector().size()

        # The expected number-density values from test_ChargeDensity.ods:test_2
        dofNumberDensityExpected = np_m.empty(dofNumberDensity_F.function.vector().size(), dtype=np_m.float64)

        # Note that these are in the vertex numbering order, not DoF order. We will
        # need to compare dofNumberDensityExpected[ivert] with dofNumberDensityCalc[idof],
        # where idof is the DoF number corresponding to ivert (see conversion below).
        dofNumberDensityExpected[0] = 0.0
        dofNumberDensityExpected[1] = 1.0e10
        dofNumberDensityExpected[2] = 0.0
        dofNumberDensityExpected[3] = 4.2e10
        dofNumberDensityExpected[4] = 2.8e10
        dofNumberDensityExpected[5] = 3.0e9
        dofNumberDensityExpected[6] = 4.0e9
        dofNumberDensityExpected[7] = 3.0e9
        dofNumberDensityExpected[8] = 0.0

        # Compute the DoF density vector using the DnT functions we're testing here.
        for p in particles:
            dofNumberDensity_F.interpolate_delta_function_to_dofs(p)

        # Get a np_m array of the calculated values
        dofNumberDensityCalc = dofNumberDensity_F.function.vector().get_local()
#        print fncName, dofNumberDensityCalc

        functionSpace = dofNumberDensity_F.function_space
        gdim = dofNumberDensity_F.mesh_gdim
        # Reshape the DoF coordinates to get (x), or (x,y), or (x,y,z) tuples.
        if df_m.__version__ > "1.5.0":
            dofcoords = functionSpace.tabulate_dof_coordinates().reshape((-1, gdim))
        else:
            print('\n!!!WARNING!!!: ', fncName, ": DOLFIN too old.  Skipping rest of test")
            return

#        print "dofcoords=", dofcoords

        # Plot the result

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": number density"
        df_m.plot(dofNumberDensity_F.function, title=plotTitle)
        mplot_m.show()
#        yesno = raw_input("Just called show() in test_2_interpolate_particle_density_to_2D_mesh")

        ## Check the values vs. those computed in test_ChargeDensity.ods:test_2

        # Convert vertex indices to DoF indices.  For CG1 elements, there are
        # the same number of vertices as DoFs.
        v2d=df_m.vertex_to_dof_map(functionSpace)

        for ivert in range(len(dofNumberDensityExpected)):
            idof = v2d[ivert] # The DoF number for this vertex.
#            print "vertex", ivert, "is DoF index", idof
#            print "dofNumberDensityCalc, dofNumberDensityExpected =", dofNumberDensityCalc[idof], dofNumberDensityExpected[ivert]
            # "places" means numbers after the decimal point:
            self.assertAlmostEqual(dofNumberDensityCalc[idof], dofNumberDensityExpected[ivert], places=4, msg="Wrong value of dofNumberDensity")


        ###### Cell density vector

        # Allocate storage for the cell number-density values.

        cellNumberDensityElementType = 'DG'
        cellNumberDensityElementDegree = 0
        cellNumberDensityFieldType = 'scalar'
        cellNumberDensity_F = Field_C(pmesh2d_M,
                                      element_type=cellNumberDensityElementType,
                                      element_degree=cellNumberDensityElementDegree,
                                      field_type=cellNumberDensityFieldType)

#        print "size of cellNumberDensity:", cellNumberDensity_F.function.vector().size()

        # Compute the cell density using the DnT functions that we're testing
        for p in particles:
            cellNumberDensity_F.add_weight_to_cell(p)
        cellNumberDensity_F.divide_by_cell_volumes()

        # Get an array of the density values
        cellNumberDensityCalc = cellNumberDensity_F.function.vector().get_local()
#        print fncName, "cellNumberDensityCalc =", cellNumberDensityCalc

        # The expected number-density values from test_ChargeDensity.ods:test_2
        cellNumberDensityExpected = np_m.empty(cellNumberDensity_F.function.vector().size(), dtype=np_m.float64)

        # Note that these are in the cell numbering order, not DoF order. We will
        # need to compare cellNumberDensityExpected[cellIndex] with
        # cellNumberDensityCalc[dofIndex], where dofIndex is the DoF number
        # corresponding to cellIndex.
        cellNumberDensityExpected[0] = 0.0
        cellNumberDensityExpected[1] = 4.0e8
        cellNumberDensityExpected[2] = 0.0
        cellNumberDensityExpected[3] = 0.0
        cellNumberDensityExpected[4] = 8.0e8
        cellNumberDensityExpected[5] = 0.0
        cellNumberDensityExpected[6] = 6.0e8
        cellNumberDensityExpected[7] = 0.0

        # Compare results

        for cell in df_m.cells(pmesh2d_M.mesh):
            cellIndex = cell.index()
#            cellVol = mesh_M.cell_volume_dict[cellIndex]
            # There's only 1 DoF in the cell, since this function uses constant DG elements
            dofIndex = cellNumberDensity_F.function_space.dofmap().cell_dofs(cellIndex) # return type: np_m.ndarray
            self.assertAlmostEqual(cellNumberDensityCalc[dofIndex], cellNumberDensityExpected[cellIndex], places=4, msg="Wrong value of cellNumberDensity")

        return
#    def test_2_particle_densities_on_2D_mesh(self):ENDDEF

#class TestCPP(unittest.TestCase):
    def test_3_particle_densities_on_2D_mesh(self):
        """Compute the number-densities generated by line particles on a 2D mesh, using the
           C++ functions interpolate_weights_to_dofs3D() and add_weights_to_cells3D().

           This is a copy of
           test_ChargeDensity.py:test_2_interpolate_particle_density_to_2D_mesh(self),
           using pseg_cpp.so instead of pure Python.

           Macroparticles with 3D coordinates are created within a 2D meshed region
           and are weighted to nodal points on the mesh to get the nodal
           particle-density (times a nodal volume).  To get the cell
           particle-density, the macroparticle weights are summed in each cell and
           then divided by the cell volume (area in this case).
           
           See test_ChargeDensity.ods:test_2 for calculated values of the
           density.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')
        
        # Describe a 2D mesh from (-10,-10) to (10,10) with 2 cells on a side.
        # The mesh is triangular, so there's a total of 4x2 = 8 cells.
        umi2d_I = UserMeshInput_C()
        umi2d_I.pmin = df_m.Point(-10.0, -10.0)
        umi2d_I.pmax = df_m.Point(10.0, 10.0)
        umi2d_I.cells_on_side = (2, 2)
        umi2d_I.diagonal = 'left'

        # UserMesh_FE_XYZ_Module can make the mesh from the above input.
        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": XY mesh"

        pmesh2d_M = UserMesh_C(umi2d_I, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle)
#        pmesh2d_M.compute_cell_vertex_dict()
#        pmesh2d_M.compute_cell_dict()

        ########## Kinetic particles ##########

        # Create an instance of the DTparticleInput class
        pin = self.pin = ParticleInput_C()
        # Set up particle variables
        pin.precision = np_m.float64

        pin.particle_integration_loop = 'loop-on-particles'
        pin.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pin.force_components = ['x', 'y', 'z']
        pin.force_precision = np_m.float64

        ### Particle species input

        # Give the properties of the particle species.  The charges and masses
        # are normally those of the physical particles, and not the
        # computational macroparticles.

        speciesName = 'plasma_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'implicit'
        plasmaElectron_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # speciesName = 'H_plus'
        # charge = 2.0*MyPlasmaUnits_C.elem_charge
        # mass = 1.0*MyPlasmaUnits_C.AMU
        # dynamics = 'implicit'
        # Hplus_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these two species to particle input
#        pin.particle_species = (plasmaElectron_S, Hplus_S,)
        pin.particle_species = (plasmaElectron_S,)

        # Make the particle object from pin...
        particles_P = Particle_C(pin, print_flag=False)

        # ...and attach the particle mesh
        particles_P.pmesh_M = pmesh2d_M

        # Put 3 particles inside the meshed region

        x0 = -5.0; y0 = -5.0; z0 = 0.0
        ux0 = 0.0; uy0 = 0.0; uz0 = 0.0
        weight0 = 2.0e10 # number of electrons per macroparticle
        bitflags0 = 0b0
        cell_index0 = 1
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p0 = (x0,y0,z0, x0,y0,z0, ux0,uy0,uz0, weight0, bitflags0, cell_index0, unique_ID, crossings)

        # 2nd particle
        x1 = 1.0; y1 = 1.0; z1 = 1.0
        ux1 = uy1 = 0.0; uz1 = -uz0
        weight1 = 3.0e10
        bitflags1 = 0b0
        cell_index1 = 6
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p1 = (x1,y1,z1, x1,y1,z1, ux1,uy1,uz1, weight1, bitflags1, cell_index1, unique_ID, crossings)

        # 3nd particle
        x2 = -9.0; y2 = 1.0; z2 = 1.0
        ux2 = uy2 = 0.0; uz2 = -uz0
        weight2 = 4.0e10
        bitflags2 = 0b0
        cell_index2 = 4 # Particle lies on boundary between 0 and 1
        unique_ID = Particle_C.UNIQUE_ID_COUNTER; Particle_C.UNIQUE_ID_COUNTER += 1
        crossings = 0

        p2 = (x2,y2,z2, x2,y2,z2, ux2,uy2,uz2, weight2, bitflags2, cell_index2, unique_ID, crossings)

        # Put the particles into an ndarray with the above type
        nparticles = 3
        particles = (p0, p1, p2)

        ### Put these particles into storage

        # Electrons
        species_name = 'plasma_electrons'

        number_of_macroparticles = len(particles)

        pseg_arr = particles_P.pseg_arr[species_name] # The SegmentedArray_C object for this species

        for i in range(number_of_macroparticles):
#            print ('species_name, particles[i] = ', species_name, particles[i])
            p, pindex = pseg_arr.put(particles[i])

        # Ions
#         species_name = 'H_plus'

#         number_of_macroparticles = 3

#         pseg_arr = particles_P.pseg_arr[species_name] # The SegmentedArray_C object for this species

#         for i in range(number_of_macroparticles):
# #            print 'species_name, particles[i] = ', species_name, particles[i]
#             p, pindex = pseg_arr.put(particles[i])

        ###### DoF density vector
        
        # Allocate storage for the number-density values.  The number-density array
        # stores the integral of the physical density-distribution times the element
        # basis functions.

        dofNumberDensityElementType = 'Lagrange'
        dofNumberDensityElementDegree = 1
        dofNumberDensityFieldType = 'scalar'
        dofNumberDensity_F = Field_C(pmesh2d_M,
                                  element_type=dofNumberDensityElementType,
                                  element_degree=dofNumberDensityElementDegree,
                                  field_type=dofNumberDensityFieldType)

#        print "size of dofNumberDensity:", dofNumberDensity_F.function.vector().size()

        # The expected number-density values from test_ChargeDensity.ods:test_2
        dofNumberDensityExpected = np_m.empty(dofNumberDensity_F.function.vector().size(), dtype=np_m.float64)

        # Note that these are in the vertex numbering order, not DoF order. We will
        # need to compare dofNumberDensityExpected[ivert] with dofNumberDensityCalc[idof],
        # where idof is the DoF number corresponding to ivert (see conversion below).
        dofNumberDensityExpected[0] = 0.0
        dofNumberDensityExpected[1] = 1.0e10
        dofNumberDensityExpected[2] = 0.0
        dofNumberDensityExpected[3] = 4.2e10
        dofNumberDensityExpected[4] = 2.8e10
        dofNumberDensityExpected[5] = 3.0e9
        dofNumberDensityExpected[6] = 4.0e9
        dofNumberDensityExpected[7] = 3.0e9
        dofNumberDensityExpected[8] = 0.0

        # Compute the DoF density vector using the C++ functions we're testing here.

        psa = particles_P.pseg_arr[species_name] # segmented array for this species
        (npSeg, pseg) = psa.init_out_loop()

        while isinstance(pseg, np_m.ndarray):

            # Call the (x, y, z) version:
            pseg_cpp.interpolate_weights_to_dofs3D(pseg, dofNumberDensity_F.function._cpp_object)
            (npSeg, pseg) = psa.get_next_segment('out')

        # Get a np_m array of the calculated values
        dofNumberDensityCalc = dofNumberDensity_F.function.vector().get_local()
#        print fncName, dofNumberDensityCalc

        functionSpace = dofNumberDensity_F.function_space
        gdim = dofNumberDensity_F.mesh_gdim
        # Reshape the coordinates to get (x), or (x,y), or (x,y,z) tuples.
        ndofcoords = functionSpace.tabulate_dof_coordinates().reshape((-1, gdim))

#        print "dofcoords=", dofcoords

        # Plot the result

        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": number density"
        df_m.plot(dofNumberDensity_F.function, title=plotTitle)
        mplot_m.show()
#        yesno = raw_input("Just called show() in test_2_interpolate_particle_density_to_2D_mesh")

        ## Check the values vs. those computed in test_ChargeDensity.ods:test_2

        # Convert vertex indices to DoF indices.  For CG1 elements, there are
        # the same number of vertices as DoFs.
        v2d=df_m.vertex_to_dof_map(functionSpace)

        for ivert in range(len(dofNumberDensityExpected)):
            idof = v2d[ivert] # The DoF number for this vertex.
#            print "vertex", ivert, "is DoF index", idof
#            print "dofNumberDensityCalc, dofNumberDensityExpected =", dofNumberDensityCalc[idof], dofNumberDensityExpected[ivert]
            # "places" means numbers after the decimal point:
            self.assertAlmostEqual(dofNumberDensityCalc[idof], dofNumberDensityExpected[ivert], places=4, msg="Wrong value of dofNumberDensity")


        ###### Cell density vector

        # Allocate storage for the cell number-density values.

        cellNumberDensityElementType = 'DG'
        cellNumberDensityElementDegree = 0
        cellNumberDensityFieldType = 'scalar'
        cellNumberDensity_F = Field_C(pmesh2d_M,
                                      element_type=cellNumberDensityElementType,
                                      element_degree=cellNumberDensityElementDegree,
                                      field_type=cellNumberDensityFieldType)

#        print "size of cellNumberDensity:", cellNumberDensity_F.function.vector().size()

        psa = particles_P.pseg_arr[species_name] # segmented array for this species
        (npSeg, pseg) = psa.init_out_loop()

        while isinstance(pseg, np_m.ndarray):

            # Call the 3D (x, y, z) version:
            pseg_cpp.add_weights_to_cells3D(pseg, cellNumberDensity_F.function._cpp_object)
            
            (npSeg, pseg) = psa.get_next_segment('out')

        # dolfin_cpp.print_dict(particles_P.pmesh_M.cell_volume_dict)

        # Convert the cell values to a cell density in C++
        if cellNumberDensity_F is not None:
            dolfin_cpp.divide_by_cell_volumes(cellNumberDensity_F.function._cpp_object, particles_P.pmesh_M.cell_volume_dict)
            
        # Get an array of the density values
        cellNumberDensityCalc = cellNumberDensity_F.function.vector().get_local()
        #print(fncName, "cellNumberDensityCalc =", cellNumberDensityCalc)

        # Make a numpy array to hole the expected number-density values from
        # test_ChargeDensity.ods:test_2
        cellNumberDensityExpected = np_m.empty(cellNumberDensity_F.function.vector().size(), dtype=np_m.float64)

        # Note that these are in the cell numbering order, not DoF order. We will
        # need to compare cellNumberDensityExpected[cellIndex] with
        # cellNumberDensityCalc[dofIndex], where dofIndex is the DoF number
        # corresponding to cellIndex.
        cellNumberDensityExpected[0] = 0.0
        cellNumberDensityExpected[1] = 4.0e8
        cellNumberDensityExpected[2] = 0.0
        cellNumberDensityExpected[3] = 0.0
        cellNumberDensityExpected[4] = 8.0e8
        cellNumberDensityExpected[5] = 0.0
        cellNumberDensityExpected[6] = 6.0e8
        cellNumberDensityExpected[7] = 0.0

        # Compare results

        for cell in df_m.cells(pmesh2d_M.mesh):
            cellIndex = cell.index()
            # There's only 1 DoF in the cell, since this function uses constant DG elements
            dofIndex = cellNumberDensity_F.function_space.dofmap().cell_dofs(cellIndex) # return type: np_m.ndarray
#            print("cellNumberDensityCalc[dofIndex][0]=", cellNumberDensityCalc[dofIndex][0], "cellNumberDensityExpected[cellIndex]=",cellNumberDensityExpected[cellIndex]) 
            self.assertAlmostEqual(cellNumberDensityCalc[dofIndex][0], cellNumberDensityExpected[cellIndex], places=4, msg="Wrong value of cellNumberDensity")

        return
#    def test_3_particle_densities_on_2D_mesh(self):ENDDEF    


#class TestCPP(unittest.TestCase): ENDCLASS



if __name__ == '__main__':
    unittest.main()
