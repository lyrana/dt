#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_M
import importlib as im_M
import unittest

import dolfin as df_M

from UserUnits_Module import MyPlasmaUnits_C
from Dolfin_Module import Mesh_C
from Particle_Module import *

from UserMesh_FE_XYZ_Module import *

#STARTCLASS
class TestParticleCellIndex(unittest.TestCase):
    """Test cell functions Particle_Module"""
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = ParticleInput_C()

        # Set up particle variables
        pinCI.precision = np_M.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y', 'z'] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = np_M.float64

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particle_module below) giving the distribution
        # functions, and can vary from particle to particle.

#         pinCI.particle_species = (('plasmaelectrons',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 12,
#                               }
#                              ),
#                             ('Hplus', 
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : 1.0*MyPlasmaUnits_C.elem_charge,
#                               'mass' : 1.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 6,
#                               }
#                              ),
#                             ('He', 
#                              {'initial_distribution_type' : None,
#                               'charge' : 0.0,
#                               'mass' : 4.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'implicit',
# #                              'number_per_cell' : 1,
#                               }
#                              ),
#                             )

        speciesName = 'plasma_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        plasmaElectronsCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'H_plus'
        charge = 1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.proton_mass
        dynamics = 'explicit'
        HPlusCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        speciesName = 'He'
        charge = 0.0
        mass = 4.0*MyPlasmaUnits_C.AMU
        dynamics = 'explicit'
        HeCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pinCI.particle_species = (plasmaElectronsCI, HPlusCI, HeCI,
                                 )
        # Make the particle object from pinCI
        self.particleCI = Particle_C(pinCI, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_H_He_e"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        self.particleCI.user_particle_module = userParticleModule
        self.particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### plasma_electrons are present at t=0

        ## plasma electrons are present at t=0
        # Name the initialized species (it should be in species_names above)
        speciesName = 'plasma_electrons'
        # Check that this species has been defined above
        if speciesName not in self.particleCI.species_names:
            print fncName + "The species", speciesName, "has not been defined"
            sys.exit()

        # Specify how the species will be initialized
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticleClass, speciesName):
            if printFlag: print fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticleClass
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in %s for species %s " % (speciesName, userParticleModule, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        plasmaElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }


        ### H+ ions are present at t=0
        speciesName = 'H_plus'
        # Check that this species has been defined above
        if speciesName not in self.particleCI.species_names:
            print fncName + "The species", speciesName, "has not been defined"
            sys.exit()

        # Specify how the species will be initialized
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticleClass, speciesName):
            if printFlag: print fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticleClass
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in UserParticle.py for species %s " % (speciesName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        HPlusParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                       }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_plasma_electrons': (plasmaElectronParams,),
                                'initial_H_plus': (HPlusParams,),
                                }

        self.particleCI.initial_particles_dict = initialParticlesDict

        # Create the initial particles
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']
            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(s, False)

        plotFlag = False
        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create meshes that the particles can be tested against

        # 1d mesh input
        umi1DCI = UserMeshInput_C()
        umi1DCI.pmin = df_M.Point(-10.0)
        umi1DCI.pmax = df_M.Point(10.0)
        umi1DCI.cells_on_side = (4,) # Need the comma to indicate a tuple
        # Create mesh
        self.pmesh1DCI = UserMesh_C(umi1DCI, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag)

        # 2D mesh input
        umi2DCI = UserMeshInput_C()
        umi2DCI.pmin = df_M.Point(-0.03, -0.03)
        umi2DCI.pmax = df_M.Point(0.03, 0.03)
        umi2DCI.cells_on_side = (4, 4)
        # Create mesh
        self.pmesh2DCI = UserMesh_C(umi2DCI, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag)
#        self.mesh2DCI.compute_cell_vertex_dict()
#        self.mesh2DCI.compute_cell_dict()

        # 3D mesh input
        umi3DCI = UserMeshInput_C()
        umi3DCI.pmin = df_M.Point(-0.03, -0.03, -0.03)
        umi3DCI.pmax = df_M.Point(0.03, 0.03, 0.03)
        umi3DCI.cells_on_side = (4, 4, 4)
        # Create mesh
        self.pmesh3DCI = UserMesh_C(umi3DCI, compute_tree=True, plot_flag=plotFlag)
        self.pmesh3DCI.compute_cell_entity_index_dict('vertex')
        self.pmesh3DCI.compute_cell_dict()

        # pmesh is the owner of the compute_index function?

        return
#    def setUp(self):END


    def test_1_cell_index(self):
        """Test the mesh's vertex-to-cell dict and the
           bounding-box-tree search for a point.

           Check the mesh's cell-index list for some test vertex
           indices in 2D and 3D.  The vertex-cell dictionary of a mesh
           gives the list of cells that share a given vertex.

           Check that a search for the index of a cell, using the
           mesh's bounding-box tree and the coordinates of the cell's
           midpoint, returns the correct cell index. Tests 2D and 3D
           mesh.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest:', fncName, '('+__file__+')'

        # 1D mesh

        # Test vertex dictionary for vertex 3
#        print 'vertex 3 is in cells', self.pmesh1DCI.vertex_cell_dict[3]
        self.assertEqual([2, 3], self.pmesh1DCI.vertex_cell_dict[3], msg = "cell list for vertex 3 is not correct")

        # Compute a cell index by locating the cell that contains
        # the cell midpoint using the BB tree.  Compare that to the
        # index value stored in the cell object.
        # Make a particle-like point from the midpoint
        print "1D test"
        for cell in df_M.cells(self.pmesh1DCI.mesh):
            midpt = cell.midpoint()
            # Have to convert the Point to a 1-element pseg-type array
            # in order to call compute_cell_index()
            p1_dtype = {'names' : ['x'], 'formats': [np_M.float64]}
            p1 = np_M.array([(midpt.x(),)], dtype=p1_dtype)

            cell_index = self.pmesh1DCI.compute_cell_index(p1[0])
#            print "cell midpoint:", midpt.x()
#            print "cell index:", cell.index(), "computed index:", cell_index
            self.assertEqual(cell_index, cell.index(), msg = "1D cell index is not correct")

        # 2D mesh

        # Test vertex dictionary for vertex 8 and vertex 20
#        print 'vertex 20 is in cells', self.pmesh2DCI.vertex_cell_dict[20]
        self.assertEqual([25], self.pmesh2DCI.vertex_cell_dict[20], msg = "cell list for vertex 20 is not correct")
#        print 'vertex 8 is in cells', self.pmesh2DCI.vertex_cell_dict[8]
        self.assertEqual([4, 5, 7, 12, 14, 15], self.pmesh2DCI.vertex_cell_dict[8], msg = "cell list for vertex 8 is not correct")

        # Compute a cell index by locating the cell that contains
        # the cell midpoint using the BB tree.  Compare that to the
        # index value stored in the cell object.
        # Make a particle-like point from the midpoint
        print "2D test"
        for cell in df_M.cells(self.pmesh2DCI.mesh):
            midpt = cell.midpoint()
            # Have to convert the Point to a 1-element pseg-type array
            # in order to call compute_cell_index()
            p2_dtype = {'names' : ['x', 'y'], 'formats': [np_M.float64]*2}
            p2 = np_M.array([(midpt.x(), midpt.y())], dtype=p2_dtype)

            cell_index = self.pmesh2DCI.compute_cell_index(p2[0])
#            print "cell midpoint:", midpt.x(), midpt.y()
#            print "cell index:", cell.index(), "computed index:", cell_index
            self.assertEqual(cell_index, cell.index(), msg = "2D cell index is not correct")

        # 3D mesh
        print "3D test"
        # Make a particle-like point from the midpoint
        for cell in df_M.cells(self.pmesh3DCI.mesh):
            midpt = cell.midpoint()
            p3_dtype = {'names' : ['x', 'y', 'z'], 'formats': [np_M.float64]*3}
            p3 = np_M.array([(midpt.x(), midpt.y(), midpt.z())], dtype=p3_dtype)

            cell_index = self.pmesh3DCI.compute_cell_index(p3[0])
#            print "cell midpoint:", midpt.x(), midpt.y(), midpt.z()
#            print "cell index:", cell.index(), "computed index:", cell_index
            self.assertEqual(cell_index, cell.index(), msg = "3D cell index is not correct")

        return
#    def test_1_cell_index(self):END


    def test_2_cell_index(self):
        """Test the function compute_cell_index().

           Check that the particles stored in a Particle_C object are
           in fact inside the cell computed by the function
           compute_cell_index(). Uses the function is_inside() to do
           the test.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest:', fncName, '('+__file__+')'

        # List all the possible spatial coordinates
        spatial_coordinates = ('x','y','z')

        print "1D, 2D, 3D tests"
        # Loop on the species
        isp = 0
        for sp in self.particleCI.species_names:
            if self.particleCI.get_species_particle_count(sp) == 0: continue

            psaCI = self.particleCI.pseg_arr[sp] # segmented array for this species

            # Loop on the particles for this species
            (np_seg, pseg) = psaCI.init_out_loop()
#            pseg = psaCI.get_next_segment()

            while isinstance(pseg, np_M.ndarray):
                # Check if the particle is still on the meshed region
#                for ip in xrange(pseg.size):
                for ip in xrange(np_seg):

                    # pseg[i] is 'x', 'y', 'z', 'vx', 'vy',... values of ith item
                    # So pseg[i][0:3] is 'x', 'y', 'z'.
                    # Can't use slice syntax here, because the dtype is inhomogeneous.

                    # Compute and store the cell index
                    for dim in range(1,4):
                        if dim == 1:
                            pmeshCI = self.pmesh1DCI
#                        spatial_components = spatial_coordinates[0:dim]
#                        p = np_M.array([pseg[ip][comp] for comp in spatial_components])
                        if dim == 2:
                            pmeshCI = self.pmesh2DCI
                        elif dim == 3:
                            pmeshCI = self.pmesh3DCI
# Put a particle outside the mesh to see what cell index is returned:
#                            p[0] = p[0] + 100.0
#                            pseg[ip]['cell_index'] = self.pmesh2DCI.compute_cell_index(df_M.Point(p))
                        # Compute the cell index containing the particle, and save it.
                        pseg[ip]['cell_index'] = pmeshCI.compute_cell_index(pseg[ip])
                        # print "Coordinates", p, "are in cell", pseg[ip]['cell_index']
                        if pseg[ip]['cell_index'] != Mesh_C.NO_CELL:
#                            cell_vertices = pmeshCI.cell_entity_index_dict['vertex'][pseg[ip]['cell_index']]
                            # Look up the cell index in the particle datalist
                            c = pmeshCI.cell_dict[pseg[ip]['cell_index']]
                        else:
                            c = None

                        if c is not None:
#                            self.assertTrue(c.contains(df_M.Point(p)), msg = "The computed cell does not contain the particle")
                            # Verify that this cell does actually contain the particle.
                            self.assertTrue(pmeshCI.is_inside(pseg[ip], pseg[ip]['cell_index']), msg = "The computed cell does not contain the particle")
                        else:
                            self.assertTrue(False, msg = "A particle is outside the mesh")

                # Done with this segment.
                # Get the next one, if it exists.
                (np_seg, pseg) = psaCI.get_next_segment('out')

            # Move on to the next species
            isp += 1
        return
#    def test_2_cell_index(self):END

#class TestParticleCellIndex(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
