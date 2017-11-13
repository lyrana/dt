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

from DT_Module import DTcontrol_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import *

# Here's the mesh definition for this test
from UserMesh_y_Fields_FE_XYZ_Module import *

#STARTCLASS
class TestParticleMigration(unittest.TestCase):
    """Test of tracking neutral particles that drift across the mesh.

       Meshes on 1D, 2D and 3D rectangles are used.

    """
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pin = ParticleInput_C()

        # Set up particle variables
        pin.precision = np_M.float64
        pin.particle_integration_loop = 'loop-on-particles'
        # Use the 3 position coordinates, since we're doing 1, 2, and 3D particle motion
        pin.position_coordinates = ['x', 'y', 'z'] # Determines particle storage dimension
        # Neutral particles: No forces.
        """
        pin.force_components = ['x', 'y',]
        """
        pin.force_precision = np_M.float64

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particle_module below) giving the distribution
        # functions, and can vary from particle to particle.

#         pin.particle_species = (('neutral_H',
#                              {'initial_distribution_type' : 'listed',
#                               'charge' : 0.0,
#                               'mass' : 1.0*MyPlasmaUnits_C.AMU,
#                               'dynamics' : 'explicit',
# #                              'number_per_cell' : 12,
#                               }
#                              ),
#                             )

        speciesName = 'neutral_H'
        charge = 0.0
        mass = 1.0*MyPlasmaUnits_C.AMU
        dynamics = 'explicit'
        neutralH_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pin.particle_species = (neutralH_S,
                                 )
        # Make the particle object from pin
        self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary-condition callbacks, source regions, etc.)
        # Particles have a 3D/3D phase-space even though mesh is just 2D
        userParticleModule = "UserParticles_3D"
#        print "setUp: UserParticle module is", userParticleModule

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)

        self.particle_P.user_particle_module = userParticleModule
        self.particle_P.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### neutral H atoms are present at t=0
        speciesName = 'neutral_H'
        # Check that this species has been defined above
        if speciesName not in self.particle_P.species_names:
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
        neutralHParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                       }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {
                                'initial_neutral_H': (neutralHParams,),
                                }

        self.particle_P.initial_particles_dict = initialParticlesDict

        # Create the initial particles
        for ip in initialParticlesDict:
            ipList = initialParticlesDict[ip]
            ipParams = ipList[0]
            s = ipParams['species_name']
            initialDistributionType = ipParams['initial_distribution_type']
            if initialDistributionType == 'listed':
                # Put user-listed particles into the storage array
                self.particle_P.create_from_list(s, False)

        plotFlag = False
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name

        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create 1D, 2D and 3D meshes that the particles can be tested against

        # 1d mesh input
        umi1D = UserMeshInput_C()
        umi1D.pmin = df_M.Point(-10.0)
        umi1D.pmax = df_M.Point(10.0)
        umi1D.cells_on_side = (4,)

        # Identify where particle boundary-conditions will be imposed
        xminIndx = 1
        xmaxIndx = 2
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                }

        umi1D.particle_boundary_dict = particleBoundaryDict

        # Create a 1D particle mesh
        self.pmesh1DCI = UserMesh_C(umi1D, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 1D")
#        self.pmesh1DCI.compute_cell_vertex_dict()
#        self.pmesh1DCI.compute_cell_dict()

        # 2D mesh input
        umi2D = UserMeshInput_C()
        umi2D.pmin = df_M.Point(-10.0, -10.0)
        umi2D.pmax = df_M.Point(10.0, 10.0)
        umi2D.cells_on_side = (4, 2)

        # Identify where particle boundary-conditions will be imposed
        xminIndx = 1
        xmaxIndx = 2
        yminIndx = 4
        ymaxIndx = 8
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                'ymin': yminIndx,
                                'ymax': ymaxIndx,
                                }

        umi2D.particle_boundary_dict = particleBoundaryDict

        # Create a 2D particle mesh
        self.pmesh2DCI = UserMesh_C(umi2D, compute_dictionaries=True, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 2D")
#        self.pmesh2DCI.compute_cell_vertex_dict()
#        self.pmesh2DCI.compute_cell_dict()

        # 3D mesh input
        umi3D = UserMeshInput_C()
        umi3D.pmin = df_M.Point(-10.0, -10.0, -10.0)
        umi3D.pmax = df_M.Point(10.0, 10.0, 10.0)
        umi3D.cells_on_side = (4, 4, 4)

        # Identify where particle boundary-conditions will be imposed
        xminIndx = 1
        xmaxIndx = 2
        yminIndx = 4
        ymaxIndx = 8
        zminIndx = 16
        zmaxIndx = 32
        particleBoundaryDict = {'xmin': xminIndx,
                                'xmax': xmaxIndx,
                                'ymin': yminIndx,
                                'ymax': ymaxIndx,
                                'zmin': zminIndx,
                                'zmax': zmaxIndx,
                                }

        umi3D.particle_boundary_dict = particleBoundaryDict


        # Create a 3D particle mesh
        self.pmesh3D_M = UserMesh_C(umi3D, compute_tree=True, plot_flag=plotFlag, plot_title=plotTitle + ": 3D")
        # Explicitly compute dictionaries needed
        self.pmesh3D_M.compute_cell_entity_index_dict('vertex')
        self.pmesh3D_M.compute_cell_entity_index_dict('facet')
        self.pmesh3D_M.compute_cell_dict()
        self.pmesh3D_M.compute_cell_facet_normals_dict()
        self.pmesh3D_M.compute_cell_neighbor_dict()

        # pmesh is the owner of the compute_index function?

        return
#    def setUp(self):ENDDEF

    def test_1D_particle_migration(self):
        """Test particle migration across cells on a 1D mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        # List all the possible spatial coordinates
#        spatial_coordinates = ('x','y','z')

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.n_timesteps = 19

        # Run for 1D mesh

        self.particle_P.pmesh_M = self.pmesh1DCI

        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()


        ### Put the expected ending results into the p_expected tuple ###

        # First particle

        xsp0 = -9.5; ysp0 =  -9.5; zsp0 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0=2.0
        bitflag0 = 2
        cell_index0 = 0

        psp0 = (xsp0,ysp0,zsp0, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0)

        # Second particle

        xsp1 = -9.5; ysp1 =  -9.5; zsp1 = -9.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1)

        p_expected = (psp0, psp1)

        # Integrate for n_timesteps
        print "Moving", self.particle_P.get_total_particle_count(), "particles for", ctrlCI.n_timesteps, "timesteps"
        for istep in xrange(ctrlCI.n_timesteps):
            self.particle_P.move_neutral_particles(ctrlCI.dt)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particle_P.pseg_arr[sp].get(ip)
#                print 'expected = ', p_expected[ip]
#                print 'calculated = ', getparticle
                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                # The cell index is in last position: [-1]
#                print fncname, "expected cell =", p_expected[ip][-1], "computed cell =", getparticle[-1]
                self.assertEqual(p_expected[ip][-1], getparticle[-1], msg="Particle is not in correct cell")

        return
#    def test_1D_particle_migration(self):ENDDEF

    def test_2D_particle_migration(self):
        """Test particle migration across cells on 2D mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.n_timesteps = 19

        # Run for 2D mesh

        self.particle_P.pmesh_M = self.pmesh2DCI

        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particle_P.neutral_species[0]
        for ip in [0, 1]:
            p = self.particle_P.pseg_arr[sp].get(ip).copy() # Have to make a copy! Otherwise you overwrite the only copy of the particle
#            p_ic.append(self.particle_P.pseg_arr[sp].get(ip)) # Don't do this: it's a reference.
            p_ic.append(p)
#            print 'ip =', ip, 'p_ic =', p_ic[ip]


        # The expected final position and cell index

        # First particle

        xsp0 = -9.5; ysp0 =  -9.5; zsp0 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0=2.0
        bitflag0 = 2
        cell_index0 = 0

        psp0 = (xsp0,ysp0,zsp0, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0)

        # Second particle

        xsp1 = -9.5; ysp1 =  -9.5; zsp1 = -9.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1)

        p_expected = (psp0, psp1)

        # Integrate for n_timesteps
        print "Moving", self.particle_P.get_total_particle_count(), "particles for", ctrlCI.n_timesteps, "timesteps"
        for istep in xrange(ctrlCI.n_timesteps):
            self.particle_P.move_neutral_particles(ctrlCI.dt)

        # Create a mesh plotter to display the trajectory
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": First & last positions"
        plotter=df_M.plot(self.particle_P.pmesh_M.mesh, title=plotTitle)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particle_P.pseg_arr[sp].get(ip)
#                print 'sp =', sp, 'expected =', p_expected[ip]
#                print 'calculated = ', getparticle
                path = np_M.array([p_ic[ip][0], p_ic[ip][1], getparticle[0], getparticle[1]])
                plotter.add_polygon(path)

                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                # The cell index is in last position: [-1]
#                print fncname, "expected cell =", p_expected[ip][-1], "computed cell =", getparticle[-1]
                self.assertEqual(p_expected[ip][-1], getparticle[-1], msg="Particle is not in correct cell")

        plotter.plot()
#        df_M.interactive() # Stops the plot from disappearing

        return
#    def test_2D_particle_migration(self):ENDDEF

    def test_3D_particle_migration(self):
        """Test particle migration across cells on a 3D mesh.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.n_timesteps = 19

        # Run for 3D mesh

        self.particle_P.pmesh_M = self.pmesh3D_M

        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particle_P.neutral_species[0]
        for ip in [0, 1]:
            p = self.particle_P.pseg_arr[sp].get(ip).copy() # Have to make a copy!
            p_ic.append(p)

        # The expected final position and cell

        # First particle

        xsp0 = -9.5; ysp0 =  -9.5; zsp0 = 0.0
        vxsp0 = -2.0; vysp0 = 0.0; vzsp0 = 0.0

        weight0=2.0
        bitflag0 = 2
        cell_index0 = 98

        psp0 = (xsp0,ysp0,zsp0, vxsp0,vysp0,vzsp0, weight0, bitflag0, cell_index0)

        # Second particle

        xsp1 = -9.5; ysp1 =  -9.5; zsp1 = -9.5
        vxsp1 = -2.0; vysp1 = -2.0; vzsp1 = -2.0

        weight1 = 3.0
        bitflag1 = 2
        cell_index1 = 0

        psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1, cell_index1)

        p_expected = (psp0, psp1)

        # Integrate for n_timesteps
        print "Moving", self.particle_P.get_total_particle_count(), "particles for", ctrlCI.n_timesteps, "steps"
        for istep in xrange(ctrlCI.n_timesteps):
            self.particle_P.move_neutral_particles(ctrlCI.dt)

        # Create a mesh plotter to display the trajectory (just the
        # first and last positions)
        plotTitle = os.path.basename(__file__) + ": " + sys._getframe().f_code.co_name + ": First & last positions"
        plotter=df_M.plot(self.particle_P.pmesh_M.mesh, title=plotTitle)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        for sp in self.particle_P.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particle_P.pseg_arr[sp].get(ip)
#                print 'expected = ', p_expected[ip]
#                print 'calculated = ', getparticle
                path = np_M.array([p_ic[ip][0], p_ic[ip][1], p_ic[ip][2], getparticle[0], getparticle[1], getparticle[2]])
                plotter.add_polygon(path)

                for ic in range(ncoords):
                    self.assertAlmostEqual(p_expected[ip][ic], getparticle[ic], places=6, msg="Particle is not in correct position")
                # The index of the cell containing the particle is in the
                # last position of the stored particle attributes,
                # i.e., at location [-1]
#                print fncname, "expected cell =", p_expected[ip][-1], "computed cell =", getparticle[-1]
                self.assertEqual(p_expected[ip][-1], getparticle[-1], msg="Particle is not in correct cell")

        plotter.plot()
        df_M.interactive() # Stops the plot from disappearing

        return
#    def test_3D_particle_migration(self):ENDDEF

    '''
    def test_1_particle_migration(self):
        """
           Test particle migration to neighboring cell on 1D, 2D, 3D
           meshes.

        """

        return

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        # List all the possible spatial coordinates
        spatial_coordinates = ('x','y','z')

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.n_timesteps = 18

        # Run for 1, 2, and 3D meshes

        for dim in range(1,4):
#        for dim in range(2,3):   # 2D only
            spatial_components = spatial_coordinates[0:dim]
            if dim == 1:
                meshCI = self.mesh1DCI
            elif dim == 2:
                meshCI = self.mesh2DCI
            elif dim == 3:
                meshCI = self.mesh3DCI

            # Initialize or reinitialize the particles
            for sp in self.particle_P.species_names:
                if self.particle_P.initial_distribution_type[sp] == 'listed':
                    # Put user-listed particles into the storage array
                    self.particle_P.create_from_list(sp, print_flag=False, resetCounters=True)

            # Get the initial cell index of each particle.
            self.particle_P.compute_mesh_cell_indices(meshCI)


            # The expected results
        		
            # First species

            xsp1 = 0.2032648918; ysp1 =  0.0165297836; zsp1 = -0.1702053246
            vxsp1 = 1241.1798510396; vysp1 = -1517.6402979208; vzsp1 = -4276.4604468811

            weight1=1.0
            bitflag1 = 1

            psp1 = (xsp1,ysp1,zsp1, vxsp1,vysp1,vzsp1, weight1, bitflag1)

            # Second species

            xsp2 = 8.26835865540104E-005; ysp2 = 0.0001353672; zsp2 = 0.0001880508
            vxsp2 = 1.2578833919; vysp2 = 2.1157667838; vzsp2 = 2.9736501757

            weight2 = 3.0
            bitflag2 = 1

            psp2 = (xsp2,ysp2,zsp2, vxsp2,vysp2,vzsp2, weight2, bitflag2)

            p_expected = (psp1, psp2)

            # Integrate for n_timesteps
            for istep in xrange(ctrlCI.n_timesteps):
                self.particle_P.move_particles_without_fields(ctrlCI.dt, meshCI=meshCI)

                    #self.assertTrue(meshCI.is_inside(pseg[ip], pseg[ip]['cell_index']), msg = "The computed cell does not contain the particle")

        return
#    def test_1_particle_migration(self):ENDDEF
    '''

#class TestParticleMigration(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
