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
from DT_Module import DTparticleInput_C
from UserUnits_Module import MyPlasmaUnits_C
from Particle_Module import Particle_C

# Here's the mesh definition for this test
from UserMesh_FE_XYZ_Module import UserMesh_C

class DTmeshInput_C(object):
    """Input for the field mesh.  List the variables that can be
       modified by the user.
    """

    def __init__(self):
        """ List the mesh variables that the user can set in MAIN.py
        """
        self.pmin = None
        self.pmax = None
        self.cells_on_side = None

        self.field_boundary_dict = None
        self.particle_boundary_dict = None

        return

#class DTmeshInput_C(object): ENDCLASS


class TestParticleMigration(unittest.TestCase):
    """Test of tracking neutral particles that drift across the mesh.

       Meshes on 1D, 2D and 3D rectangles are used.

    """
    
    def setUp(self):

        # Initializations performed before each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()

        # Set up particle variables
        pinCI.precision = np_M.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        # Use the 3 position coordinates, since we're doing 1, 2, and 3D particle motion
        pinCI.position_coordinates = ['x', 'y', 'z'] # Determines particle storage dimension
        # Neutral particles: No forces.
        """
        pinCI.force_components = ['x', 'y',]
        """
        pinCI.force_precision = np_M.float64

        # Give the properties of the particle species.  The charges
        # and masses are normally those of the physical particles, and
        # not the computational macroparticles.  Macroparticle weights
        # are specified or computed in a separate file (see
        # user_particles_module below) giving the distribution
        # functions, and can vary from particle to particle.

        pinCI.particle_species = (('neutral_H',
                             {'initial_distribution_type' : 'listed',
                              'charge' : 0.0,
                              'mass' : 1.0*MyPlasmaUnits_C.AMU,
                              'dynamics' : 'explicit',
#                              'number_per_cell' : 12,
                              }
                             ),
                            )

        # Provide the name of the .py file containing the user-provided particle distributions
        pinCI.user_particles_module = "UserParticles_3D"
#        print "setUp: UserParticle module is", pinCI.user_particles_module
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_C = UPrt_M.UserParticleDistributions_C

        self.pinCI = pinCI

        # Make the particle object using pinCI
        self.particleCI = Particle_C(pinCI, printFlag=False)

        # Store the particles
        for sp in self.particleCI.species_names:
            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(sp, False)

        plotFlag = False
        # Turn off plotting if there's no DISPLAY

        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create 1D, 2D and 3D meshes that the particles can be tested against

        # 1d mesh input
        mi1DCI = DTmeshInput_C()
        mi1DCI.pmin = df_M.Point(-10.0)
        mi1DCI.pmax = df_M.Point(10.0)
        mi1DCI.cells_on_side = (4)
        # Create a 1D particle mesh
        self.pmesh1DCI = UserMesh_C(mi1DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)
#        self.pmesh1DCI.compute_cell_vertex_dict()
#        self.pmesh1DCI.compute_cell_dict()

        # 2D mesh input
        mi2DCI = DTmeshInput_C()
        mi2DCI.pmin = df_M.Point(-10.0, -10.0)
        mi2DCI.pmax = df_M.Point(10.0, 10.0)
        mi2DCI.cells_on_side = (4, 2)
        # Create a 2D particle mesh
        self.pmesh2DCI = UserMesh_C(mi2DCI, computeDictionaries=True, computeTree=True, plotFlag=plotFlag)
#        self.pmesh2DCI.compute_cell_vertex_dict()
#        self.pmesh2DCI.compute_cell_dict()

        # 3D mesh input
        mi3DCI = DTmeshInput_C()
        mi3DCI.pmin = df_M.Point(-10.0, -10.0, -10.0)
        mi3DCI.pmax = df_M.Point(10.0, 10.0, 10.0)
        mi3DCI.cells_on_side = (4, 4, 4)

        # Create a 3D particle mesh
        self.pmesh3DCI = UserMesh_C(mi3DCI, computeTree=True, plotFlag=plotFlag)
        # Explicitly compute dictionaries needed
        self.pmesh3DCI.compute_cell_entity_index_dict('vertex')
        self.pmesh3DCI.compute_cell_entity_index_dict('facet')
        self.pmesh3DCI.compute_cell_dict()
        self.pmesh3DCI.compute_cell_facet_normals_dict()
        self.pmesh3DCI.compute_cell_neighbor_dict()

        # pmesh is the owner of the compute_index function?

        return
#    def setUp(self):ENDDEF

    def test_1D_particle_migration(self):
        """Test particle migration across cells on a 1D mesh.

        """

        fncname = sys._getframe().f_code.co_name + '():'
        print '\ntest: ', fncname, '('+__file__+')'

        # List all the possible spatial coordinates
#        spatial_coordinates = ('x','y','z')

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.nsteps = 19

        # Run for 1, 2, and 3D meshes

        self.particleCI.pmeshCI = self.pmesh1DCI

        # Initialize or reinitialize the particles
#        for sp in self.particleCI.species_names:
#            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
#                self.particleCI.create_from_list(sp, printFlag=False, resetCounters=True)
#                self.particleCI.create_from_list(sp, printFlag=False)

        # Get the initial cell index of each particle.
        self.particleCI.compute_mesh_cell_indices()


        # Put the expected ending results into the p_expected tuple

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

        # Integrate for nsteps
        print "Moving", self.particleCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"
        for istep in xrange(ctrlCI.nsteps):
            self.particleCI.move_neutral_particles(ctrlCI.dt)

        # Check the results
        ncoords = self.particleCI.particle_dimension # number of particle coordinates to check
        for sp in self.particleCI.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particleCI.pseg_arr[sp].get(ip)
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

        fncname = sys._getframe().f_code.co_name + '():'
        print '\ntest: ', fncname, '('+__file__+')'

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.nsteps = 19

        # Run for 2D mesh

        self.particleCI.pmeshCI = self.pmesh2DCI

        # Get the initial cell index of each particle.
        self.particleCI.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particleCI.neutral_species[0]
        for ip in [0, 1]:
            p = self.particleCI.pseg_arr[sp].get(ip).copy() # Have to make a copy! Otherwise you overwrite the only copy of the particle
#            p_ic.append(self.particleCI.pseg_arr[sp].get(ip)) # Don't do this: it's a reference.
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

        # Integrate for nsteps
        print "Moving", self.particleCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "timesteps"
        for istep in xrange(ctrlCI.nsteps):
            self.particleCI.move_neutral_particles(ctrlCI.dt)

        # Create a mesh plotter to display the trajectory
        plotter=df_M.plot(self.particleCI.pmeshCI.mesh, title="First & last positions")

        # Check the results
        ncoords = self.particleCI.particle_dimension # number of particle coordinates to check
        for sp in self.particleCI.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particleCI.pseg_arr[sp].get(ip)
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

        fncname = sys._getframe().f_code.co_name + '():'
        print '\ntest:', fncname, '('+__file__+')'

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.nsteps = 19

        # Run for 3D mesh

        self.particleCI.pmeshCI = self.pmesh3DCI

        # Get the initial cell index of each particle.
        self.particleCI.compute_mesh_cell_indices()

        # Save the initial conditions (p_ic) for plotting
        p_ic = []
        sp = self.particleCI.neutral_species[0]
        for ip in [0, 1]:
            p = self.particleCI.pseg_arr[sp].get(ip).copy() # Have to make a copy!
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

        # Integrate for nsteps
        print "Moving", self.particleCI.get_total_particle_count(), "particles for", ctrlCI.nsteps, "steps"
        for istep in xrange(ctrlCI.nsteps):
            self.particleCI.move_neutral_particles(ctrlCI.dt)

        # Create a mesh plotter to display the trajectory (just the
        # first and last positions)
        plotter=df_M.plot(self.particleCI.pmeshCI.mesh, title="First & last positions")

        # Check the results
        ncoords = self.particleCI.particle_dimension # number of particle coordinates to check
        for sp in self.particleCI.neutral_species:
            for ip in [0, 1]:
                getparticle = self.particleCI.pseg_arr[sp].get(ip)
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
#        df_M.interactive() # Stops the plot from disappearing

        return
#    def test_3D_particle_migration(self):ENDDEF

    '''
    def test_1_particle_migration(self):
        """
           Test particle migration to neighboring cell on 1D, 2D, 3D
           meshes.

        """

        return

        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname, '('+__file__+')'

        # List all the possible spatial coordinates
        spatial_coordinates = ('x','y','z')

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 0.5
        ctrlCI.nsteps = 18

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
            for sp in self.particleCI.species_names:
                if self.particleCI.initial_distribution_type[sp] == 'listed':
                    # Put user-listed particles into the storage array
                    self.particleCI.create_from_list(sp, printFlag=False, resetCounters=True)

            # Get the initial cell index of each particle.
            self.particleCI.compute_mesh_cell_indices(meshCI)


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

            # Integrate for nsteps
            for istep in xrange(ctrlCI.nsteps):
                self.particleCI.move_particles_without_fields(ctrlCI.dt, meshCI=meshCI)

                    #self.assertTrue(meshCI.is_inside(pseg[ip], pseg[ip]['cell_index']), msg = "The computed cell does not contain the particle")

        return
#    def test_1_particle_migration(self):ENDDEF
    '''

#class TestParticleMigration(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
