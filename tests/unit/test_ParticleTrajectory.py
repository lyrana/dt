#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTmeshInput_C
from DT_Module import DTparticleInput_C
from DT_Module import DTcontrol_C
from DT_Module import DTtrajectoryInput_C

#from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

from Particle_Module import Particle_C
from Particle_Module import ParticleMeshBoundaryConditions_C

from Trajectory_Module import Trajectory_C

from SegmentedArrayPair_Module import SegmentedArray_C

from UserMesh_y_Fields_FE2D_Module import UserMesh_C

from UserUnits_Module import MyPlasmaUnits_C

class TestParticleTrajectory(unittest.TestCase):
    """Test classes in Trajectory_Module.
    """
    
    def setUp(self):

        # initializations for each test go here...

        self.ctrlCI = DTcontrol_C()

        self.ctrlCI.dt = 1.0e-6

        ### Particle species input

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle-species properties
        # 1. electrons
        pinCI.particle_species = (('trajelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
#                              'dynamics' : 'implicit',
                              }
                             ),
                            )

        # Provide the particle distributions for the above species
        pinCI.user_particles_module = "UserParticles_2D_e"
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_M.UserParticleDistributions_C

        # Make the particle object from pinCI
        self.particleCI = Particle_C(pinCI, printFlag=False)

        ###  Mesh creation

        miCI = DTmeshInput_C()

        # Make the mesh & fields from saved files
        miCI.mesh_file = 'quarter_circle_mesh_crossed.xml'
        miCI.particle_boundary_file='Pbcs_quarter_circle_mesh_crossed.xml'

        # Create the particle mesh object
# 'crossed' diagonals
#        self.pmesh2DCI = Mesh_C(meshFile="quarter_circle_mesh_crossed.xml", computeDictionaries=True, computeTree=True, plotFlag=False)
# 'left' diagonal
#        self.pmesh2DCI = Mesh_C(meshFile="quarter_circle_mesh_left.xml", computeDictionaries=True, computeTree=True, plotFlag=False)

        # These are the boundary-name -> int pairs used to mark mesh facets:
        rminIndx = 1
        rmaxIndx = 2
        thminIndx = 4
        thmaxIndx = 8
        particleBoundaryDict = {'rmin': rminIndx,
                                'rmax': rmaxIndx,
                                'thmin': thminIndx,
                                'thmax': thmaxIndx,
                                }

        miCI.particle_boundary_dict = particleBoundaryDict

        # Read the mesh from an existing file

#        pinCI.pmeshCI = UserMesh_C(meshFile='quarter_circle_mesh_crossed.xml', particleBoundaryFile='Pbcs_quarter_circle_mesh_crossed.xml', computeDictionaries=True, computeTree=True, plotFlag=False)

# Can this be attached to Particle_C after Particle_C construction? YES
        pmeshCI = UserMesh_C(miCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)

        self.particleCI.pmeshCI = pmeshCI

        ### Field creation

        # The following value should correspond to the element degree
        # used in the potential from which negE was obtained
        phi_element_degree = 1

        if phi_element_degree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous Galerkin
            # elements.
            electric_field_element_type = "DG"
        else:
            electric_field_element_type = "Lagrange"

        ## Create the electric field from an existing file

        # Put the negative electric field directly on the particle mesh
        self.neg_electric_field = Field_C(meshCI=pmeshCI,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_M.File("negE2D_crossed.xml")
        file >> self.neg_electric_field.function

        ## Particle boundary-conditions

        # UserParticleMeshFunctions_C is where the facet-crossing callback
        # functions are defined.
        userPMeshFnsClass = UPrt_M.UserParticleMeshFunctions_C # abbreviation

        # Make the particle-mesh boundary-conditions object and add it
        # to the particle object.
        spNames = self.particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmeshCI, userPMeshFnsClass, printFlag=False)
        self.particleCI.pmesh_bcCI = pmeshBCCI

        ### Create a particle trajectory object

        # Use an input object to collect initialization data for the trajectory object
        self.trajinCI = DTtrajectoryInput_C()

        self.trajinCI.maxpoints = None # Set to None to get every point

        # Specify which particle variables to save.  This has the
        # form of a numpy dtype specification.
        self.trajinCI.explicitDict = {'names': ['x', 'ux', 'y', 'uy', 'Ex', 'Ey'], 'formats': [numpy.float32]*6}
        self.trajinCI.implicitDict = {'names': ['x', 'ux', 'phi'], 'formats': [numpy.float32]*3}
        self.trajinCI.neutralDict = {'names': ['x', 'ux', 'y', 'uy'], 'formats': [numpy.float32]*4}

        return

#class TestParticleTrajectory(unittest.TestCase):
    def test_1_trajectory_init(self):
        """ Check that the trajectory variables are saved correctly in
            the trajectory DataList.

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName

        # Abbreviations
        pCI = self.particleCI

        self.ctrlCI.nsteps = 13

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        trajCI = Trajectory_C(self.trajinCI, self.ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)
        self.particleCI.trajCI = trajCI


        # if os.environ.get('DISPLAY') is None:
        #     plotFlag=False
        # else:
        #     plotFlag=True

        # Create particles that are selected for trajectories.
# Could call particleCI.initialize_distributions() instead?

        pCI.create_from_list('trajelectrons', printFlag=True)

        # Check the results

#        print "The explicit species with trajectory particles are", pCI.trajCI.explicit_species
        expectedList = self.trajinCI.explicitDict['names']
        for sp in pCI.trajCI.explicit_species:
#            print "  with particle indices", pCI.trajCI.ParticleIdList[sp]
#            print "  and values", pCI.trajCI.DataList[sp][0].dtype.names
            gotList = pCI.trajCI.DataList[sp][0].dtype.names
            for i in range(len(expectedList)):
                self.assertEqual(gotList[i], expectedList[i], msg = "Trajectory variable is not correct for and explicit particle")

#        print "The implicit species are", pCI.trajCI.implicit_species
        expectedList = self.trajinCI.implicitDict['names']
        for sp in pCI.trajCI.implicit_species:
#            print "  with indices", pCI.trajCI.ParticleIdList[sp]
#            print "  and values", pCI.trajCI.DataList[sp][0].dtype.names
            gotList = pCI.trajCI.DataList[sp][0].dtype.names
            for i in range(len(expectedList)):
                self.assertEqual(gotList[i], expectedList[i], msg = "Trajectory variable is not correct for an implicit particle")
#                self.assertAlmostEqual(getparticle[ix], p_expected[isp][ix], msg="Particle is not in correct position")
        return
#    def test_1_trajectory_init(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):
    def test_2_record_trajectory(self):
        """ Record and plot the requested trajectory data.
            Checks the final particle position.
        """

        testName = sys._getframe().f_code.co_name
        fncName = '('+__file__+') ' + testName + '():\n'
        print '\ntest: ', fncName

        pCI = self.particleCI

        self.ctrlCI.nsteps = 13
        dt = self.ctrlCI.dt

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        trajCI = Trajectory_C(self.trajinCI, self.ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)
        self.particleCI.trajCI = trajCI

        # Create particles that are selected for trajectories.
        pCI.create_from_list('trajelectrons', printFlag=True)

        # Get the initial cell index of each particle.
        pCI.compute_mesh_cell_indices()

        print "Moving", pCI.get_total_particle_count(), "particles for", self.ctrlCI.nsteps, "timesteps"
        for istep in xrange(self.ctrlCI.nsteps):
            # needs dt; doesn't need nsteps

            # Gather particle trajectory data
            # First point is the initial particle condition at istep = 0
            if pCI.trajCI is not None:
#                print 'pCI.trajCI.skip:', pCI.trajCI.skip
                if istep % pCI.trajCI.skip == 0:
                    pCI.record_trajectory_data(neg_E_field=self.neg_electric_field)

            # Do the implicit species first
            if len(pCI.implicit_species) != 0:
                self.iterate_implicit_electrostatic_particles(dt, pCI.implicit_species)

            # Then move the explicit species
            if len(pCI.explicit_species) != 0:
                pCI.move_particles_in_electrostatic_field(dt, self.neg_electric_field)

        # Record the LAST points on the particle trajectories
        if pCI.trajCI is not None:
                pCI.record_trajectory_data(neg_E_field=self.neg_electric_field)

        # Check the results

#        print 'trajelectrons data x:', pCI.trajCI.DataList['trajelectrons'][0][:]['x']
#        print 'trajelectrons data ux:', pCI.trajCI.DataList['trajelectrons'][0][:]['ux']
#        print 'trajelectrons data Ex:', pCI.trajCI.DataList['trajelectrons'][0][:]['Ex']

        # The expected results from ParticleNonuniformE.ods
        		
        # First electron
        xp1 = 0.77759792; yp1 = 0.78651935; zp1 = 0.0
        p1 = (xp1,yp1,)

        # Second electron
        xp2 = 1.02059224; yp2 = 0.42618424; zp2 = 0.0
        p2 = (xp2,yp2,)

        p_expected = (p1, p2)

        # Check the results
        ncoords = pCI.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in pCI.species_names:
#            print 'pCI.get_species_particle_count(sp)', pCI.get_species_particle_count(sp)
            if pCI.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = pCI.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        # Plot the trajectory onto the particle mesh

        mesh = pCI.pmeshCI.mesh
        holdPlot = True # Set to True to stop the plot from disappearing.
        pCI.trajCI.plot_trajectories_on_mesh(mesh, testName, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        # Plot the trajectory in phase-space

        plotFlag = False
        if os.environ.get('DISPLAY') is not None and plotFlag is True:
            pCI.trajCI.plot_trajectories() # Phase-space plot of trajectory

        return
#    def test_2_record_trajectory(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):
    def test_3_out_of_bounds(self):
        """ Record and plot the requested trajectory data.
            Make the particles leave the mesh.
        """

        testName = sys._getframe().f_code.co_name
        fncName = '('+__file__+') ' + testName + '():\n'
        print '\ntest: ', fncName

        pCI = self.particleCI

        self.ctrlCI.nsteps = 14 # This puts the particle out-of-bounds
        dt = self.ctrlCI.dt

        ## Create the trajectory object and attach it to the particle object.
        # No trajectory storage is created until particles
        # with TRAJECTORY_FLAG on are encountered.
        trajCI = Trajectory_C(self.trajinCI, self.ctrlCI, pCI.explicit_species, pCI.implicit_species, pCI.neutral_species)
        self.particleCI.trajCI = trajCI

        # Create particles that are selected for trajectories.
        pCI.create_from_list('trajelectrons', printFlag=True)

        # Get the initial cell index of each particle.
        pCI.compute_mesh_cell_indices()

        print "Moving", pCI.get_total_particle_count(), "particles for", self.ctrlCI.nsteps, "timesteps"
        for istep in xrange(self.ctrlCI.nsteps):
            # needs dt; doesn't need nsteps

            # Gather particle trajectory data
            # First point is the initial particle condition at istep = 0
            if pCI.trajCI is not None:
#                print 'pCI.trajCI.skip:', pCI.trajCI.skip
                if istep % pCI.trajCI.skip == 0:
                    pCI.record_trajectory_data(neg_E_field=self.neg_electric_field)

            # Do the implicit species first
            if len(pCI.implicit_species) != 0:
                self.iterate_implicit_electrostatic_particles(dt, pCI.implicit_species)

            # Then move the explicit species
            if len(pCI.explicit_species) != 0:
                pCI.move_particles_in_electrostatic_field(dt, self.neg_electric_field)

        # Record the LAST points on the particle trajectory
        if pCI.trajCI is not None:
                pCI.record_trajectory_data(neg_E_field=self.neg_electric_field)

        # Check the results

#        print 'trajelectron 0 data x:', pCI.trajCI.DataList['trajelectrons'][0][:]['x']
#        print 'trajelectron 1 data x:', pCI.trajCI.DataList['trajelectrons'][1][:]['x']
#        print 'trajelectron 0 data ux:', pCI.trajCI.DataList['trajelectrons'][0][:]['ux']
#        print 'trajelectron 1 data Ex:', pCI.trajCI.DataList['trajelectrons'][0][:]['Ex']

        # The expected results from ParticleNonuniformE.ods
        		
        # First electron
        xp1 = 0.77759792; yp1 = 0.78651935; zp1 = 0.0
        p1 = (xp1,yp1,)

        # Second electron
        xp2 = 1.02059224; yp2 = 0.42618424; zp2 = 0.0
        p2 = (xp2,yp2,)

        p_expected = (p1, p2)

        # Check the results
        ncoords = pCI.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in pCI.species_names:
            if pCI.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = pCI.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        # Plot the trajectory onto the particle mesh

        mesh = pCI.pmeshCI.mesh
        holdPlot = True # Set to True to stop the plot from disappearing.
        pCI.trajCI.plot_trajectories_on_mesh(mesh, testName, hold_plot=holdPlot) # Plots trajectory spatial coordinates on top of the particle mesh

        # Plot the trajectory in phase-space

#        import matplotlib.pyplot as plot_M
#        tvals = np_M.linspace(tmin, tmax, 

        plotFlag = False
        if os.environ.get('DISPLAY') is not None and plotFlag is True:
            pCI.trajCI.plot_trajectories()

        return
#    def test_3_out_of_bounds(self):ENDDEF

#class TestParticleTrajectory(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
