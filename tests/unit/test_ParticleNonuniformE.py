#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTcontrol_C

from UserMesh_y_Fields_FE2D_Module import *

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

from Particle_Module import *

from SegmentedArrayPair_Module import SegmentedArray_C

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestParticleNonuniformE(unittest.TestCase):
    """Test classes in Particle_Module that push particles in a
       uniform E field
    """
    
    def setUp(self):

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'

        # initializations for each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = ParticleInput_C()

        # Initialize particles
        pinCI.precision = numpy.float64
#        pinCI.copy_field_mesh = False
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle properties
        # 1. electrons
        # pinCI.particle_species = (('testelectrons',
        #                      {'initial_distribution_type' : 'listed',
        #                       'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
        #                       'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
        #                       'dynamics' : 'explicit',
        #                       }
        #                      ),
        #                     )

        speciesName = 'test_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        testElectronsCI = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pinCI.particle_species = (testElectronsCI,
                                 )
        # Make the particle object from pinCI
        self.particleCI = Particle_C(pinCI, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticleModule = "UserParticles_2D_e"

        # Import this module
        UPrt_M = im_M.import_module(userParticleModule)
        self.particleCI.user_particle_class = userParticleClass = UPrt_M.UserParticleDistributions_C

        ### test_electrons are present at t=0

        # Name the initialized species (it should be in species_names above)
        speciesName = 'test_electrons'
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
        testElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_test_electrons': (testElectronParams,),
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

        ### Mesh creation

        umiCI = UserMeshInput_C()

        # Create mesh from a file
        umiCI.mesh_file = 'quarter_circle_mesh_crossed.xml'
        umiCI.particle_boundary_file='Pbcs_quarter_circle_mesh_crossed.xml'

        # These are the (int boundary-name) pairs used to mark mesh
        # facets. The string value of the int is used as the index.
        rmin_indx = 1
        rmax_indx = 2
        thmin_indx = 4
        thmax_indx = 8
        particleBoundaryDict = {'rmin': rmin_indx,
                                'rmax': rmax_indx,
                                'thmin': thmin_indx,
                                'thmax': thmax_indx,
                                }

        umiCI.particle_boundary_dict = particleBoundaryDict

        pmesh2DCI = UserMesh_C(umiCI, compute_dictionaries=True, compute_tree=True, plot_flag=False)

        self.particleCI.pmeshCI = pmesh2DCI

        ### Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
        userPBndFnsClass = UPrt_M.UserParticleBoundaryFunctions_C # abbreviation

        spNames = self.particleCI.species_names
        pmeshBCCI = ParticleMeshBoundaryConditions_C(spNames, pmesh2DCI, userPBndFnsClass, print_flag=False)
        self.particleCI.pmesh_bcCI = pmeshBCCI

        # The following value should correspond to the element degree
        # used in the potential from which negE was obtained
        phi_element_degree = 1

        # Create a Field_C object to store the field, and to provide
        # interpolation methods.

        if phi_element_degree == 1:
            # For linear elements, grad(phi) is discontinuous across
            # elements. To represent this field, we need Discontinuous
            # Galerkin elements.
            electric_field_element_type = "DG"
        else:
            electric_field_element_type = "Lagrange"

        # Create the negative electric field directly on the particle mesh
        self.neg_electric_field = Field_C(meshCI=pmesh2DCI,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_M.File("negE2D_crossed.xml")
        file >> self.neg_electric_field.function

        # Get the initial cell index of each particle.
        self.particleCI.compute_mesh_cell_indices()

        return
#    def setUp(self):ENDDEF

    def test_electric_field_push_1step(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 1 step.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print '\ntest: ', fncName, '('+__file__+')'

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 1.0e-5
        ctrlCI.n_timesteps = 1

        dt = ctrlCI.dt

        # The expected results from ParticleNonuniformE.ods
        		
        # First electron
        xp1 = 1.005609924; yp1 = 0.8131367827; zp1 = 0.0
        vxp1 = -235314.728666394; vyp1 = -254562.042790713; vzp1 = 0.0
        weight1=1.0
        p1 = (xp1,yp1, vxp1,vyp1, weight1)

        # Second electron
        xp2 = -9.3684645671; yp2 = -10.1983686032; zp2 = 0.0
        vxp2 = -1014628.2027; vyp2 = -1097618.606318; vzp2 = 0.0
        weight2 = 3.0
        p2 = (xp2,yp2, vxp2,vyp2, weight2)

        p_expected = (p1, p2)

#        species_names = self.particleCI.names

# first particle:        
#        getp = self.particleCI.pseg_arr['testelectrons'].get(0)
#        print 'getp is a', type(getp)
#        print '1st electron:', getp
# second particle:
#        getparticle = self.particleCI.pseg_arr['testelectrons'].get(1)
#        print '2nd electron:', getparticle

        # Advance the particles one timestep
        print "Moving", self.particleCI.get_total_particle_count(), "particles for one timestep"
        ctrlCI.time_step = 0
        ctrlCI.time = 0.0

        self.particleCI.move_particles_in_electrostatic_field(dt, self.neg_electric_field)

        # Check the results
        ncoords = self.particleCI.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in self.particleCI.species_names:
#            print 'species count = ', self.particleCI.get_species_particle_count(sp)
            if self.particleCI.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = self.particleCI.pseg_arr[sp].get(ip)
#                print 'calculated = ', getparticle
#                print 'expected = ', p_expected[ip]
                for ic in range(ncoords):
#                    print "result:", getparticle[ic]/p_expected[ip][ic]
# Note: for different field solver, may have to reduce the places:
                    self.assertAlmostEqual(getparticle[ic]/p_expected[ip][ic], 1.0, places=6, msg="Particle is not in correct position")
#                    print "ic", ic, "is OK"
            isp += 1

        return
#    def test_electric_field_push_1step(self):ENDDEF

    # def test_electric_field_push_10steps(self):
    #     """ Check that the electric field push is correct.  Push
    #         sample particles for 10 steps.
    #     """
    #     fncname = sys._getframe().f_code.co_name
    #     print '\ntest: ', fncname

    #     ctrlCI = DTcontrol_C()

    #     ctrlCI.dt = 1.0e-5
    #     ctrlCI.n_timesteps = 10

    #     return

    # def test_magnetic_field_push(self):            
    #     """ 1. Check that the magnetic field push is correct.
    #         2. 
    #     """
    #     fncname = sys._getframe().f_code.co_name
    #     print '\ntest: ', fncname
    #     pass

    #     return

    def test_something(self):
        """ Check the number of particles in each species.
        """
        fncName = sys._getframe().f_code.co_name
        print '\ntest: ', fncName, '('+__file__+')'
        pass

        return

#class TestParticleNonuniformE(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
