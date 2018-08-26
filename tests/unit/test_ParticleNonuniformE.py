#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_m
import unittest

import dolfin as df_m

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
        pin = ParticleInput_C()

        # Initialize particles
        pin.precision = numpy.float64
#        pin.copy_field_mesh = False
        pin.particle_integration_loop = 'loop-on-particles'
        pin.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pin.force_components = ['x', 'y',]
        pin.force_precision = numpy.float64

        # Specify the particle properties
        # 1. electrons
        # pin.particle_species = (('testelectrons',
        #                      {'initial_distribution_type' : 'listed',
        #                       'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
        #                       'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
        #                       'dynamics' : 'explicit',
        #                       }
        #                      ),
        #                     )

        speciesName = 'two_electrons'
        charge = -1.0*MyPlasmaUnits_C.elem_charge
        mass = 1.0*MyPlasmaUnits_C.electron_mass
        dynamics = 'explicit'
        testElectrons_S = ParticleSpecies_C(speciesName, charge, mass, dynamics)

        # Add these species to particle input
        pin.particle_species = (testElectrons_S,
                                 )
        # Make the particle object from pin
        self.particle_P = Particle_C(pin, print_flag=False)

        # Give the name of the .py file containing additional particle data (lists of
        # particles, boundary conditions, source regions, etc.)
        userParticlesModuleName = "UserParticles_2D_e"

        # Import this module
        userParticlesModule = im_m.import_module(userParticlesModuleName)
        self.particle_P.user_particles_class = userParticlesClass = userParticlesModule.UserParticleDistributions_C

        ### two_electrons are present at t=0

        # Name the initialized species (it should be in species_names above)
        speciesName = 'two_electrons'
        # Check that this species has been defined above
        if speciesName not in self.particle_P.species_names:
            print(fncName + "The species", speciesName, "has not been defined")
            sys.exit()

        # Specify how the species will be initialized
        initialDistributionType = 'listed'
        # Check that there's a function listing the particles particles
        printFlag = True
        if hasattr(userParticlesClass, speciesName):
            if printFlag: print(fncName + "(DnT INFO) Initial distribution for", speciesName, "is the function of that name in", userParticlesClass)
        # Write error message and exit if no distribution function exists
        else:
            errorMsg = fncName + "(DnT ERROR) Need to define a particle distribution function %s in %s for species %s " % (speciesName, userParticlesModuleName, speciesName)
            sys.exit(errorMsg)

        # Collect the parameters into a dictionary
        # The 'listed' type will expect a function with the same name as the species.
        testElectronParams = {'species_name': speciesName,
                              'initial_distribution_type': initialDistributionType,
                              }

        # The dictionary keys are mnemonics for the initialized particles
        initialParticlesDict = {'initial_two_electrons': (testElectronParams,),
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

        ### Mesh creation

        umi = UserMeshInput2DCirc_C()

        # Create mesh from a file
        umi.mesh_file = 'mesh_quarter_circle_crossed.xml'
        umi.particle_boundary_file='mesh_quarter_circle_crossed_Pbcs.xml'

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

        umi.particle_boundary_dict = particleBoundaryDict

        pmesh2D_M = UserMesh2DCirc_C(umi, compute_dictionaries=True, compute_tree=True, plot_flag=False)

        self.particle_P.pmesh_M = pmesh2D_M

        ### Particle boundary-conditions

        # UserParticleBoundaryFunctions_C is where the facet-crossing callback
        # functions are defined.
        # userPBndFnsClass = userParticlesModule.UserParticleBoundaryFunctions_C # Now need an instantiation, not just a class name:
        userPBndFns = userParticlesModule.UserParticleBoundaryFunctions_C(self.particle_P.position_coordinates, self.particle_P.dx)

        spNames = self.particle_P.species_names
        pmeshBCs = ParticleMeshBoundaryConditions_C(spNames, pmesh2D_M, userPBndFns, print_flag=False)
        self.particle_P.pmesh_bcs = pmeshBCs

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
        self.neg_electric_field = Field_C(pmesh2D_M,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_m.File("negE_test_2_2D.xml")
        file >> self.neg_electric_field.function

        # Get the initial cell index of each particle.
        self.particle_P.compute_mesh_cell_indices()

        return
#    def setUp(self):ENDDEF

    def test_electric_field_push_1step(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 1 step.
        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        ctrl = DTcontrol_C()

        ctrl.dt = 1.0e-5
        ctrl.n_timesteps = 1
                
        dt = ctrl.dt

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

#        species_names = self.particle_P.names

# first particle:        
#        getp = self.particle_P.pseg_arr['testelectrons'].get(0)
#        print 'getp is a', type(getp)
#        print '1st electron:', getp
# second particle:
#        getparticle = self.particle_P.pseg_arr['testelectrons'].get(1)
#        print '2nd electron:', getparticle

        # Advance the particles one timestep
        print("Moving", self.particle_P.get_total_particle_count(), "particles for one timestep")
        ctrl.time_step = 0
        ctrl.time = 0.0

        self.particle_P.move_particles_in_electrostatic_field(ctrl, self.neg_electric_field)

        # Check the results
        ncoords = self.particle_P.particle_dimension # number of particle coordinates to check
        isp = 0
        for sp in self.particle_P.species_names:
#            print 'species count = ', self.particle_P.get_species_particle_count(sp)
            if self.particle_P.get_species_particle_count(sp) == 0: continue

            # Check that the first two particles in the array reaches the correct values
            for ip in [0, 1]:
                getparticle = self.particle_P.pseg_arr[sp].get(ip)
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

    #     ctrl = DTcontrol_C()

    #     ctrl.dt = 1.0e-5
    #     ctrl.n_timesteps = 10

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
        print('\ntest: ', fncName, '('+__file__+')')
        pass

        return

#class TestParticleNonuniformE(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
