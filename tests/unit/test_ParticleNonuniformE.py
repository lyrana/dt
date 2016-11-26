#!/usr/bin/env python

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import importlib as im_M
import unittest

import dolfin as df_M

from DT_Module import DTparticleInput_C
from DT_Module import DTcontrol_C

from UserMesh_y_Fields_FE2D_Module import UserMesh_C

from Dolfin_Module import Mesh_C
from Dolfin_Module import Field_C

from Particle_Module import Particle_C

from SegmentedArrayPair_Module import SegmentedArray_C

from UserUnits_Module import MyPlasmaUnits_C

class TestParticleNonuniformE(unittest.TestCase):
    """Test classes in Particle_Module that push particles in a
       uniform E field
    """
    
    def setUp(self):

        # initializations for each test go here...

        # Create an instance of the DTparticleInput class
        pinCI = DTparticleInput_C()
        # Initialize particles
        pinCI.precision = numpy.float64

#        pinCI.copy_field_mesh = False
        pinCI.particle_integration_loop = 'loop-on-particles'
        pinCI.position_coordinates = ['x', 'y',] # determines the particle-storage dimensions
        pinCI.force_components = ['x', 'y',]
        pinCI.force_precision = numpy.float64

        # Specify the particle properties
        # 1. electrons
        pinCI.particle_species = (('testelectrons',
                             {'initial_distribution_type' : 'listed',
                              'charge' : -1.0*MyPlasmaUnits_C.elem_charge,
                              'mass' : 1.0*MyPlasmaUnits_C.electron_mass,
                              'dynamics' : 'explicit',
                              }
                             ),
                            )

        # Provide the particle distributions (lists of particles, functions, etc.)
        pinCI.user_particles_module = "UserParticles_2D_e"
        UPrt_M = im_M.import_module(pinCI.user_particles_module)
        pinCI.user_particles_class = UPrt_C = UPrt_M.ParticleDistributions_C

        self.pinCI = pinCI

        # Make the particle storage array for all species
        self.particleCI = Particle_C(pinCI, printFlag=False)

        # Store the particles
        for sp in self.particleCI.species_names:
            if self.particleCI.initial_distribution_type[sp] == 'listed':
                # Put user-listed particles into the storage array
                self.particleCI.create_from_list(sp, False)

        # Make the mesh & fields from saved files
        # The is needed to construct the function space for E

        # Create mesh from a file
        mesh2DCI = Mesh_C(meshFile="quarter_circle_mesh_crossed.xml", computeDictionaries=True, computeTree=True, plotFlag=False)

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

        self.neg_electric_field = Field_C(meshCI=mesh2DCI,
                                          element_type=electric_field_element_type,
                                          element_degree=phi_element_degree-1,
                                          field_type='vector')

        file = df_M.File("negE2D_crossed.xml")
        file >> self.neg_electric_field.function

        # Get the initial cell index of each particle.
        self.particleCI.compute_mesh_cell_indices(mesh2DCI)

        return
#    def setUp(self):ENDDEF

    def test_electric_field_push_1step(self):
        """ Check that the electric field push is correct.  Push
            sample particles for 1 step.
        """
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname

        ctrlCI = DTcontrol_C()

        ctrlCI.dt = 1.0e-5
        ctrlCI.nsteps = 1

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
        self.particleCI.move_particles_in_electrostatic_field(dt, self.neg_electric_field)


        # Check the results
        ncoords = self.particleCI.dimension # number of particle coordinates to check
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
    #     ctrlCI.nsteps = 10

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
        fncname = sys._getframe().f_code.co_name
        print '\ntest: ', fncname
        pass

        return

#class TestParticleNonuniformE(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
