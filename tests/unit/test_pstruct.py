#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2018, 2019 L. D. Hughes'
#__all__ = []

import sys
import os
import numpy as np_m
import unittest

# import dolfin as df_m
# import matplotlib.pyplot as mplot_m

# from Dolfin_Module import Mesh_C
# from Dolfin_Module import Field_C
# from Particle_Module import *

#from SegmentedArrayPair_Module import SegmentedArrayPair_C

# Here's the user's mesh definition for test_2
# from UserMesh_y_Fields_FE_XYZ_Module import *

# from UserUnits_Module import MyPlasmaUnits_C

# Use the C++ functions in the dnt_pstruct.so library
import dnt_pstruct

class TestPstruct(unittest.TestCase):
    """Test the C++ functions in dnt_pstruct.so"""
    
    def setUp(self):

        # initializations for each test go here...

        self.segment_length = 100
        self.delete_flag = 0b01 # the lowest bit is 1
        self.trajectory_flag = 0b10 # the second lowest bit is 1
        self.precision = np_m.float64

        return

#class TestPstruct(unittest.TestCase):
    def test_1_Cpp_cartesian_x(self):
        """Create a SegmentedArray for the "cartesian_x" particle type.

           Make two 1D particle tuples and put them into a SegmentedArrayPair_Cpp
           object of type "cartesian_x".

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        # Create C++ version of a SegmentedArray for cartesian_x particles
        seg_array_obj_Cpp_cartesian_x = dnt_pstruct.SegmentedArrayPair_Cpp_cartesian_x(self.segment_length)
        
        particle_dimension = 1
        x=1.5; x0=1.0; ux=3.0; weight = 101.1

        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        NO_CELL = -1
        cell_index = NO_CELL

        unique_ID = 7
        crossings = 5

        # Put two particles into a SegmentedArray store
        
        seg_arr = seg_array_obj_Cpp_cartesian_x

        # particle #1
        putparticle = (x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings)
#        print("#1 particle is:", putparticle)
        (seg_index, full_index) = seg_arr.push_back(putparticle)
#        print("#1 seg_index", seg_index, "full_index", full_index)

        # Retrieve the particle using the returned Numpy array it's in.
        (parray, offset) = seg_arr.get_segment_and_offset(full_index)
        getparticle = parray[offset]
#        print("#1 particle at offset", offset, "is", getparticle)
        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

        # Retrieve the particle using get_as_tuple(full_index)
        getparticle = seg_arr.get_as_tuple(full_index)
#        print("#1 particle at full_index", full_index, "is", getparticle)
        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")
        
        # particle #2
        putparticle = (x+0.5, x0+0.5, ux+0.5, weight, bitflags, cell_index, unique_ID+1, crossings)
#        print("#2 particle is:", putparticle)
        (seg_index, full_index) = seg_arr.push_back(putparticle)
#        print("#2 seg_index", seg_index, "full_index", full_index)

        # Retrieve the particle using the returned Numpy array it's in.
        (parray, offset) = seg_arr.get_segment_and_offset(full_index)
        getparticle = parray[offset]
#        print("#2 particle at offset", offset, "is", getparticle)
        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

        # Retrieve the particle using get_as_tuple(full_index)
        getparticle = seg_arr.get_as_tuple(full_index)
#        print("#2 particle at full_index", full_index, "is", getparticle)
        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")
        
        return
#    def test_1_Cpp_cartesian_x(self):ENDDEF

    def test_2_Cpp_cartesian_x_y(self):
        """Create a SegmentedArray for the "cartesian_x_y" particle type.

           Make a 3D particle tuple and put it into a SegmentedArrayPair_Cpp object
           of type "cartesian_x_y".

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        # Create C++ version of a SegmentedArray object for cartesian_x_y particles
        seg_array_obj_Cpp_cartesian_x_y = dnt_pstruct.SegmentedArrayPair_Cpp_cartesian_x_y(self.segment_length)        
        
        # Create a cartesian_x_y particle and put it into the SegmentedArray
        x=0.0; x0=x; y=1.0; y0=y; ux=3.0; uy=4; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        NO_CELL = -1
        cell_index = NO_CELL

        unique_ID = 7
        crossings = 5

        # Make a tuple
        putparticle = (x,y, x0,y0, ux,uy, weight, bitflags, cell_index, unique_ID, crossings)

        # Put this particle into the SegmentedArray store
        
        seg_arr = seg_array_obj_Cpp_cartesian_x_y
        (seg_index, full_index) = seg_arr.push_back(putparticle)

        # Retrieve the particle using the returned Numpy array it's in.
        (parray, offset) = seg_arr.get_segment_and_offset(full_index)
#        print("#1 particle at offset", offset, "is", parray[offset])
        # Retrieve the particle using get_as_tuple(full_index)
        getparticle = seg_arr.get_as_tuple(full_index)
#        print("#1 particle at full_index", full_index, "is", getparticle)

        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")
        
        return
#    def test_2_Cpp_cartesian_x_y(self):ENDDEF


    def test_3_Cpp_cartesian_x_y_z(self):
        """Create a SegmentedArray for the "cartesian_x_y_z" particle type.

           Make a 3D particle tuple and put it into a SegmentedArrayPair_Cpp object
           of type "cartesian_x_y_z".

        """

        fncName = '('+__file__+') ' + sys._getframe().f_code.co_name + '():\n'
        print('\ntest: ', fncName, '('+__file__+')')

        # Create C++ version of a SegmentedArray object for cartesian_x_y_z particles
        seg_array_obj_Cpp_cartesian_x_y_z = dnt_pstruct.SegmentedArrayPair_Cpp_cartesian_x_y_z(self.segment_length)        
        
        # Create a cartesian_x_y_z particle and put it into the SegmentedArray
        x=0.0; x0=x; y=1.0; y0=y; z=2.0; z0=z; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        NO_CELL = -1
        cell_index = NO_CELL

        unique_ID = 7
        crossings = 5

        # Make a tuple
        putparticle = (x,y,z, x0,y0,z0, ux,uy,uz, weight, bitflags, cell_index, unique_ID, crossings)

        # Put this particle into the SegmentedArray store
        
        seg_arr = seg_array_obj_Cpp_cartesian_x_y_z
        (seg_index, full_index) = seg_arr.push_back(putparticle)

        # Retrieve the particle using the returned Numpy array it's in.
        (parray, offset) = seg_arr.get_segment_and_offset(full_index)
#        print("#1 particle at offset", offset, "is", parray[offset])
        # Retrieve the particle using get_as_tuple(full_index)
        getparticle = seg_arr.get_as_tuple(full_index)
#        print("#1 particle at full_index", full_index, "is", getparticle)

        # Check the returned values against the input values
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")
        
        return
#    def test_3_Cpp_cartesian_x_y_z(self):ENDDEF

    def test_4_several_segments(self):
        """ Test the creation of several new segments.
        """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        # Create C++ version of a SegmentedArray object for cartesian_x_y_z particles
        seg_array_obj_Cpp_cartesian_x_y_z = dnt_pstruct.SegmentedArrayPair_Cpp_cartesian_x_y_z(self.segment_length)        

        x=0.0; x0=x; y=1.0; y0=y; z=2.0; z0=z; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        NO_CELL = -1
        cell_index = NO_CELL

        unique_ID = 1
        crossings = 5

        dx = 0.2

        seg_arr = seg_array_obj_Cpp_cartesian_x_y_z
        
        # Put in more particles than one segment can hold
        for i in range(self.segment_length+1):
            putparticle = (x,y,z, x0,y0,z0, ux,uy,uz, weight, bitflags, cell_index, unique_ID, crossings)
            seg_arr.push_back(putparticle)
            x += dx
            unique_ID += 1

        # Add the same number again, so there are 2 particles in the 3rd segment
        for i in range(self.segment_length+1):
            putparticle = (x,y,z, x0,y0,z0, ux,uy,uz, weight, bitflags, cell_index, unique_ID, crossings)
            seg_arr.push_back(putparticle)
            x += dx
            unique_ID += 1

        # Check that there are now 3 segments
        nseg_in, nseg_out = seg_arr.get_number_of_segments()
        self.assertEqual(nseg_out, 3, msg="Should have 3 segments now.")

        # Check the current capacity
        nmax_expected = nseg_out*self.segment_length
        nmax_in, nmax_out = seg_arr.get_capacity()
        self.assertEqual(nmax_expected, nmax_out, msg="Capacity returned is not correct")

        # Check the current number of megabytes allocated for the particle arrays
        n_double64 = 10 # Number of doubles in the structure
        n_bytes_per_double = 8 # Bytes per double
        n_int32 = 4 # Number of ints in the structure
        n_bytes_per_int = 4 # Bytes per int

        mb_expected = (n_double64*n_bytes_per_double+n_int32*n_bytes_per_int)*nmax_out/1.0e6
        mb_in, mb_out = seg_arr.get_number_of_mbytes()
        self.assertAlmostEqual(mb_out, mb_expected, msg="MB allocated is not correct")

        # Check the number of items currently stored
        n_expected = (nseg_out-1)*self.segment_length+2
        n = seg_arr.get_number_of_items()
        self.assertEqual(n_expected, n, msg="Count of stored items is not currect")

        # Get the last particle out using a full index and check it
        getparticle = seg_arr.get_as_tuple(2*self.segment_length+1)
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Last particle in Segment 3 is not correct")
        return
#    def test_4_several_segments(self): ENDDEF

#class TestPstruct(unittest.TestCase):ENDCLASS

if __name__ == '__main__':
    unittest.main()
