#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import unittest

from Dolfin_Module import Mesh_C
from SegmentedArrayPair_Module import SegmentedArrayPair_C

class TestSegmentedArray(unittest.TestCase):
    """Test the classes in SegmentedArray_Module.py"""
    
    def setUp(self):

        # initializations for each test go here...

        self.segment_length = 100
        self.delete_flag = 0b01 # the lowest bit is 1
        self.trajectory_flag = 0b10 # the second lowest bit is 1

        precision = numpy.float64
        position_coordinates = ('x', 'y', 'z')
        momentum_coordinates = ('ux', 'uy', 'uz')
        phase_coordinates = position_coordinates + momentum_coordinates
        pvars = [coord for coord in phase_coordinates]

        pvars.append('weight')
        self.pvars = pvars
#        pvartypes = [precision for coord in phase_coordinates]
        pvartypes = [precision for var in pvars]

        pvars.append('bitflags')
        pvartypes.append(numpy.int32)

        pvars.append('cell_index')
        pvartypes.append(numpy.int32) # The size determines how many local cells you can have.

        self.particle_dtype = {'names' : pvars, 'formats': pvartypes}

        metadata = {'names' : pvars, 'formats': pvartypes}
        
        self.seg_array_obj = SegmentedArrayPair_C(self.segment_length, metadata)

        # put a particle into the array


        return

    def test_1_nSeg(self):
        """ Test segment count at initialization. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        for iSA in (0, 1):
            nSeg = self.seg_array_obj.nSeg[iSA]
            number_of_segments = 1
            self.assertEqual(nSeg, number_of_segments, msg="nSeg should be 1")

    def test_2_field_names(self):
        """ Test field names. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

#        names = self.seg_array_obj.SegList[0].dtype.names
        for iSA in (0, 1):
            names = self.seg_array_obj.SegListPair[iSA][0].dtype.names
            for i in range(len(names)):
                self.assertEqual(names[i], self.pvars[i], msg="Field names are not correct")

    def test_3_put_and_getitem(self):
        """ Test putting particle data into the array using 'put'. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
#        cell_index = -1
        cell_index = Mesh_C.NO_CELL

        # Trim the number of coordinates here to match "position_coordinates" variable above
        putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cell_index)

        self.seg_array_obj.put(putparticle)

        # Get the item back out by subscripting and check it
        getparticle = self.seg_array_obj[0]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

        return
#    def test_3_put_and_getitem(self):ENDDEF
    
    def test_4_setitem_and_getitem(self):
        """ Test putting particle data into the array using a subscript. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL
        putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cell_index)
        self.seg_array_obj[0] = putparticle

        # Get the item back out by subscripting and check it
        getparticle = self.seg_array_obj[0]

        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

    def test_5_one_new_segment(self):
        """ Test the creation of a new segment when needed. Also test
        the get_number_of_segments() function.
        """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL
        dx = 0.1

        # Put in more particles than one segment can hold
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cell_index,)
            self.seg_array_obj.put(putparticle)
            x += dx

        # Check that there are now 2 segments
        nseg_in, nseg_out = self.seg_array_obj.get_number_of_segments()
        self.assertEqual(nseg_out, 2, msg="Should have 2 segments now.")

        # Get the last particle out by subscripting and check it
        getparticle = self.seg_array_obj[self.segment_length]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle in Segment 2 is not correct")

    def test_6_several_segments(self):
        """ Test the creation of several new segments.
        """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cell_index = Mesh_C.NO_CELL
        dx = 0.2

        # Put in more particles than one segment can hold
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cell_index,)
            self.seg_array_obj.put(putparticle)
            x += dx

        # Add more particles, so there are 2 particles in the 3rd segment
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cell_index,)
            self.seg_array_obj.put(putparticle)
            x += dx

        # Check that there are now 3 segments
        nseg_in, nseg_out = self.seg_array_obj.get_number_of_segments()
        self.assertEqual(nseg_out, 3, msg="Should have 3 segments now.")

        # Check the current capacity
        nmax_expected = nseg_out*self.segment_length
        nmax_in, nmax_out = self.seg_array_obj.get_item_capacity()
        self.assertEqual(nmax_expected, nmax_out, msg="Capacity returned is not correct")

        # Check the current number of megabytes allocated for the particle arrays
        n_double64 = 7
        n_bytes_per_double = 8
        n_int32 = 2
        n_bytes_per_int = 4
        mb_expected = (n_double64*n_bytes_per_double+n_int32*n_bytes_per_int)*nmax_out/1.0e6
        mb_in, mb_out = self.seg_array_obj.get_number_of_mbytes()
        self.assertAlmostEqual(mb_out, mb_expected, msg="MB allocated is not correct")

        # Check the number of items currently stored
        n_expected = (nseg_out-1)*self.segment_length+2
        n = self.seg_array_obj.get_number_of_items()
        self.assertEqual(n_expected, n, msg="Count of stored items is not currect")

        # Get the last particle out using [index] and check it
        getparticle = self.seg_array_obj[2*self.segment_length+1]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Last particle in Segment 3 is not correct")
        return
#    def test_6_several_segments(self): ENDDEF




#class TestSegmentedArray(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
