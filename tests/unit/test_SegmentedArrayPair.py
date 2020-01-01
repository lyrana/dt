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
        positionCoordinates = ('x', 'y', 'z')
        momentumCoordinates = ('ux', 'uy', 'uz')
        phaseCoordinates = positionCoordinates + momentumCoordinates
        pvars = [coord for coord in phaseCoordinates]

        pvars.append('weight')
        self.pvars = pvars
        pvartypes = [precision for var in pvars]

        pvars.append('bitflags')
        pvartypes.append(numpy.int32)

        pvars.append('cell_index')
        pvartypes.append(numpy.int32) # The size determines how many local cells you can have.

        self.particle_dtype = {'names' : pvars, 'formats': pvartypes}

        metadata = {'names' : pvars, 'formats': pvartypes}
        
        self.sap = SegmentedArrayPair_C(self.segment_length, metadata)

        return

    def test_1_nSeg(self):
        """ Test segment count at initialization. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        for iSA in (0, 1):
            nSeg = self.sap.nseg[iSA]
            number_of_segments = 1
            self.assertEqual(nSeg, number_of_segments, msg="nSeg should be 1")

        return
#    def test_1_nSeg(self):ENDDEF

    def test_2_field_names(self):
        """ Test field names. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        for iSA in (0, 1):
            names = self.sap.seg_list_pair[iSA][0].dtype.names
            for i in range(len(names)):
                self.assertEqual(names[i], self.pvars[i], msg="Field names are not correct")

        return
#    def test_2_field_names(self):ENDDEF

    def test_3_push_and_getitem(self):
        """ Test putting particle data into the array using 'put'. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
#        cellIndex = -1
        cellIndex = Mesh_C.NO_CELL

        # Trim the number of coordinates here to match "positionCoordinates" variable above
        putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cellIndex)

        self.sap.push_back(putparticle)

        # Get the item back out by subscripting and check it
        getparticle = self.sap[0]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

        return
#    def test_3_push_and_getitem(self):ENDDEF
    
    def test_4_setitem_and_getitem(self):
        """ Test putting particle data into the array using a subscript. """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cellIndex = Mesh_C.NO_CELL
        putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cellIndex)
        self.sap[0] = putparticle

        # Get the item back out by subscripting and check it
        getparticle = self.sap[0]

        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")

        return
#    def test_4_setitem_and_getitem(self):ENDDEF


    def test_5_one_new_segment(self):
        """ Test the creation of a new segment when needed. Also test
        the get_number_of_segments() function.
        """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cellIndex = Mesh_C.NO_CELL
        dx = 0.1

        # Put in more particles than one segment can hold
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cellIndex,)
            self.sap.push_back(putparticle)
            x += dx

        # Check that there are now 2 segments
        nsegIn, nsegOut = self.sap.get_number_of_segments()
        self.assertEqual(nsegOut, 2, msg="Should have 2 segments now.")

        # Get the last particle out by subscripting and check it
        getparticle = self.sap[self.segment_length]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Particle in Segment 2 is not correct")

        return
#    def test_5_one_new_segment(self):ENDDEF
            
    def test_6_several_segments(self):
        """ Test the creation of several new segments.
        """

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
        cellIndex = Mesh_C.NO_CELL
        dx = 0.2

        # Put in more particles than one segment can hold
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cellIndex,)
            self.sap.push_back(putparticle)
            x += dx

        # Add more particles, so there are 2 particles in the 3rd segment
        for i in range(self.segment_length+1):
            putparticle = (x, y, z, ux, uy, uz, weight, bitflags, cellIndex,)
            self.sap.push_back(putparticle)
            x += dx

        # Check that there are now 3 segments
        nsegIn, nsegOut = self.sap.get_number_of_segments()
        self.assertEqual(nsegOut, 3, msg="Should have 3 segments now.")

        # Check the current capacity
        nmaxExpected = nsegOut*self.segment_length
        nmaxIn, nmaxOut = self.sap.get_capacity()
        self.assertEqual(nmaxExpected, nmaxOut, msg="Capacity returned is not correct")

        # Check the current number of megabytes allocated for the particle arrays
        #     These sizes need to be set for each case separately
        nDouble64 = 7
        nBytesPerDouble = 8
        nInt32 = 2
        nBytesPerInt = 4
        mbExpected = (nDouble64*nBytesPerDouble+nInt32*nBytesPerInt)*nmaxOut/1.0e6
        mbIn, mbOut = self.sap.get_number_of_mbytes()
        self.assertAlmostEqual(mbOut, mbExpected, msg="MB allocated is not correct")

        # Check the number of items currently stored
        nExpected = (nsegOut-1)*self.segment_length+2
        n = self.sap.get_number_of_items()
        self.assertEqual(nExpected, n, msg="Count of stored items is not currect")

        # Get the last particle out using [index] and check it
        getparticle = self.sap[2*self.segment_length+1]
        for i in range(len(getparticle)):
            self.assertEqual(getparticle[i], putparticle[i], msg="Last particle in Segment 3 is not correct")

        return
#    def test_6_several_segments(self): ENDDEF


#class TestSegmentedArray(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
