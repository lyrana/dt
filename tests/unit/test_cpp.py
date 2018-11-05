#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

import sys
import numpy
import unittest

from Dolfin_Module import Mesh_C
from SegmentedArrayPair_Module import SegmentedArray_C

# The C++ functions in a .so library
import pseg_cpp

class TestCPP(unittest.TestCase):
    """Test the C++ functions in DnT_cpp.so"""
    
    def setUp(self):

        # initializations for each test go here...

        self.segment_length = 100
        self.delete_flag = 0b01 # the lowest bit is 1
        self.trajectory_flag = 0b10 # the second lowest bit is 1

        precision = numpy.float64
        pvars = ['x', 'x0', 'ux',]

        pvars.append('weight')
        
#        pvartypes = [precision for coord in phase_coordinates]
        pvartypes = [precision for var in pvars]

        pvars.append('bitflags')
        pvartypes.append(numpy.int32)

        pvars.append('cell_index')
        pvartypes.append(numpy.int32) # The size determines how many local cells you can have.

        pvars.append('unique_ID')
        pvartypes.append(numpy.int32)
        
        pvars.append('crossings')
        pvartypes.append(numpy.int32)

#         pvars = ['x_', 'x0_', 'ux_',]

#         pvars.append('weight_')
        
# #        pvartypes = [precision for coord in phase_coordinates]
#         pvartypes = [precision for var in pvars]

#         pvars.append('bitflags_')
#         pvartypes.append(numpy.int32)

#         pvars.append('cell_index_')
#         pvartypes.append(numpy.int32) # The size determines how many local cells you can have.

#         pvars.append('unique_ID_')
#         pvartypes.append(numpy.int32)
        
#         pvars.append('crossings_')
#         pvartypes.append(numpy.int32)
        
        self.pvars = pvars
        
        self.particle_dtype = {'names' : pvars, 'formats': pvartypes}

        metadata = {'names' : pvars, 'formats': pvartypes}
        
        self.seg_array_obj = SegmentedArray_C(self.segment_length, metadata)

        return

    def test_1_print_pseg(self):
        """ Test print_pseg"""

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        x=1.5; x0=1.0; ux=3.0; weight = 101.1
        bitflags = 0b00 # initialize all bits to 0
        bitflags = bitflags | self.trajectory_flag # turn on trajectory flag
#        cell_index = -1
        cell_index = Mesh_C.NO_CELL

        unique_ID = 7
        crossings = 5

        # Trim the number of coordinates here to match "position_coordinates" variable above

        putparticle = (x, x0, ux, weight, bitflags, cell_index, unique_ID, crossings)
        self.seg_array_obj.put(putparticle)

        putparticle = (x+0.5, x0+0.5, ux+0.5, weight, bitflags, cell_index, unique_ID+1, crossings)
        self.seg_array_obj.put(putparticle)

        # Get the item back out by subscripting and check it
        getparticle = self.seg_array_obj[0]
#        for i in range(len(getparticle)):
#            self.assertEqual(getparticle[i], putparticle[i], msg="Particle variables are not correct")
        print("getparticle =", getparticle)

        pseg_arr = self.seg_array_obj

        (npSeg, pseg) = self.seg_array_obj.init_out_loop()        
        print("npSeg=", npSeg)
        print("pseg_cpp returns", pseg_cpp.print_pseg(pseg))
        
        return
#    def test_1_print_pseg:ENDDEF
    
#class TestCPP(unittest.TestCase): ENDCLASS

if __name__ == '__main__':
    unittest.main()
