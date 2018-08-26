#!/usr/bin/env python3

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

#import math
import sys
import unittest

from UserUnits_Module import MyPlasmaUnits_C

#STARTCLASS
class TestUserUnits(unittest.TestCase):
    """Test the classes in UserUnits_Module.py"""
    
    def setUp(self):
        # initializations for each test go here...
        return

    def test_eV(self):

        fncname = sys._getframe().f_code.co_name
        print('\ntest: ', fncname, '('+__file__+')')

        # The temperature value input by the user:
        input_temp_eV = 3.0

        # Convert it to code units (Kelvin)
#        converted_temp = input_temp_eV*MyPlasmaUnits_C.eV
        # Compare the result to the correct value
#        correct_value = 3.0*11604.505

        # Convert it to code units (Joule)
        converted_temp = input_temp_eV*MyPlasmaUnits_C.elem_charge
        correct_value = 3.0*1.602176487e-19

        self.assertAlmostEqual(converted_temp, correct_value, delta = 0.1, msg = "!Temperature not converted correctly!")
        

if __name__ == '__main__':
    unittest.main()
