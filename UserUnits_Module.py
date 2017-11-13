# UserUnits_M

#__version__ = 0.1
#__author__ = 'Copyright (C) 2016 L. D. Hughes'
#__all__ = []

# To used this class, put
#
#      import UserUnits_M
#
# or an abbreviation such as
#
#      import UserUnits_M as U_M
#
# in the python script.

#STARTCLASS
class MyPlasmaUnits_C(object):
    """ MyPlasmaUnits_C is to be edited by the user.  It gives the
        conversion factors to change input-parameter units to internal
        units.  The internal units are MKS, i.e.,

        Length: meters
        Mass: kilograms
        Time: seconds
        Temperature: Kelvin
        Charge: Coulombs

        For example, if you specify temperature in eV, you should to
        define an attribute (variable) 'eV' in this class, and assign
        to it the factor needed to convert to Kelvin, i.e.,

        'eV' : 11604.0 # Kelvin per eV

        Then use the multiplier UserUnits_M.MyPlasmaUnits_C.U.eV
        when writing the value of the temperature, e.g.,

        temperature = 1.6*UserUnits_M.MyPlasmaUnits_C.eV

        Reference:http://physics.nist.gov/cuu/Constants/index.html 
    """

# These variable are static

    # Convenient-unit conversions
    ns = 1.0e-9
    microsec = 1.0e-6
    coulomb = 1.0

    # Physical constants

# Should get all these from NIST above

    elem_charge = 1.602176487e-19 # Coulomb per user's unit of electric charge
    electron_mass = 9.10938215e-31 # kg per user's electron-mass unit
    proton_mass = 1.672621637e-27  # kg
    boltzmann_constant = 1.38064852e-23 # joules per Kelvin

    epsilon_0 = 8.854187817e-12 # Farad/m

#    eV = 11604.505 # Kelvin per eV

    m_per_s = 1.0 # m/s per user's velocity unit
    number_per_m3 = 1.0 # n/m3 per user's number-density unit

    # Derived units
    AMU = 1822.9*electron_mass # kg per user's atomic mass unit

# Print the values (fix this to print attributes instead)
    @staticmethod
    def print_units():
        print "Conversion factors in ", MyPlasmaUnits_C.__name__, ":"
        U = MyPlasmaUnits_C.U
        for key in U:
            print key, " = ", U[key]

#    def __init__(self):
#        self.number_per_m3 = 1.0
#        self.eV = 1.0
        return

#STARTCLASS
class MyCGSUnits_C(object):
    """ MyCGSUnits_C is to be edited by the user.  It gives the
        conversion factors to change input-parameter units to internal
        units.  The internal units are MKS, i.e.,

        Length: meters
        Mass: kilograms
        Time: seconds
        Temperature: Kelvin
        Charge: Coulombs

        For example, if you specify temperature in eV, you need to
        define a static attribute (variable) 'eV' in this class, to convert the
        value to Kelvin, i.e.,

        eV = 11604.0 # Kelvin per eV

        and then use the multiplier UserUnits_M.MyCGSUnits_C.eV when
        giving the temperature, e.g.,

        temperature = 1.6*UserUnits_M.MyCGSUnits_C.eV
    """

# Static data

# Write the unit used on the LHS, and the factor to convert this to
# internal units on the RHS
    eV = 11604.0 # Kelvin per eV
    m = 1.0e-2 # meters per user unit
    m_per_s = 1.0e-2 # m/s per user's velocity unit
    number_per_m3 = 1.0e-6 # n/m3 per user's number-density unit
