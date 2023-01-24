# Example: Non-symmetric laminate subject to temperature change
#
# A non-symmetric laminate is subjected to a temperature change. We
# will use Tompouce to:
#
# 1. Calculate the resulting curvature
# 2. The bending moment required to flatten the laminate
#
import sys; sys.path.append('../')
from tompouce import Material, Ply, Laminate, Load
from numpy import pi

# First, a Material object is created. This time we use a dictionary
# to create the object.
mdata = {'name': "TC1225", 'manufacturer': "Toray AC",
         'fiber': "Carbon AS4", 'matrix': "Victrex LM-PAEK",
         'E1': 135E9, 'E2': 10E9, 'v12': 0.26, 'G12': 4.3E9,
         'alpha1': 4E-7, 'alpha2': 3E-5,
         'S1t': 2410E6, 'S1c': 1300E6, 'S2t': 86E6, 'S2c': 100E6, 'S6': 42E6}
TC1225 = Material(mdata)

# As a second step, we create two Ply objects: one with a fiber
# orientation of 0 degrees and one with an orientation of 90 degrees.
ply_thickness = 0.14E-3
P00 = Ply(TC1225, 0, ply_thickness)
P90 = Ply(TC1225, pi/2, ply_thickness)

# A Laminate object can be created from a list with Ply objects.
CP_laminate = Laminate([P00, P00, P00, P90, P90, P90])

# The lay-up and engineering constants can be inspected via:
# print(L_CP1225)

# The laminate is cooled down from a (stress-free) glass transition
# temperature, while no load is applied.
dT = -130
loading_conditions = {'Fx': 0.0, 'Fy': 0.0, 'Fxy': 0.0,
                      'Mx': 0.0, 'My': 0.0, 'Mxy': 0.0, 'dT': dT}
thermal_load = Load(loading_conditions)

# The resulting deformations can be calculated as:
_, d = CP_laminate.loaddef(thermal_load)

# The bending moments to flatten the laminate can be easily computed
# as well. As a first step we define a new loading condition, where we
# force the curvatures equal to zero.
loading_conditions = {'Fx': 0.0, 'Fy': 0.0, 'Fxy': 0.0,
                      'kx': 0.0, 'ky': 0.0, 'kxy': 0.0, 'dT': -130.0}
combined_load = Load(loading_conditions)

# The required bending moments can be extracted from the force vector.
force, _ = CP_laminate.loaddef(combined_load)
moments = force[3:]
