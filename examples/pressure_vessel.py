# Example: Pressure vessel
#
# Consider a pressure vessel with a radius of 0.3 mm, which has to
# withstand a pressure of 200 Bar. The layup of the pressure vessel is
# [+theta/-theta]_ns. Five winding angles `theta` are considered,
# namely 35, 45, 55, 65 and 75 degrees. Determine which of these
# vessels will fail as a result of the applied pressure.
#
from tompouce import Material, Laminate, pressure_vessel
from tompouce.failure_criteria import Tsai_Hill
import numpy as np


def layup(theta: float, n: int):
    """Function to create [+theta/-theta]_ns layup."""
    if n % 2 == 0:
        el = np.array([theta, -theta])
        half = np.repeat(el, n/2)
        return np.concatenate((half, np.flip(half)))
    else:
        raise ValueError('n should be an even number')


# Parameters
radius = 0.3
pressure = 200E5
fiber_angles = np.deg2rad([35, 45, 55, 65, 75],)
n = 32

# The material properties are imported from a json file
TC1200 = Material('../materials/TC1200_UD.json')
t = 0.15E-3  # ply thickness

# Create the Laminate objects and store these in a list
vessels = [Laminate(layup(angle, n), TC1200, t) for angle in fiber_angles]

# Load object
load = pressure_vessel(pressure, radius)

# Calculate and print results
print("Failure analysis pressure vessels")
print("---------------------------------\n")
for angle, vessel in zip(fiber_angles, vessels):
    fail = vessel.failure(load, Tsai_Hill)
    if any(fail):
        print(f"Vessel with fiber angle of {np.rad2deg(angle):3.1f} failed.")
    else:
        print(f"Vessel with fiber angle of {np.rad2deg(angle):3.1f} survived.")
