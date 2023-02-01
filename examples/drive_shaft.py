# Example: Drive shaft
#
# Consider three drive shafts with a length of 1.25 m. All shafts have
# a [45/-45]_5s layup but different radii, namely 25 mm, 50 mm and 75
# mm. The shaft is loaded with a torque of 500 Nm. Calculate the
# maximum angular deformation for each shaft.
#
import sys; sys.path.append('../')
from tompouce import Material, Laminate, torsion_shaft
from numpy import rad2deg, deg2rad

# Parameters
radii = [25E-3, 50E-3, 75E-3]
torque = 500.0
length = 1.25

# The material properties are imported from a json file
TC1200 = Material('../materials/TC1200_UD.json')

# Create the Laminate object
ply_thickness = 0.15E-3
layup = [45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0,
         -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0]
bias_laminate = Laminate(deg2rad(layup), TC1200, ply_thickness)

# Create the three Load objects, and store these in a list
loads = [torsion_shaft(torque, radius) for radius in radii]

# Calculate and print deformations
print("Resulting deformations\n----------------------")
for load, radius in zip(loads, radii):
    _, d = bias_laminate.loaddef(load)
    angle = d[2]*length/radius
    print(f"Radius: {radius*1E3:4.0f} mm")
    print(f"Angle: {rad2deg(angle):4.2f} deg.\n")
