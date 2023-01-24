# Example: Flexural rigidity of a sandwich beam
#
# Calculate the flexural rigidity of sandwhich beams. All beams have a
# width of 2 cm. Three different face sheet thicknesses are used,
# namely 0.25, 0.5 and 1.0 mm, while the core thickness is veried
# between 1 and 20 mm.
#
import sys; sys.path.append('../')
from tompouce import Material, Laminate, QI_layup, CP_layup
import numpy as np
import matplotlib.pyplot as plt

# List of cutting angles
theta = np.linspace(0, 2*np.pi, 50)

# The skin material properties are imported from a json file
TC1200 = Material('../materials/TC1200_UD.json')
ply_thickness = 0.15E-3

# Create the Laminate objects and store in array
QI_laminate = Laminate(QI_layup(24), TC1200, ply_thickness)
CP_laminate = Laminate(CP_layup(24), TC1200, ply_thickness)

# Calculate results
modulus = np.zeros((2, len(core_thickness)))
modulus[0, 0] = QI_laminate.engineering_constants()['Ex']
modulus[1, 0] = CP_laminate.engineering_constants()['Ex']

dtheta = theta[1] - theta[0]
for i in range(1, len(theta)):
    QI_laminate.layup_rotate(-dtheta)
    CP_laminate.layup_rotate(-dtheta)
    modulus[0, i] = QI_laminate.engineering_constants()['Ex']
    modulus[1, i] = CP_laminate.engineering_constants()['Ex']

# Plot results
with plt.style.context('ggplot'):
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.plot(theta, modulus[0], label="Quasi-Isotropic")
    ax.plot(theta, modulus[1], label="Cross-Ply")
    ax.legend()
    fig.tight_layout()
