# Example: Quasi-isotropy
#
# Calculate the modulus of a quasi-isotropic and cross-ply laminate as
# a function of in-plane orientation angle.
#
from tompouce import Material, Laminate, QI_layup, CP_layup
import numpy as np
import matplotlib.pyplot as plt

# List of orientation angles
theta = np.linspace(0, 2*np.pi, 200)

# The skin material properties are imported from a json file
TC1200 = Material('../materials/TC1200_UD.json')
ply_thickness = 0.15E-3

# Create the Laminate objects and store in array
QI_laminate = Laminate(QI_layup(24), TC1200, ply_thickness)
CP_laminate = Laminate(CP_layup(24), TC1200, ply_thickness)

# Calculate results
modulus = np.zeros((2, len(theta)))
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
    fig, ax = plt.subplots(figsize=(4.5, 3),
                           subplot_kw={'projection': 'polar'})
    ax.plot(theta, modulus[0]/1E9, label="Quasi-Isotropic")
    ax.plot(theta, modulus[1]/1E9, label="Cross-Ply")
    ax.legend(bbox_to_anchor=(1, 1), bbox_transform=fig.transFigure)
    ax.set_title("Tensile modulus [GPa]", fontsize=10)
    fig.tight_layout()
    plt.savefig('../img/quasi_isotropy.png', bbox_inches='tight')
