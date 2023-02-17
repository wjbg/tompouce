# Example: 24-ply QI laminate subjected to uniaxial load
#
# A quasi-isotropic fiber reinforced composite test specimen with a
# width of 25 mm that is subjected to a tensile test. We will:
#
# 1. Determine the Young's modulus in X-direction
# 2. Check for failure in case a load of 50 kN is applied
#
from tompouce import Material, Laminate, Load, QI_layup
from tompouce.failure_criteria import Tsai_Hill
import matplotlib.pyplot as plt

# First, we will need an instance of the Material class, which holds
# all relevant material properties. There are two ways to create a
# Material object. We can:
#
# 1. Import data from a json file
# 2. Provide a dictionary with the relevant data
#
# Here, we will import data from a json file.
TC1200 = Material('../materials/TC1200_UD.json')

# We can simply print the object to learn more about the material
# properties:
print(TC1200)

# Next we need to create a Laminate object. Again, there are two ways
# to do so. We can:
#
# 1. Provide a list of ply instances, which is demonstrated in
#    the next example
# 2. Provide a Material object, a list with fiber orientation angles,
#    and a ply thickness (equal for all plies).
#
# Here we use the second option, making use of the function `QI_layup()`
# to generate the orientation angles for a quasi-isotropic layup.
ply_thickness = 0.15E-3
layup = QI_layup(24)
QI_laminate = Laminate(layup, TC1200, ply_thickness)

# In case we would like information about our Laminate objects, such
# as the Young's modulus, we can again simply use the print statement:
print(QI_laminate)

# Next, a Load object is needed. There are two ways to create such an
# object. We can:
#
# 1. Import data from a json file
# 2. Provide a dictionary with the relevant data
#
# Here, we will create a Load object from a dictionary. The dictionary
# should have the following keys:
# - Fx OR ex      : normal load or strain in x-direction
# - Fy OR ey      : normal load or strain in y-direction
# - Fxy OR exy    : shear load or shear strain
# - Mx OR kx      : bending moment or curvature in x-direction
# - My OR ky      : bending moment or curvature in y-direction
# - Mxy OR kxy    : twisting mometn or curvature
# - dT (optional) : temperature difference
F, width = 50E3, 25E-3
loading_conditions = {'Fx': F/width, 'Fy': 0.0, 'Fxy': 0.0,
                      'Mx': 0.0, 'My': 0.0, 'Mxy': 0.0, 'dT': 0.0}
uniaxial_tension = Load(loading_conditions)

# A Load object serves as argument to several Laminate methods. For
# example, we can determine the deformation of the laminate due to the
# applied load.
F, d = QI_laminate.loaddef(uniaxial_tension)

# Alternatively, we can also print this information in a pretty way:
QI_laminate*uniaxial_tension

# We can also plot the stress distribution in the laminate. Here, we
# plot the 1st normal stress in ply CS.
with plt.style.context('ggplot'):
    fig, ax = plt.subplots(figsize=(4.5, 3))
    QI_laminate.plot_stress(uniaxial_tension, comp=0, CS='ply', ax=ax)
    ax.set_xlabel("Stress [MPa]", fontsize=10)
    ax.set_ylabel("Z-location [mm]", fontsize=10)
    ax.set_title("Stress in 1* direction", fontsize=10)
    ax.get_legend().remove()
    fig.tight_layout()
    plt.savefig('../img/uniaxial.png', bbox_inches='tight')

# Lastly, we can check for 1st ply failure according to specified
# failure criterion.
QI_laminate.print_failure(uniaxial_tension, Tsai_Hill)
