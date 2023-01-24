# Example: Flexural rigidity of a sandwich beam
#
# Calculate the flexural rigidity of sandwhich beams. All beams have a
# width of 2 cm. Three different face sheet thicknesses are used,
# namely 0.25, 0.5 and 1.0 mm, while the core thickness is veried
# between 1 and 20 mm.
#
import sys; sys.path.append('../')
from tompouce import Material, Ply, Laminate
import numpy as np
import matplotlib.pyplot as plt

# Core thicknesses and specimen width
core_thickness = np.linspace(1E-3, 2E-2, 20, endpoint=True)
b = 0.02

# The skin material properties are imported from a json file
TC1100 = Material('../materials/TC1100_5HS.json')
ply_thickness = 0.25E-3
ply = Ply(TC1100, 0.0, ply_thickness)

# The core material is generated from a dict
foam = Material({'name': 'Core',
                 'E1': 250E6, 'E2': 400E6,
                 'v12': 0.32, 'G12': 850E6})
core = Ply(foam, 0.0, core_thickness[0])

# Create the Laminate objects and store in array
sandwich = [Laminate([ply, core, ply]),
            Laminate([ply, ply, core, ply, ply]),
            Laminate([ply, ply, ply, ply, core, ply, ply, ply, ply])]

# Calculate results
rigidity = np.zeros((3, len(core_thickness)))
for i, thickness in enumerate(core_thickness):
    core.t = thickness
    rigidity[0, i] = b/sandwich[0].abd()[3, 3]
    rigidity[1, i] = b/sandwich[1].abd()[3, 3]
    rigidity[2, i] = b/sandwich[2].abd()[3, 3]

# Plot results
with plt.style.context('ggplot'):
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.plot(core_thickness*1E3, rigidity[0], label="t = 0.25 mm")
    ax.plot(core_thickness*1E3, rigidity[1], label="t = 0.50 mm")
    ax.plot(core_thickness*1E3, rigidity[2], label="t = 1.00 mm")
    ax.legend()
    ax.set_title('Rigidity as fuction of face and core thickness',
                 fontsize=10)
    ax.set_xlabel('core thickness [mm]', fontsize=10)
    ax.set_ylabel('flexural rigidity [Nm$^2$]', fontsize=10)
    fig.tight_layout()
