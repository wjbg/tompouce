<h1 align="center">
<img src="img/tompouce.svg" width="600">
</h1><br>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

An object-oriented implementation of the Classical Lamination Theory
in Python.

## Example

Consider a C/PEEK drive shaft with a length of 1.25 m, a wall
thickness of 2 mm, and a radius of 50 mm. The shaft has a [45/-45]s
layup and is subjected to a torque of 500 Nm. We will determine the
maximum angular deformation of the shaft.

```python
from tompouce import Material, Laminate, torsion_shaft
from numpy import rad2deg, deg2rad

# Parameters
length = 1.25
radius = 50E-3
torque = 500.0

# Create a Material object created using data from a json file
TC1200 = Material('materials/TC1200_UD.json')

# Create a Laminate object
ply_thickness = 0.5E-3
layup = [45.0, -45.0, 45.0, -45.0]
shaft = Laminate(deg2rad(layup), TC1200, ply_thickness)

# Create a Load object
load = torsion_shaft(torque, radius)

# Print the loads and deformations
shaft*load
```

Which provides a summary of the loads and the resulting deformations:

```
Deformation                               Load        Thermal Load

{ 0.00E+00}     [ a a a | b b b ]    ( { 0.00E+00}    { 0.00E+00} )
{ 0.00E+00}     [ a a a | b b b ]    ( { 0.00E+00}    { 0.00E+00} )
{ 5.28E-04}     [ a a a | b b b ]    ( { 3.18E+04}    { 0.00E+00} )
{---------}  =  [---------------]    ( {---------} +  {---------} )
{ 1.64E-01}     [ b b b | d d d ]    ( { 0.00E+00}    { 0.00E+00} )
{ 1.64E-01}     [ b b b | d d d ]    ( { 0.00E+00}    { 0.00E+00} )
{ 0.00E+00}     [ b b b | d d d ]    ( { 0.00E+00}    { 0.00E+00} )
```

The deformation can also be calculated and then used to determine the
angular deformation:

```python
# Calculate and print the deformation
_, d = shaft.loaddef(load)
angle = d[2]*length/radius
print(f"The angular deformation of the shaft is {rad2deg(angle):.2f} degrees.")
```

```
The angular deformation of the shaft is 0.76 degrees.
```

[Here](examples.md) you can find a small example gallery, which should
give you a pretty good idea what tompouce can do and how things work.

## Install

You can clone the repository to your folder of choice using
[git](https://git-scm.com/downloads):

```
git clone https://github.com/wjbg/tompouce.git
```

Alternatively, you can also download this repository as zip file by
clicking the green button at the top of this page. Tompouce is
reasonably well-documented - if I may say so - and the annotated
examples should be sufficient to get you started.

## Alternatives

There are many alternative implementations of the CLT in Python. One
of these is my own [wjbg/py.CLT](https://github.com/wjbg/py.CLT),
which follows a functional programming paradigm. For more
alternatives, please check out this (non-exhaustive) list:
- [Eacean-CLT](https://github.com/Eacaen/CLT-material-properties)
- [CLamPY](https://github.com/e-dub/CLamPy)
- [pyPLY](https://github.com/Rlee13/pyPLY)
- [Lagerweij-CCLTC](https://github.com/AJJLagerweij/Classical-Composite-Laminate-Theory-Calculator)
- [Compysite](https://github.com/echaffey/Compysite)


## License

Free as defined in the [MIT](https://choosealicense.com/licenses/mit/)
license.
