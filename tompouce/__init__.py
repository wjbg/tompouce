"""
tompouce: An object-oriented implementation of the CLT
======================================================

Provides a set of objects and functions to analyze the
thermomechanical behavior of a thin laminated plate.


Documentation
-------------

Documentation is available in the docstrings. In addition, annotated
examples are available at: `https://github.com/wjbg/tompouce`.


Available objects
-----------------

Material
  Class used to represent a material.
Ply
  Class used to represent a ply. It takes a Material, a fiber angle
  and a thickness as an input.
Laminate
  Class used to represent a laminate, which consists of several plies.
Load
  Class to represent a load case, which can be applied to a laminate.


Available subpackage
--------------------

failure_criteria
  Collection of failure criteria for lamina.

"""

from tompouce import *
