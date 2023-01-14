"""Failure criteria for lamina in plane stress.

Set of functions with failure criteria for composite lamina. The
functions take (at least) two arguments:

1. The in-plane stress state which should be a NDarray of length 3
   with the longitudinal and the transverse normal stresses and the
   in-plane shear stress.
2. An instance of the Material class or any other class with the
   lamina strength data stored in attributes. Most functions use the
   following data:
   - S1t : Longitudinal tensile strength
   - S1c : Longitudinal compressive strength
   - S2t : Transverse tensile strength
   - S2c : Transverse compressive strength
   - S6  : Shear strength
   Some criteria may require different or additional strength data.

All functions provide a boolean as output with True in case the lamina
failed and False in case the lamina did not fail.

So far, the following failure criteria have been implemented:

- Tsai-Hill
- Norris
- Tsai-Wu
- Hoffman
- Hashin-Rotem
- Hashin
- Maximum stress criterion
- Maximum strain criterion

"""


import numpy as np
import numpy.typing as npt
from typeguard import typechecked
from tompouce import Material


vector = matrix = npt.NDArray[np.float_]


@typechecked
def Tsai_Hill(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to Tsai-Hill criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength

    Note: the compressive and shear strengths should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    S1 = mat.S1c if stress[0] < 0.0 else mat.S1t
    S2 = mat.S2c if stress[1] < 0.0 else mat.S2t
    S6 = mat.S6
    TH = (stress[0]**2/S1**2 - stress[0]*stress[1]/S1**2 +
          stress[1]**2/S2**2 + stress[2]**2/S6**2)
    failed = True if TH >= 1.0 else False
    return failed


@typechecked
def Norris(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to Norris criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength

    Note: the compressive and shear strengths should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    S1 = mat.S1c if stress[0] < 0.0 else mat.S1t
    S2 = mat.S2c if stress[1] < 0.0 else mat.S2t
    S6 = mat.S6
    TH = (stress[0]**2/S1**2 - (stress[0]/S1)*(stress[1]/S2) +
          stress[1]**2/S2**2 + stress[2]**2/S6**2)
    failed = True if TH >= 1.0 else False
    return failed


@typechecked
def Tsai_Wu(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to Tsai-Wu criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength (it is assumed that S6c = S6t)

    Note: the compressive and shear strengths should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    F1 = 1/mat.S1t - 1/mat.S1c
    F2 = 1/mat.S2t - 1/mat.S2c
    F11 = 1/(mat.S1t*mat.S1c)
    F22 = 1/(mat.S2t*mat.S2c)
    F12 = -np.sqrt(F11 * F22)/2
    F66 = 1/(mat.S6**2)
    TW = (F1*stress[0] + F2*stress[1] + 2*F12*stress[0]*stress[1] +
          F11*stress[0]**2 + F22*stress[1]**2 + F66*stress[2]**2)
    failed = True if TW >= 1.0 else False
    return failed


@typechecked
def Hoffman(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to Hoffman criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength

    Note: the compressive and shear strengths should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    F1 = 1/mat.S1t - 1/mat.S1c
    F2 = 1/mat.S2t - 1/mat.S2c
    F11 = 1/(mat.S1t*mat.S1c)
    F22 = 1/(mat.S2t*mat.S2c)
    F66 = 1/(mat.S6**2)
    F12 = -1/(mat.S1t*mat.S1c)
    H = (F1*stress[0] + F2*stress[1] + F11*stress[0]**2 + F22*stress[1]**2 +
         F66*stress[2]**2 + F12*stress[0]*stress[1])
    failed = True if H >= 1.0 else False
    return failed


@typechecked
def Hashin_Rotem(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to Hashin-Rotem criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength

    Note: the compressive and shear strengths should be defined with a
    positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    if stress[0] >= 0.0:
        HR1 = (stress[0]/mat.S1t)**2
    elif stress[0] < 0.0:
        HR1 = (stress[0]/mat.S1c)**2
    if stress[1] >= 0.0:
        HR2 = (stress[1]/mat.S2t)**2 + (stress[2]/mat.S6)**2
    elif stress[1] < 0.0:
        HR2 = (stress[1]/mat.S2c)**2 + (stress[2]/mat.S6)**2
    failed = True if HR1 >= 1.0 | HR2 >= 1.0 else False
    return failed


@typechecked
def Hashin(stress: vector, mat: Material, alpha: float = 1.0) -> bool:
    """Checks if stress causes failure according to Hashin criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A material object with the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S4  : transverse (2-3) shear strength
        - S6  : in-plane shear strength
    alpha : float (defaults to 1.0)
        Coefficient that determines the contribution of the in-plane shear
        stress to fiber tensile failure. Allowable values: 0 <= alpha <= 1.0.

    Note: the compressive and shear strengths should be defined with a
    positive sign.

    Note 2: The Hashin failure criterion requires the transverse shear
    strength (S4) of the ply.

    Returns
    -------
    failed : bool
        True if failed.

    """
    if stress[0] >= 0.0:
        H1 = (stress[0]/mat.S1t)**2 + alpha*(stress[2]/mat.S6)**2
    elif stress[0] < 0.0:
        H1 = (stress[0]/mat.S1c)**2
    if stress[1] >= 0.0:
        H2 = (stress[1]/mat.S2t)**2 + (stress[2]/mat.S6)**2
    elif stress[1] < 0.0:
        H2 = ((stress[2]/mat.S6)**2 + (stress[1]/(2*mat.S4))**2 +
              ((mat.S2c/(2*mat.S4))**2 - 1) * (stress[1]/mat.S2c))
    failed = True if H1 >= 1.0 | H2 >= 1.0 else False
    return failed


@typechecked
def max_stress(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to max stress criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A Material object the following strength data:
        - S1t : tensile strength in 1-direction
        - S1c : compressive strength in 1-direction
        - S2t : tensile strength in 2-direction
        - S2c : compressive strength in 2-direction
        - S6  : in-plane shear strength

    Note: the compressive and shear strengths should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    failed = False
    if stress[0] <= -abs(mat.S1c) | stress[0] >= mat.S1t:
        failed = True
    if stress[1] <= -abs(mat.S2c) | stress[1] >= mat.S2t:
        failed = True
    if abs(stress[2]) >= abs(mat.S6):
        failed = True
    return failed


@typechecked
def max_strain(stress: vector, mat: Material) -> bool:
    """Checks if stress causes failure according to max stress criterion.

    Arguments
    ---------
    stress : NDArray(dtype=float, dim=1)
        Array with the in-plane stress components (1, 2, 6) in material CS.
    mat : Material class
        A Material object the following strength data:
        - E1t : ultimate tensile strain in 1-direction
        - E1c : ultimate compressive strain in 1-direction
        - E2t : ultimate tensile strain in 2-direction
        - E2c : ultimate compressive strain in 2-direction
        - E6  : ultimate in-plane shear strain

    Note: the compressive and shear strain limits should be defined
    with a positive sign.

    Returns
    -------
    failed : bool
        True if failed.

    """
    strain = mat.C() @ stress
    failed = False
    if strain[0] <= -abs(mat.E1c) | strain[0] >= mat.E1t:
        failed = True
    if strain[1] <= -abs(mat.E2c) | strain[1] >= mat.E2t:
        failed = True
    if abs(strain[2]) >= abs(mat.E6):
        failed = True
    return failed
