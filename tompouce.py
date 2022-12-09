# Object-Oriented implementation of the CLT in Python

from __future__ import annotations
import numpy as np
from numpy.linalg import inv
import json
from typing import Union, Optional, Callable
import numpy.typing as npt
from typeguard import typechecked

vector = matrix = npt.NDArray[np.float_]


@typechecked
class Material():
    """Class used to represent a material.

    Attributes
    ==========

    The class attributes can be subdivided into three categories,
    namely i) meta data, 2) thermoelastic properties and 3)
    strengths. The numerical data should be provided in SI units.

    Meta data
    ---------
    name : str
    manufacturer : str
    matrix : str
    fiber : str

    Thermoelastic properties
    ------------------------
    E1 : float
    E2 : float
    G12 : float
    v12 : float
    alpha1 : float
    alpha2 : float

    Strengths
    ---------
    S1c : float
    S1t : float
    S2c : float
    S2t : float
    S6 : float (shear strength)

    Methods
    =======
    load_from_dict(data)
        Loads data from dictionary and sets keys as ojbect's attributes.
    load_json(fname)
        Loads data from json file.
    save_json(fname)
        Saves data to json file.
    C()
        Returns stiffness matrix in material CS.
    S()
        Returns compliance matrix in material CS.
    alpha()
        Returns CTE vector in material CS.

    """

    def __init__(self, inp: Union[str, dict]):
        """Initializes Material object.

        Parameters
        ----------
        inp : str or dict
            Data is read either from a json file or from a dictionary.
            In the former case, inp should be a string with the
            filename, while in the latter case inp should be a
            dictionary whose keys will be the object's attributes.

        """
        if isinstance(inp, dict):
            self.load_from_dict(inp)
        if isinstance(inp, str):
            self.load_json(inp)

    def load_from_dict(self, data: dict[str, str | float]):
        """Loads data from dictionary."""
        for key, value in data.items():
            setattr(self, key, value)

    def load_json(self, fname: str):
        """Loads data from json file."""
        with open(fname, 'r') as fn:
            self.load_from_dict(json.load(fn))

    def save_json(self, fname: str):
        """Saves data to json file."""
        data = vars(self)
        with open(fname, 'w') as fn:
            json.dump(data, fn, sort_keys=True, indent=4)

    def C(self) -> vector:
        """Returns stiffness matrix."""
        E1, E2, v12, G12 = self.E1, self.E2, self.v12, self.G12
        v21 = E2*v12/E1
        C = np.array([[E1/(1-v12*v21), v21*E1/(1-v12*v21), 0],
                      [v21*E1/(1-v12*v21), E2/(1-v12*v21), 0],
                      [0, 0, G12]])
        return C

    def alpha(self) -> vector:
        """Returns CTE vector."""
        return np.array([self.alpha1, self.alpha2, 0.0])

    def S(self) -> matrix:
        """Returns compliance matrix."""
        C = self.C()
        return inv(C)

    def failure(self, stress: vector, criterion: Callable) -> bool:
        """Check whether failure occurs for given stress and criterion.

        Arguments
        ---------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in material CS.
        criterion : function(stress, Material)
            Stress criterion to use.

        Return
        ------
        failure : bool
            True if material fails.

        """
        pass

    def __str__(self) -> str:
        s = (f"\n{self.name}\n" +
             "-----------------------------\n" +
             f"Manufacturer:   {self.manufacturer}\n" +
             f"Matrix:         {self.matrix}\n" +
             f"Fiber:          {self.matrix}\n" +
             "\n" +
             "Thermoelastic properties\n" +
             "-----------------------------\n" +
             f"E1:             {self.E1/1E9:7.1f} GPa\n" +
             f"E2:             {self.E2/1E9:7.1f} GPa\n" +
             f"G12:            {self.G12/1E9:7.1f} GPa\n" +
             f"v12:            {self.v12:7.2f}\n" +
             f"alpha1:         {self.alpha1*1E6:7.2f} um/m\n" +
             f"alpha2:         {self.alpha2*1E6:7.2f} um/m\n" +
             "\n" +
             "Strength\n" +
             "-----------------------------\n" +
             f"S1t:            {self.S1t/1E6:7.1f} MPa\n" +
             f"S1c:            {self.S1c/1E6:7.1f} MPa\n" +
             f"S2t:            {self.S2t/1E6:7.1f} MPa\n" +
             f"S2c:            {self.S2c/1E6:7.1f} MPa\n" +
             f"S6:             {self.S6/1E6:7.1f} MPa\n")
        return s


TC1200 = Material('materials/TC1200UD.json')


@typechecked
class Ply():
    """Class used to represent a ply.

    Attributes
    ==========
    mat : Material
        Ply material.
    theta : float
        Fiber orientation angle.
    t : float
        Ply thickness.

    Methods
    =======
    C()
        Returns stiffness matrix in ply CS.
    S()
        Returns compliance matrix in ply CS.
    alpha()
        Returns CTE vector in ply CS.

    """
    _R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])  # Reuter matrix

    def __init__(self, material: Material, theta: float, t: float):
        """Initializes Ply instance.

        Parameters
        ----------
        material : Material object
            Instance from Material class.
        theta : float
            Fiber orientation angle.
        t : float
            Ply thickness.

        """
        self.mat = material
        self.theta = theta
        self.t = t

    def _T(self) -> matrix:
        """Returns transformation matrix."""
        n = np.sin(self.theta)
        m = np.cos(self.theta)
        T = np.array([[m**2, n**2,     2*n*m],
                      [n**2, m**2,    -2*n*m],
                      [-m*n,  m*n, m**2-n**2]])
        return T

    def C(self) -> matrix:
        """Returns stiffness matrix."""
        T = self._T()
        return inv(T) @ self.mat.C() @ Ply._R @ T @ inv(Ply._R)

    def S(self) -> matrix:
        """Returns compliance matrix."""
        return inv(self.C())

    def alpha(self) -> vector:
        """Returns CTE vector."""
        T = self._T()
        return Ply._R @ inv(T) @ inv(Ply._R) @ self.mat.alpha()

    def stress(self, strain: vector, CS: str = 'ply') -> vector:
        """Returns stress in ply or material CS.

        Arguments
        ---------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in ply CS.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in selected coordinate system.

        """
        if CS == 'ply':
            return self.C() @ strain
        elif CS == 'mat':
            return self.T() @ self.C() @ strain
        else:
            raise ValueError("expected 'ply' or 'mat' as second argument")

    def strain(self, stress: vector, CS: str = 'ply') -> vector:
        """Returns strain in ply or material CS.

        Arguments
        ---------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in ply CS.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in selected coordinate system.

        """
        if CS == 'ply':
            return self.S() @ stress
        elif CS == 'mat':
            return Ply._R @ self.T() @ inv(Ply._R) @ self.S() @ stress
        else:
            raise ValueError("expected 'ply' or 'mat' as second argument")

    def __str__(self) -> str:
        s = f"{self.mat.name:20s} {self.t*1E3:10.2f} mm" + \
            f"{self.theta*180/np.pi:10.1f} deg"
        return s


P = Ply(TC1200, np.pi/2, 0.15E-3)


@typechecked
class Laminate():
    """Class used to represent a laminate.

    Attributes
    ==========
    layup : list[Ply]
        List of Ply instances.

    Methods
    =======

    Classical Lamination Theory
    ---------------------------
    ABD()
        Returns ABD matrix.
    abd()
        Return abd matrix.
    H()
        Returns laminate thickness.
    engineering_constants()
        Return engineering constants for laminate.

    Layup manipulation
    -----------------
    append(Ply)
        Append Ply to bottom of layup.
    insert(Ply, int)
        Insert Ply at ith position.
    remove(int)
        Remove i-th Ply from layup (defaults to last Ply)
    split(int)
        Split layup at i-th position and return two new Laminates.
    rotate(int)
        Rotate layup with provided angle.

    """

    def __init__(self, layup: list[float] | list[Ply],
                 material: Optional[Material] = None,
                 thickness: Optional[float] = None):
        """Initializes Laminate object.

        Parameters
        ----------
        layup : list[Ply] or list[float]
            List of Ply objects or a list with orientation angles.
        material : Material (required when layup is a list[float])
            Material used in the layup.
        thickness : float (required when layup is a list[float])
            Ply thickness.

        """
        if isinstance(layup[0], Ply):
            self.layup = layup
        elif isinstance(layup[0], float):
            self.layup = [Ply(material, phi, thickness) for phi in layup]

    def ABD(self) -> matrix:
        """Returns ABD matrix."""
        z = self._ply_edges()
        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))
        for i, ply in enumerate(self.layup):
            A = A + ply.C() * (z[i+1] - z[i])
            B = B + ply.C() * (z[i+1]**2 - z[i]**2)/2
            D = D + ply.C() * (z[i+1]**3 - z[i]**3)/3
        ABD = np.block([[A, B], [B, D]])
        return ABD

    def abd(self) -> matrix:
        """Returns abd matrix."""
        return inv(self.ABD())

    def H(self) -> float:
        """Returns laminate thickness."""
        return np.sum([ply.t for ply in self.layup])

    def _ply_edges(self) -> vector:
        """Returns location of ply edges"""
        H, N = self.H(), len(self.layup)
        z = np.linspace(-H/2, H/2, N+1)
        return z

    def engineering_constants(self) -> dict[str, float]:
        """Returns laminate engineering constants"""
        abd = self.abd()
        H = self.H()
        engcon = {'Ex': 1/abd[0, 0]/H,
                  'Ey': 1/abd[1, 1]/H,
                  'Gxy': 1/abd[2, 2]/H,
                  'vxy': -abd[1, 0]/abd[0, 0],
                  'vyx': -abd[0, 1]/abd[1, 1],
                  'Efx': 12/abd[3, 3]/H**3,
                  'Efy': 12/abd[4, 4]/H**3}
        return engcon

    def append(self, ply: Ply):
        """Appends a Ply object to the layup."""
        if isinstance(ply, Ply):
            self.layup.append(ply)
        else:
            raise TypeError("expected a Ply object")

    def insert(self, ply: Ply, i: int):
        """Inserts a Ply object at the i-th position in the layup."""
        if isinstance(ply, Ply):
            self.layup.insert(i)
        else:
            raise TypeError("expected a Ply object")

    def remove(self, i: int = -1):
        """Removes i-th ply from layup, defaults to last ply."""
        self.layup.pop(i)

    def split(self, i: int) -> tuple[Laminate, Laminate]:
        """Splits laminate before i-th ply and returns both halves."""
        L1, L2 = self.layup[:i], self.layup[i:]
        return Laminate(layup=L1), Laminate(layup=L2)

    def rotate(self, phi: float):
        """Rotates all plies in laminates by angle phi."""
        for ply in self.layup:
            ply.theta += phi


class Loads():
    pass


@typechecked
def CP_layup(N: int) -> list:
    """Returns layup for symmetric cross-ply laminate.

    Parameters
    ----------
    N : int
        Number of plies.

    Returns
    -------
    layup : list
        List with layup angles in radians.

    """
    if N % 4 == 0:
        unit = [0.0, np.pi/2]
        k = int(N/4)
        half = unit*k
        layup = half[:] + half[::-1]
    else:
        raise ValueError("The number of plies should be a multiple of 4.")
    return layup


@typechecked
def QI_layup(N: int) -> list:
    """Returns layup for symmetric quasi-isotropic laminate.

    Parameters
    ----------
    N : int
        Number of plies.

    Returns
    -------
    layup : list
        List with layup angles in radians.

    """
    if N % 8 == 0:
        unit = [np.pi/4, 0.0, -np.pi/4, np.pi/2]
        k = int(N/8)
        half = unit*k
        layup = half[:] + half[::-1]
    else:
        raise ValueError("The number of plies should be a multiple of 8.")
    return layup
