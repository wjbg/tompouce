# Object-Oriented implementation of the CLT in Python

from __future__ import annotations
import numpy as np
from numpy.linalg import inv, solve
import json
from typing import Union, Optional, Callable
import numpy.typing as npt
from typeguard import typechecked
import matplotlib.pyplot as plt
from failure_criteria import Tsai_Hill


vector = matrix = npt.NDArray[np.float_]


@typechecked
class Material():
    """Class used to represent a material.

    Attributes
    ==========

    The class attributes can be subdivided into three categories,
    namely 1) meta data, 2) thermoelastic properties and 3)
    strength values. The numerical data should be provided in SI
    units.

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

    Note: the compressive and shear strengths should be defined with a
    positive sign.

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
            Check the class docstring to check for valid keys.

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
        with open(fname, 'r') as f:
            data = json.load(f)
            if data['type'] == 'Material':
                self.load_from_dict(data)
            else:
                raise KeyError("not a Load object")

    def save_json(self, fname: str):
        """Saves data to json file."""
        data = vars(self)
        data['type'] = 'Material'
        with open(fname, 'w') as f:
            json.dump(data, f, sort_keys=True, indent=4)

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
            Stress criterion to use. Should be a function that takes
            two inputs, namely the stress state and a Material or
            (data)class with strength data.

        Return
        ------
        failure : bool
            True if material fails.

        """
        return criterion(stress, self)

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

    Public methods
    ==============
    C()
        Returns stiffness matrix in ply CS.
    S()
        Returns compliance matrix in ply CS.
    alpha()
        Returns CTE vector in ply CS.
    stress(strain, dT, CS)
        Returns stress in ply or material CS.
    strain(stress, dT, CS)
        Returns strain in ply or material CS.
    failure(stress, criterion)
        Returns True in case ply failed.

    Private methods
    ===============
    _T()
        Returns transformation matrix.
    _stress_to_matCS(stress)
        Rotates stress from ply to material CS.
    _strain_to_matCS(strain)
        Rotates strain from ply to material CS.
    _wrap_pi(theta)
        Wrap angles between -pi/2 and pi/2

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
        self.theta = Ply._wrap_pi(theta)
        self.t = t

    def _wrap_pi(theta: float) -> float:
        """A clumsy function to wrap angles between -pi/2 < theta <= pi/2."""
        theta = (theta + np.pi) % (2 * np.pi) - np.pi
        if -np.pi/2 < theta <= np.pi/2:
            return theta
        elif np.pi/2 < theta <= np.pi:
            return theta-np.pi
        elif -np.pi <= theta <= -np.pi/2:
            return theta+np.pi

    def _T(self) -> matrix:
        """Returns transformation matrix."""
        n = np.sin(self.theta)
        m = np.cos(self.theta)
        T = np.array([[m**2, n**2,     2*n*m],
                      [n**2, m**2,    -2*n*m],
                      [-m*n,  m*n, m**2-n**2]])
        return T

    def _stress_to_matCS(self, stress: vector) -> vector:
        """Rotates stress in ply CS to material CS.

        Arguments
        ---------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in ply CS.

        Returns
        -------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in material CS.

        """
        return self.T() @ stress

    def _strain_to_matCS(self, strain: vector) -> vector:
        """Rotates stress in ply CS to material CS.

        Arguments
        ---------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in ply CS.

        Returns
        -------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in material CS.

        """
        return Ply._R @ self.T() @ inv(Ply._R) @ strain

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

    def stress(self, strain: vector, dT: float, CS: str = 'ply') -> vector:
        """Returns stress in ply or material CS.

        Arguments
        ---------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in ply CS.
        dT : float
            Temperature difference.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in selected coordinate system.

        """
        stress = self.C() @ (strain - self.alpha()*dT)
        if CS == 'ply':
            return stress
        elif CS == 'mat':
            return self._stress_to_matCS(stress)
        else:
            raise ValueError("expected 'ply' or 'mat' as second argument")

    def strain(self, stress: vector, dT: float, CS: str = 'ply') -> vector:
        """Returns strain in ply or material CS.

        Arguments
        ---------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in ply CS.
        dT : float
            Temperature difference.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        strain : np.ndarray(dim=1, dtype=float)
            Strain vector in selected coordinate system.

        """
        strain = self.S() @ stress + self.alpha()*dT
        if CS == 'ply':
            return strain
        elif CS == 'mat':
            return self._strain_to_matCS(strain)
        else:
            raise ValueError("expected 'ply' or 'mat' as second argument")

    def failure(self, stress: vector, criterion: Callable) -> bool:
        """Check whether failure occurs for given stress and criterion.

        Arguments
        ---------
        stress : np.ndarray(dim=1, dtype=float)
            Stress vector in ply CS.
        criterion : function(stress, Material)
            Stress criterion to use. Should be a function that takes
            two inputs, namely the stress state and a Material or
            (data)class with strength data.

        Return
        ------
        failure : bool
            True if material fails.

        """
        return criterion(self.T() @ stress, self.mat)

    def __str__(self) -> str:
        s = f"{self.mat.name:20s} {self.t*1E3:10.2f} mm" + \
            f"{self.theta*180/np.pi:10.1f} deg"
        return s


@typechecked
class Laminate():
    """Class used to represent a laminate.

    Attributes
    ==========
    layup : list[Ply]
        List of Ply instances.

    Public methods
    ==============

    Classical Lamination Theory
    ---------------------------
    ABD()
        Returns ABD matrix.
    abd()
        Return abd matrix.
    thickness()
        Returns laminate thickness.
    number_of_plies()
        Returns number of plies.
    engineering_constants()
        Return engineering constants for laminate.

    Load
    ----
    loaddef(load)
        Returns force and deformation vector for applied Load.
    loaddef_thermal(load)
        Returns thermal load and accompanying deformation vector.
    strain(load)
        Returns strain on ply interfaces due to load.
    stress(load, CS)
        Returns stresses due to load in material or ply CS.
    plot_stress(load, component, CS, axes)
        Plots through-thickness stress distribution.
    print_stress(load, CS)
        Prints stress state for all plies.
    failure(load, criterion)
        Checks all plies for failure.
    print_failure(load, criterion)
        Prints failure information.

    Layup manipulation
    -----------------
    layup_append(Ply)
        Append Ply to bottom of layup.
    layup_insert(Ply, int)
        Insert Ply at ith position.
    layup_remove(int)
        Remove i-th Ply from layup (defaults to last Ply)
    layup_split(int)
        Split layup at i-th position and return two new Laminates.
    layup_rotate(int)
        Rotate layup with provided angle.

    Private methods
    ===============
    _ply_interfaces()
        Returns z-location of ply interfaces (includes outer surfaces)
    _ply_top_bottom()
        Returns z-location of top and bottom surface of each ply.

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
            Ply thickness (equal for all plies).

        """
        if isinstance(layup[0], Ply):
            self.layup = layup
        elif isinstance(layup[0], float):
            self.layup = [Ply(material, phi, thickness) for phi in layup]

    def thickness(self) -> float:
        """Returns laminate thickness."""
        return sum([ply.t for ply in self.layup])

    def number_of_plies(self) -> int:
        """Returns number of plies."""
        return len(self.layup)

    def layup_append(self, ply: Ply):
        """Appends a Ply object to the layup."""
        if isinstance(ply, Ply):
            self.layup.append(ply)
        else:
            raise TypeError("expected a Ply object")

    def layup_insert(self, ply: Ply, i: int):
        """Inserts a Ply object at the i-th position in the layup."""
        if isinstance(ply, Ply):
            self.layup.insert(i)
        else:
            raise TypeError("expected a Ply object")

    def layup_remove(self, i: int = -1):
        """Removes i-th ply from layup, defaults to last ply."""
        self.layup.pop(i)

    def layup_split(self, i: int) -> tuple[Laminate, Laminate]:
        """Splits laminate before i-th ply and returns both halves."""
        L1, L2 = self.layup[:i], self.layup[i:]
        return Laminate(layup=L1), Laminate(layup=L2)

    def layup_rotate(self, phi: float):
        """Rotates all plies in laminates by angle phi."""
        for ply in self.layup:
            ply.theta += phi

    def ABD(self) -> matrix:
        """Returns ABD matrix."""
        z = self._ply_interfaces()
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

    def engineering_constants(self) -> dict[str, float]:
        """Returns dict with laminate engineering constants."""
        abd = self.abd()
        H = self.thickness()
        engcon = {'Ex': 1/abd[0, 0]/H,
                  'Ey': 1/abd[1, 1]/H,
                  'Gxy': 1/abd[2, 2]/H,
                  'vxy': -abd[1, 0]/abd[0, 0],
                  'vyx': -abd[0, 1]/abd[1, 1],
                  'Efx': 12/abd[3, 3]/H**3,
                  'Efy': 12/abd[4, 4]/H**3}
        return engcon

    def loaddef_thermal(self, load: Load) -> tuple[vector, vector]:
        """Returns thermal load and accompanying deformation vector.

        Arguments
        ---------
        load : Load class
            Load case.

        Returns
        -------
        F_th : np.ndarray(dim=1, dtype=float)
            Thermal force vector.
        d_th : np.ndarray(dim=1, dtype=float)
            Deformation vector due to thermal force.

        """
        z = self._ply_interfaces()
        N = M = np.zeros(3)
        dT = load.dT
        for i, ply in enumerate(self.layup):
            N = N + dT*ply.C() @  ply.alpha() * (z[i+1] - z[i])
            M = M + dT*ply.C() @  ply.alpha() * (z[i+1]**2 - z[i]**2)/2
        F_th = np.block([N, M])
        d_th = self.abd()@F_th
        return F_th, d_th

    def loaddef(self, load: Load) -> tuple[vector, vector]:
        """Returns force and deformation vector for applied Load.

        Arguments
        ----------
        load : Load
            Object from Load class with information about loading conditions.

        Returns
        -------
        F : np.ndarray(dim=1, dtype=float)
            Force vector.
        d : np.ndarray(dim=1, dtype=float)
            Deformation vector.

        """
        F = load.F
        d = load.d
        iF, iD = load._masks()
        _, d_thermal = self.loaddef_thermal(load)
        if all(iF):
            d = self.abd()@F + d_thermal
        elif all(iD):
            F = self.ABD()@(d - d_thermal)
        else:
            abd = self.abd()
            aap = (abd[:, iF][iD])@F[iF]
            noot = d[iD] - aap - d_thermal[iD]
            F[iD] = solve(abd[:, iD][iD], noot)
            d = abd@F + d_thermal
        return F, d

    def strain(self, load: Load, CS: str = 'ply') -> matrix:
        """Returns strain on top and bottom surface of each ply due to Load.

        Arguments
        ----------
        load : Load
            Object from Load class with information about loading conditions.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        epsilon : np.ndarray(dim=2, dtype=float)
            Strains in selected coordinate system. The returned
            matrix is of shape (2*N, 3), with N the number of plies.
            For each ply the stress state at its top and bottom are
            returned. For the ith ply:

            - column i*2 holds the stress state at its top
            - column i*2+1 holds the stress state at its bottom

            The three rows correspond to the three in-plane strains,
            i.e. the two normal strains and the shear strain.
        z_int : np.ndarray(dim=1, dtype=float)
            Array with z-locations of the corresponding strains.

        """
        _, d = self.loaddef(load)
        z_int = self._ply_top_bottom()
        strain = d[:3][..., None] + z_int*d[3:][..., None]
        if CS == 'mat':
            for i, ply in enumerate(self.layup):
                strain[:, i*2] = ply._rotate_to_matCS(strain[:, i*2])
                strain[:, i*2+1] = ply._rotate_to_matCS(strain[:, i*2+1])
        return strain, z_int

    def stress(self, load: Load, CS: str = 'ply') -> matrix:
        """Returns stress on top and bottom for each ply due to Load.

        Argument
        --------
        load : Load
            Object from Load class with information about loading conditions.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        Returns
        -------
        stress : np.ndarray(dim=2, dtype=float)
            Stresses in selected coordinate system. The returned
            matrix is of shape (2*N, 3), with N the number of plies.
            For each ply the stress state at its top and bottom are
            returned. For the ith ply:

            - column i*2 holds the stress state at its top
            - column i*2+1 holds the stress state at its bottom

            The three rows correspond to the three in-plane stresses,
            i.e. the two normal stresses and the shear stress.
        z_int : np.ndarray(dim=1, dtype=float)
            Array with z-locations of the corresponding strains.

        """
        dT = load.dT
        strain, z_int = self.strain(load)
        stress = np.zeros((3, 2*len(self.layup)))
        for i, ply in enumerate(self.layup):
            stress[:, i*2] = ply.stress(strain[:, i*2], dT, CS)
            stress[:, i*2+1] = ply.stress(strain[:, i*2+1], dT, CS)
        return stress, z_int

    def plot_stress(self, load: Load, comp: int = 0, CS: str = 'ply',
                    ax: Optional[plt.Axes] = None) -> plt.Axes:
        """Plots stress.

        Arguments
        ---------
        load : Load
            Object from Load class with information about loading conditions.
        comp : int (0, 1, 2)
            Stress component to plot:
            - 0 : 1 or 1* component
            - 1 : 2 or 2* component
            - 2 : shear stress
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.
        ax : plt.Axes (Optional)
            Axes for plotting.

        Returns
        -------
        ax : plt.Axes (Optional)
            Axes for plotting.

        """
        if ax is not None:
            fig, ax = plt.subplots()
        ax.plot(self.stress(CS)[comp], self._ply_top_bottom())
        plt.show()

    def print_stress(self, load: Load, CS: str = 'ply'):
        """Prints stress state.

        Arguments
        ---------
        load : Load
            Object from Load class with information about loading conditions.
        CS : str - 'ply' or 'mat' (defaults to 'ply')
            Coordinate system to use.

        """
        stress = self.stress(CS)/1E6
        header = f"Stresses in {CS} CS [MPa]\n" + \
            f"{'Ply #':10s}{'Top':15s}{'Bottom':15s}\n" + \
            "----------------------------------------\n"
        data = ""
        for i in range(self.number_of_plies):
            data += f"{i:10f}{stress[0,i*2]:15.2f}{stress[0,i*2+1]:15.2f}\n"
            data += f"{'':10f}{stress[1,i*2]:15.2f}{stress[1,i*2+1]:15.2f}\n"
            data += f"{'':10f}{stress[2,i*2]:15.2f}{stress[2,i*2+1]:15.2f}\n"
        return header + data

    def failure(self, load: Load,
                criterion: Callable = Tsai_Hill) -> list[bool]:
        """Checks all plies for failure.

        Arguments
        ---------
        load : Load
            Object from Load class with information about loading conditios.
        criterion : Callable
            Failure criterion to use.

        Returns
        -------
        failed : list of bools (True indicates failure)
            For the ith ply:
            - column i*2 represents its top surface
            - column i*2+1 represents its bottom surface

        """
        stress = self.stress(load)
        failed = []*len(stress)
        for i, ply in enumerate(self.layup):
            failed[i*2] = ply.failure(stress[:, i*2], criterion)
            failed[i*2+1] = ply.failure(stress[:, i*2+1], criterion)
        return failed

    def print_failure(self, load: Load, criterion: Callable = Tsai_Hill):
        """Prints failure information.

        Arguments
        ---------
        load : Load
            Object from Load class with information about loading conditios.
        criterion : Callable
            Failure criterion to use.

        """
        def txt(f): return "Failed" if f else "OK!"
        failed = self.failure(load, criterion)
        header = "Ply failure overview\n" + \
            f"{'Ply #':10s}{'Top':15s}{'Bottom':15s}\n" + \
            "----------------------------------------\n"
        data = ""
        for i in range(self.number_of_plies):
            data += f"{i:10f}{txt(failed[i*2]):15s}{txt(failed[i*2+1]):15s}\n"
        print(header + data)

    def _ply_interfaces(self) -> vector:
        """Returns location of ply interfaces; includes outer surfaces."""
        H, N = self.thickness(), len(self.layup)
        z = np.linspace(-H/2, H/2, N+1)
        return z

    def _ply_top_bottom(self) -> vector:
        """Returns z-coordinates for the top and bottom surface of each ply."""
        z = self.laminate._ply_interfaces()
        return np.repeat(z, 2)[1:-1]

    def __str__(self) -> str:
        s = "\n" + \
           f"{'Material':24s} {'Thickness':11s} {'Orientation':10s}\n" + \
            "------------------------------------------------ Top\n"
        for ply in self.layup:
            s += ply.__str__() + "\n"
        s += "------------------------------------------------ Bottom\n" + \
            f"# plies: {self.number_of_plies():<3}" +\
            f" {1000*self.thickness():18.2f} mm\n"
        s += "\nEngineering constants\n---------------------\n"
        eng = self.engineering_constants()
        s += f"Young's Modulus X:   {eng['Ex']/1E9:4.1f} GPa\n" + \
             f"Young's Modulus Y:   {eng['Ey']/1E9:4.1f} GPa\n" + \
             f"Shear Modulus X:     {eng['Gxy']/1E9:4.1f} GPa\n" + \
             f"Poisson's ratio XY:  {eng['vxy']:4.2f}\n" + \
             f"Poisson's ratio YX:  {eng['vyx']:4.2f}\n" + \
             f"Flexural modulus X:  {eng['Efx']/1E9:4.1f} GPa\n" + \
             f"Flexural modulus Y:  {eng['Efy']/1E9:4.1f} GPa\n"
        return s


class Load():
    """Class to represent a loading condition.

    Attributes
    ==========
    F : np.ndarray(dim=1, dtype=float)
        Force vector.
    d : np.ndarray(dim=1, dtype=float)
        Deformation vector.
    dT : float
        Temperature difference.

    The force and deformation vectors have six elements:

    0) normal force (F) and strain (d) in X-direction
    1) normal force (F) and strain (d) in Y-direction
    2) shear force (F) and shear strain (d)
    3) bending moment (F) and curvature (d) in X-direction
    4) bending moment (F) and curvature (d) in Y-direction
    5) twisting moment (F) and twisting curvature (d)

    A valid load condition requires that only a force (or a moment) OR
    a strain (or a curvature) is imposed in a given direction. As an
    example, in case a normal force is applied in X-direction, the
    corresponding strain cannot be provided and must be calculated.
    The first element in the force vector F equals the imposed force,
    while the corresponding element in the deformation vector d should
    then be np.nan.

    Public methods
    ==============
    load_from_dict(data)
        Loads data from dictionary and sets keys as ojbect's attributes.
    load_json(fname)
        Loads data from json file.
    save_json(fname)
        Saves data to json file.
    valid_load()
        Returns True if load condition is valid.

    Private methods
    ===============
    _masks()
        Returns a mask for the non-NaNs in self.F and self.D.

    """

    _labels = [('Fx', 'ex'), ('Fy', 'ey'), ('Fxy', 'exy'),
               ('Mx', 'kx'), ('My', 'ky'), ('Mxy', 'kxy')]

    def __init__(self):
        self.F = np.nan * np.ones(6)
        self.d = np.nan * np.ones(6)
        self.dT = 0.0

    def load_from_dict(self, load):
        """Load data from dictionary.

        Argument
        --------
        load : dict
            Dictionary with the following keys:

            - Fx OR ex      : normal load or strain in x-direction
            - Fy OR ey      : normal load or strain in y-direction
            - Fxy OR exy    : shear load or shear strain
            - Mx OR kx      : bending moment or curvature in x-direction
            - My OR ky      : bending moment or curvature in y-direction
            - Mxy OR kxy    : twisting mometn or curvature
            - dT (optional) : temperature difference

        """
        for i, label in enumerate(Load._labels):
            if label[0] in load & label[1] in load:
                self.F[0] = load[label[0]]
            elif label[0] not in load & label[1] in load:
                self.d[0] = load[label[1]]
            elif label[0] in load & label[1] in load:
                raise KeyError("overdefined problem:" +
                               f" both {label[0]} and {label[1]} provided")
            elif label[0] not in load & label[1] not in load:
                raise KeyError("underdefined problem:" +
                               f" no {label[0]} or {label[1]} provided")
        self.dT = load['dT'] if 'dT' in load else 20.0

    def _masks(self) -> tuple[vector]:
        """Returns masks with 1 for every non-NaN in F and D"""
        iF = 1 * ~np.isnan(self.F)
        iD = 1 * ~np.isnan(self.d)
        return iF, iD

    def save_json(self, fname):
        """Saves data to json file."""
        data = vars(self)
        data['type'] = 'Load'
        with open(fname, 'w') as f:
            json.dump(data, f, sort_keys=True, indent=4)

    def load_json(self, fname):
        """Loads data from json file."""
        with open(fname, 'r') as f:
            data = json.load(f)
            if data['type'] == 'Load':
                self.load_from_dict(data)
            else:
                raise KeyError("not a Load object")

    def valid_load(self):
        """Returns True if loading condition is valid."""
        iF, iD = self._masks()
        if all(iF + iD == 1):
            return True
        else:
            return False

    def __str__(self):
        s = f"{'Force':>16s}      {'Deformation':12s}\n" + \
            "--------------------------------------\n"
        for i in range(3):
            s += f"{self.F[i]/1E3:>11.2f} kN/m {self.d[i]:8.2f} m/m\n"
        for i in range(3, 6):
            s += f"{self.F[i]/1E3:>10.2f} kNm/m  {self.d[i]:7.2f} 1/m\n"
        s += f"\nDelta T: {self.dT:5.1f} C\n"
        return s.replace("nan", "???")


@typechecked
def pressure_vessel(P: float, R: float) -> Load:
    """Returns Load object for a pressure vessel.

    Arguments
    ---------
    P : float
        Internal pressure.
    R : float
        Radius of the pressure vessel.

    Returns
    -------
    load : Load object
        Loading condition.

    """
    load = Load()
    load.F = np.zeros(6)
    load.F[0] = P*R/2
    load.F[1] = R*R
    return load


@typechecked
def torsion_shaft(T: float, R: float) -> Load:
    """Returns Load object for a torsion shaft.

    Arguments
    ---------
    T : float
        Torque.
    R : float
        Radius.

    Returns
    -------
    load : Load object
        Loading condition.

    """
    load = Load()
    load.F = np.zeros(6)
    load.F[2] = T/(2*np.pi*R**2)
    return load


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


TC1200 = Material('materials/TC1200UD.json')
layup = QI_layup(8)
L = Laminate(layup=layup, material=TC1200, thickness=0.15E-3)
F = pressure_vessel(1E6, 0.4)
