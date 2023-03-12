"""Tests for tompouce.

Based on `pytest-7.2.1.` Tests can be run using:

```
$ pytest -vv
```

from the prompt.

"""

from pytest import approx
import numpy as np
from tompouce import Material, Ply, Laminate, Load
from tompouce import torsion_shaft, pressure_vessel
from tompouce.failure_criteria import max_stress

ZERO_13 = np.zeros(3)
ZERO_33 = np.zeros((3, 3))


class TestMaterial:
    """Tests for tompouce Material object."""
    data = {'name': "Test Material - Do not change!",
            'E1': 100E9, 'E2': 10E9, 'v12': 0.25, 'G12': 5E9,
            'alpha1': 1E-6, 'alpha2': 1E-5}

    def test_C(self):
        material = Material(self.data)
        C = np.array([[1.0063E11, 2.5157E9, 0.0],
                      [2.5157E9, 1.0063E10, 0.0],
                      [0.0, 0.00, 5E9]])
        assert material.C() - C == approx(ZERO_33, abs=1E7)

    def test_S(self):
        material = Material(self.data)
        S = np.array([[1E-11, -2.5E-12, 0.0],
                      [-2.5E-12, 1E-10, 0.0],
                      [0.0, 0.0, 2E-10]])
        assert material.S() - S == approx(ZERO_33)

    def test_alpha(self):
        material = Material(self.data)
        alpha = np.array([1E-6, 1E-5, 0.0])
        assert material.alpha() - alpha == approx(ZERO_13)

    def test_C_from_json(self):
        material = Material('test_file.json')
        C = np.array([[1.0063E11, 2.5157E9, 0.0],
                      [2.5157E9, 1.0063E10, 0.0],
                      [0.0, 0.00, 5E9]])
        assert material.C() - C == approx(ZERO_33, abs=1E7)

    def test_failure(self):
        material = Material('test_file.json')
        s0 = np.array([0.0, 0.0, 0.0])
        s1c, s1t = np.array([-600E6, 0.0, 0.0]), np.array([1.1E9, 0.0, 0.0])
        s2c, s2t = np.array([0.0, -110E6, 0.0]), np.array([0.0, 60E6, 0.0])
        s6 = np.array([0.0, 0.0, -60E6])
        assert ((not material.failure(s0, max_stress)) and
                material.failure(s1c, max_stress) and
                material.failure(s1t, max_stress) and
                material.failure(s2c, max_stress) and
                material.failure(s2t, max_stress) and
                material.failure(s6, max_stress))


class TestPly:
    """Tests for tompouce Ply object."""

    data = {'name': "Test Material - Do not change!",
            'E1': 100E9, 'E2': 10E9, 'v12': 0.25, 'G12': 5E9,
            'alpha1': 1E-6, 'alpha2': 1E-5}
    mat = Material(data)
    ply = Ply(mat, np.pi/4, 0.1E-3)

    def test_T(self):
        T = np.array([[0.5, 0.5, 1.0],
                      [0.5, 0.5, -1.0],
                      [-0.5, 0.5, 0.0]])
        assert self.ply._T() - T == approx(ZERO_33)

    def test_R(self):
        R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
        assert self.ply._R - R == approx(ZERO_33)

    def test_C(self):
        C = np.array([[3.393E10, 2.393E10, 2.264E10],
                      [2.393E10, 3.393E10, 2.264E10],
                      [2.264E10, 2.264E10, 2.642E10]])
        assert self.ply.C() - C == approx(ZERO_33, abs=1E7)

    def test_S(self):
        S = np.array([[7.625E-11, -2.375E-11, -4.500E-11],
                      [-2.375E-11, 7.625E-11, -4.500E-11],
                      [-4.500E-11, -4.500E-11, 1.150E-10]])
        assert self.ply.S() - S == approx(ZERO_33, abs=1E-14)

    def test_alpha(self):
        alpha = np.array([5.5E-6, 5.5E-6, -9E-6])
        assert self.ply.alpha() - alpha == approx(ZERO_13, abs=1E-8)

    def test_stress_to_matCS(self):
        ply_stress = np.array([100E6, 10E6, -30E6])
        mat_stress = np.array([25E6, 85E6, -45E6])
        assert (self.ply._stress_to_matCS(ply_stress) - mat_stress ==
                approx(ZERO_13, abs=1))

    def test_strain_to_matCS(self):
        ply_strain = np.array([1E-4, -3E-5, 7E-5])
        mat_strain = np.array([7E-5, 5.781E-21, -1.3E-4])
        assert (self.ply._strain_to_matCS(ply_strain) - mat_strain ==
                approx(ZERO_13, abs=1E-10))

    def test_ply_stress(self):
        ply_strain = np.array([9E-3, -2.5E-4, -8E-3])
        ply_stress = np.array([118.263E6, 25.763E6, -13.208E6])
        assert (self.ply.stress(ply_strain, 0.0) - ply_stress ==
                approx(ZERO_13, abs=1E3))

    def test_ply_strain(self):
        ply_stress = np.array([100E6, 10E6, -30E6])
        ply_strain = np.array([8.738E-3, -2.625E-4, -8.4E-3])
        assert (self.ply.strain(ply_stress, 0.0) - ply_strain ==
                approx(ZERO_13, abs=1E-6))

    def test_ply_stress_T(self):
        T = 100.0
        ply_strain = np.array([9E-3, -2.5E-4, -8E-3])
        ply_stress = np.array([106.816E6, 14.3160E6, -14.3396E6])
        assert (self.ply.stress(ply_strain, T) - ply_stress ==
                approx(ZERO_13, abs=1E3))

    def test_ply_strain_T(self):
        T = 100.0
        ply_stress = np.zeros(3)
        ply_strain = np.array([0.55E-3, 0.55E-3, -0.9E-3])
        assert (self.ply.strain(ply_stress, T) - ply_strain ==
                approx(ZERO_13, abs=1E-6))

    def test_mat_stress(self):
        ply_strain = np.array([9E-3, -2.5E-4, -8E-3])
        mat_stress = np.array([58.805E6, 85.220E6, -46.25E6])
        assert (self.ply.stress(ply_strain, 0.0, 'mat') - mat_stress ==
                approx(ZERO_13, abs=1E3))

    def test_mat_strain(self):
        ply_stress = np.array([100E6, 10E6, -30E6])
        mat_strain = np.array([3.75E-5, 8.438E-3, -9.0E-3])
        assert (self.ply.strain(ply_stress, 0.0, 'mat') - mat_strain ==
                approx(ZERO_13, abs=1E-6))


class TestLaminate:
    """Tests for tompouce Laminate object."""

    data = {'name': "Test Material - Do not change!",
            'E1': 100E9, 'E2': 10E9, 'v12': 0.25, 'G12': 5E9,
            'alpha1': 1E-6, 'alpha2': 1E-5}
    mat = Material(data)
    ply_00 = Ply(mat, 0.0, 0.1E-3)
    ply_90 = Ply(mat, np.pi/2, 0.1E-3)
    ply_45 = Ply(mat, np.pi/4, 0.1E-3)

    def test_init_by_angles(self):
        lam = Laminate([0.0, np.pi/2], self.mat, 0.1E-3)
        assert (lam.layup[0].theta == 0.0 and lam.layup[1].theta == np.pi/2)

    def test_init_by_plies(self):
        lam = Laminate([self.ply_00, self.ply_90])
        assert (lam.layup[0].theta == 0.0 and lam.layup[1].theta == np.pi/2)

    def test_layup_append(self):
        lam = Laminate([self.ply_00, self.ply_90])
        lam.layup_append(self.ply_45)
        assert lam.layup[2].theta == np.pi/4

    def test_layup_insert(self):
        lam = Laminate([self.ply_00, self.ply_90])
        lam.layup_insert(self.ply_45, 1)
        assert lam.layup[1].theta == np.pi/4

    def test_layup_remove(self):
        lam = Laminate([self.ply_00, self.ply_90, self.ply_45])
        lam.layup_remove(1)
        assert lam.layup[1].theta == np.pi/4

    def test_number_of_plies(self):
        lam = Laminate([self.ply_00, self.ply_90, self.ply_45])
        assert lam.number_of_plies() == 3

    def test_thickness(self):
        lam = Laminate([self.ply_00, self.ply_90, self.ply_45])
        assert lam.thickness() - 0.3E-3 == approx(0, abs=1E-10)

    def test_thickness_sandwich(self):
        core = Ply(self.mat, 0.0, 2E-3)
        lam = Laminate([self.ply_00, core, self.ply_00])
        assert lam.thickness() - 2.2E-3 == approx(0, abs=1E-10)

    def test_ply_interfaces(self):
        core = Ply(self.mat, 0.0, 2E-3)
        lam = Laminate([self.ply_00, core, self.ply_00])
        assert (lam._ply_interfaces() - 1E-3*np.array([-1.1, -1, 1, 1.1]) ==
                approx(np.zeros(4), abs=1E-10))

    def test_ABD_A(self):
        lam = Laminate([self.ply_00, self.ply_90])
        A = np.array([[1.1069E7, 5.0315E5, 0.0],
                      [5.0315E5, 1.1069E7, 0.0],
                      [0.0, 0.0, 1E6]])
        assert A - lam.ABD()[:3, :3] == approx(ZERO_33, abs=1E3)

    def test_ABD_B(self):
        lam = Laminate([self.ply_00, self.ply_90])
        B = np.array([[-4.528E2, 0.0, 0.0],
                      [0.0, 4.528E2, 0.0],
                      [0.0, 0.0, 0.0]])
        assert B - lam.ABD()[3:, :3] == approx(ZERO_33, abs=1E-1)

    def test_ABD_D(self):
        lam = Laminate([self.ply_00, self.ply_90])
        D = np.array([[3.6898E-2, 1.6771E-3, 0.0],
                      [1.6771E-3, 3.6898E-2, 0.0],
                      [0.0, 0.0, 3.3333E-3]])
        assert D - lam.ABD()[3:, 3:] == approx(ZERO_33, abs=1E-6)

    def test_abd_a(self):
        lam = Laminate([self.ply_00, self.ply_90])
        a = np.array([[1.8219E-7, -8.2813E-9, 0.0],
                      [-8.2813E-9, 1.8219E-7, 0.0],
                      [0.0, 0.0, 1.0E-6]])
        assert a - lam.abd()[:3, :3] == approx(ZERO_33, abs=1E-11)

    def test_engineering_constants(self):
        lam = Laminate([0.0, np.pi/4, np.pi/4, 0.0], self.mat, 0.1E-3)
        EC = lam.engineering_constants()
        assert (EC['Ex'] - 57.267E9 == approx(0.0, abs=1E6) and
                EC['Ey'] - 13.404E9 == approx(0.0, abs=1E6) and
                EC['Gxy'] - 9.538E9 == approx(0.0, abs=1E6) and
                EC['vxy'] - 0.3660 == approx(0.0, abs=1E-4) and
                EC['vyx'] - 0.0857 == approx(0.0, abs=1E-4) and
                EC['Efx'] - 89.814E9 == approx(0.0, abs=1E6) and
                EC['Efy'] - 11.814E9 == approx(0.0, abs=1E6))

    def test_loaddef_F(self):
        lam = Laminate([0.0, np.pi/2, np.pi/2, 0.0], self.mat, 0.1E-3)
        load = Load({'Fx': 50E3, 'Fy': -3E3, 'Fxy': 14E4,
                     'Mx': 200, 'My': -12, 'Mxy': -2.5, 'dT': 0.0})
        _, d = lam.loaddef(load)
        assert (d[:3] - np.array([2.269E-3, -2.387E-4, 7E-2]) ==
                approx(ZERO_13, abs=1E-6) and
                d[3:] - np.array([4.243E2, -1.551E2, -9.375E1]) ==
                approx(ZERO_13, abs=1E-1))

    def test_loaddef_d(self):
        lam = Laminate([0.0, np.pi/2, np.pi/2, 0.0], self.mat, 0.1E-3)
        load = Load({'ex': 1E-3, 'ey': -2E-3, 'exy': 1E-4,
                     'kx': 10.0, 'ky': -5.2, 'kxy': -2.5, 'dT': 0.0})
        F, _ = lam.loaddef(load)
        assert (F[:3] - np.array([2.013E4, -4.327E4, 2E2]) ==
                approx(ZERO_13, abs=1E1) and
                F[3:] - np.array([4.693, -4.589E-1, -6.667E-2]) ==
                approx(ZERO_13, abs=1E-2))

    def test_loaddef_T(self):
        lam = Laminate([0.0, np.pi/2, np.pi/2, 0.0], self.mat, 0.1E-3)
        load = Load({'Fx': 0.0, 'Fy': 0.0, 'Fxy': 0.0,
                     'Mx': 0.0, 'My': 0.0, 'Mxy': 0.0, 'dT': 50.0})
        _, d = lam.loaddef(load)
        assert (d[:3] - np.array([9.891E-5, 9.891E-5, 0.0]) ==
                approx(ZERO_13, abs=1E-8) and
                d[3:] - np.array([0.0, 0.0, 0.0]) ==
                approx(ZERO_13, abs=1E-12))

    def test_strain(self):
        lam = Laminate([0.0, np.pi/2, np.pi/2, 0.0], self.mat, 0.1E-3)
        load = Load({'Fx': 50E3, 'Fy': -3E3, 'Fxy': 14E4,
                     'Mx': 200, 'My': -12, 'Mxy': -2.5, 'dT': 0.0})
        strain, _ = lam.strain(load)
        assert (strain[0, 0] + 8.258E-2) == approx(0, abs=1E-4)

    def test_stress(self):
        lam = Laminate([0.0, np.pi/2, np.pi/2, 0.0], self.mat, 0.1E-3)
        load = Load({'Fx': 50E3, 'Fy': -3E3, 'Fxy': 14E4,
                     'Mx': 200, 'My': -12, 'Mxy': -2.5, 'dT': 0.0})
        stress, _ = lam.stress(load)
        assert (stress[0, 0] + 8.232E9) == approx(0, abs=1E6)


class TestLoad:
    """Tests for tompouce Load object."""

    def test_shaft(self):
        L = torsion_shaft(100.0, 0.1)
        F = np.array([0.0, 0.0, 1591.549, 0.0, 0.0, 0.0])
        assert (L.F - F) == approx(np.zeros(6), abs=1E-2)

    def test_tank(self):
        L = pressure_vessel(100.0, 0.1)
        F = np.array([5.0, 10.0, 0.0, 0.0, 0.0, 0.0])
        assert (L.F - F) == approx(np.zeros(6), abs=1E-1)
