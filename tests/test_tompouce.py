from pytest import approx
import numpy as np
from tompouce import Material, Ply

ZERO_13 = np.zeros(3)
ZERO_33 = np.zeros((3, 3))
ZERO_66 = np.zeros((6, 6))


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
    pass


class TestLoad:
    pass
