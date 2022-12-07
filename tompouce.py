# Object-Oriented implementation of the CLT in Python

import numpy as np
from numpy.linalg import inv
import json


class Material():

    def __init__(self, **kw):
        if 'data' in kw:
            self.load_from_dict(kw['data'])
        if 'fname' in kw:
            self.load_json(kw['fname'])

    def load_from_dict(self, props):
        for key, value in props.items():
            setattr(self, key, value)

    def save_json(self, fname):
        data = vars(self)
        with open(fname, 'w') as fn:
            json.dump(data, fn, sort_keys=True, indent=4)

    def load_json(self, fname):
        with open(fname, 'r') as fn:
            self.load_from_dict(json.load(fn))

    def C(self):
        """Returns stiffness matrix."""
        try:
            E1, E2, v12, G12 = self.E1, self.E2, self.v12, self.G12
            v21 = E2*v12/E1
            C = np.array([[E1/(1-v12*v21), v21*E1/(1-v12*v21), 0],
                          [v21*E1/(1-v12*v21), E2/(1-v12*v21), 0],
                          [0, 0, G12]])
        except Exception:
            print("Check whether the elastic constants are provided.")
            C = None
        return C

    def alpha(self):
        """Returns CTE vector."""
        try:
            alpha = np.array([self.alpha1, self.alpha2, 0.0])
        except Exception:
            print("Check whether the CTEs are provided.")
            alpha = None
        return alpha

    def S(self):
        """Returns compliance matrix."""
        C = self.C()
        return inv(C)

    def __str__(self):
        s = (f"\n{self.name}\n" +
             "-------------------------------------\n" +
             f"Manufacturer:   {self.manufacturer}\n" +
             f"Matrix:         {self.matrix}\n" +
             f"Fiber:          {self.matrix}\n" +
             "\n" +
             "Thermoelastic properties\n" +
             "-------------------------------------\n" +
             f"E1:             {self.E1/1E9:.1f} GPa\n" +
             f"E2:             {self.E2/1E9:.1f} GPa\n" +
             f"G12:            {self.G12/1E9:.1f} GPa\n" +
             f"v12:            {self.v12:.2f}\n" +
             f"alpha1:         {self.alpha1*1E6:.2f} um/m\n" +
             f"alpha2:         {self.alpha2*1E6:.2f} um/m\n" +
             "\n" +
             "Strengths\n" +
             "-------------------------------------\n" +
             f"S1t:            {self.S1t/1E6:.1f} MPa\n" +
             f"S1c:            {self.S1c/1E6:.1f} MPa\n" +
             f"S2t:            {self.S2t/1E6:.1f} MPa\n" +
             f"S2c:            {self.S2c/1E6:.1f} MPa\n" +
             f"S6:             {self.S6/1E6:.1f} MPa\n")
        return s


TC1200 = Material(fname='materials/TC1200UD.json')


class Ply():

    def __init__(self, material, theta, t):
        self.mat = material
        self.theta = theta
        self.t = t

    def transformation_matrix(self):
        n = np.sin(self.theta)
        m = np.cos(self.theta)
        T = np.array([[m**2, n**2,     2*n*m],
                      [n**2, m**2,    -2*n*m],
                      [-m*n,  m*n, m**2-n**2]])
        return T

    def C(self):
        T = self.transformation_matrix()
        R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
        C = self.mat.C()
        return inv(T)@C@R@T@inv(R)

    def alpha(self):
        T = self.transformation_matrix()
        R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
        alpha = self.mat.alpha()
        return R@inv(T)@inv(R)@alpha

    def __str__(self):
        pass


P = Ply(TC1200, np.pi/2, 0.15E-3)


class Laminate():
    pass


class Loads():
    pass
