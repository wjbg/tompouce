"""Tests for failure criteria.

Based on `pytest-7.2.1.` Tests can be run using:

```
$ pytest -vv
```

from the prompt.

"""

import numpy as np
from tompouce import Material
import tompouce.failure_criteria as fc


class TestFailureCriteria:
    """Tests for tompouce Material object."""
    material = Material('test_file.json')
    s0 = np.array([0.0, 0.0, 0.0])
    s1c, s1t = np.array([-600E6, 0.0, 0.0]), np.array([1.1E9, 0.0, 0.0])
    s2c, s2t = np.array([0.0, -110E6, 0.0]), np.array([0.0, 60E6, 0.0])
    s6 = np.array([0.0, 0.0, -60E6])

    def test_max_stress(self):
        assert ((not self.material.failure(self.s0, fc.max_stress)) and
                self.material.failure(self.s1c, fc.max_stress) and
                self.material.failure(self.s1t, fc.max_stress) and
                self.material.failure(self.s2c, fc.max_stress) and
                self.material.failure(self.s2t, fc.max_stress) and
                self.material.failure(self.s6, fc.max_stress))

    def test_Tsai_Hill(self):
        assert ((not self.material.failure(self.s0, fc.Tsai_Hill)) and
                self.material.failure(self.s1c, fc.Tsai_Hill) and
                self.material.failure(self.s1t, fc.Tsai_Hill) and
                self.material.failure(self.s2c, fc.Tsai_Hill) and
                self.material.failure(self.s2t, fc.Tsai_Hill) and
                self.material.failure(self.s6, fc.Tsai_Hill))

    def test_Norris(self):
        assert ((not self.material.failure(self.s0, fc.Norris)) and
                self.material.failure(self.s1c, fc.Norris) and
                self.material.failure(self.s1t, fc.Norris) and
                self.material.failure(self.s2c, fc.Norris) and
                self.material.failure(self.s2t, fc.Norris) and
                self.material.failure(self.s6, fc.Norris))

    def test_Tsai_Wu(self):
        assert ((not self.material.failure(self.s0, fc.Tsai_Wu)) and
                self.material.failure(self.s1c, fc.Tsai_Wu) and
                self.material.failure(self.s1t, fc.Tsai_Wu) and
                self.material.failure(self.s2c, fc.Tsai_Wu) and
                self.material.failure(self.s2t, fc.Tsai_Wu) and
                self.material.failure(self.s6, fc.Tsai_Wu))

    def test_Hoffman(self):
        assert ((not self.material.failure(self.s0, fc.Hoffman)) and
                self.material.failure(self.s1c, fc.Hoffman) and
                self.material.failure(self.s1t, fc.Hoffman) and
                self.material.failure(self.s2c, fc.Hoffman) and
                self.material.failure(self.s2t, fc.Hoffman) and
                self.material.failure(self.s6, fc.Hoffman))

    def test_Hashin_Rotem(self):
        assert ((not self.material.failure(self.s0, fc.Hashin_Rotem)) and
                self.material.failure(self.s1c, fc.Hashin_Rotem) and
                self.material.failure(self.s1t, fc.Hashin_Rotem) and
                self.material.failure(self.s2c, fc.Hashin_Rotem) and
                self.material.failure(self.s2t, fc.Hashin_Rotem) and
                self.material.failure(self.s6, fc.Hashin_Rotem))
