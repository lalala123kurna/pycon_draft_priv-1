from .. import ellastic_collision_2d_angle as ec_2d_angle
from .. import ellastic_collision_2d as ec_2d

import pytest

#something wrong
def test_collision_1():
    v1_f, v2_f, beta = ec_2d_angle.collision_2d_rest(v1_i=1, alpha=0)
    assert v1_f == 0
    assert v2_f == 1

#idea for tests collision_2d
# - checking for alpha=0
# - checking for alpha 90
# - checking alpha's unit
