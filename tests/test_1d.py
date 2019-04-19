from .. import ellastic_collision_1d as ec_1d
import pytest

def test_position_1():
    x_f = ec_1d.position(x_i=0, v=2)
    assert x_f == 2

def test_position_2():
    x_f = ec_1d.position(x_i=0, v=2, dt=2)
    assert x_f == 4


def test_collision_1():
    v1_f, v2_f = ec_1d.collision_1d(v1_i=1, v2_i=-2)
    assert v1_f == -2
    assert v2_f == 1


def test_collision_2():
    v1_f, v2_f = ec_1d.collision_1d(v1_i=1, v2_i=-2, m1=2, m2=2)
    assert v1_f == -2
    assert v2_f == 1

def test_collision_3():
    v1_f, v2_f = ec_1d.collision_1d(v1_i=1, v2_i=-2, m1=1, m2=1e6)
    assert v2_f == pytest.approx(-2, rel=1e-3)


def test_simulation_collision_1():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=1, v2=-1)

    assert (x1_f, x2_f) == (5, 5)
    assert time_f == 5
    assert v1_f == -1
    assert v2_f == 1


@pytest.mark.skip(reason=" this will be infinite loop, del_x should be better defined")
def test_simulation_collision_2():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=3, v2=-1)


@pytest.mark.skip(reason=" this will be infinite loop, we're not checking if they can meet")
def test_simulation_collision_3():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=1, v2=1)
