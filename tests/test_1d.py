from .. import ellastic_collision_1d as ec_1d
import pytest

def test_advection_1d_1():
    x_f = ec_1d.advection_1d(x_i=0, v=2, dt=1)
    assert x_f == 2

def test_advection_1d_2():
    x_f = ec_1d.advection_1d(x_i=0, v=2, dt=2)
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

@pytest.mark.parametrize("dt", [1, 0.01])
def test_simulation_1(dt):
    # initial condition and simulation parameters
    domain_x = [-2,12]
    dt = dt
    t_max = 5
    t = 0
    loc_0 = [0, 10]
    vel_0 = [1, -1]
    radius = 1

    # create movie
    movie = ec_1d.Movie(dt, t_max - dt, loc_0, vel_0, domain_x, radius)                             
    movie.animate("pytest_movie_dt_"+str(dt)) 

    # run the simulation
    loc = loc_0
    vel = vel_0
    while(t<t_max):
        loc, vel = ec_1d.simulation_step(dt, loc[0], loc[1], vel[0], vel[1], domain_x, radius)
        t += dt

    # test location and velocities after colision
    if dt == 1:
        assert (loc[0], loc[1]) == (5, 5)
    assert vel[0] == -1
    assert vel[1] == 1

@pytest.mark.skip(reason=" this will be infinite loop, del_x should be better defined")
def test_simulation_collision_2():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=3, v2=-1)


@pytest.mark.skip(reason=" this will be infinite loop, we're not checking if they can meet")
def test_simulation_collision_3():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=1, v2=1)
