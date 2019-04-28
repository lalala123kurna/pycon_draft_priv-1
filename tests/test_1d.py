#from .. import ellastic_collision_1d as ec_1d
from .. import collisions as cl
from .. import ellastic_collision_2d_many as ec
from .. import movies as mv

import pytest
import numpy as np

def test_collision_1d_1():
    v1_f, v2_f = cl.collision_1d(v1_i=1, v2_i=-2)
    assert v1_f == -2
    assert v2_f == 1

def test_collision_1d_2():
    v1_f, v2_f = cl.collision_1d(v1_i=1, v2_i=-2, m1=2, m2=2)
    assert v1_f == -2
    assert v2_f == 1

def test_collision_1d_3():
    v1_f, v2_f = cl.collision_1d(v1_i=1, v2_i=-2, m1=1, m2=1e6)
    assert v2_f == pytest.approx(-2, rel=1e-3)

@pytest.mark.parametrize("dt", [1])#, 0.01])
def test_simulation_1d(dt):
    # initial condition and simulation parameters
    domain = ([-2, 12], [0, 3])
    dt = dt
    t_max = 6
    t = 0
    loc_0 = np.array([[0, 1.5],[10, 1.5]])
    vel_0 = np.array([[1, 0], [-1, 0]])
    radius = 1
    mass = [1, 1]

    loc = np.copy(loc_0)
    vel = np.copy(vel_0)
    # create movie
    movie = mv.Movie_2d(ec.simulation_step, dt, t_max - dt, loc, vel, domain, mass, radius)                             
    movie.animate("pytest_movie_1d_dt_"+str(dt)) 

    # run the simulation
    loc = np.copy(loc_0)
    vel = np.copy(vel_0)
    while(t<t_max):
        print("loc", loc)
        print("vel", vel)
        loc, vel = ec.simulation_step(dt, mass, radius, loc, vel, domain)
        print("loc", loc)
        print("vel", vel)
        t += dt

    # test location and velocities after colision
    if dt == 1:
        assert (loc[0][0], loc[1][0]) == (5, 5)
    assert vel[0][0] == -1
    assert vel[1][0] == 1

@pytest.mark.skip(reason=" todo ")
def test_energy_concervation():

    domain_x = [-2,12]
    dt = 1
    t_max = 5
    t = 0
    loc_0 = [0, 10]
    vel_0 = [1, -1]
    radius = 1

    E_ini = 0.5 * (pow(vel_0[0],2) + pow(vel_0[1],2))

    # run the simulation
    loc = loc_0
    vel = vel_0
    while(t<t_max):
        loc, vel = ec_1d.simulation_step(dt, loc[0], loc[1], vel[0], vel[1], domain_x, radius)
        t += dt

    E_fin = 0.5 * (pow(vel[0],2) + pow(vel[1],2))

    assert E_ini == E_fin
 

@pytest.mark.skip(reason=" this will be infinite loop, del_x should be better defined")
def test_simulation_collision_2():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=3, v2=-1)


@pytest.mark.skip(reason=" this will be infinite loop, we're not checking if they can meet")
def test_simulation_collision_3():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=1, v2=1)
